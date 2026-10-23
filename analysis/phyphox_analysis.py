"""
phyphox_analysis.py — turn phone accelerometer recordings from the Bumpy Road kit into
natural frequency, damping ratio and transmissibility T(f), and compare with theory.

Record with the free phyphox app, experiment "Acceleration (without g)", phone flat in the
tray. Export: menu > Export Data > CSV (comma, decimal point) or Excel-style (semicolon,
decimal comma). The .zip that phyphox makes can be passed directly.

    python phyphox_analysis.py pluck  pluck.csv  [--mass 0.25]          # Task 2: fn, zeta
    python phyphox_analysis.py sweep  sweep.csv  [--fn 6.6 --zeta 0.08] # Task 3: T(f)
    python phyphox_analysis.py base   base.csv                          # calibrate cam amplitude Y
    python phyphox_analysis.py compare g15.csv g3.csv g1.csv --labels "15 mm" "3 mm" "1 mm"
    python phyphox_analysis.py analyse data/my_runs --mass 0.27       # EASIEST: whole folder, see below

`analyse` matches files by name: base.csv (optional), pluck_<tag>.csv, sweep_<tag>.csv
(tag = your set-up, e.g. 15mm, 3mm, 1mm). Each sweep is plotted against the theory from the
pluck with the same tag. Results (figures, summary.csv, notes.txt) go to <folder>/results.
    python phyphox_analysis.py demo                                     # synthetic data, no rig needed

THE RIG MODEL (not the textbook one). The lower springs (k_l) join the cam carriage to
the platform; the upper springs (k_u) and the eddy-current damper join the platform to the
FIXED frame (bridge, damper upright). With platform motion x and road motion y:

    m x'' + c x' + (k_l + k_u) x = k_l y
    T = X/Y = alpha / sqrt((1 - r^2)^2 + (2 zeta r)^2),   alpha = k_l / (k_l + k_u)

with fn = (1/2pi) sqrt(k/m), k = k_l + k_u, zeta = c / (2 sqrt(k m)), r = f/fn.
Four springs above and four below (as built) give alpha = 0.5: at low frequency the platform moves half
as much as the road. Because the damper is grounded ("skyhook"), more damping lowers T at
EVERY frequency, unlike the textbook base-excited system (damper to the road), where
damping worsens isolation above sqrt(2) fn. The textbook curve is plotted for comparison.

T comes from the phone alone: the road input is the cam, y = Y sin(2 pi f t) with
Y = 2 mm, so the road acceleration is (2 pi f)^2 Y. Check Y with the `base` command
(phone taped to the carriage) if in doubt.
"""

import argparse
import csv
import io
import math
import os
import re
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DESIGNATHON = os.path.abspath(os.path.join(HERE, ".."))
END_STOP_MM = 21.0                   # platform travel to the end stops (manufacturing/params.py)
PREDICTED = {"fn": 9.44, "zeta": 0.07, "alpha": 0.5}   # hand_calcs.py (k = 880 N/m: 8 springs, m = 250 g)


# --- models ---------------------------------------------------------------------
def transmissibility_rig(f, fn, zeta, alpha=0.5):
    """Platform/road displacement ratio for THIS rig (upper springs and damper grounded)."""
    r = np.asarray(f, dtype=float) / fn
    return alpha / np.sqrt((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2)


def transmissibility_textbook(f, fn, zeta):
    """Textbook base excitation: spring AND damper between road and mass."""
    r = np.asarray(f, dtype=float) / fn
    return np.sqrt((1 + (2 * zeta * r) ** 2) / ((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2))


# --- loading ----------------------------------------------------------------------
def load(path):
    """Read a phyphox export (.csv or .zip). Returns (t, {axis: values}, header)."""
    if path.lower().endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.lower().endswith(".csv")]
            name = next((n for n in names if "raw data" in n.lower()), names[0])
            text = z.read(name).decode("utf-8-sig")
    else:
        with open(path, encoding="utf-8-sig") as f:
            text = f.read()
    first = text.splitlines()[0]
    delim = "\t" if "\t" in first else (";" if ";" in first else ",")
    rows = list(csv.reader(io.StringIO(text), delimiter=delim))
    header = [h.strip() for h in rows[0]]
    data = []
    for r in rows[1:]:
        if not r or all(not c.strip() for c in r):
            continue
        data.append([float(c.strip().replace(",", ".")) if c.strip() else np.nan for c in r])
    data = np.array(data, dtype=float)
    axes = {}
    for j, name in enumerate(header[1:], start=1):
        low = name.lower()
        if "absolute" in low:
            continue
        m = re.search(r"\b([xyz])\b", low)
        if m:
            axes[m.group(1)] = data[:, j]
    if not axes:
        raise ValueError(f"no x/y/z acceleration columns in {path}: {header}")
    return data[:, 0], axes, header


def choose_axis(axes, axis="auto"):
    if axis != "auto":
        return axis, axes[axis]
    best = max(axes, key=lambda k: np.nanvar(axes[k]))
    return best, axes[best]


def resample(t, a):
    """phyphox timestamps are not perfectly uniform: interpolate onto a uniform grid."""
    ok = np.isfinite(t) & np.isfinite(a)
    t, a = t[ok], a[ok]
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, a), 1.0 / dt


def bandpass(a, fs, lo, hi):
    from scipy.signal import butter, sosfiltfilt
    hi = min(hi, 0.45 * fs)
    sos = butter(4, [lo, hi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, a - np.mean(a))


def dominant_freq(x, fs, band):
    """Peak of the Hann-windowed, zero-padded spectrum inside band, parabolic-refined."""
    n = len(x)
    nfft = 8 * 2 ** int(math.ceil(math.log2(max(n, 2))))
    spec = np.abs(np.fft.rfft((x - np.mean(x)) * np.hanning(n), nfft))
    f = np.fft.rfftfreq(nfft, 1.0 / fs)
    sel = np.where((f >= band[0]) & (f <= band[1]))[0]
    i = sel[np.argmax(spec[sel])]
    if 0 < i < len(spec) - 1:
        y0, y1, y2 = spec[i - 1], spec[i], spec[i + 1]
        den = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        return float(f[i] + shift * (f[1] - f[0]))
    return float(f[i])


def _r2(y, yhat):
    ss = np.sum((y - np.mean(y)) ** 2)
    return float(1 - np.sum((y - yhat) ** 2) / ss) if ss > 0 else 0.0


# --- Task 2: pluck test -------------------------------------------------------------
def analyse_pluck(t, a, mass=0.25, band=(1.5, 30.0), skip=1):
    """Free decay after a pluck. Log decrement (viscous) and linear decrement (Coulomb
    friction) are both fitted to the positive peaks; the better fit is reported."""
    from scipy.signal import find_peaks
    tu, au, fs = resample(t, a)
    af = bandpass(au, fs, *band)
    i0 = int(np.argmax(np.abs(af)))
    quiet = [af[int(0.85 * len(af)):]]
    if i0 > int(0.5 * fs):
        quiet.append(af[: i0 - int(0.2 * fs)])
    noise = min(float(np.std(q)) for q in quiet if len(q) > 10)
    seg, tseg = af[i0:], tu[i0:]
    fd0 = dominant_freq(seg[: min(len(seg), int(3 * fs))], fs, band)
    pk, _ = find_peaks(seg, distance=max(1, int(0.6 * fs / fd0)))
    thr = max(4 * noise, 0.01 * float(np.max(np.abs(seg))))
    def above(pks):
        out = []
        for p in pks:
            if seg[p] < thr:
                break
            out.append(p)
        return out
    kept = above(pk[skip:])
    if len(kept) < 3:                       # heavy damping: keep the first peak too
        kept = above(pk)
    if len(kept) < 2:
        raise ValueError(f"only {len(kept)} usable peaks: pluck harder, record longer, or open the magnet gap")
    kept = np.array(kept)
    # sub-sample peak times (parabola through 3 samples): phones sample at only 100-400 Hz
    dtp = np.zeros(len(kept))
    for j, p in enumerate(kept):
        if 0 < p < len(seg) - 1:
            y0, y1, y2 = seg[p - 1], seg[p], seg[p + 1]
            den = y0 - 2 * y1 + y2
            dtp[j] = 0.5 * (y0 - y2) / den / fs if den != 0 else 0.0
    tp, A = tseg[kept] + dtp, seg[kept]
    n = np.arange(len(A))
    sv, iv = np.polyfit(n, np.log(A), 1)
    sc, ic = np.polyfit(n, A, 1)
    if len(A) >= 3:
        r2v, r2c = _r2(np.log(A), sv * n + iv), _r2(A, sc * n + ic)
    else:
        r2v, r2c = math.nan, math.nan
    delta = -sv
    zeta_ld = delta / math.sqrt(4 * math.pi ** 2 + delta ** 2)
    fd_ld = (len(tp) - 1) / (tp[-1] - tp[0])
    # time-domain fit of a damped cosine: robust when heavy damping leaves few peaks
    from scipy.optimize import curve_fit
    sel = (tseg >= tp[0]) & (tseg <= tp[-1] + 1.0 / fd_ld)
    tt, yy = tseg[sel] - tp[0], seg[sel]
    model = lambda t, A0, s, f, ph: A0 * np.exp(-s * t) * np.cos(2 * np.pi * f * t + ph)   # noqa: E731
    try:
        (A0, s, fdf, _), _ = curve_fit(model, tt, yy, p0=[A[0], delta * fd_ld, fd_ld, 0.0], maxfev=20000)
        wdf = 2 * np.pi * abs(fdf)
        zeta_fit = s / math.hypot(s, wdf)
        fd_fit = abs(fdf)
    except RuntimeError:
        zeta_fit, fd_fit = math.nan, math.nan
    use_fit = len(A) < 5 and np.isfinite(zeta_fit)
    zeta = zeta_fit if use_fit else zeta_ld
    fd = fd_fit if use_fit else fd_ld
    fn = fd / math.sqrt(1 - zeta ** 2)
    wn = 2 * math.pi * fn
    k = mass * wn ** 2
    friction = k * (-sc / wn ** 2) / 4                     # amplitude falls 4F/k per cycle
    return {"fn_Hz": fn, "fd_Hz": fd, "zeta": zeta, "log_decrement": delta,
            "zeta_log_decrement": zeta_ld, "zeta_curve_fit": zeta_fit,
            "zeta_method": "curve fit (few peaks)" if use_fit else "log decrement",
            "k_total_N_per_m": k, "c_N_s_per_m": 2 * zeta * math.sqrt(k * mass),
            "friction_N_if_coulomb": friction, "r2_viscous": r2v, "r2_coulomb": r2c,
            "damping_model": ("undetermined (too few peaks)" if not np.isfinite(r2v) else
                              "viscous (exponential decay)" if r2v >= r2c else "Coulomb friction (linear decay)"),
            "n_peaks": int(len(A)), "first_amplitude_mm": 1e3 * A[0] / wn ** 2,
            "noise_m_s2": noise, "mass_kg": mass,
            "_plot": {"t": tseg, "a": seg, "tp": tp, "A": A, "fit_v": (sv, iv), "fit_c": (sc, ic)}}


# --- Task 3: frequency sweep ----------------------------------------------------------
def _windows(t, a, win, step, fband):
    """Per window: dominant road frequency f, platform DISPLACEMENT amplitude and purity.
    Displacement comes from the acceleration spectrum divided by w^2 over 0.5f-1.5f, so the
    sidebands from uneven hand cranking (+-crank rate) are counted, not lost. Purity is the
    fraction of acceleration power in that band (knocks, silence, fast speed changes: low)."""
    tu, au, fs = resample(t, a)
    af = bandpass(au, fs, 1.0, 40.0)
    n, hop = int(win * fs), max(1, int(step * fs))
    w = np.hanning(n)
    w2 = np.mean(w ** 2)
    fk = np.fft.rfftfreq(n, 1.0 / fs)
    total_band = (fk >= 1.5) & (fk <= 25.0)
    out = []
    for s in range(0, len(af) - n + 1, hop):
        seg = af[s:s + n]
        f = dominant_freq(seg, fs, fband)
        X = np.fft.rfft((seg - seg.mean()) * w)
        near = (fk >= 0.5 * f) & (fk <= 1.5 * f)
        p_near = np.sum(np.abs(X[near]) ** 2)
        p_tot = np.sum(np.abs(X[total_band]) ** 2)
        D = X[near] / (2 * np.pi * fk[near]) ** 2                  # displacement spectrum
        rms = math.sqrt(2 * np.sum(np.abs(D) ** 2) / n ** 2 / w2)
        out.append((tu[s] + win / 2, f, math.sqrt(2) * rms, p_near / p_tot if p_tot > 0 else 0.0))
    return np.array(out)


def analyse_sweep(t, a, Y_mm=2.0, win=2.0, step=0.5, fband=(2.0, 16.0), min_fit=0.7, bin_hz=0.25):
    """Short windows: road frequency f and platform displacement amplitude X; T = X / Y
    (the cam fixes the road DISPLACEMENT amplitude Y whatever the crank speed). Windows with
    purity < min_fit (knocks, silence, fast speed changes) are discarded. Fits fn, zeta and
    alpha of the rig model to the binned T(f)."""
    from scipy.optimize import curve_fit
    w = _windows(t, a, win, step, fband)
    tc, f, X, fit = w.T
    ok = fit >= min_fit
    T = X / (Y_mm * 1e-3)
    disp_mm = 1e3 * X
    edges = np.arange(fband[0], fband[1] + bin_hz, bin_hz)
    fb, Tb, nb = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = ok & (f >= lo) & (f < hi)
        if sel.any():
            fb.append(float(np.median(f[sel])))
            Tb.append(float(np.median(T[sel])))
            nb.append(int(sel.sum()))
    fb, Tb = np.array(fb), np.array(Tb)
    res = {"Y_mm": Y_mm, "n_windows": int(len(f)), "n_used": int(ok.sum()),
           "f_range_Hz": (float(f[ok].min()), float(f[ok].max())) if ok.any() else (math.nan, math.nan),
           "max_disp_mm": float(disp_mm[ok].max()) if ok.any() else math.nan,
           "bins": {"f_Hz": fb, "T": Tb, "n": np.array(nb)},
           "_plot": {"t": tc, "f": f, "T": T, "ok": ok}}
    res["end_stop_warning"] = bool(res["max_disp_mm"] > 0.8 * END_STOP_MM)
    if len(fb) >= 5:
        p0 = [float(fb[np.argmax(Tb)]), 0.1, 0.5]
        try:
            p, cov = curve_fit(transmissibility_rig, fb, Tb, p0=p0,
                               bounds=([1.0, 0.005, 0.05], [30.0, 1.0, 1.5]))
            sd = np.sqrt(np.diag(cov))
            res["fit"] = {"fn_Hz": p[0], "zeta": p[1], "alpha": p[2],
                          "sd": {"fn_Hz": sd[0], "zeta": sd[1], "alpha": sd[2]}}
        except RuntimeError:
            res["fit"] = None
    else:
        res["fit"] = None
    res["peak_T"] = float(Tb.max()) if len(Tb) else math.nan
    res["f_at_peak_Hz"] = float(fb[np.argmax(Tb)]) if len(Tb) else math.nan
    return res


def analyse_base(t, a, win=2.0, step=0.5, fband=(2.0, 16.0), min_fit=0.7):
    """Phone on the CARRIAGE: measured road displacement amplitude Y (expect 2 mm)."""
    w = _windows(t, a, win, step, fband)
    ok = w[:, 3] >= min_fit
    Y = 1e3 * w[ok, 2]
    return {"Y_mm_median": float(np.median(Y)), "Y_mm_iqr": tuple(np.percentile(Y, [25, 75])),
            "n_used": int(ok.sum())}


# --- synthetic recordings (no rig needed) -------------------------------------------
def simulate(kind, fn=6.68, zeta=0.08, alpha=0.5, friction_N=0.0, mass=0.25, Y_mm=2.0,
             fs=200.0, noise=0.03, seed=0, lobes=4, rps=(0.8, 3.0), duration=None, X0_mm=6.0,
             jitter=0.0):
    """Integrate the rig model and return (t, acc_z) as a phone would record it.
    kind = 'pluck' (released from X0 at t = 1 s) or 'sweep' (crank speed ramps over rps).
    jitter: crank speed varies by +-jitter within each revolution (hand cranking is uneven)."""
    from scipy.integrate import solve_ivp
    wn = 2 * math.pi * fn
    k = mass * wn ** 2
    c = 2 * zeta * math.sqrt(k * mass)
    Y = Y_mm * 1e-3
    if kind == "pluck":
        T_end = duration or 8.0
        t0 = 1.0
        y = lambda t: 0.0                                           # noqa: E731
        x_init = [X0_mm * 1e-3, 0.0]
    elif kind == "sweep":
        T_end = duration or 60.0
        t0 = 0.0
        r0, r1 = rps
        th0 = lambda t: 2 * math.pi * (r0 * t + 0.5 * (r1 - r0) * t ** 2 / T_end)     # noqa: E731
        theta = lambda t: th0(t) + jitter * math.sin(th0(t))                          # noqa: E731
        y = lambda t: Y * math.sin(lobes * theta(t))                                 # noqa: E731
        x_init = [0.0, 0.0]
    else:
        raise ValueError(kind)

    def rhs(t, s):
        x, v = s
        fr = friction_N * math.tanh(v / 1e-3)
        return [v, (-c * v - k * x + alpha * k * y(t) - fr) / mass]

    t_eval = np.arange(t0, T_end, 1.0 / fs)
    sol = solve_ivp(rhs, (t0, T_end), x_init, t_eval=t_eval, max_step=0.5 / fs, rtol=1e-7, atol=1e-10)
    acc = np.array([rhs(ti, si)[1] for ti, si in zip(sol.t, sol.y.T)])
    t = sol.t
    if kind == "pluck":                                             # held still before release
        pre = np.arange(0.0, t0, 1.0 / fs)
        t, acc = np.concatenate([pre, t]), np.concatenate([np.zeros_like(pre), acc])
    rng = np.random.default_rng(seed)
    return t, acc + noise * rng.standard_normal(len(acc))


def write_phyphox_csv(path, t, az, noise=0.03, seed=1):
    """Same columns as phyphox 'Acceleration (without g)' export."""
    rng = np.random.default_rng(seed)
    ax, ay = noise * rng.standard_normal(len(t)), noise * rng.standard_normal(len(t))
    with open(path, "w", newline="") as f:
        w = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        w.writerow(["Time (s)", "Linear Acceleration x (m/s^2)", "Linear Acceleration y (m/s^2)",
                    "Linear Acceleration z (m/s^2)", "Absolute acceleration (m/s^2)"])
        for row in zip(t, ax, ay, az, np.sqrt(ax ** 2 + ay ** 2 + az ** 2)):
            w.writerow([float(f"{v:.6g}") for v in row])


# --- plots ------------------------------------------------------------------------------
def plot_pluck(r, out, title="Pluck test"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p = r["_plot"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [2, 1]})
    a1.plot(p["t"], p["a"], lw=0.8, color="0.3", label="platform acceleration (filtered)")
    a1.plot(p["tp"], p["A"], "o", color="C3", ms=4, label="peaks used")
    n = np.arange(len(p["A"]))
    sv, iv = p["fit_v"]
    sc, ic = p["fit_c"]
    a1.plot(p["tp"], np.exp(sv * n + iv), color="C0", label=f"viscous fit, R² = {r['r2_viscous']:.3f}")
    a1.plot(p["tp"], sc * n + ic, "--", color="C1", label=f"Coulomb fit, R² = {r['r2_coulomb']:.3f}")
    a1.set_xlabel("time (s)")
    a1.set_ylabel("acceleration (m/s²)")
    a1.set_xlim(p["tp"][0] - 0.3, p["tp"][-1] + 0.5)
    a1.legend(fontsize=8)
    a2.plot(n, np.log(p["A"]), "o", color="C3")
    a2.plot(n, sv * n + iv, color="C0")
    a2.set_xlabel("cycle number n")
    a2.set_ylabel("ln(peak)")
    a2.set_title(f"slope = −δ = {sv:.3f}", fontsize=9)
    fig.suptitle(f"{title}:  fn = {r['fn_Hz']:.2f} Hz,  ζ = {r['zeta']:.3f},  "
                 f"best fit: {r['damping_model']}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_sweep(r, out, theory=None, title="Frequency sweep"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p, b = r["_plot"], r["bins"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), gridspec_kw={"height_ratios": [3, 1]})
    ff = np.linspace(1.0, 16.0, 400)
    a1.plot(p["f"][p["ok"]], p["T"][p["ok"]], ".", color="0.65", ms=4, label="2 s windows")
    a1.plot(b["f_Hz"], b["T"], "o", color="C3", label="measured (median per 0.25 Hz)")
    th = theory or PREDICTED
    src = "from pluck test" if theory else "hand-calc prediction"
    a1.plot(ff, transmissibility_rig(ff, th["fn"], th["zeta"], th.get("alpha", 0.5)), color="C0",
            label=f"rig model {src}: fn {th['fn']:.2f} Hz, ζ {th['zeta']:.2f}, α {th.get('alpha', 0.5):.2f}")
    if r.get("fit"):
        fp = r["fit"]
        a1.plot(ff, transmissibility_rig(ff, fp["fn_Hz"], fp["zeta"], fp["alpha"]), "--", color="C2",
                label=f"rig model fitted to sweep: fn {fp['fn_Hz']:.2f} Hz, ζ {fp['zeta']:.2f}, α {fp['alpha']:.2f}")
    a1.plot(ff, transmissibility_textbook(ff, th["fn"], th["zeta"]), ":", color="0.5",
            label="textbook (damper to road), NOT this rig")
    a1.axvline(th["fn"], color="0.8", lw=0.8)
    a1.axvline(math.sqrt(2) * th["fn"], color="0.8", lw=0.8, ls="--")
    a1.text(math.sqrt(2) * th["fn"], a1.get_ylim()[1] * 0.9 if a1.get_ylim()[1] < 20 else 5, " √2·fn", fontsize=8, color="0.5")
    a1.axhline(1.0, color="0.85", lw=0.8)
    a1.set_yscale("log")
    a1.set_xlim(1, 16)
    a1.set_ylabel("transmissibility T = X / Y")
    a1.legend(fontsize=7, loc="lower left")
    a1.set_title(title, fontsize=10)
    a2.plot(p["t"][~p["ok"]], p["f"][~p["ok"]], "x", color="0.7", ms=4, label="discarded")
    a2.plot(p["t"][p["ok"]], p["f"][p["ok"]], ".", color="C3", ms=4, label="used")
    a2.set_xlabel("time (s)")
    a2.set_ylabel("road freq. (Hz)")
    a2.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)



def compare_sweeps(files, labels=None, Y_mm=2.0, axis="auto"):
    """Analyse several sweeps (e.g. three magnet gaps) and return [(label, result), ...]."""
    out = []
    for i, f in enumerate(files):
        t, axes, _ = load(f)
        r = analyse_sweep(t, choose_axis(axes, axis)[1], Y_mm=Y_mm)
        out.append(((labels[i] if labels and i < len(labels) else os.path.basename(f)), r))
    return out


def plot_compare(results, out, title="Transmissibility vs road frequency", theories=None):
    """Measured T(f) for each run (points) with its fitted rig model (line), plus a summary table.
    theories: optional {label: {"fn", "zeta", "alpha"}} drawn dotted (prediction from the pluck test)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, axt) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [2.2, 1]})
    ff = np.linspace(1.0, 16.0, 400)
    rows = []
    for i, (lab, r) in enumerate(results):
        col = f"C{i}"
        b = r["bins"]
        ax.plot(b["f_Hz"], b["T"], "o", color=col, ms=4, label=f"{lab}: measured")
        fp = r.get("fit")
        if fp:
            ax.plot(ff, transmissibility_rig(ff, fp["fn_Hz"], fp["zeta"], fp["alpha"]), "-", color=col, lw=1.4,
                    label=f"{lab}: model fit")
            rows.append([lab, f"{fp['fn_Hz']:.2f}", f"{fp['zeta']:.3f}", f"{fp['alpha']:.2f}", f"{r['peak_T']:.2f}"])
        else:
            rows.append([lab, "-", "-", "-", f"{r['peak_T']:.2f}"])
        th = (theories or {}).get(lab)
        if th:
            ax.plot(ff, transmissibility_rig(ff, th["fn"], th["zeta"], th.get("alpha", 0.5)), ":", color=col, lw=1.4,
                    label=f"{lab}: predicted from pluck")
    ax.axhline(1.0, color="0.85", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlim(1, 16)
    ax.set_xlabel("road frequency f (Hz)  [= 4 x crank rev/s]")
    ax.set_ylabel("transmissibility T = X / Y")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7)
    axt.axis("off")
    tab = axt.table(cellText=rows, colLabels=["run", "fn (Hz)", "zeta", "alpha", "peak T"], loc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1, 1.6)
    axt.set_title("fitted to each sweep", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return rows



# --- one command for a folder of recordings ------------------------------------------
def _tag(path, prefix):
    stem = os.path.splitext(os.path.basename(path))[0]
    t = stem[len(prefix):].lstrip("_- ")
    return t or "run"


def _tag_key(tag):
    """Set-ups in a natural order: numeric tags largest first (15mm, 3mm, 1mm), then names."""
    m = re.search(r"[-+]?\d*\.?\d+", tag)
    return (0, -float(m.group())) if m else (1, tag)


def analyse_folder(folder, mass=0.25, Y_mm=None, alpha=0.5, out=None, axis="auto"):
    """Analyse every recording in `folder`, matched by file name:
        base*.csv            phone taped to the carriage -> measured road amplitude Y (optional)
        pluck_<tag>.csv      e.g. pluck_15mm.csv -> fn, zeta for that set-up
        sweep_<tag>.csv      e.g. sweep_15mm.csv -> T(f); theory curve uses pluck_<tag>'s fn, zeta
    (.zip exports from phyphox work too.) Writes figures and summary.csv to out
    (default <folder>/results) and returns the summary rows."""
    import csv as _csv
    import glob
    out = out or os.path.join(folder, "results")
    os.makedirs(out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(folder, "*.csv")) + glob.glob(os.path.join(folder, "*.zip")))
    pick = lambda pre: [f for f in files if os.path.basename(f).lower().startswith(pre)]   # noqa: E731
    notes = []
    # 1. road amplitude
    Y_used, Y_src = (Y_mm, "given") if Y_mm else (2.0, "nominal cam (no base file)")
    if not Y_mm and pick("base"):
        t, axes, _ = load(pick("base")[0])
        Yb = analyse_base(t, choose_axis(axes, axis)[1])
        if Yb["n_used"] >= 3:
            Y_used, Y_src = Yb["Y_mm_median"], f"measured ({os.path.basename(pick('base')[0])})"
        else:
            notes.append("base file had too few clean windows: using the nominal 2 mm")
    # 2. plucks
    plucks = {}
    for f in pick("pluck"):
        tag = _tag(f, "pluck")
        try:
            t, axes, _ = load(f)
            r = analyse_pluck(t, choose_axis(axes, axis)[1], mass=mass)
            plot_pluck(r, os.path.join(out, f"pluck_{tag}.png"), f"Pluck test ({tag})")
            plucks[tag] = r
        except Exception as e:                                       # noqa: BLE001  keep going, report it
            notes.append(f"pluck {tag}: {e}")
    # 3. sweeps, each against the pluck with the same tag (else the hand-calc prediction)
    sweeps, theories = [], {}
    rows = []
    for f in pick("sweep"):
        tag = _tag(f, "sweep")
        try:
            t, axes, _ = load(f)
            r = analyse_sweep(t, choose_axis(axes, axis)[1], Y_mm=Y_used)
        except Exception as e:                                       # noqa: BLE001
            notes.append(f"sweep {tag}: {e}")
            continue
        p = plucks.get(tag)
        th = {"fn": p["fn_Hz"], "zeta": p["zeta"], "alpha": alpha} if p else None
        plot_sweep(r, os.path.join(out, f"sweep_{tag}.png"), th, f"Frequency sweep ({tag})")
        if th:
            theories[tag] = th
        else:
            notes.append(f"sweep {tag}: no pluck_{tag} file, theory curve = hand-calc prediction")
        if r["end_stop_warning"]:
            notes.append(f"sweep {tag}: platform near the end stops, peak T clipped")
        sweeps.append((tag, r))
    # 4. summary (one row per tag)
    sweeps.sort(key=lambda tr: _tag_key(tr[0]))
    tags = sorted(set(plucks) | {t for t, _ in sweeps}, key=_tag_key)
    sw = dict(sweeps)
    for tag in tags:
        p, r = plucks.get(tag), sw.get(tag)
        fit = (r or {}).get("fit") or {}
        rows.append({
            "setup": tag,
            "fn_pluck_Hz": round(p["fn_Hz"], 3) if p else "",
            "zeta_pluck": round(p["zeta"], 4) if p else "",
            "decay": p["damping_model"] if p else "",
            "fn_sweep_fit_Hz": round(fit["fn_Hz"], 3) if fit else "",
            "zeta_sweep_fit": round(fit["zeta"], 4) if fit else "",
            "alpha_fit": round(fit["alpha"], 3) if fit else "",
            "peak_T": round(r["peak_T"], 3) if r else "",
            "f_at_peak_Hz": round(r["f_at_peak_Hz"], 2) if r else "",
            "fn_predicted_Hz": PREDICTED["fn"], "Y_mm": round(Y_used, 3), "mass_kg": mass,
        })
    if len(sweeps) >= 2:
        plot_compare(sweeps, os.path.join(out, "compare.png"), "Measured T(f) for each set-up", theories)
    with open(os.path.join(out, "summary.csv"), "w", newline="") as fh:
        if rows:
            w = _csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    with open(os.path.join(out, "notes.txt"), "w") as fh:
        fh.write(f"Y = {Y_used:.3f} mm ({Y_src}); mass = {mass} kg; alpha = {alpha}\n")
        fh.write("\n".join(notes) + ("\n" if notes else ""))
    return rows, notes, Y_used, Y_src

# --- CLI --------------------------------------------------------------------------------
def _print(d, indent="  "):
    for k, v in d.items():
        if k.startswith("_") or k == "bins":
            continue
        if isinstance(v, dict):
            print(f"{indent}{k}:")
            _print(v, indent + "  ")
        elif isinstance(v, float):
            print(f"{indent}{k:24s} {v:.4g}")
        else:
            print(f"{indent}{k:24s} {v}")


def demo(outdir=None, figdir=None):
    outdir = outdir or os.path.join(DESIGNATHON, "data", "synthetic")
    figdir = figdir or os.path.join(DESIGNATHON, "figures")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)
    truth = {"fn": 9.4, "zeta": 0.07, "alpha": 0.5, "friction_N": 0.02, "crank_unevenness": 0.1}
    t, a = simulate("pluck", fn=truth["fn"], zeta=truth["zeta"], friction_N=truth["friction_N"], fs=100)
    pp = os.path.join(outdir, "pluck_synthetic.csv")
    write_phyphox_csv(pp, t, a)
    t, a = simulate("sweep", fn=truth["fn"], zeta=truth["zeta"], friction_N=truth["friction_N"], seed=2,
                    fs=100, duration=45, jitter=0.1)
    sp = os.path.join(outdir, "sweep_synthetic.csv")
    write_phyphox_csv(sp, t, a)
    print(f"SYNTHETIC data (truth: {truth}) written to {outdir}")
    tt, axes, _ = load(pp)
    rp = analyse_pluck(tt, choose_axis(axes)[1])
    print("\npluck:")
    _print(rp)
    plot_pluck(rp, os.path.join(figdir, "phyphox_demo_pluck.png"), "Example output (synthetic data): pluck test")
    tt, axes, _ = load(sp)
    rs = analyse_sweep(tt, choose_axis(axes)[1])
    print("\nsweep:")
    _print(rs)
    plot_sweep(rs, os.path.join(figdir, "phyphox_demo_sweep.png"),
               theory={"fn": rp["fn_Hz"], "zeta": rp["zeta"], "alpha": 0.5}, title="Example output (synthetic data): frequency sweep")
    files = []
    for lab, z in (("gap 15 mm", 0.05), ("gap 3 mm", 0.10), ("gap 1 mm", 0.20)):
        t, a = simulate("sweep", fn=truth["fn"], zeta=z, friction_N=truth["friction_N"], seed=3, fs=100,
                        duration=45, jitter=0.1)
        fp = os.path.join(outdir, f"sweep_synthetic_{lab.split()[1]}mm.csv")
        write_phyphox_csv(fp, t, a)
        files.append(fp)
    plot_compare(compare_sweeps(files, ["15 mm", "3 mm", "1 mm"]),
                 os.path.join(figdir, "phyphox_demo_compare.png"), "Example output (synthetic data): three damper gaps")
    print(f"\nfigures: {figdir}/phyphox_demo_pluck.png, phyphox_demo_sweep.png, phyphox_demo_compare.png")
    return rp, rs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("pluck", "sweep", "base"):
        s = sub.add_parser(name)
        s.add_argument("file")
        s.add_argument("--axis", default="auto", choices=["auto", "x", "y", "z"])
        s.add_argument("--out", help="figure path (PNG)")
        if name == "pluck":
            s.add_argument("--mass", type=float, default=0.25, help="moving mass in kg (phone + platform)")
        if name == "sweep":
            s.add_argument("--Y", type=float, default=2.0, help="cam lobe amplitude in mm")
            s.add_argument("--fn", type=float, help="fn from your pluck test (Hz)")
            s.add_argument("--zeta", type=float, help="zeta from your pluck test")
            s.add_argument("--alpha", type=float, default=0.5, help="k_lower / (k_lower + k_upper)")
            s.add_argument("--csv", help="write the binned T(f) here")
    s = sub.add_parser("compare", help="overlay several sweeps (e.g. magnet gaps) with fitted models + table")
    s.add_argument("files", nargs="+")
    s.add_argument("--labels", nargs="+", help="one label per file, e.g. '15 mm' '3 mm' '1 mm'")
    s.add_argument("--Y", type=float, default=2.0, help="cam lobe amplitude in mm (measure with `base`)")
    s.add_argument("--axis", default="auto", choices=["auto", "x", "y", "z"])
    s.add_argument("--out", default="compare.png")
    s = sub.add_parser("analyse", help="analyse a whole folder: base*.csv, pluck_<tag>.csv, sweep_<tag>.csv")
    s.add_argument("folder")
    s.add_argument("--mass", type=float, default=0.25, help="moving mass in kg (phone + platform): WEIGH IT")
    s.add_argument("--Y", type=float, help="cam amplitude in mm (default: measured from base*.csv, else 2.0)")
    s.add_argument("--alpha", type=float, default=0.5, help="k_lower / (k_lower + k_upper)")
    s.add_argument("--axis", default="auto", choices=["auto", "x", "y", "z"])
    s.add_argument("--out", help="output folder (default <folder>/results)")
    sub.add_parser("demo")
    args = ap.parse_args(argv)
    if args.cmd == "analyse":
        rows, notes, Y, src = analyse_folder(args.folder, args.mass, args.Y, args.alpha, args.out, args.axis)
        out = args.out or os.path.join(args.folder, "results")
        print(f"Y = {Y:.2f} mm ({src}), mass = {args.mass} kg")
        if rows:
            print(f"{'set-up':>10s} {'fn pluck':>9s} {'zeta pluck':>11s} {'fn sweep':>9s} {'zeta sweep':>11s} {'peak T':>7s}")
            for r in rows:
                print(f"{r['setup']:>10s} {str(r['fn_pluck_Hz']):>9s} {str(r['zeta_pluck']):>11s} "
                      f"{str(r['fn_sweep_fit_Hz']):>9s} {str(r['zeta_sweep_fit']):>11s} {str(r['peak_T']):>7s}")
        else:
            print("no pluck_*.csv or sweep_*.csv files found")
        for n in notes:
            print("  note:", n)
        print(f"figures + summary.csv: {out}")
        return
    if args.cmd == "demo":
        demo()
        return
    if args.cmd == "compare":
        res = compare_sweeps(args.files, args.labels, args.Y, args.axis)
        rows = plot_compare(res, args.out)
        print(f"{'run':>12s} {'fn Hz':>7s} {'zeta':>7s} {'alpha':>6s} {'peak T':>7s}")
        for r in rows:
            print(f"{r[0]:>12s} {r[1]:>7s} {r[2]:>7s} {r[3]:>6s} {r[4]:>7s}")
        print(f"  figure: {args.out}")
        return
    t, axes, _ = load(args.file)
    axis, a = choose_axis(axes, args.axis)
    print(f"{args.file}: {len(t)} samples, {t[-1] - t[0]:.1f} s, axis {axis}")
    stem = os.path.splitext(args.file)[0]
    if args.cmd == "pluck":
        r = analyse_pluck(t, a, mass=args.mass)
        _print(r)
        out = args.out or stem + "_pluck.png"
        plot_pluck(r, out)
        if r["first_amplitude_mm"] > 0.8 * END_STOP_MM:
            print("  WARNING: first cycle near the end stops; pluck more gently")
    elif args.cmd == "sweep":
        r = analyse_sweep(t, a, Y_mm=args.Y)
        _print(r)
        theory = {"fn": args.fn, "zeta": args.zeta, "alpha": args.alpha} if args.fn and args.zeta else None
        out = args.out or stem + "_sweep.png"
        plot_sweep(r, out, theory)
        if r["end_stop_warning"]:
            print("  WARNING: platform displacement near the end stops; T near resonance is clipped")
        if args.csv:
            b = r["bins"]
            np.savetxt(args.csv, np.column_stack([b["f_Hz"], b["T"], b["n"]]), delimiter=",",
                       header="f_Hz,T,n_windows", comments="")
    else:
        r = analyse_base(t, a)
        _print(r)
        out = None
    if out:
        print(f"  figure: {out}")


if __name__ == "__main__":
    main()
