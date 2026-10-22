"""Designathon phyphox analysis: model limits, parameter recovery from synthetic phone
recordings (pluck and sweep), Coulomb-friction detection, file formats and the CLI."""
import os
import sys
import zipfile

import numpy as np
import pytest

pytest.importorskip("scipy")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import phyphox_analysis as P                                     # noqa: E402


def test_rig_model_limits():
    f = np.array([0.01, 100.0])
    assert P.transmissibility_rig(f, 6.68, 0.1, 0.5)[0] == pytest.approx(0.5, rel=1e-3)
    assert P.transmissibility_textbook(f, 6.68, 0.1)[0] == pytest.approx(1.0, rel=1e-3)
    ff = np.linspace(1, 16, 50)
    # grounded ("skyhook") damper: more damping lowers T at every frequency
    assert np.all(P.transmissibility_rig(ff, 6.68, 0.3) < P.transmissibility_rig(ff, 6.68, 0.1))
    # textbook damper: more damping RAISES T well above sqrt(2) fn
    assert P.transmissibility_textbook(15.0, 6.68, 0.3) > P.transmissibility_textbook(15.0, 6.68, 0.1)


@pytest.mark.parametrize("zeta", [0.05, 0.2])
def test_pluck_recovers_fn_and_zeta(zeta):
    t, a = P.simulate("pluck", fn=6.68, zeta=zeta, fs=200, noise=0.05)
    r = P.analyse_pluck(t, a)
    assert r["fn_Hz"] == pytest.approx(6.68, rel=0.03)
    assert r["zeta"] == pytest.approx(zeta, rel=0.10)


def test_pluck_flags_coulomb_friction():
    t, a = P.simulate("pluck", fn=6.68, zeta=0.005, friction_N=0.05, noise=0.05)
    r = P.analyse_pluck(t, a, mass=0.25)
    assert r["damping_model"].startswith("Coulomb")
    assert r["friction_N_if_coulomb"] == pytest.approx(0.05, rel=0.25)


def test_sweep_fit_recovers_rig_parameters():
    t, a = P.simulate("sweep", fn=6.68, zeta=0.1, fs=100, noise=0.05, duration=45, jitter=0.1)
    r = P.analyse_sweep(t, a)
    fit = r["fit"]
    assert fit["fn_Hz"] == pytest.approx(6.68, rel=0.03)
    assert fit["zeta"] == pytest.approx(0.1, rel=0.25)
    assert fit["alpha"] == pytest.approx(0.5, rel=0.15)
    assert not r["end_stop_warning"]


def test_base_calibration_recovers_cam_amplitude():
    t = np.arange(0, 30, 0.01)
    th = 2 * np.pi * (0.8 * t + 0.5 * 2.2 * t ** 2 / 30)
    y = 2e-3 * np.sin(4 * th)
    a = np.gradient(np.gradient(y, t), t) + 0.05 * np.random.default_rng(0).standard_normal(len(t))
    assert P.analyse_base(t, a)["Y_mm_median"] == pytest.approx(2.0, rel=0.08)


def test_loader_reads_comma_semicolon_and_zip(tmp_path):
    t, a = P.simulate("pluck", duration=3.0, fs=100)
    p = tmp_path / "raw.csv"
    P.write_phyphox_csv(str(p), t, a)
    t1, ax1, _ = P.load(str(p))
    assert P.choose_axis(ax1)[0] == "z"
    text = p.read_text().replace(",", ";")
    text = "\n".join(text.splitlines()[:1] + [ln.replace(".", ",") for ln in text.splitlines()[1:]])
    q = tmp_path / "excel.csv"
    q.write_text(text)
    t2, ax2, _ = P.load(str(q))
    z = tmp_path / "export.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.write(p, "Raw Data.csv")
    t3, ax3, _ = P.load(str(z))
    for tt, ax in ((t2, ax2), (t3, ax3)):
        assert np.allclose(tt, t1) and np.allclose(ax["z"], ax1["z"])


def test_cli_pluck_writes_figure(tmp_path):
    t, a = P.simulate("pluck", fs=200)
    p = tmp_path / "pluck.csv"
    P.write_phyphox_csv(str(p), t, a)
    P.main(["pluck", str(p), "--out", str(tmp_path / "fig.png")])
    assert (tmp_path / "fig.png").stat().st_size > 10_000


def test_compare_orders_damping_and_writes_figure(tmp_path):
    files = []
    for i, z in enumerate((0.05, 0.2)):
        t, a = P.simulate("sweep", fn=6.68, zeta=z, fs=100, noise=0.05, duration=40, jitter=0.1, seed=i)
        p = tmp_path / f"s{i}.csv"
        P.write_phyphox_csv(str(p), t, a)
        files.append(str(p))
    res = P.compare_sweeps(files, ["low", "high"])
    zl, zh = res[0][1]["fit"]["zeta"], res[1][1]["fit"]["zeta"]
    assert zl < zh and res[0][1]["peak_T"] > res[1][1]["peak_T"]
    P.main(["compare", *files, "--labels", "low", "high", "--out", str(tmp_path / "c.png")])
    assert (tmp_path / "c.png").stat().st_size > 10_000


def test_analyse_folder_matches_files_and_uses_pluck_theory(tmp_path):
    d = tmp_path / "runs"
    d.mkdir()
    t = np.arange(0, 30, 0.01)
    th = 2 * np.pi * (0.8 * t + 0.5 * 2.2 * t ** 2 / 30)
    y = 1.9e-3 * np.sin(4 * th)
    P.write_phyphox_csv(str(d / "base.csv"), t, np.gradient(np.gradient(y, t), t))
    for tag, z in (("15mm", 0.05), ("1mm", 0.2)):
        tp, ap_ = P.simulate("pluck", fn=6.68, zeta=z, fs=200, noise=0.05)
        P.write_phyphox_csv(str(d / f"pluck_{tag}.csv"), tp, ap_)
        ts, as_ = P.simulate("sweep", fn=6.68, zeta=z, Y_mm=1.9, fs=100, noise=0.05, duration=40, jitter=0.1)
        P.write_phyphox_csv(str(d / f"sweep_{tag}.csv"), ts, as_)
    P.main(["analyse", str(d), "--mass", "0.25"])
    out = d / "results"
    rows = {r["setup"]: r for r in __import__("csv").DictReader(open(out / "summary.csv"))}
    assert set(rows) == {"15mm", "1mm"}
    assert float(rows["15mm"]["Y_mm"]) == pytest.approx(1.9, rel=0.08)          # Y measured from base.csv
    for tag, z in (("15mm", 0.05), ("1mm", 0.2)):
        assert float(rows[tag]["fn_pluck_Hz"]) == pytest.approx(6.68, rel=0.04)
        assert float(rows[tag]["zeta_pluck"]) == pytest.approx(z, rel=0.15)
    assert float(rows["15mm"]["peak_T"]) > float(rows["1mm"]["peak_T"])
    for f in ("pluck_15mm.png", "sweep_1mm.png", "compare.png", "notes.txt"):
        assert (out / f).exists()
