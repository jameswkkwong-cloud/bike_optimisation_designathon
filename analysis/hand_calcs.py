"""
hand_calcs.py — design calculations for the "Bumpy Road" quarter-bike vibration kit.

Single-degree-of-freedom base excitation: road (cam) displacement y(t) = Y sin(wt)
drives a sprung "rider" mass m (phone + platform) through stiffness k and damping c.

    fn  = (1/2pi) sqrt(k/m)                         natural frequency
    T   = sqrt((1+(2 zeta r)^2) / ((1-r^2)^2 + (2 zeta r)^2)),  r = f/fn   (TEXTBOOK)
    sag = m g / k                                   static deflection
    isolation (T < 1) for r > sqrt(2)

THIS RIG differs (corrected 20 Oct): only the two LOWER springs connect the road (carriage)
to the platform; the two upper springs and the eddy-current damper connect it to the FIXED
frame. So m x'' + c x' + (k_l + k_u) x = k_l y and

    T_rig = alpha / sqrt((1-r^2)^2 + (2 zeta r)^2),   alpha = k_l / (k_l + k_u) = 0.5

fn is unchanged (k = k_l + k_u), but T -> 0.5 at low frequency, and damping lowers T at all
frequencies (grounded "skyhook" damper). See phyphox_analysis.py and lab_sheet.md.

Run:  python designathon/scripts/hand_calcs.py
"""

import math

G = 9.81


def natural_frequency(k, m):
    return math.sqrt(k / m) / (2 * math.pi)


def transmissibility(f, fn, zeta):
    r = f / fn
    return math.sqrt((1 + (2 * zeta * r) ** 2) / ((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2))


def transmissibility_rig(f, fn, zeta, alpha=0.5):
    r = f / fn
    return alpha / math.sqrt((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2)


def stiffness_for(fn, m):
    return m * (2 * math.pi * fn) ** 2


def main():
    m_phone, m_platform = 0.180, 0.070          # kg (phones vary ~0.15-0.23 kg)
    m = m_phone + m_platform
    k = 440.0                                     # N/m: 4 x RS 821-352 (0.11 N/mm); 560 with RS 751-477 (0.14)
    Y = 2.0e-3                                    # m, cam lobe amplitude
    lobes, crank_rps = 3, (1.0, 3.0)

    fn = natural_frequency(k, m)
    sag = m * G / k
    print("BUMPY ROAD — design calculations")
    print(f"  moving mass m          = {m*1000:.0f} g")
    print(f"  spring rate k          = {k:.0f} N/m  (k for fn = 6 Hz: {stiffness_for(6.0, m):.0f} N/m)")
    print(f"  natural frequency fn   = {fn:.2f} Hz")
    print(f"  static sag mg/k        = {sag*1000:.1f} mm")
    print(f"  road frequency range   = {lobes*crank_rps[0]:.0f}–{lobes*crank_rps[1]:.0f} Hz "
          f"({lobes}-lobe cam at {crank_rps[0]:.0f}–{crank_rps[1]:.0f} rev/s)")
    print(f"  isolation starts at    = {math.sqrt(2)*fn:.2f} Hz (r = sqrt 2)")
    for lo, hi in ((0.150, 0.230),):
        print(f"  fn range for phones {lo*1000:.0f}–{hi*1000:.0f} g: "
              f"{natural_frequency(k, lo+m_platform):.2f}–{natural_frequency(k, hi+m_platform):.2f} Hz")
    print("\n  THIS RIG (upper springs + damper grounded, alpha = 0.5):")
    print("  zeta   peak T   peak amplitude   amplitude > sag?   peak accel   T(9 Hz)  [textbook T(9 Hz)]")
    for z in (0.05, 0.10, 0.15, 0.40):
        Tpk = max(transmissibility_rig(fn * (0.5 + i / 1000), fn, z) for i in range(1001))
        X = Tpk * Y
        acc = X * (2 * math.pi * fn) ** 2 / G
        print(f"  {z:4.2f}   {Tpk:6.2f}   {X*1000:9.1f} mm      {'YES -> preload/stops' if X > sag else 'no':20s} "
              f"{acc:5.2f} g      {transmissibility_rig(9.0, fn, z):5.2f}    [{transmissibility(9.0, fn, z):5.2f}]")
    # eddy-current damper, order of magnitude: c ~ sigma * t * B^2 * A (thin plate, magnet face A)
    c_crit = 2 * math.sqrt(k * m)
    for label, sigma, B, A, n in (("1 magnet, 3 mm Al, B~0.3 T", 3.5e7, 0.3, 2e-4, 1),
                                  ("magnet pair, 3 mm Cu, B~0.3 T", 5.9e7, 0.3, 2e-4, 2)):
        c = n * sigma * 3e-3 * B ** 2 * A
        print(f"  eddy damper ({label}): c ~ {c:.1f} N s/m  -> zeta ~ {c/c_crit:.2f} (c_crit = {c_crit:.1f})")
    print("    (order of magnitude only: real B and effective area depend on gap; MEASURE zeta in the pluck test)")
    print("\n  Design consequences:")
    print("   - at low damping the platform moves further than its static sag, so free springs")
    print("     would unload and the platform would hop: springs are preloaded (top bridge) and")
    print("     travel is limited by end stops (+-25 mm)")
    print("   - peak acceleration < 3 g, inside a phone accelerometer's range (check +-8 g)")
    print("   - a 4-lobe cam reaches ~10 Hz at 2.5 rev/s (3 rev/s is fast by hand)")


if __name__ == "__main__":
    main()
