# Bumpy Road: lecture demo and student follow-up

**Module:** MECH0023 Dynamics and Control · **Side A:** 10-minute lecture demo ·
**Side B:** optional 20-minute student follow-up (kit on loan or drop-in)

> ⚠️ **Safety**
> - **Magnets:** strong; they pinch, so keep them away from cards and pacemakers.
> - **Crank by the handle only**, keeping fingers clear of the cam and the platform, at no more than about 3 rev/s.
> - **Strap the phone in** before cranking.
> - Pluck gently (about 5 mm).
> - phyphox needs no account; data stays on the phone.

## The rig and its model

A hand-cranked **cam** moves the **carriage** (the road) up and down by Y = 2 mm, 4 bumps per turn.

- **Four lower springs** carry the road motion up to the **platform** (the rider: platform plus your phone).
- **Four upper springs** and a **magnet** moving past a fixed **aluminium plate** (an eddy-current damper) hold it to the fixed frame.

$$m\ddot x + c\dot x + (k_l+k_u)\,x = k_l\,y \qquad
T=\frac{X}{Y}=\frac{\alpha}{\sqrt{(1-r^2)^2+(2\zeta r)^2}}$$

where:
- $f_n=\tfrac{1}{2\pi}\sqrt{k/m}$, with $k = k_l + k_u$;
- $r=f/f_n$;
- $\zeta=c/2\sqrt{km}$;
- $\alpha=k_l/k=0.5$.

Road frequency f = 4 × crank turns per second.

---

## Side A: lecture demo (about 10 min)

**Before the lecture (2 min):**

1. Strap the phone in the tray and open phyphox → **Acceleration (without g)**.
2. **⋮ → Allow remote access**, then open the URL shown on the lecture laptop: the live graph goes on the projector.
   The phone and laptop must be on the same network. **eduroam often blocks this, so use a phone hotspot**, and test it beforehand.
3. Set the magnet gap to **15 mm** and start a metronome app.

Each demo is **predict → observe → explain**. Ask the room to commit to a prediction first.

### Demo 1: Pluck (natural frequency and damping)

- **Predict:** k = 8 springs × 0.11 N/mm = 880 N/m, m ≈ 0.25 kg. How many bounces per second? *(fₙ ≈ 9.4 Hz)*
- **Observe:** push the platform down about 5 mm and let go. A decaying wave appears on screen.
- **Explain:**
  - The bounce rate is $f_n=\tfrac{1}{2\pi}\sqrt{k/m}$, set by stiffness and mass only.
  - How fast it dies away is the damping ratio ζ: each peak is a fixed fraction of the last, so $\zeta \approx \tfrac{1}{2\pi}\ln(a_0/a_1)$.

### Demo 2: Crank slow → fast (resonance)

- **Predict:** at what crank speed will the rider shake most? *(fₙ/4 ≈ 2.4 rev/s ≈ 140 rpm)*
- **Observe:** crank in metronome steps, holding each:
  - **60 bpm (slow):** the platform moves *less* than the road;
  - **about 140 bpm:** the platform moves *far more* than the road. **Resonance.**
  - **about 170 bpm:** the motion starts to fall again.
- **Explain:**
  - At low speed T ≈ α = 0.5, because only the lower springs connect the rider to the road.
  - T peaks at f = fₙ.
  - The rider would be isolated above √2·fₙ (≈ 13 Hz, ≈ 200 rpm): too fast to crank by hand with this stiff, 8-spring
    "seatpost". Ask the class how to bring isolation within reach: softer springs, more mass or more cam lobes.

### Demo 3: Close the damper

- **Predict:** close the magnet gap to 1 mm. What happens to the resonance peak?
- **Observe:** crank through resonance again. The peak is much smaller.
- **Explain:**
  - A moving magnet induces eddy currents in the plate, which oppose the motion: $c \propto \sigma\, t\, B^2 A$, and B falls fast with distance.
  - Here the damper is fixed to the frame ("skyhook"), so more damping helps at **every** frequency.
  - A car's damper connects to the road and can't do this. Only **active control** can imitate it, which links to the control half of the module.

---

## Side B: student follow-up (optional, about 20 min)

1. **Record** four files in phyphox (⋮ → Export data → CSV):

   | File | How |
   |---|---|
   | `pluck_15mm.csv` | Gap 15 mm: pluck about 5 mm, record 5 s |
   | `sweep_15mm.csv` | Gap 15 mm: crank in metronome steps, 60, 80, 100, 115, 125, 135, 140, 145, 155, 170 bpm, about 6 s each |
   | `pluck_1mm.csv`, `sweep_1mm.csv` | The same with a 1 mm gap |

2. **Analyse.** Put the files in `data/my_runs/` of the GitHub repo and run:

   ```bash
   python analysis/phyphox_analysis.py analyse data/my_runs --mass <your phone + platform in kg>
   ```

   You get each sweep plotted against the theory from its pluck, plus `compare.png` and `summary.csv`.
   **No laptop?** Read fₙ off the phyphox graph, and take $\zeta \approx \tfrac{1}{2\pi n}\ln(a_0/a_n)$ from two peaks n cycles apart.

3. **Answer:**
   1. How close is your measured fₙ to $\tfrac{1}{2\pi}\sqrt{k/m}$? Which input is most likely off?
   2. Why does T start at about **0.5** and not 1? Which springs connect the rider to the road?
   3. Did the pluck decay look **exponential** (viscous damping) or like a **straight line** (friction)? Where could friction come from?
   4. *(Our capstone link)* This rig resonates near **9 Hz**. Is that a good 'seatpost' for a road that shakes at 9 Hz? To isolate the rider, would you make the springs softer or stiffer, and add or remove mass? What would that cost a bike in handling stiffness?

## What this supports in MECH0023

| Learning outcome | Where |
|---|---|
| Measure the dynamic response of a simple system and identify its parameters | Demo 1, Side B |
| Predict frequency-domain behaviour: resonance, isolation | Demos 2–3, Q2, Q4 |
| Theoretical model vs practical system, and its limits | Q1, Q3 |

<sub>Instructor notes.
- **Expected values:** fₙ ≈ 8.5–10 Hz with 8 springs (phones of 150–230 g); T(low) ≈ 0.5; ζ ≈ 0.05 at a 15 mm gap, rising at 1 mm (measure it).
- **Printed sleeve bearings add friction:** expect some straight-line decay. A drop of light oil on the rods helps.
- **Uneven cranking smears the peak:** use the metronome.
- **Full version** (model derivation, spring-rate task, seatpost design challenge, simulation and sample-rate extensions): `lab_sheet_full.md`.</sub>
