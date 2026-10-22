# Bumpy Road: the quarter-bike vibration lab

**Module:** MECH0023 Dynamics and Control (Year 3, Term 1)  
**Time:** 45–60 min core, plus about 20 min for the recommended extension E1 · **Who:** on your own or in pairs, in the lab or at home  
**You need:**
- the kit;
- **your phone** with [phyphox](https://phyphox.org) (free, iOS/Android, no account);
- a laptop with Python 3 (`pip install numpy scipy matplotlib`);
- a ruler;
- about 100 g of coins or 4 AA batteries;
- a metronome app (optional).

> **Why this matters.** On a road bike, every bump passes through the tyres and frame into
> the rider. Designers trade **comfort** (isolating the rider) against **stiffness** (handling,
> power transfer). Your phone is the rider, the cam is the road, and the springs and magnet are
> the "seatpost". You will measure how much road vibration reaches the rider, then design a
> seatpost that meets a comfort target.

| MECH0023 learning outcome | Where you practise it |
|---|---|
| Measure the dynamic response of a simple system and identify the system parameters | Tasks 1–3 |
| Create a model of a dynamic system and implement it in a numerical simulation | Pre-lab, Extension E2 |
| Predict dynamic performance in the frequency domain (resonance, bandwidth) | Pre-lab, Task 3, Extension E1 |
| Understand the difference between a theoretical model and a practical system, and its limits | Discussion questions in Tasks 2–3, E1 |
| Topics: base excitation and vibration isolators; measuring frequency response; sample rate | Task 3, E1, E3 |

> ⚠️ **Safety**
> - **Magnets** (N42 neodymium): they snap together hard and can pinch. Keep them away from
>   bank cards and pacemakers. Small parts: not for children under 14.
> - **Strap the phone in** with the velcro strap before cranking.
> - **Keep your fingers clear** of the cam and roller. Crank by the handle only, and no faster
>   than about 3 rev/s.
> - End stops limit the platform's travel. If you hear knocking, slow down.
> - **Data:** phyphox stores recordings on your phone only and needs no account.

---

## 0. The rig and its model

```
   fixed bridge ────┬──────────────┬──── fixed damper plate (aluminium)
                    │ k_u (2 springs) ┆  c  (eddy-current: magnet on platform, plate fixed)
   platform + phone ┷━━━━ m ━━━━━━━━┷      x(t)
                    │ k_l (2 springs)
   carriage ────────┷──── roller on cam    y(t) = Y sin(2πft),  Y = 2 mm,  f = lobes × rev/s
```

The lower springs carry the road into the platform. The **upper springs and the damper are
attached to the fixed frame**, not to the road. This is *not* the textbook base-excitation
system, where the spring and damper both connect to the road.

---

## Pre-lab (15 min, before you touch the kit)

**P1. Model.** Draw a free-body diagram of the platform and show that

$$m\ddot{x} + c\dot{x} + (k_l + k_u)\,x = k_l\,y
\quad\Rightarrow\quad
T(f) = \frac{X}{Y} = \frac{\alpha}{\sqrt{(1-r^2)^2 + (2\zeta r)^2}},$$

where:
- $\alpha = k_l/(k_l + k_u)$;
- $r = f/f_n$;
- $f_n = \tfrac{1}{2\pi}\sqrt{k/m}$ with $k = k_l + k_u$;
- $\zeta = c/(2\sqrt{km})$.

What is T at very low frequency, and why is it not 1?

**P2. Predict.** Each spring has a datasheet rate of 0.11 N/mm (set A) or 0.14 N/mm (set B);
there are 2 above and 2 below. Take m = m_phone + m_platform, using the platform mass printed
on the kit label. Calculate for set A:

- $f_n$;
- $\sqrt2\,f_n$;
- the crank speed (rev/s) that hits resonance with the 4-lobe cam.

**P3. Compare.** Plot this rig's T(f) and the textbook
$T = \sqrt{(1+(2\zeta r)^2)/((1-r^2)^2+(2\zeta r)^2)}$ for ζ = 0.1 and 0.3. In which one does
**more damping always help**? Why can't a car suspension simply copy this rig ("skyhook"
damping)? Hint: what is the damper attached to?

---

## Task 1. Spring rate (10 min)

1. Hang one spare spring from its end, with a 100 g bag of coins on it.
2. Measure the extension and calculate k_each = mg/δ.
3. Then k = 4·k_each.

**Q:** How close is your k to the datasheet? What tolerance would you expect on a cheap spring?

## Task 2. Pluck test: natural frequency and damping (15 min)

1. Phone flat in the tray, screen up, strapped. In phyphox, open **"Acceleration (without g)"**.
2. Open the magnet gap to 15 mm with the sled wing nuts.
3. Press ▶, wait 1 s, push the platform down about 5 mm, and **let go**. Record 5 s. Do this 3 times.
4. Repeat with the gap at **3 mm** and at **1 mm**.
5. Export (⋮ → Export data → CSV) and run:

   `python phyphox_analysis.py pluck my_pluck.csv --mass 0.25`

   Use your own m in kg. The script prints $f_n$ and ζ and saves a figure.

| Gap (mm) | $f_n$ (Hz) | ζ | Best-fitting decay | k = m(2πf_n)² (N/m) |
|---|---|---|---|---|
| 15 | | | | |
| 3 | | | | |
| 1 | | | | |

**No laptop?** Read the peak spacing off the phyphox graph to get $f_d$. Then take the log
decrement by hand from two peaks n cycles apart: $\delta = \tfrac1n\ln(a_0/a_n)$ and
$\zeta = \delta/\sqrt{4\pi^2+\delta^2}$.

**Discuss:**
- Does $f_n$ match P2? If not, which input (k or m) is most likely wrong?
- The script fits an **exponential** decay (viscous damping) and a **straight-line** decay
  (Coulomb friction). Which fits your data? Where could friction come from?
- How does ζ change with the gap? Eddy-current force ∝ B², and B falls quickly with distance.

## Task 3. Frequency sweep: transmissibility T(f) (20 min)

1. Fit the **4-lobe cam**. Set the gap to 3 mm.
2. Start recording, then crank in **steps**, holding each speed for about 6 s. Set the
   metronome to one beat per crank turn: 60, 75, 85, 95, 100, 105, 115, 130, 145, 160 bpm.
   That gives road frequencies of 4–11 Hz. Steady cranking matters: uneven speed smears the
   peak.
3. Run:

   `python phyphox_analysis.py sweep my_sweep.csv --fn <Task 2> --zeta <Task 2>`

   The script computes T = X/Y for each 2 s window, overlays the model using your Task 2
   values, and fits $f_n$, ζ and α directly to the sweep.

**Discuss:**
- Where is the peak? What is T at the lowest frequencies? Compare with α from P1.
- Above what frequency does the rider move less than at low speed? Compare with $\sqrt2 f_n$.
- Is the ζ fitted to the sweep the same as the pluck-test ζ? Give **two** reasons it might
  differ. Think about amplitude, friction and your cranking.
- Estimate the half-power bandwidth and compare it with $2\zeta f_n$.
- **Design question:** if the platform were a saddle on a rough road at about 9 Hz, would you make
  the springs softer or stiffer, and add or remove mass? What does that cost in handling
  (stiffness)? Test your answer in E1.

---

## Extensions

- **E1. Design the seatpost (recommended, 20 min).**

  **Specification:**
  - **Comfort:** at a road frequency of **9 Hz**, T ≤ **0.4**.
  - **Safety:** when passing through resonance, T_peak ≤ **3**. That keeps the platform's travel
    under about 6 mm.

  **Choices:**
  - spring set **A or B**;
  - added mass **0 or about 100 g** (coins or batteries, taped under the strap);
  - magnet gap **1–15 mm**.

  1. **Predict first.** Use the model and *your* measured k, m and ζ(gap) to fill the table.
  2. Choose a design, build it, and **verify** with a sweep (as in Task 3).

  | Springs | Added mass | Gap | Predicted $f_n$ | Predicted T(9 Hz) | Predicted T_peak | Meets spec? | Measured T(9 Hz) |
  |---|---|---|---|---|---|---|---|
  | | | | | | | | |

  **Discuss:**
  - Which design is the **stiffest** that still meets the comfort target? On a bike, why would
    you want that one?
  - What did the model not capture?

- **E2. Simulation.** Write the equation of motion in Python or MATLAB/Simulink, adding
  friction $F\,\mathrm{sign}(\dot x)$. Drive it with your sweep's crank speeds. What value of F
  reproduces your pluck decay?
- **E3. Sample rate.** phyphox shows the sensor rate. What would a 9 Hz road look like if you
  sampled at 12 Hz? Sketch it and compute the alias frequency.
- **E4. Configuration.** If the damper plate were mounted on the carriage (damper to the road),
  how would T(f) change above $\sqrt2 f_n$? Use P3.

## Hand in

1. The T(f) plot with measured points, your model from Task 2, and the fitted model.
2. The prediction-vs-measurement table from Task 2 (and E1 if done).
3. Three sources of discrepancy, **each with a test that would confirm it**.

---

<sub>Instructor notes. Expected values:
- $f_n$ ≈ 6–7 Hz (phones 150–230 g, set A);
- T(f → 0) ≈ 0.5;
- ζ ≈ 0.05 at a wide gap, rising to about 0.1–0.3 at 1 mm with the stacked magnets (to be confirmed on the prototype);
- some Coulomb friction from the linear bearings is expected.

The E1 spec is met by set A + 100 g with ζ ≳ 0.1; set A alone and set B fail on comfort.

Analysis script and synthetic example data: `analysis/phyphox_analysis.py`
(`python phyphox_analysis.py demo`).</sub>
