# Bumpy Road: quarter-bike vibration kit (analysis software)

UCL MECH0073 Designathon 2026. An educational kit for **MECH0023 Dynamics and Control**:
a hand-cranked cam makes a "road", springs carry it to a platform, and the student's **own phone**
(free [phyphox](https://phyphox.org) app) is both the rider and the sensor. Students measure
natural frequency, damping and transmissibility, and compare them with a model they derive.

This repository holds the **laptop software** and the **student lab sheet**. The CAD, laser and
print files are submitted on Moodle.

| File | What |
|---|---|
| [`lab_sheet.md`](lab_sheet.md) | Lecture demo (pluck, crank through resonance, close the damper) + optional 20-min student follow-up |
| [`lab_sheet_full.md`](lab_sheet_full.md) | Full lab version: model derivation, spring rate, seatpost design challenge, extensions |
| [`analysis/phyphox_analysis.py`](analysis/phyphox_analysis.py) | Phone data → fₙ, ζ, T(f), with the rig model fitted; figures |
| [`analysis/hand_calcs.py`](analysis/hand_calcs.py) | Design calculations (fₙ, sag, peak amplitude, damper estimate) |
| `data/synthetic/` | **Synthetic** example recordings (made by the model, not measured) to try the tools |
| `figures/` | Example outputs from the synthetic data |
| `tests/` | Automated checks: the analysis recovers known fₙ and ζ from synthetic data |

## Quick start (any laptop, Python 3.9+)

```bash
pip install -r requirements.txt
python analysis/phyphox_analysis.py demo          # synthetic data -> figures/, no rig needed
python -m pytest -q                                # 10 tests
```

## Using it with the kit

1. In phyphox open **"Acceleration (without g)"**, phone flat in the tray, strapped.
2. Record, then **⋮ → Export data → CSV**. The `.zip` that phyphox makes works too.
3. **Easiest:** put the files in `data/my_runs/` named `base.csv` (optional), `pluck_<setup>.csv`,
   `sweep_<setup>.csv` (e.g. `pluck_15mm.csv`, `sweep_15mm.csv`) and run **one command**:

   ```bash
   python analysis/phyphox_analysis.py analyse data/my_runs --mass 0.27
   ```

   It measures Y from `base.csv`, analyses each pluck, plots each sweep against the theory from
   the pluck with the same name, and writes figures + `summary.csv` to `data/my_runs/results/`.
4. Or run the steps one at a time:

| Test | What you do | Command | You get |
|---|---|---|---|
| **Pluck** | Push the platform down ~5 mm, release, record 5 s | `python analysis/phyphox_analysis.py pluck pluck.csv --mass 0.25` | fₙ, ζ (log decrement / curve fit), viscous vs friction decay, figure |
| **Base** (once) | Phone taped to the **carriage**, short crank | `... base base.csv` | Real road amplitude Y (nominal 2 mm) |
| **Sweep** | Crank in steady steps (metronome), slow → fast | `... sweep sweep.csv --Y 2.0 --fn 9.4 --zeta 0.07` | T(f) with model overlay and fitted fₙ, ζ, α |
| **Compare** | One sweep per magnet gap | `... compare g15.csv g3.csv g1.csv --labels "15 mm" "3 mm" "1 mm"` | All T(f) curves + fitted table on one figure |

Use your measured moving mass for `--mass` (phone + platform, kg) and the Y from `base`.

## The model (derived in the lab sheet)

The lower springs k_l connect the road (cam carriage, motion y) to the platform (mass m, motion x);
the upper springs k_u and the eddy-current damper c connect the platform to the **fixed** frame:

$$m\ddot x + c\dot x + (k_l + k_u)\,x = k_l\,y \quad\Rightarrow\quad
T = \frac{X}{Y} = \frac{\alpha}{\sqrt{(1-r^2)^2 + (2\zeta r)^2}}$$

with α = k_l/(k_l+k_u) = 0.5, r = f/fₙ, fₙ = (1/2π)√(k/m), ζ = c/(2√(km)). Because the damper is
grounded ("skyhook"), more damping lowers T at every frequency, unlike the textbook base-excited
system (damper to the road). The textbook curve is plotted for comparison.

## Limits (read before trusting a number)

- **Validated on synthetic data only so far:** fₙ within ~2 %, ζ within ~10 %. Real-rig results
  will be added after testing.
- Hand cranking must be steady (metronome): speed varying more than ±20 % within a turn biases ζ up.
- Bearing friction (Coulomb) also raises ζ at small amplitude; the pluck output reports whether
  the decay looks viscous or frictional.
- T uses the cam amplitude Y: measure it with `base` rather than assuming 2 mm.

## Acknowledgements

AI use (UCL Category 2) is declared in [`AI_ACKNOWLEDGEMENT.md`](AI_ACKNOWLEDGEMENT.md).
Licence: MIT ([`LICENSE`](LICENSE)).
