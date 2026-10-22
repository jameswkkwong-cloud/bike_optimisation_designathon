# Put your phyphox exports here

Name the files so the script can match them (`.csv` or the phyphox `.zip`):

| File | Recording |
|---|---|
| `base.csv` (optional) | Phone taped to the **carriage**, crank a few speeds: measures the real road amplitude Y |
| `pluck_<setup>.csv` | Pluck test for one set-up, e.g. `pluck_15mm.csv` (magnet gap 15 mm) |
| `sweep_<setup>.csv` | Frequency sweep for the same set-up, e.g. `sweep_15mm.csv` |

Then, from the repository folder:

```bash
python analysis/phyphox_analysis.py analyse data/my_runs --mass 0.27
```

(use your weighed phone + platform mass in kg). Results appear in `data/my_runs/results/`:
one figure per pluck and sweep, `compare.png` (all set-ups), `summary.csv` and `notes.txt`.
Each sweep is plotted against the theory predicted from the pluck with the same `<setup>` name.
