# SB02 vessel-proxy validation study

**Question:** does satellite-detected vessel presence predict *acoustically*
detected ship presence at a hydrophone? If it does, a satellite proxy can carry
some of the load a hydrophone carries. If it does not, that is worth knowing
before anyone builds a layer on the assumption.

**Status: scoped and prepared. Not run.** Nothing here has touched Earth Engine.
No result is claimed, and none should be quoted from this directory.

## Why SB02

| | |
|---|---|
| Site | SanctSound SB02, Stellwagen Bank NMS, 42.470793 N, −70.24294 W |
| Ground truth | 15,255 timestamped ship-detection intervals, **45.7% duty cycle** |
| Window | 2018-11-12 → 2021-12-15 (the *ships* product, ~37 months) |
| SAR | 278 Sentinel-1 IW scenes, all dual-pol, 3 relative orbits (142/40/62), 93 ASC / 185 DESC |

The 45.7% base rate makes this a naturally balanced classification problem.
Sentinel-1B failed 2021-12-23, eight days after the ships record ends, so the
whole usable window has two-satellite revisit.

**Note:** the ships product is ~37 months, *not* the 44 months of the TOL
sound-level record. Exporting 44 months would pull scenes with no ground truth.

## Files

| File | Role |
|---|---|
| `sar_params.md` | Detection parameter specification, with derivations |
| `sb02_sar_export.js` | GEE Code Editor export (ES5 only — see below) |
| `fetch_ships.py` | Pulls the SanctSound ground truth; data is **not** committed |
| `acoustic_label.py` | Was a ship audible at time *t* ± tolerance |
| `flag_ambiguities.py` | Offline azimuth-ambiguity flagging on the detection table |

## Run order — not optional

1. `python3 fetch_ships.py` — writes `sb02_ships.json` (~715 KB, gitignored)
2. In GEE: read **CHECK 1–10** before exporting anything.
   **The 278 / 93 ASC / 185 DESC figures are PROVISIONAL.** They came from an
   exploratory query whose window and filters were never recorded beside
   them, and a live run of the shipped configuration
   (`2018-11-12 .. 2021-12-16`, end exclusive) returned **184 scenes across
   2 relative orbits**. One is wrong and the file cannot say which. CHECK 1b
   prints MATCH or MISMATCH; CHECK 3b names any configured orbit the window
   does not contain. Settle the window you want, re-measure, and update
   `EXPECTED_*` together with the query that produced it.
3. **CHECK 10 is a decision point.** The persistence histogram must be
   *bimodal* for the 0.20 threshold to be defensible. See the warning below.
   It reduces at `analysisScaleM` (10 m), the scale the mask is actually
   applied at. It used to reduce at `MASK_SCALE_M` (30 m), which measured a
   different detector — `annulusKernel` is a fixed PIXEL kernel with no
   `.reproject()`, so a 30 m request inflated the 400/150 m annulus to
   ~1200/450 m and the 300 m² minimum target to ~2,700 m². That detector
   finds fewer small fixed objects, suppressing the near-1.0 mode the test
   looks for. **If CHECK 10 times out at 10 m, shrink the AOI and say which
   one you used — do not coarsen the scale back.** The scale is printed in
   each CHECK 10 label.
4. Run STAGE A (persistence assets) with `PERSISTENCE_MODE: 'compute'`.
5. Switch to `PERSISTENCE_MODE: 'asset'`, then run the per-scene and
   detection exports. `PERSISTENCE_ASSET_PREFIX` already points at
   `projects/ee-lathavis/assets/`; change it if you export elsewhere. Both
   the write and the read path append `_<PARAM_SET_ID>`, so bumping
   `PARAM_SET_ID` for a sweep requires re-running STAGE A for that set.
6. Offline: join detections to acoustic labels, sweep radius, flag ambiguities.

Running step 5 in `'compute'` mode is legal and gives identical numbers, but it
re-evaluates the detector over every scene of every orbit inside every scene's
own computation. On 278 scenes at 10 m it is very unlikely to finish.

## Known-bad versions

Four defects have made this script fail or mislead. All four passed `node --check`.

**Died before CHECK 1 could print:**

1. `buildPersistence()` read a collection named `SCENES`, which exists nowhere.
2. The block calling it ran at module level ~115 lines *above* the `s1Joined`
   assignment. Fixing (1) alone would have swapped the `ReferenceError` for
   `Cannot read properties of undefined (reading 'filter')`.
3. `Math.log10` — an ES2015 *library* addition the Code Editor sandbox does not
   implement, the same class as `String.repeat` and `Object.assign`, both of
   which broke the main tool historically. On the default path, reached
   synchronously by `map()` at graph-build time.

**Failed hours later, after the tasks were queued:**

4. `ee.Image('ESA/WorldCover/v200')` — the STAC record says
   `gee:type: image_collection`. A lazy node, so it fails only on *evaluation*.
   Since the water mask feeds CHECK 7, CHECK 10/10b, the detector and all three
   exports, the symptom was CHECK 1–6 printing normally, CHECK 7 throwing, and
   every export task failing with an opaque asset-load error.

**And one that produced a confident wrong answer rather than an error:**
`flag_ambiguities.py` computed the azimuth axis from `asin()`, which returns only
the ascending branch — 168.88° at SB02, where the descending branch is 11.12°,
**22.24° away against a 20° tolerance**. 185 of the 278 scenes are descending, so
ghosts on two thirds of the data could never be flagged. The damage is not the
missing flags: `separation_histogram()` is the mandatory empirical check, and an
absent 4.5–6.5 km spike is documented as meaning ambiguities are not a material
contaminant. The old self-test could not catch it — it sited its synthetic ghosts
with the same function it then tested, so a wrong axis moved the ghost too.

Run `node tests/ci.js` (114 checks) and `python3 flag_ambiguities.py` before
pasting anything into the Code Editor. Section 15 *executes* every script under
`studies/` against the Earth Engine stub; section 16 tests the ES5 gate itself.
Neither can catch a type error like (4) — the stub cannot tell `ee.Image` from
`ee.ImageCollection`.

The copy before commit `HEAD` died in the Code Editor on line 778 with
`SCENES is not defined`, before CHECK 1 could print. Two separate defects:

1. `buildPersistence()` read a collection named `SCENES`, which appears nowhere
   else in the file.
2. The block that calls it ran at module level, ~115 lines *above* the line
   that assigns `s1Joined`. Fixing the name alone would have swapped the
   `ReferenceError` for `Cannot read properties of undefined (reading 'filter')`.

Both passed `node --check` and the ES5 gate. `tests/ci.js` section 15 now
*executes* every script under `studies/` against the Earth Engine stub, which
catches an undefined identifier and a use-before-assignment alike. Run
`node tests/ci.js` before pasting this file into the Code Editor.

## What a live run found (2026-09, config as shipped)

| | |
|---|---|
| CHECK 1 | **184**, not the recorded 278. CHECK 5a = 5b, so nothing was lost to the dual-pol filter |
| CHECK 3 | **2** of the 3 configured relative orbits present — one contributes zero scenes |
| CHECK 10, orbit 142 | `Image.divide: ... Got 0 and 1` — the empty-orbit crash, now fixed |
| CHECK 10, orbits 40 / 62 | `User memory limit exceeded` at 10 m over the 40 km disc |

The empty-orbit crash was real: `ee.ImageCollection([]).sum()` has zero bands, so
dividing it by the scene count threw, and the first empty orbit killed the whole
diagnostic. An empty orbit now produces an all-zero mask — the correct reading of
no evidence — and **CHECK 3b names it**, so a zero is never read as "nothing is
persistent here" when it means "nothing was looked at".

The memory limit is the cost the previous change warned about. `CHECK10_RADIUS_M`
defaults to **10 km** so CHECK 10 fits; both labels print the radius used. Raise
it once `PERSISTENCE_MODE` is `'asset'` and the detector is no longer re-evaluated
per scene. **Shrink the region, never coarsen the scale.**

## Still open, deliberately

**Persistence is built at production `k = 8`, not the spec's `k = 5`.**
`sar_params.md` §4.4 step 1 calls for a deliberately loose pre-pass; `buildPersistence()`
passes the same `PARAMS` as the detector, so `thresholdK` is 8.0. The consequence is
one-directional: at `k = 8` a marginal fixed object (a buoy, a small wreck, the
SanctSound mooring itself) is detected in fewer scenes, its `persistence_frac` falls
below 0.20, and it survives the mask — contributing a detection to every scene of its
orbit. That is a per-scene constant a count-vs-acoustic regression absorbs into the
intercept and never surfaces.

This is left as-is on purpose. It changes what the mask *is*, which is a decision about
the study rather than a defect, and it should be made after reading a CHECK 10 histogram
that can be trusted — which is what the scale fix above provides.

## Things that will bite you

**The persistence mask can delete real signal.** SB02 sits *in* a shipping lane.
If traffic is channelled tightly enough that one 10 m pixel is occupied in ≥20%
of scenes, the mask removes it as if it were a buoy — at exactly the site under
study. Adopt 0.20 only if CHECK 10's histogram is bimodal.

**The persistence mask is circular.** It is built *by* the detector using the
current `PARAMS`. Change the threshold, polarisation or analysis scale and the
mask is invalid. Asset names carry `param_set_id`; reusing a mask across
parameter sets is a silent error, not a caught one.

**Wind is a shared confounder, and it is the dangerous one.** High wind raises
SAR sea clutter (hiding vessels) *and* raises the acoustic noise floor (masking
detection). Predictor and ground truth move together for a reason that has
nothing to do with vessels. ERA5 wind is emitted per scene; a benign-sea-state
subset (3–10 m/s) must be reported alongside the pooled result or a positive
finding is not trustworthy.

**Tolerance above ±15 min is degenerate.** Measured on the real record:
±0 → 45.4% present, ±15 → 67.3%, ±30 → 78.7%, ±60 → **88.9%**. At ±60 a
classifier that always answers "present" scores 89%. The sweep is capped at
{0, ±15}; ±30 is a flagged low-power sensitivity run only.

**The two smallest radii are blurred by SAR physics.** A moving vessel is imaged
displaced along-track by `Δx = (Vr/Vg)·R_slant` — up to ~1.4 km at 19 kn. That
is ~65% of a 2 km radius. Do not read structure below ~2 km as detector quality.

**Power.** At n=278 and a 47% scene-weighted base rate, the test detects a
moderate proxy (Δ≥0.20 in detection rate) at 93% power but cannot resolve a weak
one (Δ≤0.10, 42%). A null here means *no strong proxy effect*, not *no effect*.

**Multiplicity.** 5 radii × 2 tolerances = 10 cells. Apply Benjamini–Hochberg
across the grid; report the skill-versus-radius curve as the primary output and
treat individual cell p-values as exploratory.

**Scene coverage is unverified.** The 278 count came from filtering on a *point*,
not a disc. The export emits `cov_frac_*` and `full_cov_*` so partial-coverage
scenes are excludable — check them before trusting any R=40 km result.

## ES5

`sb02_sar_export.js` runs in the Earth Engine Code Editor, which does **not**
support ES6. `node --check` will not catch that: ES6 *methods* are syntactically
valid ES5 and fail only at runtime in the browser. CI runs `tests/es5.js` over
every `.js` in this directory for that reason.
