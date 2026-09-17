# SB02 vessel-proxy validation study

**Question:** does satellite-detected vessel presence predict *acoustically*
detected ship presence at a hydrophone? If it does, a satellite proxy can carry
some of the load a hydrophone carries. If it does not, that is worth knowing
before anyone builds a layer on the assumption.

**Status: prototype profile (p002). No result is claimed**, and none should be
quoted from this directory. The pre-flight CHECKs have been run live; no export
has produced a table yet.

## The p002 prototype profile, and why persistence was dropped

STAGE A did not deliver. Two orbit tasks ran 8 and 10 hours and left the assets
empty, so persistence masking is unavailable — and the study had been blocked on
it for a full cycle.

**It should not have been.** A persistent fixed object — a navigation buoy, a
wreck, the SanctSound mooring itself — contributes a roughly *constant* count to
every scene of a relative orbit. A constant offset does not change a rank
correlation or an AUC at all; it moves a logistic fit's intercept, not its
discrimination. For the question this study actually asks — *does SAR vessel
presence predict acoustic ship presence?* — persistence masking is close to
irrelevant. It matters for absolute counts and per-km² densities, which the
prototype does not claim.

So `p002` is: `PERSISTENCE_MODE: 'off'`, `RADII_M` to 10 km, everything else
identical to `p001`. That makes the table exports affordable — roughly 1.3 h
rather than the ~20 h the 40 km version needs — and unblocks the whole question
today. The 10 km disc also covers the plausible acoustic detection range for
shipping-lane traffic, and `range_m` in the detection table gives 2 and 5 km
offline for free.

**What this profile does NOT let you off.** A fixed object's *detectability*
varies with sea state, and wind drives acoustic noise too, so there is a
confound path from wind to both sides of the comparison. That is what the ERA5
columns are for, and `proxy_skill.py` conditions on them rather than trusting a
raw correlation — see below.

Every row now carries `p_persistence_mode` and `p_persistence_threshold`, so a
persistence-off table can never be mistaken for a persistence-on one. Before
that, the two were distinguishable only by their numbers.

## The 278 is reconstructed, and the file's own orbit list is incomplete

CHECK 3c–3f ran. The verdict line said **NOT CONFIRMED** — a 40 km bounds
filter yields 523, not 278 — and that is correct as far as it goes. But the
detail prints reconstruct the recorded figures exactly.

At 40 km the window holds **four** relative orbits, not three:

| orbit | scenes @ 40 km | scenes @ point | pass |
|---|---|---|---|
| 135 | 154 | 0 | ASCENDING |
| 142 | 93 | 0 | DESCENDING |
| 40 | 92 | 92 | DESCENDING |
| 62 | 184 | 92 | ASCENDING |

The pass assignment is **proven, not guessed**: `DESC = 185` has exactly one
subset solution, `{40, 142}` = 92 + 93, leaving `ASC = {62, 135}` = 338. That
agrees with the point query's 92 ASC / 92 DESC over orbits 62 and 40.

And then the recorded figures fall out:

```
orbit 40  (DESC)   92
orbit 142 (DESC)   93   ->  DESC = 185   exact match
orbit 62  (ASC)    93   ->  ASC  =  93   exact match
                  278                    exact match
```

So the orbit-142 idea was right — it *is* the missing descending orbit with 93
scenes, and it closes DESC exactly. What was wrong was the radius. The recorded
query used bounds **between the point and 40 km**: wide enough to take in all of
orbit 142, narrow enough to exclude orbit 135 entirely and to catch only 93 of
orbit 62 rather than 184. `CHECK3C_SWEEP_RADII_M` now sweeps six radii so the
next run measures that bound instead of inferring it.

**Orbit 135 is not in `RELATIVE_ORBITS` and never was.** Under the shipped
`'point'` filter it contributes nothing, so nothing is wrong today. But the
config comment has always described 142/40/62 as *the* orbits, and at 40 km that
is incomplete. Anyone widening `FILTER_BOUNDS` without adding 135 would get a
persistence mask that silently ignores a quarter of the scenes.

## The 80 MiB tile limit: a ladder, and how to read the error

Three exports and a CHECK have now died on it. The error is more informative
than any estimate, and it reads the same way every time:

```
Output of image computation is too large (1 bands for N pixels = M MiB > 80.0 MiB)
```

Two numbers name the node:

- **M / N gives bytes per pixel.** 8.0 means a float64 (or int64) band.
- **√N − 4096, halved, gives the halo.** EE's tile is 4096 px; the halo is
  whatever neighbourhood the node needs.

Worked, for the CHECK 8 failure: 132 MiB / 17,305,600 px = **8.0 bytes/px**, and
√17,305,600 = 4160 = 4096 + 2×**32**. A 32 px halo at the 10 m analysis scale is
`maxComponentSizePx` **1024 = 32²**, exactly — so the node is one of the bands
carrying `connectedPixelCount`'s halo, and the two that are a single float64
band are the area accumulators, both built on `ee.Image.pixelArea()` (a double).

**Casting is not one decision, it is two, and they are not equally free:**

| band | cast is | why |
|---|---|---|
| `areaImg`, `tAreaImg` | **free** | outputs; nothing compares them to a threshold, so no detection can change. And exact — the product tops out near 1.0e5, float32 is exact to 2^24 |
| land distance | **free** | compared in squared units; equivalent by construction |
| scene valid mask | **free** | `gt(0)` before `focal_min` says the same thing in a byte |
| detector `work`, `mu`, `thr` | **not free** | these *are* the threshold comparison; casting changes decisions in their last bits, so it must be reported alongside any result |

`tileScale` does **not** help here. It partitions the reducer's aggregation, and
these failures are in the image graph, whose neighbourhood ops tile themselves —
which is why the tile stayed 4096-based at `tileScale: 16`.

## Read the BUILD line first, every time

The Code Editor holds a **copy** of `sb02_sar_export.js`. Every fix merged here
has to be re-pasted by hand, and nothing in the console used to say which
version was in the tab. That cost a full cycle: two exports failed the 80 MiB
tile limit, the fix raising `tileScale` to 16 merged **27 minutes later**, and
the next look at the task list showed the same two failures — same task IDs,
same timestamps — which read as if the fix had not worked. It had not been run.

So the pre-flight now opens with a line like:

```
BUILD (read from live config - if this disagrees with the repo, re-paste the
script): param_set=p002  radii_km=2/5/10  persistence=off  tileScale=16
analysisScale=10m  landDilate=distance  medianBg=on
```

Every value is read from the live `CONFIG`/`PARAMS`, never typed into a string,
so it cannot go stale the way a hand-maintained version stamp does — which this
project got wrong three times (see the CHANGELOG). CI asserts the printed line
**agrees with the config parsed from source**, which is the check that catches a
stamp drifting from what it claims to describe.

If a value disagrees with what you expect, the tab is out of date. Re-paste
before reading anything below it.

## The p002 export failed a THIRD way, and it was mine

Not memory, not a timeout — a per-tile output limit:

> Output of image computation is too large (1 bands for 18939904 pixels =
> 144.5 MiB > 80.0 MiB). If this is a reduction, try specifying a larger
> `tileScale` parameter.

The arithmetic names the culprit. 144.5 MiB / 18,939,904 px is **8.0 bytes per
pixel exactly** — a float64 tile — and 18,939,904 is 4352², which is a 4096 px
tile plus a 128 px halo. The distance-transform land mask built that double: it
did `.sqrt().multiply(pixelArea().sqrt())`, and both are doubles.

Two fixes, in the order they should be reached for:

1. **`tileScale` 4 → 16** (the maximum), which is the remedy the error message
   itself names. Larger tileScale means *smaller* tiles. It is a pure
   compute-partitioning knob — it changes how the work is divided and never a
   number in the output. Not free: more, smaller tiles means more per-tile
   overhead, so a job that already fits runs somewhat slower.
2. **The land buffer is compared in squared units, in float32.** This is
   *exactly* equivalent, not an approximation: `fastDistanceTransform` returns
   d_px², `pixelArea` is side², so d_px²·side² = d_m², and for non-negative
   quantities `d_m <= B` says the same thing as `d_m² <= B²`. It drops two sqrt
   nodes, and float32 halves the tile to 72.2 MiB. The cast must come *before*
   the multiply — casting after still builds the double-typed product first.

**If it fails this way again**, `LAND_DILATE_METHOD: 'focal'` reverts to the
byte-typed `focal_max` path in one line. It is the slower mask, but it cannot
produce a wide-typed tile. The next lever after that is the detector's own
`work` image, which is also a double — casting it changes the numbers in their
last bits, so it is a change to make deliberately and report, not quietly.

## The two ways this prototype could produce a confident wrong answer

**Wind confounds both sides.** Higher wind means more sea clutter (fewer or
noisier SAR detections at a fixed threshold) *and* more ambient noise (ship
detection harder acoustically). A healthy-looking raw correlation can be wind
moving both, with no causal path between them. `proxy_skill.py` therefore always
prints the pooled AUC *and* the AUC within wind terciles, together. A pooled
figure well above the three stratified ones is the signature of the confound,
not of a working proxy. A self-test constructs exactly that case — pooled 0.767,
all three strata 0.500 — so the control is one that has been seen to fire, and a
mirror test confirms it does not flatten a real effect.

**A recording gap is not a silent ocean.** `acoustic_label.ship_present()`
answers True/False, which is the right contract for it and the wrong one for
scoring: between SanctSound deployments nobody is listening, and a SAR scene in
such a gap would be scored as "no ship". Those scenes keep their SAR counts, so
they land in the negative class carrying whatever the satellite saw, dragging
the AUC toward 0.5. `proxy_skill.py` adds an explicit **unknown** and drops
those rows, using a gap-threshold heuristic that it sweeps (6 h and 24 h) rather
than fixes, because recorder-on periods cannot be recovered exactly from
detection intervals alone.

It also excludes partial-coverage scenes (`full_cov_*`), which undercount by
construction, and reports every exclusion count next to every statistic.

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
| `proxy_skill.py` | **Prototype** skill evaluation: joins the per-scene table to the acoustic record |

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
3. **CHECK 10 is a decision point, but not at this step.** In `'compute'` mode
   it prints a **random-sample spot check** — 10,000 px, binned to 0.05, over
   `CHECK10_SPOT_RADIUS_M` (2 km) and the first `CHECK10_MAX_SCENES` (20)
   acquisitions of the orbit. An exhaustive reduction here is ~1.2e12 kernel
   operations and failed at 40 km and again at 10 km; **sampling alone did not
   fix it**, because `sample()` bounds how many pixels come *back*, not how
   many are computed. The input is now bounded on both axes — region *and*
   scene count. Neither lever is `analysisScaleM`: coarsening that measures a
   different detector. **CHECK 10b is skipped in `'compute'` mode** — the
   masked area is a property of the mask the export actually applies, and
   there is no cheap honest form of it until STAGE A has run. Read the spot
   check; if it is already clearly unimodal, stop and reconsider before
   spending a STAGE A run.
4. Run STAGE A (persistence assets) with `PERSISTENCE_MODE: 'compute'`.
   **Check the asset root first** — it is the one pre-flight check no CHECK in
   the script can make. The export dialog validates the root only when you
   press RUN, so a root that does not exist fails there and nowhere earlier
   (`projects/ee-lathavis` was a guess and returned *Resource
   projects/ee-lathavis could not be found*). The Code Editor's **Assets** tab
   lists every root the account can write to.
   **Run one orbit first and let it finish** before starting the next. Each is
   ~50 M pixels × ~92 scenes × a 4,316-weight annulus median; batch tasks do
   not share the interactive timeout but can still time out after hours, and
   one orbit tells you whether the budget holds for half the cost of finding
   out the hard way.
   **Run the task for every orbit in `RELATIVE_ORBITS`, including an empty
   one.** Step 5 loads an asset per *configured* orbit, not per orbit that has
   scenes — `ee.Image(prefix + orbit + '_' + PARAM_SET_ID)` for each — and
   CHECK 10 then reduces every one of them. A missing asset is a lazy node: it
   builds fine and throws `Image.load: asset not found` on evaluation, the same
   failure class as the WorldCover `ee.Image`-vs-`ImageCollection` bug below.
   So skipping the empty orbit does not yield an all-zero histogram, it yields
   an error. Its task is cheap — no scenes means no detector runs, just a
   50 M-pixel write of zeros — so there is nothing to save by skipping it.
   Do **not** run the two table exports yet: in `'compute'` mode they
   re-evaluate the detector over every scene of the orbit inside every scene's
   own computation. They are step 5. Running one costs a task: a live run of
   `SB02_S1_vessel_detections_p001` in `'compute'` mode failed after 12 minutes
   and 147 EECU-seconds with `User memory limit exceeded`, which is the
   documented behaviour and not a defect.
5. Switch to `PERSISTENCE_MODE: 'asset'` and **read CHECK 10 again** — now
   exhaustive over the full 40 km disc, cheap because the detector is no longer
   being re-evaluated. **This is the histogram to base the threshold on.** Then
   run the per-scene and detection exports. `PERSISTENCE_ASSET_PREFIX` points at
   `projects/stem-marine-engine/assets/`; change it if you export elsewhere —
   **in the config, never in the export dialog.** The dialog sets where that one
   task writes, while the read path builds its id from the config string, so
   fixing only the dialog leaves STAGE A writing where step 2 will not look.
   Both the write and the read path append `_<PARAM_SET_ID>`, so bumping
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

Run `node tests/ci.js` (131 checks) and `python3 flag_ambiguities.py` before
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
| CHECK 6 | `User memory limit exceeded` — on `s1Joined.size()`, which is only a count |
| CHECK 10, **sampled**, all three orbits | `User memory limit exceeded` again at 10 km, **including the empty orbit 142** |
| CHECK 10b | `Earth Engine memory capacity exceeded` |

## What the next live run found (after the memory fixes)

Every memory error is gone. CHECK 6 returns **184**, equal to CHECK 1, so every
scene found an ERA5 hour and no row will carry `era5_matched = 0`. CHECK 7
returns, and reads healthy — see below. One failure remains, and it is a
different one.

| | |
|---|---|
| CHECK 2 | **92 ASCENDING, 92 DESCENDING** — not the recorded 93 / 185 |
| CHECK 3 | orbit 40 → **92**, orbit 62 → **92**; orbit 142 absent |
| CHECK 7 | 1.244e7, 7.773e7, 3.109e8, 1.244e9, 4.882e9 m² |
| CHECK 10, orbits 40 / 62 | **`Computation timed out`** — no longer a memory error |

**CHECK 7 is healthy, and the way to see that is the ratios.** Against the true
circle areas the inner four radii come in at 0.98964, 0.98966, 0.98966, 0.98967
— a *constant* fractional deficit. Land cannot do that: the nearest coast is
~25 km away, so a land deficit would be zero at 2 km and grow with radius. A
constant ratio is a shape property, and 0.9896 is the area of a ~25-gon
inscribed in a circle — the polygon `Geometry.buffer()` approximates the
geodesic disc with. At 40 km a further 1.85% of the polygon is gone, ≈ 92 km²,
which is the Cape Ann sliver plus the harbour boxes. So the mask is **not**
eating ocean, which is exactly what CHECK 7 exists to rule out.

**The 278 figure now has a candidate explanation, not yet a confirmed one.**
92 ASC + 92 DESC = 184 here. If a third, *descending* orbit of ~93 scenes were
included, DESC would be 92 + 93 = **185** and the total 278 — the recorded
numbers exactly, with ASC 92 → 93 from one extra scene. That is consistent with
a query whose bounds took in orbit 142, which covers part of the 40 km disc
without covering the hydrophone. `FILTER_BOUNDS: 'maxdisc'` is the one-line
test; until someone runs it this stays a hypothesis and `EXPECTED_*` stays as
recorded, with CHECK 1b reporting MISMATCH.

**The CHECK 10 timeout is the interactive wall clock, not a memory ceiling, and
it does not block STAGE A** — batch exports do not share the ~5 minute limit.
The 2 km × 20 scene spot check was an arithmetic estimate and it was wrong; an
annulus median over 4,316 weights is a per-pixel *selection*, slower per
operation than the kernel count implies. The defaults are now 1 km × 12 scenes.
If that times out too, do not tune a fourth time — go to STAGE A, and read the
histogram in `'asset'` mode where it is a raster read.

The empty-orbit crash was real: `ee.ImageCollection([]).sum()` has zero bands, so
dividing it by the scene count threw, and the first empty orbit killed the whole
diagnostic. An empty orbit now produces an all-zero mask — the correct reading of
no evidence — and **CHECK 3b names it**, so a zero is never read as "nothing is
persistent here" when it means "nothing was looked at".

The memory limit is the cost the previous change warned about. **Shrink the
region, never coarsen the scale** still holds — but the run above shows that
shrinking the region was not sufficient, and points at what else was wrong.

**Orbit 142 is the diagnostic.** CHECK 3b reports it as empty, so its `frac`
image is a *constant* — it cannot cost 1.2e12 kernel operations, and yet its
sampled histogram failed identically to the other two. So the detector was
never the only cost. Two things were:

1. **`sample()` bounds the output, not the input.** `numPixels` caps how many
   pixels come back; the image is still evaluated tile by tile across the whole
   region, so every pixel still paid for its 81×81 annulus on every scene. The
   claim that sampling was "~300× cheaper" was wrong. The spot check now caps
   the **scene count** (`CHECK10_MAX_SCENES`) as well as the region
   (`CHECK10_SPOT_RADIUS_M`), which is what actually bounds the input.
2. **The land mask was a 1000 m `focal_max`.** That kernel is sized in *pixels
   at the request scale*: radius 100 px, ~31,400 weights per output pixel, at
   the 10 m analysis scale — and `WATER_MASK` is read by the detector on every
   scene and by CHECK 7, CHECK 10 and CHECK 10b. It is what made the *empty*
   orbit expensive. `LAND_DILATE_METHOD: 'distance'` now grows the buffer with
   `fastDistanceTransform`, which selects **the identical set of pixels** — a
   dilation by a disc of radius *r* is exactly the set within distance *r* —
   in linear time. `'focal'` reproduces the old path.

**CHECK 6 was not a cost to tune; it was a join to remove.** It printed
`s1Joined.size()` over an `ee.Join.saveBest` matching 184 scenes against
~27,000 hourly ERA5 images, and the join had to be materialised before the
result could be counted. That collection also fed the per-scene table, the
detection table and `buildPersistence()`, so the cost sat under every export,
not just the print that exposed it. `era5At()` now looks the nearest hour up
per scene — `filterDate` a ±1 h window (≈3 candidates), sort by |Δt|, take the
first. Same value, same single server-side pass, ~3 images per scene instead of
27,000 for the collection. `era5_dt_min` is now genuinely **signed** (ERA5 hour
minus acquisition); `saveBest`'s `measureKey` was an absolute difference, so the
old comment promising a signed value was wrong about its own output. CHECK 6
now counts scenes that *found* an ERA5 hour, which is the question it was for.

`CHECK10_RADIUS_M` (10 km) remains the **asset-mode** region for CHECK 10b;
every label prints the radius it used.

## Why exporting before settling the threshold is safe

The old instruction — read CHECK 10 *before* exporting — was circular: the full
histogram is not affordable in `'compute'` mode. The circle breaks on a fact
about the export, not the compute budget.

**STAGE A writes both bands, `persist` AND `frac`.** `persistenceThreshold` is
applied only inside `buildPersistence()`, so `frac` is the raw per-pixel
detection frequency with the threshold nowhere in it. Choosing a different
threshold is a **re-read of `frac`, not a re-export**. Only a change to the
*detector* — `k`, `analysisScaleM`, `detectionPols` — invalidates the asset.

So committing to STAGE A before the threshold is settled costs nothing that
cannot be recovered, and it is the only point at which the evidence for the
threshold becomes affordable.

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
