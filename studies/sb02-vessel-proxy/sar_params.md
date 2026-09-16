# Sentinel-1 SAR Vessel-Detection Parameter Specification
## SanctSound SB02, Stellwagen Bank NMS — satellite vs. acoustic ship presence

**Site:** 42.470793 N, -70.242940 W (Massachusetts Bay)
**Collection:** `ee.ImageCollection('COPERNICUS/S1_GRD')`, `instrumentMode == 'IW'`, 2018-11-12 → 2021-12-16
**Scene inventory (verified by query):** 278 scenes, all `transmitterReceiverPolarisation == ['VV','VH']`, `resolution == 'H'`.
Relative orbits: 142 (93), 40 (92), 62 (93). Pass: ASCENDING 93 / DESCENDING 185.
**Ground truth:** 15,255 SanctSound ship-detection intervals, 45.7% duty cycle.
**Radius sweep:** R ∈ {2, 5, 10, 20, 40} km.

---

## 0. Asset facts (verified) and their consequences

Verified against the Earth Engine STAC record for `COPERNICUS/S1_GRD`
(`https://storage.googleapis.com/earthengine-stac/catalog/COPERNICUS/COPERNICUS_S1_GRD.json`)
and ESA Sentinel-1 user documentation. `developers.google.com` was blocked; nothing below
rests on that domain.

| Fact | Value | Consequence for this spec |
|---|---|---|
| Bands | `VV`, `VH` (this collection), plus `angle` | Dual-pol logic available on all 278 scenes |
| Band units | **decibels (dB)**, `10*log10(σ⁰)` | **All CFAR statistics must be done in linear power, not dB.** Convert first. |
| Preprocessing applied by GEE | thermal-noise removal → radiometric calibration → **terrain correction (SRTM 30 m / ASTER)** | Already orthorectified; no extra geocoding needed. Over water the DEM is flat, so TC adds little error offshore. |
| Band `gsd` | 10 m for VV/VH; **20,000 m for `angle`** | `angle` is a coarse interpolation of the geolocation grid — fine for a per-scene mean incidence, useless for per-pixel correction |
| IW incidence range | 29.1° – 46.0° | Real inter-orbit backscatter differences; see §5 |
| IW max NESZ | −22 dB | Sets the low-wind noise floor; see §9 |
| IW GRDH pixel spacing / resolution | 10 m × 10 m spacing; ~20 m × 22 m resolution | A point target occupies ≥ 2×2 pixels; sets the minimum blob size in §6 |
| IW GRDH ENL | ≈ **4.4** (nominal; ESA) | Sets the CFAR constant; **measure it empirically anyway** (§2.4) |
| Linear-power sibling collection | `COPERNICUS/S1_GRD_FLOAT` | Documented in the S1_GRD collection description. Preferred input if verified available — avoids a dB→linear round trip. **Treat availability as unverified until you `.first()` it.** |
| Orbit local time | LTAN 18:00 → ASCENDING ≈ 18:00 local, DESCENDING ≈ 06:00 local | **Severe diel aliasing.** See §9.6. |

**Property names to use verbatim:** `relativeOrbitNumber_start`, `orbitProperties_pass`,
`platformHeading`, `transmitterReceiverPolarisation`, `instrumentMode`, `resolution_meters`,
`system:time_start`, `system:index`.

### 0.1 Geometry correction — please re-check the brief

The brief states the R = 40 km disc "reaches Cape Ann and Boston Harbor". Haversine distances
from 42.470793, -70.242940:

| Feature | Distance | Inside R=40 km? |
|---|---|---|
| Gloucester Harbor mouth / Eastern Point Light | 36.7 km | **yes, just** |
| Cape Ann (Rockport tip) | 38.5 km | **yes, just** |
| Race Point, Provincetown | 43.5 km | no |
| Minots Ledge Light | 48.0 km | no |
| Salem Sound | 48.6 km | no |
| Scituate Harbor | 49.4 km | no |
| Marblehead Neck | 50.5 km | no |
| Boston Light (Little Brewster I.) | 55.4 km | no |
| Deer Island, Boston Harbor | 60.3 km | no |
| **Boston Inner Harbor** | **65.9 km** | **no** |

So Cape Ann and the outer approaches to Gloucester **are** in the R=40 disc; Boston Harbor is
**not**, and neither is Provincetown. The land-masking requirement in §4 is real but its target
is the Cape Ann / Gloucester shoreline, not Boston. Masking is still specified generously enough
that a later extension to R=60–80 km would not need re-derivation.

### 0.2 Disc areas / pixel budgets

| R (km) | Disc area (km²) | 10 m pixels |
|---|---|---|
| 2 | 12.6 | 125,664 |
| 5 | 78.5 | 785,398 |
| 10 | 314.2 | 3,141,593 |
| 20 | 1,256.6 | 12,566,371 |
| 40 | 5,026.5 | 50,265,482 |

### 0.3 Run the detector ONCE per scene, not once per (scene, radius)

Detect over a single disc of radius **40.5 km** (40 km + 0.5 km so the background annulus in §2
has valid data at the disc edge). Emit a **per-detection long table** carrying each detection's
great-circle range to the hydrophone. Every R in the sweep is then a `filter(range_km <= R)` on
that table. Only the masked-water-area denominators (§8) must be computed per R. This cuts
compute ~5× and, more importantly, guarantees the five radii are derived from one identical
detection set, so differences across R are pure geometry and not detector noise.

---

## 1. Polarisation

**Primary: VH. Secondary/diagnostic: VV. Do not average or combine them into a single index.**

**Why VH is primary.** Over open water, the dominant sea-surface return is Bragg scattering,
which is strongly co-polarised. Ocean clutter is therefore much weaker in cross-pol (VH) than in
co-pol (VV), while ships — dominated by dihedral and multiple-bounce structures — depolarise and
retain substantial cross-pol return. The result is a higher **target-to-clutter ratio** in VH,
which is what a detector actually cares about. Absolute ship σ⁰ is higher in VV, but so is the
background, and the background rises much faster with wind (§9.1). Using cross-pol as the primary
ship-detection channel for dual-pol C-band over water is standard operational practice.

**Roles of VV:**
1. **Wind / sea-state proxy.** Mean VV σ⁰ over the masked-water disc is a monotone proxy for
   local wind speed at C-band. Emit it (`bg_mean_VV_dB`, §8) and use it as a covariate; it is
   more directly relevant than reanalysis wind because it is co-located and co-timed.
2. **Optional confirmation channel.** Run the identical CFAR on VV and record, per detection,
   whether it also passes in VV (`vv_confirmed`). **Do not require it by default** — VV
   confirmation raises precision but cuts recall exactly on the small/low-RCS vessels that are
   already marginal. Emit it as a flag so the strict variant is a post-hoc filter.
3. **VH noise-floor diagnosis.** When `bg_mean_VH_dB` approaches −22 dB (NESZ) but
   `bg_mean_VV_dB` is well above it, the VH channel is noise-limited on that scene and the VH
   CFAR result should be down-weighted or the scene flagged (§9.1).

**Do not** use a VH/VV ratio or a polarimetric span as the detection image. Both re-import VV's
wind-driven variance into the detection statistic and destroy the CFAR's clean single-channel
speckle model.

---

## 2. Detection algorithm — two-parameter CFAR

### 2.1 Choice: adaptive (CFAR), not fixed threshold

A **fixed** backscatter threshold is rejected. Ocean σ⁰ over this collection varies by roughly
15–20 dB with wind speed alone, plus several dB with incidence angle across the three relative
orbits (§5), plus sub-swath-to-sub-swath steps. Any fixed threshold is simultaneously too loose
on windy scenes and too tight on calm ones, and the resulting count time series would be
dominated by weather rather than by ships — fatal for a study whose whole point is correlating
counts against an acoustic signal that *also* responds to weather.

A local-statistics (CFAR) threshold is **self-normalising**: it absorbs wind, incidence angle,
sub-swath steps and orbit differences into the local background estimate. This is the single
most important design decision in the spec.

Keep a fixed-threshold run (**VH σ⁰ > −15.0 dB**, absolute) as a *secondary diagnostic only*,
to sanity-check that the CFAR is not chasing noise on calm scenes.

### 2.2 Work in linear power

```
lin = ee.Image(10.0).pow(img.select('VH').divide(10.0))   // σ⁰ linear power, unitless
```
(Or take `COPERNICUS/S1_GRD_FLOAT` directly if verified available.) Everything below is in
linear power. Report thresholds in dB only for human inspection.

### 2.3 Kernel geometry — exact sizes in metres

Square kernels, `ee.Kernel.square(radius, 'meters')`, evaluated at the native 10 m scale.

| Element | Radius (m) | Box size (m) | Box size (px) | Pixel count |
|---|---|---|---|---|
| Test pixel | 0 | 10 × 10 | 1 × 1 | 1 |
| **Guard window** (excluded from background) | **150** | 310 × 310 | 31 × 31 | 961 |
| **Background outer window** | **400** | 810 × 810 | 81 × 81 | 6,561 |
| **Background annulus** (outer − guard) | — | — | — | **5,600** |

**Why 150 m guard radius.** The guard must exclude the target's own energy from the background
estimate ("self-screening"). It must cover: the largest plausible vessel in the Boston approaches
(container ships / tankers up to ~300 m LOA → ±150 m if centred), plus azimuth smearing from
target motion (§7.1, up to ~1.4 km in the worst case, but only for the *displaced* image, not for
smear extent), plus sidelobe energy. A 310 m guard box accommodates a 300 m vessel exactly. This
is the parameter most likely to need adjustment; see §2.6.

**Why 400 m outer radius.** The annulus must contain enough independent looks that the background
mean is precise relative to the threshold offset. With 5,600 pixels at ENL 4.4, the fractional
standard error of the mean is ≈ 1/√(5600 × 4.4) ≈ **0.6%**, i.e. ≈ 0.03 dB — negligible against a
~7 dB threshold offset. It is also small enough that the wind field, incidence angle and sub-swath
gain are effectively constant across it (810 m spans ~0.15° of incidence angle).

**If GEE kernel cost is prohibitive** (81×81 `reduceNeighborhood` over 50 M pixels × 278 scenes
is the expensive step), reduce the outer radius to **300 m** (61×61 = 3,721 px, annulus 2,760 px,
standard error 0.9%). The statistical loss is negligible. Do **not** shrink the guard radius to
save compute.

### 2.4 Threshold rule

Because the background is speckle-limited and the intensity of an L-look product is
Gamma-distributed with σ = μ/√L, the two-parameter rule

```
T = μ_bg + k · σ_bg        with   σ_bg = μ_bg / √L
⇔  T = μ_bg · (1 + k/√L)                                   [linear power]
⇔  T_dB = μ_bg_dB + 10·log10(1 + k/√L)                     [dB, offset is constant]
```

**Do not estimate σ_bg as the sample standard deviation of the annulus.** A ship, a wake, a buoy
or an azimuth ambiguity inside the annulus inflates the sample σ and *raises* the threshold
exactly where you most want it low — the classic CFAR self-masking failure, and a real risk here
because vessels in a traffic lane travel in convoy. Instead:

- Estimate **only the background mean**, robustly, and derive σ from the known ENL.
- Robust mean: take the annulus **median** and de-bias it, `μ_bg = median / m(L)`, where m(L) is
  the median of a unit-mean Gamma(L). Computed values:

  | ENL L | median/mean m(L) | μ̂ = median × |
  |---|---|---|
  | 3.0 | 0.8914 | 1.1219 |
  | 4.0 | 0.9180 | 1.0893 |
  | **4.4** | **0.9253** | **1.0807** |
  | 5.0 | 0.9342 | 1.0705 |
  | 6.0 | 0.9450 | 1.0582 |

  In GEE: `lin.reduceNeighborhood(ee.Reducer.median(), kernel)`. If the median reducer over an
  81×81 kernel proves too slow, use the **arithmetic mean over the annulus with the top 2% of
  samples trimmed**, implemented as a mean-of-sums difference:
  ```
  S_out = lin.reduceNeighborhood(ee.Reducer.sum(), K400)
  S_in  = lin.reduceNeighborhood(ee.Reducer.sum(), K150)
  N_out = valid.reduceNeighborhood(ee.Reducer.sum(), K400)   // valid = water-mask as 0/1
  N_in  = valid.reduceNeighborhood(ee.Reducer.sum(), K150)
  mu_bg = S_out.subtract(S_in).divide(N_out.subtract(N_in))
  ```
  **The `N_out`/`N_in` counts must come from the water mask, not from the kernel size**, or every
  pixel within 400 m of masked land gets a biased background.

**ENL must be measured, not assumed.** 4.4 is the ESA nominal figure for IW GRDH and is a
reasonable starting value, but GEE's thermal-noise removal and terrain correction alter it, and
it differs between sub-swaths. **Tuning procedure:** pick 5 low-wind scenes (ERA5 10 m wind
< 5 m/s), draw three 5 km × 5 km boxes in open water ≥ 20 km from shore and away from the traffic
lanes, one in each IW sub-swath, and compute `L̂ = mean(lin)² / variance(lin)`. Use the median
of the resulting L̂ values. If L̂ differs from 4.4 by more than ~15%, use L̂ in the formula above
and re-derive k from §2.5.

### 2.5 Choice of k, with a false-alarm budget

For an L-look Gamma background, PFA = Q(L, L·(1 + k/√L)) (regularised upper incomplete gamma).
Expected false-alarm **pixels** per scene in the R = 40 km disc (50,265,482 pixels), at L = 4.4:

| k | T/μ | Threshold above background (dB) | PFA per pixel | FA pixels / scene @ R=40 |
|---|---|---|---|---|
| 4 | 2.91 | 4.63 | 2.1e-03 | 105,923 |
| 5 | 3.38 | 5.29 | 4.2e-04 | 20,922 |
| 6 | 3.86 | 5.87 | 7.8e-05 | 3,900 |
| 7 | 4.34 | 6.37 | 1.4e-05 | 695 |
| **8** | **4.81** | **6.82** | **2.4e-06** | **119** |
| 10 | 5.77 | 7.61 | 6.4e-08 | 3.2 |
| 12 | 6.72 | 8.27 | 1.6e-09 | 0.1 |

**Adopt k = 8.0 as the primary value.** Justification:

- It puts the threshold **6.8 dB above the local background**. Vessels in cross-pol typically sit
  10–25 dB above ocean clutter, so 6.8 dB retains large margin for detection — k = 8 is
  conservative on false alarms without being anywhere near the signal level of a real ship.
- It yields ~119 isolated false-alarm *pixels* per scene at R = 40 km. Combined with the
  ≥ 3-connected-pixel rule of §6, the surviving false **detections** are effectively zero under
  an independence assumption. (Speckle is spatially correlated over ~2 pixels, so real
  performance is worse than independent — which is precisely why k = 8 is preferred over k = 6,
  and why the min-size rule is not optional.)
- It is not so high that it depends on a precise ENL: at L = 3.0 the same k = 8 gives 384 FA
  pixels, at L = 6.0 it gives 43. The result is stable across the plausible ENL range, which
  matters because ENL is measured, not known.

**Sweep k ∈ {5, 6, 8, 10, 12} and report detection counts for each.** If the k-sweep curve of
count vs k has a plateau, the plateau is the physical vessel population and the rising limb below
it is false alarms — that plateau, not a theoretical PFA, is the real justification for the final
k. **This is an explicitly empirical choice and should be reported as one.** k = 8 is the
recommended prior, not a verified optimum for this site.

### 2.6 Parameters that must be tuned empirically (flagged)

| Parameter | Prior | How to tune |
|---|---|---|
| ENL L | 4.4 | Measure on calm scenes, §2.4 |
| k | 8.0 | k-sweep plateau, §2.5; then validate against AIS |
| Guard radius | 150 m | Inspect 20 large-vessel chips; confirm target energy is fully inside the guard box. Increase to 250 m if bright smears leak into the annulus. |
| Min blob size | 3 px | AIS matching, §6 |
| Land buffer | 1,000 m | Visual inspection at Cape Ann, §4 |

**Recommended external validation:** the NOAA Marine Cadastre national AIS archive covers US
waters including Massachusetts Bay over 2018–2021. Match detections to AIS positions interpolated
to the scene timestamp, within a 1.5 km tolerance (to absorb the moving-target azimuth
displacement of §7.1). This gives an empirical recall-vs-vessel-length curve and a false-alarm
rate, and is the only defensible way to fix k, the min blob size, and the claimed minimum
detectable vessel length. Report those as *measured* if this is done, and as *assumed* if not.

---

## 3. Speckle handling

**Recommendation: apply no speckle filter to the detection (test) image. Do not pre-smooth.**

The tension in the brief is real and resolves in favour of not filtering, for three reasons:

1. **The CFAR already does the speckle averaging, in the right place.** The background estimate
   averages ~5,600 pixels. Smoothing the *test* pixel as well is double-counting: it reduces the
   test statistic's variance but also reduces the target's peak, and for small targets (2–4
   pixels) a 3×3 or 5×5 filter mixes target energy with surrounding clutter and can cut the peak
   by 5–10 dB. Small vessels are exactly the marginal population; a pre-filter removes them.
2. **The min-size rule is a better speckle suppressor than a filter.** Isolated speckle spikes
   are removed by requiring ≥ 3 connected pixels (§6) at zero cost to peak amplitude.
3. **Filters that preserve edges (Lee, Refined Lee, Frost) preserve bright point targets poorly**
   and, worse, are nonlinear — they change the effective ENL in a spatially varying way, which
   invalidates the σ = μ/√L relation the threshold rests on.

**Where filtering IS allowed:**
- **Background path only.** If the median reducer is too slow and you use the trimmed-mean path,
  the annulus mean is itself a 5,600-pixel average — no further filtering needed.
- **Sensitivity run only.** As a robustness check, re-run the full pipeline with
  `focal_median(radius=15, units='meters')` (3×3 px) applied to the test image and report the
  delta in vessel counts. Expect counts to *fall*. If counts fall by more than ~20%, the count
  series is dominated by marginal 1–3-pixel targets and that fact must be reported prominently,
  because it means the count is sensitive to sea state rather than to traffic.

**Do not** apply a multi-temporal speckle filter. The scenes are 4+ days apart and vessels move;
temporal averaging would erase every target of interest.

---

## 4. Land and clutter masking

Four masks, applied in sequence. All are computed once and cached as assets.

### 4.1 Land mask — `ESA/WorldCover/v200`

```
wc   = ee.Image('ESA/WorldCover/v200').select('Map')      // 10 m, verified
// class 80 = 'Permanent water bodies' (verified from STAC class table)
land = wc.neq(80).unmask(0)    // unmapped open ocean -> 0 (= water)
```
**Verified:** asset ID `ESA/WorldCover/v200`, band `Map`, 10 m gsd, class 80 = Permanent water
bodies. The v200 epoch is 2021, inside the study period.

**Caveat, flagged:** WorldCover's footprint is land-centric; open ocean beyond its coastal tiles
is unmapped. `.unmask(0)` treats unmapped as water, which is correct for open Massachusetts Bay
but **must be verified visually** at the R = 40 km disc edge before trusting it. If the footprint
edge falls inside the disc and produces a spurious land/water boundary, substitute a manually
digitised Massachusetts Bay coastline polygon.

**Rejected alternative:** `JRC/GSW1_4/GlobalSurfaceWater`. Verified 30 m, bands
`occurrence/max_extent/...` — but the JRC product maps **inland** surface water and does not map
oceans, so it cannot define a marine coastline here. Do not use it for this purpose.

### 4.2 Land buffer — 1,000 m

```
landCoarse = land.reproject(crs, null, 100)                     // 100 m, for cheap dilation
landBuf    = landCoarse.reduceNeighborhood(
                 ee.Reducer.max(), ee.Kernel.circle(1000,'meters'))
water      = landBuf.eq(0)
```
Dilating at **100 m** scale (circle radius 1,000 m = 21×21 cells) instead of 10 m
(201×201 cells) is ~100× cheaper with ~100 m of edge granularity, which is well inside the
tolerance of a 1,000 m buffer.

**Why 1,000 m:** absorbs (a) WorldCover coastline classification error, (b) terrain-correction
geolocation error where the DEM is poor at the shore, (c) land returns spilling seaward through
range/azimuth sidelobes and near-shore ambiguities, and (d) the surf/breaker zone, which produces
bright, vessel-like returns.

**Cost, stated plainly:** a 1,000 m buffer removes all genuine vessels within 1 km of shore. At
this site that affects only the Cape Ann / Gloucester sector at R = 40 km, ~37 km from the
hydrophone — acoustically marginal anyway. **Run 500 m and 2,000 m as sensitivity variants** and
confirm the R = 40 count is insensitive to the choice (it should be; if it is not, nearshore
clutter is leaking).

### 4.3 Harbour / fixed-infrastructure exclusion polygons

Docks, piers and moored vessels are bright, static, and must not count. At R = 40 km the only
harbour inside the disc is **Gloucester** (mouth at 36.7 km; the inner harbour is just outside).
Exclude, in addition to the 1,000 m land buffer:

| Polygon | Approx. bounding box (lon_min, lat_min, lon_max, lat_max) | Rationale |
|---|---|---|
| Gloucester Harbor / Annisquam | −70.70, 42.57, −70.63, 42.64 | Piers, fishing fleet at mooring |
| Rockport / Sandy Bay | −70.63, 42.64, −70.58, 42.69 | Breakwater, moorings |

**These boxes are approximate and must be visually verified against a high-res basemap before
use.** They are a coarse net; the persistence mask in §4.4 is the precise instrument.

If the study is later extended past R = 40 km, add Salem Sound (−70.92, 42.50, −70.80, 42.58) and
Boston Harbor (−71.06, 42.30, −70.90, 42.40).

### 4.4 Persistent-target mask (data-driven) — the important one

A moving vessel has essentially zero probability of re-occupying the same 10 m pixel on many
independent overpasses. A navigational buoy, wreck, breakwater, aid-to-navigation, or permanently
moored hull does. So:

1. Run the §2 detector on all scenes **of a given relative orbit** (93, 92, 93 scenes) with a
   deliberately loose k = 5, producing a binary detection image per scene.
2. Sum across scenes → `persistence_count`; divide by scene count → `persistence_frac`.
3. **Mask any pixel with `persistence_frac ≥ 0.20`** (≈ 19 of 93 scenes).
4. **Dilate the persistence mask by 30 m** (circle kernel, radius 30 m) to absorb geolocation
   jitter between scenes.

**Compute this per relative orbit, not pooled**, because geolocation and geometry are
orbit-specific; a fixed structure lands on slightly different pixels from different orbits.

**Threshold justification:** 0.20 is conservative. With 93 scenes over 3 years and a busy traffic
lane, the busiest single 10 m pixel in a lane might see a vessel on a few percent of passes;
20% is far above that and far below the ~100% a fixed structure achieves. **This is an assumed
value — check the histogram of `persistence_frac`; it should be strongly bimodal with a clear gap.
Set the threshold in the gap.** If it is not bimodal, something is wrong and the mask should not
be applied blind.

This mask is preferable to a hand-curated list of structures because it finds unknown ones (wrecks,
seasonal aquaculture, the moorings of the SanctSound array itself) without prior knowledge.

### 4.5 Order of operations

Apply masks **before** computing the background statistics of §2, and propagate the mask into the
`N_out`/`N_in` counts, so that near-mask pixels get an unbiased (if noisier) background. Do not
apply masks only to the final detections.

---

## 5. Incidence angle / orbit normalisation

**Decision: do NOT apply an explicit incidence-angle normalisation. STRATIFY by
`relativeOrbitNumber_start`, and let the CFAR absorb the angle dependence.**

### 5.1 Why no explicit normalisation

The incidence-angle dependence of ocean σ⁰ is wind-speed and wind-direction dependent (a
cosine-power or a GMF-style model, not a fixed exponent), so any fixed correction such as
`σ⁰_norm = σ⁰ · (cos θ / cos θ_ref)^n` is wrong by a scene-dependent amount. Worse, the `angle`
band's gsd is **20,000 m** — a coarse interpolation of the geolocation grid — so per-pixel
correction is not even supported by the data. Applying a bad normalisation is strictly worse than
applying none, because it injects a wind-correlated bias into the detection statistic.

The CFAR of §2 is **already an incidence-angle normaliser**, and a better one: the local
background is estimated over 810 m, across which incidence angle changes by ~0.15°. Whatever the
true angle dependence is, it is constant over the background window, so it cancels exactly in
`T = μ_bg · (1 + k/√L)`. **This is the reason the detection statistic transfers across orbits at
all, and the main reason a fixed threshold was rejected in §2.1.**

### 5.2 Why stratify anyway

Stratification costs nothing and protects the conclusion:

- Emit `relativeOrbitNumber_start` and `orbitProperties_pass` on every record (§8).
- Compute the count-vs-acoustic relationship **separately for orbits 142, 40 and 62** (93/92/93
  scenes — well balanced, adequate n for each).
- If the three give consistent slopes, the CFAR has transferred and you may pool with an orbit
  fixed effect. If they disagree, the disagreement is the finding and pooling would have hidden it.
- Also emit the **disc-mean incidence angle** (`inc_mean_deg`, from the `angle` band reduced over
  the masked-water disc) per scene. Across IW's 29.1°–46.0° range, three orbits will land at
  three distinct mean angles; report them. Include `inc_mean_deg` as a covariate in the final
  regression so any residual angle dependence is visible rather than silent.

### 5.3 Ascending/descending is confounded, not just a covariate

ASCENDING 93 / DESCENDING 185 corresponds to ~18:00 and ~06:00 local time respectively (S1 LTAN
is 18:00). Pass direction is therefore **perfectly confounded with time of day**, and vessel
traffic has a strong diel cycle. **Never interpret an ASC/DESC difference in vessel count as a
geometry effect.** See §9.6.

---

## 6. Minimum target size

**Primary: n_min = 3 connected pixels = 300 m² = a blob of roughly 17 m × 17 m equivalent area.**
Connectivity: **8-connected** (`ee.Image.connectedPixelCount(maxSize=256, eightConnected=true)`,
then `.gte(3)`).

**Why 3, in pixel terms.** IW GRDH has 10 m pixel *spacing* but ~20 m × 22 m *resolution*. An
ideal point target is therefore spread by the impulse response over roughly 2×2 = 4 pixels at the
sampling grid before any of it is thresholded. A minimum of 3 pixels is just below that natural
footprint: it rejects isolated 1- and 2-pixel speckle spikes (which is its whole job, see §2.5)
while not requiring the target to be larger than the resolution cell. Requiring 4+ would start
rejecting genuine point-like targets whose energy happens to split unfavourably across the grid.

**Report the full size sweep.** Emit counts at n_min ∈ {1, 2, 3, 5, 10} pixels
(= 100, 200, 300, 500, 1,000 m²) from the same detection run, so the sensitivity is visible and
the choice is auditable.

### 6.1 Corresponding vessel length — stated with its uncertainty

**Blob area does not equal vessel area.** A blob is the region where returned power exceeds the
local threshold; for a strong scatterer it can exceed the vessel's physical footprint (sidelobes,
smearing) and for a weak one it is smaller (only the brightest structures pass). So converting
300 m² to a length is an approximation, not a measurement.

With that caveat: for C-band IW GRD over water at a ~7 dB CFAR offset, the practical minimum
detectable vessel is commonly around **20–30 m length overall in calm seas** and degrades to
roughly **40–60 m LOA in rough seas** (wind ≳ 12 m/s). **These figures are generic
established-practice ranges, not a measurement at this site, and must be verified against AIS
(§2.6) before being quoted in the paper.**

### 6.2 What the method will MISS — and why that is acceptable here

| Vessel class | Typical LOA | Detected? |
|---|---|---|
| Kayaks, small RIBs, jet skis | < 8 m | **Never** |
| Recreational powerboats, small sailboats | 8–20 m | **Rarely** — sailboats especially (low RCS, non-metallic hulls, no large dihedrals) |
| Lobster boats, small inshore fishing | 10–20 m | **Mostly missed**, calm seas only |
| Whale-watch vessels, larger fishing, pilot boats | 20–50 m | **Marginal** — detected in calm, missed in rough |
| Ferries, coastal freighters, tugs/barges | 50–100 m | **Reliably** |
| Cargo, tanker, container, cruise | 100–350 m | **Always** |

**This miss profile is well-matched to the research question.** Low-frequency (tens to hundreds of
Hz) radiated ship noise that propagates tens of km is dominated by large, slow-turning propellers
on commercial vessels — precisely the class the method detects reliably. The missed class
(< 20 m planing hulls) radiates mostly at higher frequencies and is far less audible at the
hydrophone at range. **However**: if the SanctSound acoustic detector is sensitive to nearby small
craft, those will appear as acoustic-positive / SAR-negative events and will depress the apparent
skill of the satellite predictor. This asymmetry must be stated in the results, and if the
acoustic detector's frequency band is known it should be reported alongside this table.

---

## 7. Azimuth ambiguities and moving-target displacement

### 7.1 Moving-target azimuth displacement — the larger risk, and it is material

A target with a line-of-sight (range-direction) velocity component `Vr` is displaced along-track
in the SAR image by `Δx ≈ (Vr / V_g) · R_slant`. With S1's ground/beam velocity ≈ 6,700 m/s and
slant range 800–950 km:

| Vr (m/s) | ≈ knots (if fully in range direction) | Δx at R_slant = 850 km |
|---|---|---|
| 1 | 2 | 127 m |
| 2 | 4 | 254 m |
| 5 | 10 | 634 m |
| 10 | 19 | **1,269 m** |

**Consequence for the radius sweep.** A vessel underway at 15–20 kn heading roughly across the
radar look direction can be imaged up to **~1.3 km from its true position**. That is 65% of the
R = 2 km disc radius and 26% of R = 5 km. **The R = 2 km and R = 5 km bins are therefore
substantially blurred by SAR physics alone, independent of any detector error**, while R = 20 and
R = 40 km are essentially unaffected. Report this explicitly: if the R-sweep shows the acoustic
correlation peaking at a small radius, part of the apparent width of that peak is imaging
displacement, not acoustic footprint. Do not over-interpret structure below ~2 km.

This also fixes the AIS matching tolerance at **1.5 km** (§2.6) rather than a few hundred metres.

### 7.2 Azimuth ambiguities (ghost targets)

`Δ_az = λ · R_slant · PRF / (2 · V_sat)`, with λ = 0.0555 m, V_sat ≈ 7,100 m/s, S1 IW PRF in the
~1,450–1,850 Hz range depending on sub-swath:

| PRF (Hz) | R_slant 800 km | R_slant 900 km |
|---|---|---|
| 1,450 | 4.53 km | 5.10 km |
| 1,650 | 5.16 km | 5.80 km |
| 1,850 | 5.78 km | 6.50 km |

So first-order ghosts appear **displaced ~4.5–6.5 km in azimuth (roughly along-track, i.e. near
N–S for S1's ~98° inclination)** from a bright parent. **The PRF values are approximate and the
exact per-sub-swath PRF was not verified for this environment — treat the 4.5–6.5 km figure as an
order-of-magnitude estimate and confirm it empirically from the data (see mitigation 3).**

**Is it material here? Yes, moderately.** The S1 IW azimuth ambiguity-to-signal ratio is
specified around −20 to −25 dB. A large ship at +25 dB above clutter therefore produces a ghost
at roughly 0 to +5 dB above clutter — usually *below* the 6.8 dB CFAR threshold, but a very bright
tanker or a corner-reflector-like superstructure at +30 dB can push a ghost to +5 to +10 dB, which
**will** pass. The Boston traffic separation scheme runs through this area, so large bright parents
are common.

**Mitigations, in order of practicality:**

1. **Flag, do not delete (default).** For each detection, search for any other detection with
   peak ≥ 15 dB brighter lying at 3.5–7.5 km along the azimuth direction (azimuth direction from
   `platformHeading`, available as a scene property). Set `ambiguity_suspect = 1`. Emit
   `n_ambiguity_suspect` per scene (§8) and report counts with and without them.
2. **Amplitude test.** Require a candidate's peak-to-background ratio to exceed the parent's by
   more than −18 dB before rejecting it as a ghost; genuine small vessels near a big one should
   not be discarded.
3. **Empirical confirmation.** Build a histogram of along-azimuth separations between every
   detection pair within 10 km. A real ambiguity population appears as a spike at the predicted
   4.5–6.5 km. If no spike exists, ambiguities are not a material contaminant at this threshold
   and mitigation 1 can be dropped. **Do this before deciding.**

**Why flagging is enough for this study.** The dependent variable is acoustic ship presence.
Ghosts are *caused by* real ships, so they inflate counts but essentially never create a
false *presence*. For a presence/absence analysis they are harmless; for a count-magnitude
analysis they add a multiplicative bias that correlates with large-vessel traffic. Report both.

---

## 8. Output schema

Three tables. All exported as CSV / `ee.FeatureCollection`.

### 8.1 `detections_long` — one row per detected blob (radius-independent)

Produced once per scene over the R = 40.5 km disc. This is the primary product; the per-radius
table is derived from it.

| Field | Type | Units | Notes |
|---|---|---|---|
| `scene_id` | string | — | `system:index` |
| `datetime_utc` | ISO8601 | UTC | from `system:time_start` |
| `det_id` | int | — | unique within scene |
| `lon`, `lat` | double | degrees | blob centroid |
| `range_km` | double | km | great-circle to 42.470793, −70.242940 |
| `bearing_deg` | double | deg true | from hydrophone to detection |
| `area_m2` | double | m² | `n_pixels × 100` |
| `n_pixels` | int | px | connected-component size |
| `peak_vh_db` | double | dB | max VH σ⁰ in blob |
| `mean_vh_db` | double | dB | mean VH σ⁰ in blob |
| `bg_local_vh_db` | double | dB | local background μ_bg at blob centroid |
| `pcr_db` | double | dB | **peak-to-clutter ratio** = `peak_vh_db − bg_local_vh_db` — the per-detection confidence proxy |
| `peak_vv_db` | double | dB | for dual-pol diagnostics |
| `vv_confirmed` | 0/1 | — | blob also passes CFAR in VV |
| `ambiguity_suspect` | 0/1 | — | §7.2 mitigation 1 |
| `dist_to_land_km` | double | km | distance to nearest masked land |
| `k_used` | double | — | so the k-sweep lives in one table |

### 8.2 `scene_radius` — one row per (scene × radius), 278 × 5 = 1,390 rows

| Field | Units | Notes |
|---|---|---|
| `scene_id`, `datetime_utc` | — | join key |
| `radius_km` | km | 2 / 5 / 10 / 20 / 40 |
| `relative_orbit` | — | `relativeOrbitNumber_start` (142/40/62) |
| `pass` | — | `orbitProperties_pass` |
| `platform_heading_deg` | deg | for §7.2 |
| `inc_mean_deg` | deg | disc-mean of `angle` band |
| **`n_vessels`** | count | **primary covariate**; min blob 3 px, k = 8 |
| `n_vessels_ge_1000m2` | count | large-vessel subset (≥ 10 px) |
| `n_vessels_vv_conf` | count | strict dual-pol variant |
| `n_vessels_nmin1/2/5/10` | count | size sweep, §6 |
| `n_vessels_k5/k6/k10/k12` | count | k sweep, §2.5 |
| `total_bright_area_m2` | m² | Σ `area_m2` |
| `max_target_area_m2` | m² | largest blob — proxy for largest vessel |
| `mean_target_area_m2` | m² | — |
| `max_pcr_db` | dB | brightest detection's peak-to-clutter ratio |
| `mean_pcr_db` | dB | **detection-confidence proxy** for the scene |
| `nearest_det_km` | km | range to closest detection — highly relevant to the acoustic question; set to `radius_km` (censored) if none |
| **`water_area_km2`** | km² | **unmasked water area actually searched** — REQUIRED, varies with R and with the masks of §4 |
| **`vessel_density_per_1000km2`** | count/1000 km² | `n_vessels / water_area_km2 × 1000` — **use this, not raw count, as the primary covariate**, because it is robust to partial scene coverage and to mask-area differences across R |
| `disc_coverage_frac` | 0–1 | fraction of the disc actually inside the scene footprint. **Filter to ≥ 0.98 or model it.** Not all 278 scenes are guaranteed to cover the full 40 km disc — the collection was filtered on a point, not a disc. **Verify this before analysis.** |
| `n_ambiguity_suspect` | count | §7.2 |
| `bg_mean_vh_db` | dB | disc-mean VH background — **wind/sea-state proxy and low-wind noise-floor flag** |
| `bg_p90_vh_db` | dB | 90th percentile of background — clutter tail |
| `bg_mean_vv_db` | dB | co-pol wind proxy |
| `enl_est` | — | ENL measured on this scene's masked water |
| `era5_wind_10m_ms` | m/s | `ECMWF/ERA5/HOURLY`, `u_component_of_wind_10m`/`v_component_of_wind_10m`, nearest hour, at site (verified: bands exist, coverage 1940–present) |
| `era5_wind_dir_deg` | deg | relative wind direction matters for clutter |
| `era5_swh_m` | m | `significant_height_of_combined_wind_waves_and_swell` (verified band) |
| `noise_floor_flag` | 0/1 | 1 if `bg_mean_vh_db < -21.0` (within 1 dB of NESZ) — §9.1 |

### 8.3 `scene_meta` — one row per scene (278 rows)

Full S1 metadata dump for provenance: all properties in §0, plus mask version IDs, k, ENL,
kernel radii, and the pipeline git hash. Without this the numbers are not reproducible.

---

## 9. Known failure modes

### 9.1 Low wind — false alarms and a noise floor (the calm-sea failure)

Below roughly **3 m/s** at 10 m, the sea surface lacks the short Bragg waves that produce C-band
return. VH σ⁰ falls toward or below the IW **NESZ of −22 dB (verified ESA spec)**. Two things then
break:

- The scene becomes **noise-limited**: the "background" is thermal noise, not clutter, and GEE's
  thermal-noise removal leaves residual sub-swath steps and near-zero/clipped values. The local μ_bg
  collapses, so `T = μ_bg(1 + k/√L)` collapses with it and *anything* passes.
- Conversely, real detections become easier (very high contrast), so counts can spike for a
  genuine reason too — the two are hard to separate.

**Mitigation (implement this):** floor the background,
```
mu_eff = mu_bg.max(NESZ_linear)   with NESZ ≈ 10^(-22/10) = 6.31e-3
```
and set `noise_floor_flag = 1` when `bg_mean_vh_db < -21.0`. **The −22 dB figure is the ESA
*maximum* NESZ specification, not the actual value at this site's incidence angles — measure the
true floor empirically from the darkest calm-scene water in each sub-swath and use that.**
Report results with and without noise-floored scenes.

### 9.2 High wind — missed vessels (the rough-sea failure)

Above roughly **12–15 m/s**, ocean clutter rises steeply, breaking waves and whitecaps produce
bright transient returns, and target-to-clutter ratio collapses. Expect the minimum detectable
vessel to degrade from ~20–30 m to ~40–60 m LOA (§6.1). The CFAR will *not* produce a false-alarm
explosion (it renormalises), but it will silently **miss** vessels. Counts fall.

### 9.3 The wind confound — the most serious threat to the study's conclusion

**Wind degrades the satellite detector and simultaneously raises the acoustic ambient noise
floor at the hydrophone.** Both the predictor and the ground truth respond to the same
meteorological variable, in the same direction (both lose sensitivity). This can:
- **manufacture a correlation** where none exists (both detectors go quiet in storms), or
- **mask a real one** (satellite misses vessels that are acoustically loud and detected).

This is not a nuisance — it is a structural confound. **Mandatory:** include
`era5_wind_10m_ms`, `era5_swh_m` and `bg_mean_vh_db` as covariates in the final model, and
additionally report the relationship restricted to a **benign-sea-state subset**
(`era5_wind_10m_ms` between 3 and 10 m/s), where both detectors are near their best. If the
satellite–acoustic relationship holds in the benign subset, it is real; if it only appears in the
pooled data, it is probably the weather.

### 9.4 Non-vessel bright targets

- **Navigational buoys, aids to navigation, wrecks, moorings, the SanctSound mooring itself** —
  handled by §4.4 persistence masking; this is why that mask matters.
- **Sea ice / brash ice**: Massachusetts Bay ices rarely and mostly inshore, but January–February
  nearshore ice has high, vessel-like cross-pol return. Flag DJF nearshore detections for review.
- **Rain cells**: C-band is less rain-sensitive than X-band, but heavy convective cells produce
  both dark (damping) and bright (scattering) patches. Expect a few summer scenes affected.
- **Internal waves, surface slicks, current fronts**: produce *dark* and *bright* banding.
  Bright bands are broad and low-contrast, so the ≥ 3-px + 6.8 dB rule usually rejects them, but
  they can bias the local background.
- **Wind streaks / gust fronts**: modulate μ_bg over hundreds of metres, i.e. at the scale of the
  background window — this is the main reason the 400 m outer radius should not be enlarged much.

### 9.5 Detector-level failures

- **CFAR self-masking in convoy.** Two large vessels within 400 m of each other each raise the
  other's background. The robust median estimator of §2.4 largely handles this; dense clusters
  (harbour anchorages) would not be, but no such anchorage is inside R = 40 km.
- **Vessel splitting.** A long ship with a bright bow and bright superstructure but a dim midships
  can split into 2+ blobs, inflating `n_vessels` and deflating `max_target_area_m2`. Mitigation:
  merge blobs whose centroids are within **100 m** before counting, and report the merged and
  unmerged counts. **This merge distance is assumed and should be checked against AIS.**
- **Anchored vs underway is invisible to SAR.** A drifting or anchored ship looks identical to one
  making 18 kn. The acoustic detector responds mostly to *underway* vessels. This is an
  irreducible asymmetry and a likely source of SAR-positive / acoustic-negative events. There is
  a partial signal — an anchored ship has Vr ≈ 0 and so shows no azimuth displacement, and no
  wake — but extracting it is out of scope; simply state the limitation.
- **Wakes are not used.** Wake detection could confirm "underway" but requires Radon/Hough
  analysis not specified here. Noted as a possible extension.

### 9.6 Sampling design — the failure mode that is not about SAR at all

- **278 scenes over 1,130 days ≈ one every 4.1 days.** Each scene is an instantaneous snapshot
  (integration time ≪ 1 s), while the acoustic ground truth is 15,255 intervals at a 45.7% duty
  cycle. The matching rule must be stated: **recommend matching the acoustic state within
  ±30 min of the scene timestamp**, and reporting sensitivity at ±15 min and ±2 h. The ±30 min
  window is a choice, not a derived value.
- **Diel aliasing is severe.** S1's LTAN is 18:00, so all 278 samples fall at ~06:00 local
  (185 DESCENDING) or ~18:00 local (93 ASCENDING). **The study therefore measures the
  satellite–acoustic relationship at two times of day only**, and both are near commercial-traffic
  transitions. Any diel structure in vessel traffic or in acoustic propagation (diurnal
  sound-speed profile changes) is entangled with pass direction (§5.3). State this as a scope
  limitation on the deliverable: the measured acoustic footprint is a 06:00/18:00 footprint.
  The exact local overpass times should be **computed from the actual `system:time_start` values**
  rather than assumed from LTAN.
- **Seasonal aliasing.** 12-day repeat × 3 orbits is regular, so the sample is roughly uniform in
  season — good — but check for gaps (S1B failed in Dec 2021, just after the study window closes;
  earlier outages exist). Plot the scene timeline before analysis.
- **Multiple-comparison risk in the R sweep.** Five radii × several detector variants is a lot of
  tests. Pre-register the primary specification (R sweep on `vessel_density_per_1000km2`,
  k = 8, n_min = 3 px, VH only, no VV confirmation) and treat everything else as sensitivity
  analysis.

---

## 10. Implementation checklist

1. Verify `COPERNICUS/S1_GRD_FLOAT` availability; if present, use it and skip the dB→linear step.
2. Verify all 278 scenes' `disc_coverage_frac ≥ 0.98` at R = 40 km. Flag or drop the rest.
3. Build and cache masks: WorldCover land + 1,000 m buffer, harbour polygons, then the
   per-relative-orbit persistence mask (requires one loose-k pre-pass over all scenes).
4. Measure ENL per sub-swath on 5 calm scenes; measure the true VH noise floor on the same.
5. Run the CFAR detector once per scene over the R = 40.5 km disc → `detections_long`.
6. Derive `scene_radius` by range-filtering; compute `water_area_km2` per R.
7. Join ERA5 wind/wave; join SanctSound acoustic state at ±30 min.
8. Analyse stratified by `relativeOrbitNumber_start` first; pool only if consistent.
9. Report the benign-sea-state subset alongside the full result (§9.3).

---

## 11. Summary of adopted values

| # | Parameter | Value |
|---|---|---|
| 1 | Primary polarisation | **VH** (VV = wind proxy + optional confirmation flag) |
| 2 | Detector | **Two-parameter CFAR in linear power**, robust local background |
| 2 | Guard window | **150 m radius** (310 × 310 m, 31 × 31 px) |
| 2 | Background window | **400 m radius** (810 × 810 m, 81 × 81 px), annulus ≈ 5,600 px |
| 2 | Threshold | **T = μ_bg · (1 + k/√L)**, k = **8.0**, L = **4.4** (measure) → **+6.8 dB** |
| 2 | Background estimator | annulus **median** ÷ 0.9253, or 2%-trimmed mean |
| 3 | Speckle filter | **none on the test image**; 3 × 3 median as a sensitivity run only |
| 4 | Land mask | `ESA/WorldCover/v200` band `Map`, class ≠ 80, **buffered 1,000 m** |
| 4 | Clutter mask | Gloucester/Rockport polygons + **persistence ≥ 20% of same-orbit scenes**, dilated 30 m |
| 5 | Incidence angle | **no normalisation**; **stratify by `relativeOrbitNumber_start`**; emit `inc_mean_deg` |
| 6 | Minimum target | **3 connected pixels = 300 m²**, 8-connected; ≈ 20–30 m LOA min vessel (calm) |
| 7 | Ambiguities | flag at **3.5–7.5 km along-azimuth** from a ≥ 15 dB brighter parent; confirm empirically |
| 7 | Moving-target displacement | up to **~1.3 km** — do not over-interpret R = 2 and R = 5 km |
| 8 | Primary covariate | **`vessel_density_per_1000km2`**, not raw count |
| 9 | Mandatory covariates | ERA5 10 m wind, significant wave height, `bg_mean_vh_db` |

**Values that are assumptions requiring empirical tuning, not derived facts:**
ENL (4.4), k (8.0), guard radius (150 m), minimum blob size (3 px), land buffer (1,000 m),
persistence threshold (0.20), blob-merge distance (100 m), acoustic matching window (±30 min),
minimum detectable vessel length (20–30 m calm / 40–60 m rough), IW PRF and hence the
4.5–6.5 km ambiguity displacement. Each has a stated tuning procedure above. None should be
reported in the paper as a fixed property of the method without that tuning being done.
