# Changelog — STEMGeoHS Marine

Version history for `GeoMarineAnalysisV10_161.js`.

The **running version** is declared once, in code, as `TOOL_VERSION` near the
top of MODULE A. It is deliberately not repeated here: three consecutive
rounds of this file shipped with a stale marker because the number was typed
in more than one place. The entries below keep their own version numbers
because they are history, not identity.

These entries were moved out of the script's header comment so the `.js`
stays under GitHub's 1 MB in-browser render limit. Nothing was edited in the
move — including the corrections and withdrawn claims, which are part of the
record.

```text
v10.172 S21 NEW: real passive-acoustic biophony at Florida Keys FK01, as a diel
  RATIO - plus a deliberate, documented REFUSAL to route it through the
  climatology or CSD/AC1 path, and a gate audit of the entire SanctSound archive.

  EVIDENCE BASIS. Every number below was derived THIS session from the live
  public bucket gs://noaa-passive-bioacoustic (SanctSound, NOAA/NPS), or from a
  Node harness executing the shipped functions. Earth Engine itself was NOT run:
  S21 adds no Earth Engine calls at all, being synchronous JS over an embedded
  table exactly as S19 is. Nothing here is asserted from memory.

  THE DATA. Site FK01 (24.43313, -81.93068), Florida Keys National Marine
    Sanctuary. 8 deployments, 16,440 hourly third-octave (TOL_1h) rows across 30
    bands TOL_25..TOL_20000, 2018-12-18T19:00Z .. 2022-06-15T12:00Z. Reduced
    outside Earth Engine - the Code Editor sandbox has no fetch/XHR - into 28
    monthly values. 2020-10 is dropped on purpose: deployment 06 started
    2020-10-30, leaving 7 crepuscular and 4 trough hours, below the >=20-hour
    floor each side of the ratio needs.

  WHY THE STATE VARIABLE IS A RATIO AND NOT A SOUND LEVEL. Deployment 11
    (2022-03..2022-06) carries a low-frequency instrument artifact. Comparing the
    same calendar month across years, 2022-05 minus 2021-05, per band:
      TOL_25  +16.95 dB   TOL_63  +8.89 dB   TOL_125 +2.90 dB
      TOL_500  -0.44 dB   TOL_2000 -1.10 dB  TOL_20000 -1.02 dB
    A monotonic rise that grows as frequency falls and vanishes above 500 Hz is
    flow noise or mooring strum - a vessel lifts 63-500 Hz broadly. An
    anthropophony variable built on the TOL_63/TOL_125 vessel indicator bands
    would have reported a ~5 dB traffic increase that did not happen, so NO
    anthropophony variable is shipped. The diel ratio is immune twice over: it
    lives in the unaffected 2-20 kHz snapping-shrimp band, and it is a
    DIFFERENCE, so any constant per-deployment calibration offset cancels.
    Confirmed empirically - per-deployment mean ratio dep01 +1.56, dep04 +1.75,
    dep06 +1.51, dep08 +1.82, dep10 +1.79, dep11 +1.74 dB. In the ratio,
    deployment 11 is unremarkable.

  DEFINITION. Power-domain mean of TOL_2000..TOL_20000 over crepuscular hours
    (10,11,12,22,23 UTC ~ local dawn 05-07 and dusk 17-18) MINUS the same over
    the night trough (06,07,08,09 UTC ~ local 01-04), in dB. Shrimp at FK01 are
    crepuscular: pooled over the record the band peaks at local 06 (110.43 dB)
    and local 18 (110.41 dB), troughing at local 03 (108.22 dB). All dB averaging
    is done in the POWER domain, 10*log10(mean(10^(x/10))) - a decibel is a
    logarithm and an arithmetic mean of dB is simply wrong.

  WHY SEASONAL (HIRSCH-SLACK) MANN-KENDALL AND NOT THE EXISTING ENGINE. Plain
    mannKendallTest() over the 29 pooled monthly shrimp levels returns
    tau=+0.227, p=0.0878 - close enough to read as an emerging trend. It is
    SAMPLING ALIASING. FK01's deployments are seasonally unbalanced: July,
    August, September, October and November each occur in exactly ONE year, and
    the shrimp band runs ~4 dB hotter in summer, so the later-weighted summer
    coverage tilts the pooled series upward by itself. Comparing like with like
    destroys it - April across 2019-2022 reads 109.32, 110.29, 109.21, 109.33 dB,
    tau exactly 0.000. One more summer deployment would have pushed that pooled
    p under 0.05 and meant nothing at all. The seasonal test ranks a value ONLY
    against the same calendar month in other years, so it cannot make that
    mistake, and it builds no climatology and imputes nothing.
    DISCLOSED: the inter-season covariance correction for serially dependent
    seasons is NOT applied - at 3-4 years per season it is estimated far too
    noisily to help, and omitting it is the conventional choice.

  FALSE-POSITIVE RATE, MEASURED. 4000 Monte Carlo draws of pure uniform noise at
    FK01's exact design (7 seasons of 3,4,4,4,3,3,3 years), run in Node against
    the shipped function: 4.05% returned p<0.05, against a nominal 5%. Slightly
    conservative, which is the safe direction. Contrast S17's SNR>=2.0 rule,
    calibrated to no false-positive rate at all (see v10.145).

  THE GATES STAY STRICT, AND THEY REFUSE THIS SITE. FK01 has 28 valid months, so
    it passes CLIM_MIN_TOTAL_SAMPLES=26 - but only 7 of 12 calendar months reach
    CLIM_MIN_SAMPLES_PER_MONTH=3, against CLIM_MIN_DISTINCT_MONTHS=9.
    computeUsableClimatology() refuses it, correctly. S21 never asks it. No gate
    constant was changed, loosened, or special-cased.

  CSD/AC1 IS UNREACHABLE FOR THE WHOLE ARCHIVE, NOT JUST THIS SITE. All 28
    SanctSound sites were gate-checked against real hourly coverage this session
    (>=72 h before a month counts). CSD_AC1_MIN_POOLED_MONTHS=48 cannot be met by
    any of them, because SanctSound spans 2018-11..2022-06 - 44 months end to
    end. Best in archive: SB02 Stellwagen Bank, 44 months, ZERO gaps, 12/12
    calendar months - still 4 short of the floor. MB01 Monterey Bay (36 months,
    11/12 calendar months, one 1-month gap) clears the climatology gates but not
    AC1. Olympic Coast, suggested as a candidate, is the THINNEST sanctuary in
    the archive: OC01 is 1 deployment and 6 months. Only MB01 and SB02 clear
    climatology at all; 26 of 28 sites fail it. A CSD-capable acoustic record
    requires a continuous archive (MBARI MARS, Ocean Networks Canada), not
    another SanctSound site.

  RESULT AT FK01. tau=+0.200, S=+6, z=+0.784, p=0.4330 across 7 usable seasons
    and 30 comparable within-season pairs; pooled seasonal Sen slope
    +0.100 dB/yr. NO DETECTABLE TREND.
    DISCLOSED LIMIT: the index spans only 1.40 dB with ~0.3 dB between-deployment
    scatter, and this record contains no reef-degradation event, so its
    sensitivity to the thing it is meant to detect is UNMEASURED. A null result
    here is NOT evidence that the reef is healthy, and the panel says so.

  NOT CO-LOCATED WITH S19, AND THE PANEL SAYS SO. An earlier plan in this session
    claimed FK gave a triple overlap with the Looe Key SeapHOx baseline. Measured,
    that is false: FK01 is 55.0 km from Looe Key (S19 radius 15 km), and the
    nearest FK hydrophone to it is FK02 at 27.6 km - still outside. The panels
    never co-fire. There is no co-located acoustic + carbonate-chemistry site in
    this archive, and S21 must not be read as ground-truthing S19.

  ALSO FIXED, found by its own test: seasonalMannKendall() returned
    "insufficient seasons" for a fully-tied series, where the real cause is zero
    pooled variance - every season identical across years, ranks carrying no
    information. Coverage and ties are different problems and the first message
    would have sent a reader hunting for more data that would not have helped.
    They are now separate branches with separate messages.

  NO SanctSound COVERAGE EXISTS FOR BOCAS DEL TORO, PANAMA, or anywhere in the
    wider Caribbean outside the Florida Keys. SanctSound is US National Marine
    Sanctuaries only. The S21 no-match branch states this explicitly rather than
    leaving a blank panel, because the site this tool is most often driven at is
    precisely the one with no coverage. Closing that gap needs a hydrophone
    deployment or a different archive - not a wider search radius.

v10.167 FIX 21: the FAI multi-node fetch, diagnosed from the live evidence
  instead of from the span hypothesis - plus S7G, a new cross-index scatter
  panel that spends no Earth Engine calls at all.

  CHANGELOG HONESTY, the standing rule: every numeric claim below either
  (a) is quoted AS OBSERVED from live browser runs at Bocas del Toro, Panama
  (9.175, -81.986) and labelled so, or (b) was re-derived THIS session in a
  Node harness that EXECUTES THIS WHOLE FILE against stubbed ee/ui/Map and
  drives the real button handlers with synthetic Earth Engine payloads, and
  states the design that produced it. Nothing is asserted from memory and no
  precision is invented around the live figures. Earth Engine itself was NOT
  run - see RESIDUAL RISK at the end of this entry.

  S1 BLOCKER - THE FAI FETCH. THE CAUSE IS NOT THE SPAN, AND NOT THE NODE
  COUNT, AND NOT A LABELLING BUG.
    OBSERVED, across real browser runs at one site:
      module     shape                              result
      S7C        3 nodes, 24 months (2023-01 +24)   WORKS, 8 valid months/node
      S7D        9 nodes, 25-month union            returns, 0 raw features
      S7D        9 nodes, 32-month union            returns, 0 raw features
      S7E/S7F    2 sites, 36-month union            NEVER RETURNS (hangs)
      (v10.162-era) 9 nodes x 12-month AFTER window returned all 108 features
    In every S7D case the FAI pre-check reported that scenes DO exist at this
    exact point - "BEFORE 9 of 24 months | AFTER 4 of 12 months" - and S7C
    proves FAI genuinely resolves at this site.

    THREE HYPOTHESES RULED OUT, each with what rules it out:
      NOT SPAN LENGTH. S7C works at 24 months and S7D fails at 25. Those two
      numbers were never comparable: S7C and S7D do not share a fetch
      function, so that pair compares two functions, not two spans.
      NOT NODE COUNT. S7E fails with TWO sites. A "9 nodes is too many"
      explanation cannot cover a 2-site failure.
      NOT A LABELLING / GROUPING BUG. groupSeriesByLabel() does silently drop
      a feature carrying no `label` - but S7D's own diagnostic counts RAW
      features BEFORE grouping (combinedFeatCount), and the live run reported
      0 there. The features never arrived; there was nothing to drop. Traced
      (not run - this path needs Earth Engine): ptsFC sets `label` on every
      feature, reduceRegions carries non-reducer properties through
      untouched, f.set('t',t) preserves them, flatten() preserves them, and
      groupSeriesByLabel() buckets on exactly that key.

    WHAT DOES SEPARATE EVERY FAILURE FROM EVERY SUCCESS, 5 observations for
    5: THE CODE PATH. Every failing FAI fetch goes through
    extractMultiNodeSeries(); every succeeding one goes through
    computeRealCSDDeseasonalized() / computeZonalSyncCSD(). Line by line
    that is the only structural difference - the two paths agree on the
    collection builder (mkMoFAIRange), the band ('fai'), the scale (20 m),
    the buffer (150 m) and the site:

      WORKING (S7C)                     FAILING (S7D / S7E / S7F)
      monthlyColl.map(img => ...)       toList() + ee.List.sequence().map()
      img.reduceRegion(ONE geometry)    img.select([b]).reduceRegions(ptsFC)
      maxPixels: 1e9 set                no maxPixels, no tileScale
      output = nMonths features         output = nMonths x nRegions features
      one region per graph              2 or 9 regions per graph
      no label (positional {t,v})       label carried through flatten()

    So the failing path asks ONE .evaluate() to build nMonths Sentinel-2
    monthly median composites AND materialize nMonths x nRegions features out
    of them, unchunked and unbounded. Ordering the observations by that
    per-call work:
      12 mo x 9 = 108 elements    -> returned, fully populated
      24 mo x 9 = 216             -> returned EMPTY
      25 mo x 9 = 225             -> returned EMPTY
      32 mo x 9 = 288             -> returned EMPTY
      36 mo x 9 = 324 (+ 36 x 2)  -> never returned
    The two failure MODES order by that number. The failure itself is the
    path.

    FIXED by chunking, not by a span rule: each series is fetched in chunks
    whose per-call work never exceeds the LARGEST CONFIGURATION EVER OBSERVED
    TO RETURN DATA ON THIS PATH - 12 months x 9 nodes = 108 features. Both
    bounds come from that one observation and nothing else; neither is a
    guess at where the real ceiling sits, because the real ceiling has not
    been measured. Chunks are fired SEQUENTIALLY (this file has hit the
    account concurrency quota before) and merged client-side, de-duplicated
    on (label, timestamp). tileScale:4 was added to the reduceRegions call -
    the standard Earth Engine remedy for the memory pressure the empty
    returns point at, free when it is not needed, and UNTESTED here.
    S7C IS NOT CHANGED. It is the one configuration observed to work.

    WHICH FIX ADDRESSES WHICH FAILURE:
      chunking + tileScale    -> the "returns 0 raw features" failure (S7D)
      chunking + the deadline -> the "never returns" failure (S7E / S7F).
      Chunking alone cannot fix a hang; it only makes each call small enough
      that hanging is less likely and names WHICH chunk stalled. The deadline
      (S2 below) is what turns a stall into a visible refusal.

    EE CALL COUNTS, re-derived this session by counting .evaluate() calls
    while driving each module's real handler to completion with fully
    populated synthetic payloads, on the 24+12 windows of the live runs
    (a 31-month joint span, so 3 chunks per series):
      S7D  3 -> 5     S7E  3 -> 5     S7F  5 -> 9
    General form, where span is the number of calendar months the two windows
    jointly cover:
      S7D = 1 pre-check + ceil(span/12) + 1 NDVI
      S7E = 1 pre-check + 1 GEBCO search + ceil(span/12)
      S7F = 1 pre-check + 1 GEBCO search + 1 NDVI + 2 x ceil(span/12)
    which is 4 / 4 / 7 whenever the span is 24 months or less and 10 / 10 /
    19 at the 95-month widest span the existing span gate accepts. Every
    on-screen string quoting a count was updated - button labels, status
    lines, progress counters, method lines, scope notes. The count went UP,
    deliberately: it buys a fetch that can complete.

  S2 BLOCKER - S7E AND S7F COULD BE LEFT WAITING FOREVER. NOW THEY REFUSE.
    v10.161 S6a stated plainly that "a callback that never fires still never
    fires. Without a timer this cannot be turned into an automatic failure."
    That was true of the GLOBAL setTimeout, which the GEE sandbox does not
    expose. The Code Editor's own ui.util namespace DOES expose setTimeout
    and clearTimeout, so every chunked call now carries a 4-minute deadline
    that names the chunk and refuses on screen. It is FEATURE-DETECTED, never
    assumed: where ui.util is absent the panels behave exactly as they did
    before and the status box SAYS SO rather than promising a deadline it
    cannot enforce.
    All three modules now refuse identically on: an errored chunk, a chunk
    that returns 0 features, and a chunk that blows the deadline. A 0-feature
    chunk is treated as a FAILED CALL, not a data gap, because reduceRegions
    emits one feature per region per image whether or not the pixels are
    masked - so a chunk that ran returns exactly chunkMonths x nRegions
    features. S7F keeps its S7D half when only the S7E half fails, and says
    the reference reef WAS found (it used to print "no usable reference site
    found" there, which would now be false).
    The v10.161 stall-explaining progress text and the per-panel run sequence
    numbers are unchanged, and a late reply from an already-refused run
    cannot overwrite the refusal.

    VERIFIED IN NODE, by executing this file under stubbed ee/ui/Map and
    driving the real button handlers:
      - 0-FEATURE CASE: first chunk returns {features: []}. S7D, S7E and S7F
        each refuse, naming the chunk, its month range and the 108 (or 24)
        features expected.
      - PARTIALLY-POPULATED CASE: chunk 3 returns 3 months x 2 of the 9
        labels instead of 7 x 9. The run completes on the merged series and
        the Console records the short chunk.
      - NEVER-RETURNS CASE: chunk 2 is issued and never answered. On SHIPPED
        v10.166 the same drive leaves S7D reading "Running: 1 / 2 batched
        calls done... Still waiting on: combined BEFORE..AFTER series" with
        no code path able to clear it. On v10.167 the deadline fires and the
        panel refuses with the chunk named; a late reply afterwards does not
        overwrite the refusal.
      - NO-TIMER CASE: with ui.util removed entirely, the file still loads
        and runs, no deadline is armed, and the status box reads "NO TIMER
        AVAILABLE in this sandbox build".
      - The chunk planner is unit-tested directly against the shipped
        functions: month grids for all three live configurations, mid-month
        start dates, every chunk within the 12-month and 108-feature bounds,
        months conserved across the plan, and merge de-duplication.

  S3 NEW - S7G CROSS-INDEX SCATTER PANEL. ZERO NEW EARTH ENGINE CALLS.
    Every algae number in this tool was tabulated and nothing was plotted
    against anything, so spatial drift across the ring was invisible. S7G
    sits below S7F, behind its own button, and reads two client-side stores
    filled by callbacks that have already paid for their data: S7B's 9-point
    scan and the current click.
    CHARTS DRAWN:
      2. NDCI vs FAI across the 9 scan points
      3. FAI vs NDVI-water (the off-reef check against the algae index)
    CHART NOT DRAWN, and why: 1. SST x Chl-a across the 9 scan points. S7B's
    single call reduces ONE image - faiImg + ndciImg + ndviWater - over its
    nine buffers. It fetches no SST and no Chl-a, so there is no 9-point SST
    or Chl-a anywhere in this tool, and the zero-new-calls constraint forbids
    creating one. The only pair in hand is the SINGLE clicked point, which is
    below the 3-point floor. That slot prints the clicked point's SST, Chl-a,
    FAI, NDCI and NDVI-water as numbers and states the reason.
    COLOUR, and which option was chosen where: NOT nine categorical hues for
    nine compass nodes. The NDCI vs FAI chart colours by the ON-REEF /
    OFF-REEF binary S7D already computes from NDVI-water (> -0.10) - two
    named series plus an explicit "unknown" bucket. The FAI vs NDVI-water
    chart uses the other permitted option, a single series with the centre
    point marked (15px diamond), because colouring THAT chart by an
    NDVI-water cut would encode only which side of x = -0.10 each point sits
    on, which its own x axis already says.
    NO DUAL-AXIS CHART ANYWHERE. One measure per axis, units on both
    (dimensionless index here; deg C and mg/m3 in the numeric slot). Grid and
    axes recessive (#e8e8e8 gridlines, 8px grey ticks, minor gridlines off);
    marks prominent (9px points). Every series named in the legend; a
    single-series chart gets no legend box because the title names it. Fewer
    than 3 plottable points after nulls are dropped renders an explicit
    message, never an empty chart frame. Each chart is followed by its own
    numeric table, because GEE chart rendering is limited and the sidebar is
    256px wide.
    ui.Chart(dataTable, 'ScatterChart', options) is used rather than
    ui.Chart.array.values() or ui.Chart.feature.byFeature(): both of those
    take ee.Array / ee.FeatureCollection arguments and evaluate them
    server-side, which is exactly the Earth Engine call this panel is
    forbidden to spend.
    VERIFIED IN NODE against the shipped functions: pairing and null/NaN
    dropping on both axes, the <3-point guard at 0 / 2 / 3 / 6 points and on
    null input, canonical series order (stable under input shuffling, so a
    series always gets the same colour), exactly one y-value per data row,
    legend suppressed on a single series, the marked centre point, and the
    text table. Driving S7B and then S7G end to end spends 0 .evaluate()
    calls; pressing S7G before any S7B scan spends 0 and says so.

  RESIDUAL RISK, stated plainly.
    THE S1 FIX CANNOT BE VERIFIED FROM HERE. The harness proves the chunk
    plan, the merge, the refusals, the deadline and the call counts. It
    cannot prove that a 12-month x 9-node reduceRegions returns data at Bocas
    del Toro, because it never talks to Earth Engine. The evidence that it
    should is that this exact shape was OBSERVED to return all 108 features
    in a live run. If the real ceiling is lower, the remedy is a smaller
    S7_FAI_MAX_ELEMENTS_PER_CALL or S7_FAI_MAX_MONTHS_PER_CALL - one constant
    each - and the failure is now a visible refusal naming the chunk instead
    of a hang.
    The server-side mechanism behind "216 elements -> empty result with no
    error" is INFERRED, not observed. Which code path fails and which does
    not is NOT inferred: that is the observation, 5 for 5.
    tileScale:4 is untested against live Earth Engine.
    ui.util.setTimeout is feature-detected but has not been exercised in a
    real Code Editor session from here; if it is missing the panels fall back
    to v10.166 behaviour and say so on screen.
    Sequential chunking makes a long span take longer in wall-clock time than
    v10.166's single call would have if that call had worked. It did not.

v10.162 FIX 20: six defects, all diagnosed from the SAME live Earth Engine
  browser run at Bocas del Toro, Panama (9.175, -81.981) that produced the
  v10.161 round. The on-screen output is quoted verbatim where it is the
  evidence.

  CHANGELOG HONESTY, the standing rule: every numeric claim below either
  (a) is quoted AS OBSERVED from that browser session and labelled so, or
  (b) was re-derived THIS session in a Node harness that EXECUTES THIS WHOLE
  FILE against stubbed ee/ui/Map and calls the shipped functions, and states
  the design that produced it. Nothing is asserted from memory and no
  precision is invented around the live figures. Earth Engine itself was NOT
  run - see RESIDUAL RISK at the end of this entry.

  S1 BLOCKER - S7E AND S7F NEVER COMPLETED, BECAUSE THE BEFORE FETCH FAILED.
    OBSERVED, S7D's own diagnostic:
      Raw features returned: BEFORE=0 | AFTER=108 | NDVI=9
         (expect BEFORE=216, AFTER=108, NDVI=9 if fully populated)
      NOTE: 1 of 3 batched calls returned no usable data.
    AFTER (12 months x 9 nodes) returned. BEFORE (24 x 9) returned nothing.
    S7D rendered on half its data; S7E and S7F sat on "Still waiting on:
    BEFORE FAI series" / "Still waiting on: S7D BEFORE series" forever,
    because both AWAIT that call and this sandbox has no timer.
    FIXED STRUCTURALLY: the separate BEFORE fetch is GONE. Each site now
    fetches ONE series spanning BEFORE-start .. AFTER-end (including the gap)
    and slices the two windows out of it CLIENT-SIDE BY REAL TIMESTAMP - the
    same fetch-once-and-slice pattern S13 STEP 4 has used since v10.123, with
    the one difference that STEP 4's windows share a start date and can slice
    by array position while these cannot.
    EE CALL COUNTS: S7D 3 -> 2, S7E 3 -> 2, S7F 6 -> 4. Every on-screen
    string quoting the old counts was updated (button labels, status lines,
    progress counters, method lines, scope notes).
    REFUSES INSTEAD OF HANGING: if the combined series does not cover both
    windows, each module refuses with the existing insufficient-data wording.
    If the two windows are further apart than beforeMonths+afterMonths+24, it
    refuses BEFORE spending any call, because past that the combined fetch
    costs more than the call it saves.
    VERIFIED IN NODE, by executing this file under stubbed ee/ui/Map and
    driving the real button handlers with synthetic Earth Engine payloads:
      - S7D happy path consumes 3 .evaluate() calls (1 FAI pre-check + 2
        module calls); S7E 3 (pre-check + GEBCO + combined); S7F 5 (pre-check
        + 3 first-stage + 1 second-stage). Under v10.161 the same drive-through
        consumed 4 / 4 / 7.
      - Feeding the SHIPPED v10.161 an empty BEFORE and a full AFTER
        reproduces the live line verbatim: "Raw features returned: BEFORE=0 |
        AFTER=108 | NDVI=9 (expect BEFORE=216, AFTER=108, NDVI=9 ...)" and it
        renders the table anyway. v10.162 on the same input refuses.
      - Combined-returns-nothing, combined-covers-AFTER-only, and
        windows-too-far-apart all end in a refusal, never in "Still waiting".
      - The timestamp slice was unit-tested: from a 42-month combined grid it
        recovers exactly 24 BEFORE months (2021-06..2023-05) and 12 AFTER
        months (2023-11..2024-10) with no overlap; it still does so when the
        input is reversed and a third of it deleted, where an array-position
        slice gives the wrong months; and it handles a day-of-month mismatch
        between the two window starts, and the two windows typed the wrong
        way round.
    NOT FIXED, stated plainly: the combined graph is HEAVIER than either
    window alone. If the live BEFORE failure was graph weight, this fetch can
    fail too - but it then fails as ONE visible refusal instead of an
    invisible half-run that hangs two other panels.

  S2 BLOCKER - S17 AND S17b DISAGREED ABOUT n, ON ONE CLICK, BOTH WAYS.
    OBSERVED, one click, salinity:
      S17 :  not yet  SNR=1.96  n=32 of 32 nominal yr (df=30)  RISING  [MEASURED r=+0.564]
      S17b:  tau=-0.333, p=0.3813 -> FALLING  not significant  [n=7 annual points]
    and at an earlier click a few hundred metres away, S17 n=12 of 32 against
    S17b n=7. Over-counting is ANTI-CONSERVATIVE: n drives Sxx=n(n^2-1)/12,
    SStot=n*sd^2 and df=n-2, which is the defect TOE-01 exists to prevent.
    CAUSE DETERMINED, not assumed. Ruled OUT by reading the code: (a) count
    counting IMAGES - every annual image is filter(year).select(band).mean()
    with 't' added as a separate band, so a masked year is masked in the value
    band and ImageCollection.reduce(count) counts unmasked values per pixel;
    (c) different collections or date ranges - S17b is handed the very same
    _ann*Coll objects S17 reduces. The REAL cause is the SAMPLING FOOTPRINT,
    different in both dimensions: S17 reduced a BARE POINT at toeScale=27750 m
    for EVERY variable, while S17b sampled a 4 km BUFFER at 4000 m (25000 m
    for pH). HYCOM salinity is ~9 km and CMEMS ocean colour is 4 km, so at
    27750 m Earth Engine serves a pyramid overview and reducing the COUNT
    image with mean() there returns the ~28 km block's MEAN valid-year count -
    fractional, then rounded by Math.floor(x+0.5) into something that looks
    like an honest integer. At a coastal pixel more masked than the water
    around it, that average is biased UP. Hypothesis (b) - the MK path
    dropping nulls more strictly - is real but secondary: the two sets
    coincide once both read the same pixel.
    FIXED: ONE shared footprint, the CONSERVATIVE one of the two already in
    use - bare point (S17's geometry, and S18's since v10.147) at the FINER
    scale (S17b's), per variable: TOE_SAMPLE_SCALE = 4000 m for
    SST/Chl/Salinity/NO2 and 25000 m for pH/DO. Both panels now reduce the
    identical pixel of the identical collection, so the count reducer and the
    Mann-Kendall value list see the same years BY CONSTRUCTION. All four of a
    variable's reductions (linearFit, stdDev, count, correlation) moved
    together, so calcToE()'s standing assumption - that they skip the same
    masked pixels - still holds. ZERO extra Earth Engine calls: the count
    already travels inside the dictionary S17 evaluates, and S17b still fires
    5 fetches.
    ALSO FIXED, found while tracing this: calcToE() treated "the count key is
    present and EXPLICITLY NULL" the same as "no count key exists", and both
    fell into the nominal-fallback - the most flattering branch there is.
    Earth Engine returns null, not 0, for a reduceRegion over a region with no
    valid pixel, so a fully-masked pixel was handed n=nominal and df=n-2.
    VERIFIED IN NODE by driving the real map-click handler: on a dictionary
    whose count key is explicitly null, SHIPPED v10.161 rendered
    "SAL: EMERGED SNR=3.20 ... RISING n=32 of 32 nominal yr" - an EMERGED
    verdict for a pixel with no data at all. v10.162 renders "n=0 of 32
    nominal yr" and cannot emerge. An UNRECOGNISED key still falls back to
    nominal with its existing upper-bound warning, which is what that branch
    was written for.
    RECONCILIATION ON SCREEN: S17 publishes its per-variable result for the
    clicked point and S17b compares against it per variable. Where the two n's
    still differ - which needs an extra EE call to resolve, so it is not
    resolved - S17b names both, uses the SMALLER (an over-counted n inflates
    Sxx and df) and says why. VERIFIED IN NODE on both live cases: the n=32/7
    case renders "[DISAGREES WITH S17: S17 n=32, here n=7 -> the SMALLER (7)
    is the one to trust]", the n=12/7 case the same with 12, and when the two
    agree it renders "[S17 agrees: n=7]".
    THE HARDCODED-SOUNDING NOTE IS GONE: S17b's "only N years had real, valid
    HYCOM data at this exact point" was salinity-specific and could contradict
    S17 on the same click. It is now variable-agnostic, reports what BOTH
    panels measured, and names the disagreement when there is one.
    DIRECTION: with the footprints unified, S17 RISING vs S17b FALLING is no
    longer two datasets - it is an ordinary-least-squares slope and a
    rank-based trend disagreeing on the SAME sample. S17b now says so and
    calls the direction UNESTABLISHED for that variable, in the detail text
    and in its headline.
    WHAT THIS COSTS: S17's slope/noise/correlation now describe THIS PIXEL
    rather than a ~28 km block, so SNR and the measured r WILL move at any
    site where the two differ - which is every site where this bug was
    visible. That is the intent, and it cannot be verified without a live
    Earth Engine session.

  S3 THE BONFERRONI CLAIM RESTED ON PERMUTATION NOISE.
    OBSERVED: a FIND SWEET SPOT run reported "Study AC1 p=0.023" for the
    48-month row against a corrected bar of 0.0250 and concluded "At least one
    window survives the STRICTER Bonferroni-corrected bar - this is real
    evidence". At 300 shuffles p=0.023 IS 7/300, and the grid step is 0.0033.
    OBSERVED independently, in an earlier run: three DUPLICATE rows -
    identical data, identical statistics (-0.70x, +0.214, -1.03x, +0.34x on
    all three) - returned control-variance p-values of 0.003 / 0.010 / 0.030.
    FIXED: 10000 shuffles on any row that can be COUNTED; 300 kept on the
    sub-floor and DUPLICATE rows, which are excluded from every tally, the
    best-window pick and the Bonferroni divisor, so no claim rests on them.
    WHY 10000, RE-DERIVED THIS SESSION against the shipped
    permutationTestDeltaFixed (BEFORE 36 valid months, AFTER 48 - the heaviest
    row this panel builds):
        B       SE at p=0.025   SE at p=0.00833   wall clock per test
       300        0.00901         0.00525            20 ms
      2000        0.00349         0.00203            49 ms
     10000        0.00156         0.00091           187 ms
     20000        0.00110         0.00064           337 ms
    The corrected bar is 0.05/k for k powered windows: 0.0500 at k=1 down to
    0.00833 at k=6. 10000 puts the Monte Carlo SE at the k=2 bar (0.0250 - the
    live case) at 0.0016, about 1/16 of the bar and 5.8x tighter than 300.
    COST: at most 6 windows x 4 tests, and only countable rows pay, so the
    realistic addition is 8-16 tests x 187 ms = 1.5-3.0 s, 4.5 s worst case -
    against this panel's 8 Earth Engine calls (20-60 s). No EE work is added
    at all; the shuffles are client-side JS on series already fetched.
    SEED-TO-SEED SPREAD, measured on ONE fixed dataset over 12 independent
    repeats of the same test: at B=300 the p-values spanned 0.6033-0.7000
    (spread 0.0967); at B=10000, 0.6493-0.6653 (spread 0.0160).
    THE PANEL NOW STATES THE SHUFFLE COUNT PER ROW, in its own table column.
    AND IT NO LONGER ASSERTS WHAT IT CANNOT: a "survives the corrected bar"
    claim is checked against its own Monte Carlo SE, and if p is within 2 SE
    of the bar the panel prints TOO CLOSE TO CALL with the numbers instead.
    RE-DERIVED: at the live 48-month figures (p=0.023, bar 0.0250, B=10000)
    SE=0.0015 and p+2SE=0.0260, which is ABOVE the bar - so that exact row
    would now be reported as too close to call, not as a survivor. It would
    take B=100000 for p=0.023 to clear 0.0250 by 2 SE. Raising B further is
    possible but the honest reading is that the window sits AT the bar.

  S4 THE STEP 3 HEADLINE CONTRADICTED ITSELF.
    OBSERVED, three lines in this order on one panel:
      STATISTICALLY SIGNIFICANT LOCAL CSD SIGNAL (AC1-confirmed, p=0.028)
      STRONG SIGNAL, LIKELY REGIONAL (AC1-confirmed, but matches control site too)
      Regional context: VERDICTS DISAGREE ... Treat NEITHER as confirmed
    v10.151 FIX 6 made the disagreement VISIBLE; it did not stop the most
    prominent banner asserting one side of it.
    ROOT CAUSE: the permutation test behind that banner runs on the STUDY SITE
    ONLY. It never looks at the control, so LOCAL - which means "different
    from the control" - was never something it could establish.
    FIXED: the banner states what IS established (the permutation p-value on
    the study site's dAC1) first, and what is NOT (local vs regional) second,
    in that order; when the two classifiers disagree it says so instead of
    picking a side. The banner is now a function and both async paths redraw
    it, so whichever lands second updates the other.
    MAGNITUDE: the numbers behind "matches control site too" were study
    dAC1=+0.512 against control dAC1=+0.134 - about 3.8x. "Matches" is a
    magnitude claim and nothing in the code checked a magnitude. A pure,
    unit-tested classifyRegionalAC1() now compares them and the headline
    carries the two real numbers. RE-DERIVED: on the live pair it returns
    "amplified", ratio 3.82x; on the Nuuk pair the old wording was written for
    (+0.234 / +0.242) it returns "comparable", 0.97x.
    THE DECISION RULE IS DELIBERATELY NOT CHANGED, and the reason is measured
    rather than asserted: on PAIRED NULL data (both sites the same AR(1)
    process, each site's dAC1 the difference of two independent lag-1 AC1
    estimates from the shipped jsLag1AC1, 20000 draws per cell) a ratio of
    2x or more occurs 14.5% of the time at n=24 phi=0.2, 16.0% at n=36
    phi=0.2, 13.9% at n=48 phi=0.5 and 13.4% at n=60 phi=0.5 - about 1 run in
    7, on pure noise. A magnitude gap is not a significance test, so it
    changes the WORDING and gates no verdict, and that 13.4-16.0% figure is
    printed on screen next to it.

  S5 A SIGNIFICANT CONTROL-SITE CHANGE WAS GOING UNREPORTED.
    OBSERVED, across several windows:
      36mo      | Study Var p=0.647 | Ctrl Var p=0.000
      48mo(40)  | Study Var p=0.553 | Ctrl Var p=0.000
      60mo(40)  | Study Var p=0.603 | Ctrl Var p=0.007
    "Local signal? no" was CORRECT - local means the study site DIVERGING from
    the control - but the panel then ended "No window shows a real local
    signal" and never said the CONTROL was changing significantly. That is a
    regional finding this tool is positioned to make: the control is the
    baseline the study site is measured against, and a baseline that is itself
    moving is a result, not a null.
    FIXED: significant control-site hits on either statistic, in either
    direction, are collected on the countable rows and reported in a REGIONAL
    FINDING block with the statistic, the direction, the p-value, the shuffle
    count, the delta, and whether the study site is significant on that
    statistic too. It states what it is NOT (a local warning) and that it
    weakens any local claim made against a moving baseline, and repeats the
    10 km / 80 km buffer caveat. Hits that appear only on UNDERPOWERED or
    DUPLICATE rows are excluded, exactly as the study-site tallies exclude
    them, and said so.
    VERIFIED IN NODE by driving the real FIND SWEET SPOT handler with a
    synthetic control-variance surge: the table reproduces the reported shape
    (Study Var not significant at 36/48/54/60mo while Ctrl Var p<=0.001) and
    the new block lists all five control hits. The SHIPPED v10.161, on the
    identical input, ends at "No window shows a real local signal, corrected
    or uncorrected." and says nothing about the control.

  S6 THE VARIANCE-ARTIFACT RULE LEAKED.
    OBSERVED in S7C:
      Center: AC1=-0.146 Var=4.14x  (1st-half=0.0000, 2nd-half=0.0000)  [NOT flagged]
      North:  AC1=0.009  Var=12.08x (1st-half=0.0000, 2nd-half=0.0000)  LIKELY ARTIFACT
    Same 0.0000 first half, one flagged and one not, because the rule was
    varFirst < 0.001 AND |ratio| > 5. The AND made the magnitude threshold
    decide something it cannot know: varTrendRatio = secondHalfVar /
    max(firstHalfVar, 1e-6), so once the first half is at the 1e-6 floor the
    denominator is a CONSTANT and the ratio's size says nothing about whether
    the ratio is trustworthy.
    FIXED: a near-zero denominator disqualifies on its own. The magnitude
    argument is still accepted (every call site passes it) but no longer
    gates the flag. S7C and S7D each had their OWN private copy of the rule
    with the same leak; both now call the single shared isVarRatioArtifact().
    RE-RAN THE FALSE-POSITIVE SWEEP against the shipped jsNodeStatsFixed:
    genuine AR(1) surges, BEFORE phi=0.2 vs AFTER phi=0.9, per-reading noise
    SD 1.0, 500 draws per cell: 0/500 flagged at 24 months and 0/500 at 36,
    BEFORE and AFTER the change (min first-half variance 0.202 and 0.234, two
    orders of magnitude above the 0.001 bar). A flat phi=0.2->0.2 control at
    36 months: 0/500 both. Synthetic rows at the live numbers: Center
    (varFirst 4e-5, 4.14x) goes unflagged -> FLAGGED; North (2e-5, 12.08x)
    stays flagged. Genuine surges on a real first half (varFirst 1.0/0.4/
    0.0012 at 12x/40x/50x) stay unflagged in both versions.
    ONE HONEST CAVEAT, found by the same sweep: the 0.001 bar is ABSOLUTE and
    therefore scale-dependent. On a series whose real variance is near it -
    noise SD 0.05, stationary variance ~2.6e-3 - flags rise from 4/500 to
    12/500. That is the pre-existing threshold, not the AND removal, and
    CSD_VAR_ARTIFACT_VARFIRST is unchanged; it is recorded here rather than
    silently absorbed.

  RESIDUAL RISK, v10.162.
    EARTH ENGINE WAS NOT RUN. Everything above was verified by executing this
    entire file in Node against stubbed ee/ui/Map and driving the real button
    and click handlers with synthetic payloads, which exercises every
    client-side path - the callback graph, the slicing, the refusals, the
    rendered strings - but NOT what Earth Engine actually returns.
    SPECIFICALLY UNVERIFIABLE WITHOUT A LIVE SESSION:
      - S1: whether the combined fetch SUCCEEDS. The slicing, the call counts
        and the refusal paths are proven in Node; whether a 42-month combined
        Sentinel-2 graph returns where a 24-month one did not is a question
        only Earth Engine can answer, and the combined graph is the heavier
        of the two. The hang is removed either way.
      - S2: whether reducing at 4000/25000 m instead of 27750 m returns what
        it should, and whether the two panels' n's then actually agree in the
        field. The reasoning is that they must, because they reduce the same
        pixel of the same collection - but pyramid behaviour at a point is not
        something this harness can reproduce. S17's SNR and measured r will
        move at coastal sites, by an amount no one can state until it is run.
      - S2 again: the exact key name the count and correlation reducers use is
        still read defensively client-side (toeNum), unchanged since v10.157,
        because it cannot be confirmed without a session.
      - S5: the regional block's wording is exercised on synthetic data; the
        real control-site p-values will differ.
    NOT FIXED, DELIBERATELY: the two STEP 3 classifiers are still two
    classifiers with different inputs (a study-only permutation test and a
    variance-based study-vs-control rule). v10.162 stops the banner asserting
    one over the other and makes the magnitudes visible; it does not merge
    them into a single local-vs-regional significance test, which would be a
    new statistical design rather than a fix, and one that could not be
    calibrated without a live session.

v10.161 FIX 19: seven defects found in a REAL Earth Engine browser run at
  Bocas del Toro, Panama (9.175, -81.981). Unlike previous rounds these are
  OBSERVED BEHAVIOURS from live satellite data, not simulated ones. The
  on-screen output is quoted verbatim where it is the evidence.

  CHANGELOG HONESTY, the standing rule in this file: every numeric claim
  below either (a) is quoted AS OBSERVED from that browser session and is
  labelled so, or (b) was re-derived THIS session in a standalone Node
  harness against the pure-JS functions extracted from this file, and states
  the design that produced it. Nothing here is asserted from memory, and no
  precision is invented around the live figures. Earth Engine itself was NOT
  run this session - see RESIDUAL RISK at the end of this entry.

  S1 VERSION MARKERS - THE MECHANISM IS FIXED, NOT THE LITERALS.
    OBSERVED: the sidebar footer read "Scroll up for measurements | v10.159
    + GEM" while the sidebar title read "STEMGeoHS Marine v10.160", in the
    same panel, on the same screen. This is the THIRD CONSECUTIVE ROUND a
    version marker has been missed (found at v10.149 in round 3; round 5
    claimed "all five reconciled"; wrong again one version later).
    Patching literals has now failed three times, so the literals are gone.
    var TOOL_VERSION near the top of MODULE A is the ONLY place
    the running version is written down. All SIX self-identifying markers -
    sidebar title, S13 section header, sidebar footer, per-click console
    banner, startup READY line, and the file header comment (which now
    carries no version at all, because a comment cannot concatenate) - are
    derived from it. VERIFIED: grep for a hardcoded self-identifying version
    literal returns exactly ONE hit, the constant.
    DELIBERATELY NOT TOUCHED, and why: well over two hundred string literals
    in this file contain a v10.1xx token (228 at the moment this entry was
    written, counted by script over non-comment lines - the figure moves with
    every changelog line added, so treat it as an order of magnitude, not a
    constant). They are (a) changelog text, (b) provenance -
    "(v10.160, Monte Carlo)" says WHEN a figure was measured, "S7D - FULL
    9-NODE ALGAE COUPLING NETWORK (v10.106)" says when a module was
    introduced. Both are WRONG if they move with the running version. Only
    markers that assert "the tool you are running is version X" were
    converted. The stale comment "// S13 - CSD EARLY WARNING TEST (v10.149)"
    above the S13 header had its version dropped rather than corrected.

  S2 THE v10.158 CHLOROPHYLL INVERSION NEVER REACHED THE MAP LEGEND.
    OBSERVED: the sidebar readout correctly said "1.024 mg/m3 (ocean)
    [enriched]" - a stressor - while the map legend on the same screen said
      < 0.1  Oligotrophic | 0.5-1.0 Moderate | 1.0-2.0 Good | > 5.0 Bloom
    i.e. 1.024 mg/m3 read as STRESS in one panel and "Good" in the other.
    v10.158 W-08 reversed the direction in computeScore()'s s2 and s5;
    v10.158/v10.159 reworded the readout; the legend was never touched.
    The legend now lists EXACTLY the six s2 bins in computeScore, including
    the 0.45 mg/m3 GBR annual-mean guideline (De'ath & Fabricius 2010) the
    readout uses, with each bin's s2 sub-score printed next to it so the two
    cannot drift again without the numbers visibly disagreeing. The old
    legend was also incomplete, not just misdirected: four bands with two
    gaps (0.1-0.5, 2.0-5.0) and a top band at >5.0 that no scoring rule has
    ever used. RE-DERIVED against the shipped computeScore this session
    (sat inputs SST 29.0, trend 0.04, NO2 1.2e-4, turb 0.15, DHW 5, Bocas
    profile): s2 = 5 / 5 / 15 / 35 / 35 / 60 / 80 / 80 / 95 / 95 at chl =
    0 / 0.05 / 0.15 / 0.3 / 0.44 / 0.46 / 1.024 / 1.5 / 2.5 / 6.0 - monotone
    non-decreasing, with the step at the 0.45 guideline.
    THE MAP PALETTE WAS ALSO WRONG, not just the labels: it ran blue ->
    GREEN -> dark green over 0.01-5.0, painting the most enriched water in
    the colour this tool uses for "safe" everywhere else (score dot: green =
    DEEP BASIN). Now a blue -> red stress ramp stretched 0-2.0, so the 0.45
    guideline and the 1.0/2.0 bin edges fall inside the ramp instead of
    being squeezed into its first fifth.
    GREPPED FOR OTHER MISSED DIRECTION TEXT: the only remaining "chlorophyll
    is good" wording is in S8 computeAquaculture() and its sidebar rows
    ("Nutrients Chl-a gate", "SST and chlorophyll both favorable"). That
    direction is CORRECT there - chlorophyll is food for a seaweed crop -
    and is deliberately left alone, as v10.158 W-08 already stated.

  S3 ECI WAS TWO DIFFERENT NUMBERS IN THE SAME PANEL.
    OBSERVED, S16, four lines apart:
      Energy risk (current):  0.159 low risk
      ... physical wave-exposure (ECI=0.89) and broad ecological risk
          (B=0.45) measure different things
    TRACED: 0.159 IS the ECI (eci_current = min(1.414, 1/sqrt(depth)); a
    ~39.5 m GEBCO depth gives 0.159). 0.89 was the local variable bEci =
    1 - ECI/1.414, a DIFFERENT quantity constructed only so the comparison
    would sit on B's 0-1 "higher = safer" scale. It was never the ECI. It is
    now named shelterIdx and labelled "shelter index = 1 - ECI/1.414", and
    the comparison line also prints the ECI itself, explicitly tagged as the
    same value the Energy risk row above displays. Two names, two numbers,
    no collision. The v10.145 point that this comparison is NOT a validation
    is unchanged and still correct.

  S4 THE FAI MODULES SPENT THEIR WHOLE CALL BUDGET BEFORE REFUSING.
    OBSERVED at Bocas del Toro, Sentinel-2 FAI:
      Study data density: 11% (0/24 BEFORE, 4/12 AFTER valid months) - INSUFFICIENT
      Reference data density: 11% (0/24 BEFORE, 4/12 AFTER valid months) - INSUFFICIENT
    ZERO valid months in the entire 24-month BEFORE window. The refusal is
    correct and well explained - but it arrived only after S7D/S7E/S7F had
    fired 3, 3 and 6 heavy reduceRegions graphs over 36 monthly composites
    and the user had waited out the run.
    NEW: faiPrecheckThenRun(), a data-density gate that runs BEFORE the
    budget is committed. ONE .evaluate() of
      S2_SR_HARMONIZED.filterBounds(point).filterDate(union of the windows)
                      .filter(CLOUDY_PIXEL_PERCENTAGE < 20)
                      .aggregate_array('system:time_start')
    - a metadata/index query that reads no pixels, runs no reducer and
    touches no imagery. The timestamps are bucketed into calendar months
    client-side and compared against the 4-valid-months-per-window floor the
    modules already refuse on.
    WHY A SCENE COUNT IS THE RIGHT PROXY: mkMoFAIRange() applies NO per-pixel
    cloud mask, so a month with no qualifying scene over the point CANNOT
    produce a valid FAI value, and a month that has one normally does. The
    scene count is an UPPER BOUND on valid months that is expected to be
    tight. It can prove insufficiency conclusively; it can never promise
    sufficiency, and the gate is written to refuse only, never to guarantee.
    S7D/S7E/S7F refuse below 4 scene-months in either window, naming the
    counts and saying "S7C-S7F cannot run here". S7C is DELIBERATELY
    advisory - its stated job is to measure whether FAI is dense enough at
    all, so refusing it for sparseness would refuse the one module that
    exists to report sparseness - and it refuses only a completely empty
    window (zero scene-months), otherwise running with the count shown up
    front and a warning below the floor. FAIL-OPEN: if the pre-check errors
    or returns nothing, the run proceeds and says the pre-check did not run.
    UNIT-TESTED this session in Node on the pure half, faiPrecheckSummary(),
    with synthetic scene-timestamp lists against windows BEFORE 2021-06-01
    +24mo and AFTER 2023-11-01 +12mo: the live-run shape (scenes only in
    4 AFTER months) gives BEFORE 0 of 24 | AFTER 4 of 12 and REFUSES; 24 and
    11 pass; exactly 4 in AFTER passes (floor is inclusive); 3 refuses; five
    scenes inside one calendar month count as ONE month; a scene at
    2023-12-31T23:30Z buckets to 2023-12, not 2024-01. Month arithmetic
    verified separately: 2021-06-01 +24 -> 2023-06-01, 2023-12-01 +1 ->
    2024-01-01, 24-month key list runs 2021-06 .. 2023-05.
    COST: NOT measured against live Earth Engine - there is no EE access
    from the Node harness. What can be said is structural: one extra
    .evaluate() of a metadata aggregation against 3-6 .evaluate()s of
    reduceRegions over 36 monthly median composites at 20 m. See RESIDUAL
    RISK.
    ASKED SEPARATELY - CAN THE REFERENCE DENSITY SILENTLY INHERIT THE
    STUDY'S? (both reported 0/24 and 4/12 at sites 20 km apart). TRACED
    THROUGH THE CODE: IT CANNOT. Study and Reference are two distinct
    labelled features of ONE FeatureCollection; extractMultiNodeSeries runs
    reduceRegions per image over both and carries each feature's own `label`
    through the flatten; groupSeriesByLabel() buckets strictly on that label,
    has no default bucket and no copy-from-sibling path, and DROPS a feature
    with no label rather than merging it. A missing Reference series would
    yield [] and 0 valid months - which would DIFFER from the study, not
    match it. They match for a real reason: with no per-pixel cloud mask,
    validity is scene availability, which is a per-granule property, and an
    S2 granule is ~110 km across - two sites 20 km apart normally share the
    same scene list month for month. Both S7E and S7F now SAY this on screen
    when the counts come back identical.

  S5 A PLACEHOLDER WAS MOVING THE HEADLINE SCORE *AND* BUYING ACCURACY.
    OBSERVED:
      F1 Urchin grazer:   0 (N=0.3/m2)
      F3 Anem density:   -6 (aN=6/m2) [ESTIMATED]
      F4 Metals:         +0 (Cd=0.006)
      F5 Recruitment:     0 (recruit=2) [ESTIMATED]
      Total correction:  -6 total
      Satellite CCS: 61/100  ->  FUSED CCS: 55/100
      Field gain: +12%  ->  Combined: 89% TOTAL
    The only non-zero correction was F3 = -6, from anem_N = 6.0, which
    getFieldProfile() itself flags anem_N_estimated:true and whose notes say
    "Densities are ESTIMATED PLACEHOLDERS". So a placeholder moved the
    headline 61 -> 55 and then earned +12 percentage points of claimed
    accuracy for having done so.
    ROUTE TAKEN, and why: the task offered "exclude estimated fields from
    accuracy_field_gain" OR "scale their contribution and label it". BOTH
    halves of the first route, and none of the second - the CORRECTION is
    left at full strength (deleting it would hide a real modelling choice),
    and it is the ACCURACY CLAIM that is scaled, because accuracy is the
    thing a placeholder cannot honestly buy. The gain is apportioned by the
    share of the correction's MAGNITUDE contributed by non-estimated fields,
    with a count-based fallback when every correction is exactly zero so a
    genuinely measured all-zero profile is not punished. Nominal and earned
    are both returned (acc_field_nominal / acc_field) so the forfeited part
    is visible rather than quietly deleted. The fused score, the total-
    correction row and the console block now all carry the split.
    ALSO FIXED HERE, a NaN swallowed by the same guards: every field guard
    tested `!== null`, which is TRUE for a MISSING KEY. The Caribbean
    profile has recruit:2.0 and NO recruit_healthy key, so
    `fp.recruit_healthy !== null` passed, 2.0/undefined gave NaN, and the
    isNaN() line turned it into a silent 0 - which is exactly why the live
    run shows "F5 Recruitment: 0 (recruit=2)" instead of taking the intended
    recruit<=5.0 fallback. Guards now use _fhas(), which rejects undefined
    and NaN.
    MEASURED IN NODE against the shipped functions, region 'Bocas del Toro,
    Panama', with the same satellite inputs throughout:
      before  F1 -0.3  F3 -6  F4 0  F5  0    fcTotal -6  acc_field 12  acc_total 89
      after   F1 -0.3  F3 -6  F4 0  F5 -3.2  fcTotal -9  acc_field  0  acc_total 77
    i.e. the correction gets BIGGER (the F5 branch now runs) and the
    accuracy claim goes to zero, because -9.2 of the -9.5 raw correction is
    placeholder and only -0.3 is a real measurement. A satellite composite
    of 61 therefore fuses to 52, not 55. STATED PLAINLY: this changes the
    headline score at Caribbean / Florida / Gulf / Atlantic-USA sites.
    REGRESSION-CHECKED on the other profiles: Great Barrier Reef keeps its
    full +7% (recruit 187/247, not flagged estimated), Mediterranean keeps
    its full +10% (urchin + metals, not flagged), Red Sea keeps its full +8%
    - its contribution is dhw_calibration, which feeds s4d INSIDE the
    satellite composite rather than any F-term, so it is counted explicitly
    as a present non-estimated contribution; without that it would have been
    silently zeroed. Pacific Coast USA and the other hasField:false regions
    stay at +0%.
    DISCLOSED LIMIT of the apportioning rule: the dhw_calibration
    contribution carries no point magnitude, so a hypothetical region with
    BOTH a real dhw_calibration AND a placeholder F-term would apportion on
    the F-terms alone. No shipped region has that combination today.
    NOT CHANGED: profiles whose notes say "NOT independently verified"
    (Mediterranean metals, GBR recruits) are UNVERIFIED, not PLACEHOLDERS,
    and carry no *_estimated flag. Only an explicit flag scales the gain.

  S6 INSUFFICIENT-DATA PATHS ABANDONED THE UI MID-RENDER (three defects).
    (a) STUCK PROGRESS LABELS. OBSERVED: long after S7F had completed with
    full results, S7D still read "Running: 2 / 3 batched calls done..." and
    S7E "Step 2/3: 1 / 2 batched calls done...", permanently.
    TRACED (not run - UI paths cannot be executed in the Node harness):
    every one of these panels updates its status inside a per-call bump()
    and only replaces it when the pending counter reaches ZERO. Both finish()
    functions DO set a final label on every path they reach, including the
    insufficient-data ones - so a frozen counter means one .evaluate()
    callback NEVER ARRIVED (quota refusal or hung request; this file has hit
    the account concurrency quota before). There is no code path that can
    clear it, and there is no setTimeout in the GEE sandbox to time it out.
    WHAT IS FIXED: the counter now NAMES the calls still outstanding; it
    carries, in the label itself, what a stalled counter means and that
    pressing RUN again is a clean retry; and a per-panel run sequence number
    stops a late callback from an abandoned run overwriting a newer run's
    label. Applied to S7C, S7D, S7E, S7F and S13 STEP 4.
    WHAT IS NOT FIXED, plainly: a callback that never fires still never
    fires. Without a timer that cannot be turned into an automatic failure.
    The S13 STEP 4 panel's claim that the counter means it "never looks
    frozen" was FALSE and is replaced with what the counter can actually
    promise - in the panel text and in the startup feature list.
    (b) AN EMPTY DENOMINATOR RENDERED AS A FRACTION. OBSERVED:
      S7D: NO SYNCHRONIZATION SIGNAL (0/7 on-reef rising) | AC1 rising at 0/0 nodes
    "0/0 nodes" is not a result: NO node had a computable AC1 in both
    windows, so the tally divided by nothing, and "tested 0 nodes" is
    indistinguishable on screen from "tested some and none rose" - opposite
    findings. New countOfTotal() renders a real fraction when the
    denominator is positive and an explicit sentence when it is zero.
    (c) A BONFERRONI CORRECTION OVER A FAMILY OF ZERO TESTS. OBSERVED:
      NOTE: 4 tests fired here; the Bonferroni bar is p<0.0125.
    printed while all four of Study AC1 / Study Var / Reference AC1 /
    Reference Var had come back NOT TESTABLE. Both numbers were HARDCODED,
    so neither could notice. The family is now the tests that actually
    returned a p-value, the bar is derived from that count, and when the
    count is zero the line says so instead of quoting a correction for tests
    that never ran. The same defect was found and fixed in S7D (its bar came
    from a fixed 17 = 8 pairs + 9 nodes regardless of how many were
    testable) and in S13 STEP 4 (with nPoweredWindows = 0 it printed
    "0 of those 6 windows are ... the CORRECTED bar below (0.05 / 0 POWERED
    windows = 0.0500)", an alpha that came from Math.max(1,0), not from
    anything measured).
    AUDIT OF EVERY OTHER INSUFFICIENT-DATA EARLY RETURN IN S7C/S7D/S7E/S7F/
    S13 - one more found: S13 STEP 5's csdSlideValidV was never cleared at
    run start and is written only far down the SUCCESS path, after the
    feats.length<4 early return. On a re-run that hit the insufficient-data
    path, the PREVIOUS run's full independent-window significance block
    stayed on screen underneath a fresh "INSUFFICIENT DATA" verdict, reading
    as current results for the new coordinates. Cleared at run start, and
    both early returns now say explicitly that the test did not run.
    The same class was introduced by S4's own pre-check gate (a refusal
    returns before the module clears its result panel) and is handled inside
    faiPrecheckThenRun rather than left for the next round to find.

  S7 SWEEP FOR OTHER x/0 DISPLAY STRINGS AND UNCLEARED PROGRESS LABELS.
    Every "a of b" / "a/b" tally that renders a denominator was checked for
    reachability of b = 0. Already guarded and left alone: S7B patchiness,
    S13 ROBUSTNESS / AC1 TALLY / FLAGGED WINDOWS, classifyToolkitConfidence,
    S7D and S7F on-reef pair counts, S17 n-of-nominal. FIXED, all reachable:
      - S13 TOOLKIT SUMMARY header, "0/0 scored indicators agree"
      - S13 COMPARE verdict, "Supporting: 0/0 agree"
      - S13 STEP 4, "0 of 0 DISTINCT windows that can support inference"
      - S7B, "WARNING SIGN COUNT: FAI elevated 0/0 | NDCI 0/0 | NDVI 0/0",
        printed directly under a PATTERN line that correctly said NO DATA
      - S17b, "0 of 0 testable variables show a REAL significant trend"
      - S17 compound, "0/0 available variables emerged"
      - S13 STEP 2 console, "Raw valid months: 0 / 0" on an empty table
      - S7D, "0 of 0 on-reef pairs tested" and "0 of 0 nodes tested"
    Progress labels: S7C/S7D/S7E/S7F/S13-STEP-4 are the only set-but-never-
    cleared ones and are covered by S6(a). Every other in-flight label
    ("Step 1/3: Finding control site...", "Fetching ... OISST data",
    "Checking for manual control site override...") was checked and has a
    path that replaces it on both its success and its error branch.

  RESIDUAL RISK - NEEDS A LIVE EARTH ENGINE SESSION.
    1. THE PRE-CHECK'S REAL EE COST IS UNMEASURED. aggregate_array over a
       point-filtered, date-filtered, property-filtered S2 collection is a
       metadata query by construction, but how long it takes and what it
       counts against the account's quota has NOT been observed. If it is
       slow at some site, narrow filterDate - do not remove the gate.
    2. THE PRE-CHECK'S TIGHTNESS IS INFERRED FROM THE CODE, NOT MEASURED.
       The claim "scene-months is an upper bound on valid months, expected
       to be tight" follows from mkMoFAIRange applying no per-pixel cloud
       mask. It has not been checked against a real run. It is used only to
       REFUSE, so a loose bound costs a wasted run, never a false refusal.
    3. ALL UI PATHS WERE TRACED, NOT EXECUTED. The S7C-S7F and S13 changes
       are inside .evaluate() callbacks and ui.Button onClick handlers that
       the Node harness cannot drive. The whole file was executed top to
       bottom under stubbed ee/ui/Map to catch ordering faults that
       `node --check` misses, and every pure function touched was unit-
       tested - but the rendered strings themselves were read, not run.
    4. THE FOUR onClick BODIES ARE NOW WRAPPED in a pre-check callback. The
       wrapping is mechanical and the file parses and executes, but a
       closure mistake inside a handler would only show at click time.
    5. THE S5 SCORE CHANGE IS REAL AND USER-VISIBLE. Caribbean-group sites
       move from a -6 field correction to -9 and from +12% claimed field
       accuracy to +0%. That is the intended, honest outcome, but any
       standing output recorded at those sites before v10.161 is not
       comparable with one recorded after it.

v10.160 FIX 18: four blockers - TWO of them regressions this fix series
  introduced, and one of them a v10.159 "fix" that turned out to be completely
  inert - plus seven should-fix defects, seven nits, and a changelog audit.

  ================================================================
  READ THIS FIRST: THE ONE SET OF NUMBERS.
  Three previous entries in this file quoted DIFFERENT figures for the same
  quantities, and several of those figures did not reproduce. Everything
  numeric below was measured THIS session in a standalone Node harness against
  the pure-JS functions extracted from THIS file, and every figure states the
  design that produced it. Where a v10.159 figure could not be reproduced it is
  WITHDRAWN by name further down rather than quietly restated. Earth Engine
  cannot be run from the harness; RESIDUAL RISK at the end says what that
  leaves unchecked.

  COMMON DESIGN for every Monte Carlo figure below unless the line says
  otherwise: monthly {t,v} series; a seasonal cycle plus white noise of
  per-reading SD 1.0; NO real change (false-positive cells) or the real change
  named in the line; random start calendar month; BEFORE->AFTER gap drawn
  uniformly from 0-3 months; 300 shuffles per permutation test; replicate count
  stated per table; nominal 5%. Monte Carlo standard error is about 0.5
  percentage points near 5% at 2000 replicates and about 0.4 at 3000, so a
  single cell moves by roughly a point between runs. No claim is made about any
  cell's third digit.

  TABLE 1 - VARIANCE-RATIO STATISTIC, THE SHIPPED WINDOW LENGTHS.
  2000 reps/cell. Ranges span seasonal amplitude 3 AND 10, three cycle shapes
  (smooth sine, sawtooth, summer spike) and two gap policies (fixed 0 months,
  uniform 0-3). Real change = AFTER window's second half at 5x the noise SD.
    windows              FPR           power against the 5x change
    18+12 (v10.159 S7D/S7E default)  5.5 - 7.0%     63.7 - 74.8%
    24+12 (v10.160 default, all 3)   4.5 - 5.9%     66.0 - 70.2%
    48+24                            4.8 - 5.4%     86.3 - 87.1%
    36+36 (STEP 3 / STEP 4)          4.3 - 4.5%     96.4 - 96.8%

  TABLE 2 - VARIANCE-RATIO STATISTIC BY POOLED TOTAL, floors lowered to 8 so
  sub-floor totals run. 3000 reps/cell, amplitude 3, smooth + sawtooth, gap 0
  and gap 0-3. This is the table CLIM_MIN_TOTAL_SAMPLES is chosen from.
    pooled    FPR            power against the 5x change
      24    1.1 - 4.1%       21.7 - 40.3%
      25    3.0 - 4.6%       34.2 - 48.8%
      26    3.4 - 5.4%       47.1 - 58.5%   <- CHOSEN
      28    5.2 - 6.2%       61.3 - 68.8%
      30    6.1 - 6.9%       70.7 - 75.3%
      36    4.9 - 5.4%       68.3 - 69.8%

  TABLE 3 - AC1 STATISTIC, same function, different statistic. 2000 reps/cell,
  amplitude 3 smooth sine, AR(1) noise, BEFORE phi=0, AFTER phi=0 (FPR) or 0.8
  (power). This is the table CSD_AC1_MIN_POOLED_MONTHS is chosen from.
    pooled  cfg      FPR     power (phi 0 -> 0.8)
      24    12+12    0.2%      0.7%
      30    18+12    2.0%      3.5%
      36    18+18    2.1%     13.1%
      48    24+24    4.3%     39.2%   <- CHOSEN floor
      72    36+36    4.8%     69.4%
      96    48+48    4.9%     89.3%

  TABLE 4 - AC1 STATISTIC, BALANCED WINDOWS, BEFORE phi=0.2 vs AFTER phi.
  2000 reps/cell, amplitude 3 smooth sine. This replaces the power table in
  CSD_POWER_TABLE_TXT, whose design was never recorded.
    windows  phi=0.2 (FPR)  phi=0.5  phi=0.7  phi=0.9
    12+12     NO p-VALUE       -        -        -
    18+18     NO p-VALUE       -        -        -
    24+24        4%           9%      17%      25%
    36+36        5%          15%      34%      52%
    48+48        5%          22%      53%      75%
  ================================================================

  BLOCKER 1 - THE v10.159 W-03 FIX WAS INERT. wActual COUNTED FEATURES, NOT
    VALID MONTHS.
    fssActualAfterMonths() returned Math.min(w, poolStudyAfterFull.length).
    mkMoSSTRange() builds nMonths images UNCONDITIONALLY - the
    ee.Algorithms.If else-branch is a fully-masked constant image, not an
    omission - so reduceRegions() emits one feature per nominal month whatever
    OISST does, and groupSeriesByLabel() pushes {t, v:null} for each masked
    one. poolStudyAfterFull.length was therefore always 60.
    REPRODUCED against an EE-shaped FeatureCollection of 60 features of which
    months 39-59 carry no band value (the exact case the v10.159 entry is
    written around), driving the real extracted code:
      v10.159: wActual [12,24,36,48,54,60], duplicateOf all null,
               nPoweredWindows 4, Bonferroni alpha 0.0125, fssShortAfter false
      v10.160: wActual [12,24,36,39,39,39], 54mo and 60mo DUPLICATE of 48mo,
               nPoweredWindows 2, alpha 0.0250, fssShortAfter true
    i.e. exactly the state v10.159 CLAIMED to have reached. Also verified at a
    36-month record (1 powered row, alpha 0.0500, three rows marked duplicate)
    and a 60-month one (4 powered rows, alpha 0.0125, no duplicates, warning
    silent) - both matching v10.159's claimed behaviour, now for real.
    fssShortAfter was false for the same reason (60 < 60), so its warning could
    never print. Both now count months that carry a value, and the two branches
    of fssActualAfterMonths - one of which counted features and the other
    nValidMonths - now measure the same quantity.

  BLOCKER 2 - THE STEP 3 HEADLINE PRINTED "NO RELIABLE CSD SIGNAL" WHEN THE
    TRUTH WAS "NOT TESTABLE", AND v10.159 MADE THAT PATH COMMON.
    classifyToolkitConfidence() returns the NOT TESTABLE label for
    permStatus:'unavailable' but with level:'preliminary', and the headline
    switch had no 'preliminary' arm - so it fell through to the final else and
    rendered "NO RELIABLE CSD SIGNAL - primary indicator (AC1) not rising" in
    the calm no-signal colour, with "AC1 (primary): ... (significance pending)"
    beneath it for a significance that was never coming. That is a claim about
    the DATA manufactured from the absence of a TEST. v10.159 made
    'unavailable' the common case (it stopped emitting raw-series p-values), so
    that release INCREASED the reachability of the exact mis-wording its own
    W-01 item 2 set out to remove. FIXED: an explicit 'preliminary' arm,
    splitting csdPermStatus==='unavailable' (NOT TESTABLE, purple, and it says
    in those words that this is not "no signal") from genuinely pending (grey),
    with ac1PrimaryTxt split on the same flag. AND the re-render: the headline
    is now a function, and csdToolkitRerender - which previously redrew only
    csdToolkitV - redraws it too, so a p-value that resolves after the spatial
    block has drawn actually reaches this box.

  BLOCKER 3 - THE AC1 PERMUTATION TEST WAS INERT AT EVERY SHIPPED WINDOW
    LENGTH, AND v10.159's CHANGELOG PRESENTED THAT AS A WIN.
    See TABLE 3. At 24-36 pooled months the measured false-positive rate of the
    AC1 delta is 0.2-3.2% against a nominal 5% and its power against a real
    AR(1) phi 0 -> 0.8 change is 0.7-13.1%: the statistic cannot produce a
    significant result whether or not one is there.
    MECHANISM, and this file documented it before re-introducing it. v10.151
    FIX 2 rejected 2-sample climatologies because subtracting a mean dominated
    by 2 samples forces that month's residual pair toward mirror images.
    v10.159 reinstated 1-2-sample climatologies everywhere by blending them
    with a fitted cycle, and the blend does not rescue it - the fitted cycle is
    itself estimated from the same readings. MEASURED against the shipped
    computeUsableClimatology()/deseasonalizeSeries() pair (12 calendar months
    each observed exactly k times, amplitude 10, noise SD 1, 4000 draws), the
    correlation between two residuals of the same calendar month is
      k=2  -0.969   k=3  -0.513   k=4  -0.326   k=5  -0.243
    v10.159's changelog reported part of this as a success - "the AC1
    statistic's false-positive rate also came down (18+12 sawtooth: 8.4% ->
    2.8%)". That is not a rate coming down. That is a statistic going dead.
    FIXED: the AC1 path gets its OWN floor, CSD_AC1_MIN_POOLED_MONTHS = 48
    pooled valid months = 12 calendar months x 4 samples, one more per month
    than v10.151's rule. CHOSEN FROM TABLE 3, not asserted: 48 is the shortest
    pooled total at which the null is calibrated (4.0-4.3%) and a balanced pair
    of windows clears the same "not missed two times in three" criterion the
    variance floor is held to (39.2% at 24+24). Below it, permutationTestDelta-
    Fixed returns NO AC1 p-value, with the measurement in the reason string;
    the VARIANCE p-value for the same two windows is unaffected and still
    reported. permutationTestAC1Fixed gets the same floor on its single window.
    WHAT THE FLOOR DOES NOT BUY, stated because it would otherwise be read in:
    power is governed by the SHORTER window, not the pooled total. At 48 pooled
    the measured power is 4.4% at 36+12, 17.3% at 30+18 and 39.2% at 24+24.
    And see TABLE 4: even at 36+36 a real phi 0.2 -> 0.7 change is missed two
    times in three. The AC1 test is now calibrated. It is still not powerful,
    and CSD_POWER_TABLE_TXT says so on screen.
    THE FALSE BLANKET CLAIM is corrected in all five places it appeared (this
    header, the permutationTestAC1Fixed comment, CSD_POWER_TABLE_TXT, STEP 4's
    on-screen block and S7E's): "the two-sided permutationTestDeltaFixed is NOT
    affected ... 4.5-6.2% across every configuration" was measured on the
    VARIANCE statistic and asserted of both.

  BLOCKER 4 - A SIGNIFICANCE DECISION COMBINED A DESEASONALIZED p-VALUE WITH A
    RAW DIRECTION (REGRESSION THIS SERIES INTRODUCED).
    S7F's _sig(t,d) took its direction from `d`, a delta computed by
    jsNodeStatsFixed with NO shared climatology, while t.pValue came from
    permutationTestDeltaFixed, which deseasonalizes against a pooled one. S7D's
    coupling test did the same with dCorr from jsPairCorrelation on RAW series.
    In v10.158 both sides were raw, so they agreed; v10.159 deseasonalized the
    test and left the direction raw. Because the p-value is TWO-SIDED, a
    significant FALL in deseasonalized AC1 paired with a RAW rise was reported
    as a rising signal. MEASURED (amplitude 3 smooth sine, noise SD 1, gap 0,
    1500 draws/cell): the displayed and tested deltas disagree in SIGN on 38.3%
    of draws at 18+12 unmasked and 40.4% at 10% masking, 35.2% / 41.0% at
    24+24, and 63.4% at 36+12; among only the draws that reached p<0.05, on
    7.4% at 24+24 (n=76) and 75.0% at 36+12 (n=76). Worst single case seen:
    displayed dAC1 -0.5729 against a tested delta of +0.7090.
    FIXED: the direction now comes from the test object itself - observedDelta
    for a delta test, corrAfter-corrBefore for the coupling test - which is the
    exact statistic each null was built around, so the two cannot disagree.
    AUDIT of every other significance decision in the file: STEP 3 COMPARE
    already used ac1Test.observedDelta (v10.139) and is correct. STEP 4's
    sAC1Sig/sVarSig/cAC1Sig/cVarSig had NO direction attached at all, so a
    SIGNIFICANT FALL counted as a "local signal" - one step further along the
    same fault, and fixed the same way. S7D's nAc1Sig and the Bonferroni hit
    count are direction-neutral by design and are worded as "changes", not
    "rises"; left alone. S7E's studySignal/refSignal are threshold-only
    verdicts with no p-value involved and are labelled as such.

  S5 - S7F's 12+12 DEFAULT WAS BAD IN THREE MEASURED WAYS. (1) Its
    false-positive rate was conservative and gap-dependent: 1.1-4.1% at a fixed
    0-month gap, 3.9-4.7% at a random 0-3 month one. (2) Power against a real
    x5 variance change was 21.7-24.8% at gap 0. (3) It pooled exactly
    CLIM_MIN_TOTAL_SAMPLES months, and that floor counts VALID months, so on
    cloud-masked Sentinel-2 FAI - the data this module serves - one masked
    month refused the whole test. MEASURED p-value emission, 3000 draws/cell,
    at 0 / 10 / 20 / 25% uniform per-month masking, under the v10.159 floor of
    24: 12+12 gives 100 / 8.7 / 0.3 / 0.1%, 18+12 gives 100 / 97.6 / 59.3 /
    33.5%, 24+12 gives 100 / 100 / 98.2 / 91.4%.
    FIXED: all three modules (S7D, S7E and S7F) now default to 24+12. S7D and
    S7E moved from 18+12, whose rate runs 5.5-7.0% (TABLE 1), about 40% above
    nominal. The per-box range stays 12-36, and a new check refuses a run whose
    POOLED total is under CLIM_MIN_TOTAL_SAMPLES before it starts, instead of
    the test discovering it later. The v10.159 sentence "the shipped defaults
    always produce a properly deseasonalized, calibrated p-value" was false on
    both halves and is replaced, not edited.

  S6 - S1's FIX NEVER REACHED THE SHIPPED PANELS. v10.159 S1 made the displayed
    and tested AC1 the same estimator, but only when both get the same
    climatology - and every shipped call site (S7D, S7E, S7F both paths) called
    jsNodeStatsFixed(series) with NO sharedClimatology. See the BLOCKER 4
    measurements: 38.3% sign disagreement at the S7D/S7E default. FIXED by a
    pooledClimFor() helper that builds one pooled climatology per BEFORE/AFTER
    pair and hands the same one to both windows, as STEP 3 COMPARE and STEP 4
    already did. AFTER: the same harness reports 0.0% sign disagreement and a
    worst |displayed - tested| of 0.0000 in every cell measured.

  S8 - HETEROSCEDASTICITY THE PERMUTATION NULL DOES NOT MODEL - DISCLOSED, NOT
    FIXED. Deseasonalizing deflates residual variance by how many own samples
    that calendar month has. MEASURED against the shipped functions (true
    per-reading noise SD 1.0, 20000 readings): residual SD by own-samples is
    1 -> 0.26, 2 -> 0.73, 3 -> 0.82, 4 -> 0.86, 5 -> 0.90, 6 -> 0.90. A
    calendar month seen once has its residual variance deflated about 13-fold.
    The test then pools all residuals and shuffles them freely, assuming an
    exchangeability the deseasonalizing destroyed. It is NOT demonstrated to
    inflate the rate at the configurations now allowed (TABLE 1), so it is
    written up in full above deseasonalizeSeries() as a latent hazard and cited
    as one reason the floors sit where they do. A correct fix is a stratified
    or restricted permutation, which is a different design.
    AND THE HEADLINE MECHANISM v10.159 LED WITH IS INERT. "A month with no
    samples is exactly the fitted value" is true of the returned object and
    irrelevant: EVERY caller builds the climatology from a SUPERSET of the
    series it deseasonalizes, so a zero-sample month contains no reading to
    look up. MEASURED by instrumenting deseasonalizeSeries() and driving every
    entry point with 20000 random ragged series: 80784 calls, 3541 of them
    against a climatology with at least one fully-imputed month, 2220827
    readings deseasonalized, and ZERO lookups of a zero-own-sample month. What
    actually removed the holes is the other half of the change - 1-2-sample
    months are now KEPT (blended) instead of REJECTED. Corrected in place.

  S9 - THE FLOOR AT 24 VIOLATED ITS OWN STATED CRITERION. The refusal string
    says "a p-value that would miss a real change two times in three is not
    reported", and at pooled 24 a real x5 variance change is missed 59.7-78.3%
    of the time (TABLE 2). CLIM_MIN_TOTAL_SAMPLES is now 26: the shortest
    pooled total at which the stated criterion is true under EVERY gap policy
    measured (power 47.1-58.5%), and the closest-to-nominal false-positive rate
    of any candidate (3.4-5.4%). 28 and 30 buy more power but run 5.2-6.9%.
    The whole 24-36 band is within a point or two on FPR; the power column is
    what separates the rows.

  S10 - S12's HEADLINE AC1 SILENTLY SWITCHED ESTIMATOR AND WAS BIMODAL.
    mkMoSST() builds exactly 24 monthly images, and v10.159 removed
    jsNodeStatsFixed's n>=48 gate, so S12 became "deseasonalized if 24 valid
    months survive, raw otherwise" - with 24 exactly on the v10.159 climatology
    floor, so ONE masked month flipped the estimator. MEASURED (24 nominal
    months, true AC1=0, amplitude 3 smooth sine, noise SD 1, 3000 draws/cell)
    under the v10.159 floor: 100.0 / 28.9 / 7.5 / 0.6% deseasonalized at 0 / 5
    / 10 / 20% masking, with mean displayed AC1 -0.0721 / +0.4534 / +0.6051 /
    +0.6478. The >0.6 and >0.3 colour bands and the realAC1>0.5 +
    varTrendRatio>1.3 trigger were all calibrated against the RAW estimator and
    were never re-cut - the same objection v10.159's own N4 raised against
    changing the variance estimator, not applied here. S9's floor of 26 pins
    S12 to ONE estimator, because 24 < 26: re-measured, 0.0% deseasonalized at
    every masking level and mean AC1 stable at +0.65 to +0.68. That is a
    coincidence of two constants, so it is CHECKED at run time instead of
    assumed - if the estimator ever changes, the panel says the thresholds are
    not calibrated for it. And the honest caveat on the raw estimator, measured
    in the same harness: on 24 months of pure seasonal cycle + white noise with
    TRUE AC1 = 0, raw AC1 exceeds 0.6 on 85.4% of draws and the "worth
    watching" trigger fires on 23.4%. That is now on screen beside the trigger.

  S11 - csdPermPAC1/csdPermPVar read .pValue directly, bypassing permUsable() -
    the exact contract v10.159 said it enforced everywhere. Harmless today, but
    they feed classifyToolkitConfidence() and the STEP 3 headline. Gated.

  NITS. (1) The STEP 4 label still read "Tests 6, 9, 12, 24, 36 and 48-month
    AFTER windows" - the pre-v10.158 list - and (2) another still hardcoded "17
    Earth Engine calls ... 1 + 6 + 6 + 4", so v10.159 N2's claim that all four
    such strings were derived was false for two of them. Both derived now.
    (3) CSD_SWEET_SPOT_NCALLS was 13 and contradicted the real budget the
    onClick handler counts (multiTotal = 17); the 4 permutation-test raw-series
    fetches were missing from it. It is now the same expression as multiTotal,
    and CSD_SWEET_SPOT_NPERMTESTS - previously referenced by nothing - is used.
    (4) nodeStatsDisclosure() pushed climatologySource unconditionally, so it
    never returned '' as its own header claimed, the "no disclosure needed"
    branch at every consumer was unreachable, and every clean run gained a
    ~180-character noise line. Gated on deseasonalized!==true || ac1PairsDropped.
    (5) The climatologyIsComplete() comment claimed the 12x3=36 rule and the
    36-month window floor "cannot drift apart". They already had: an imputed
    climatology reports complete at 26 valid months over 9 calendar months.
    Rewritten to say what "complete" now guarantees (no reading dropped) and
    what it does not (sample counts, record length). (6) S7D reported
    permSeriesNote for nodes[0] only, though nine nodes can have nine different
    permStatus values; every distinct note is now listed with the nodes it
    applies to. (7) fmtDepth's sentinel window was +-0.5, which is symmetric in
    the input but not in what is PRINTED, because fmtDepth rounds:
    fmtDepth(-9999.6,0) returned "-10000" with " very deep ocean" beside it -
    the same fabricated-10km-depth failure v10.159 S4 set out to remove, a
    tenth of a metre outside the guard. Window widened to +-1.0 and the ROUNDED
    value re-tested. (8) One line of trailing whitespace removed.

  CHANGELOG AUDIT - WHAT DID NOT REPRODUCE, WITHDRAWN BY NAME.
    Every headline figure below appeared in more than one place in this file
    with a different value. Re-measured this session with the design stated:
    - "18+12 false-positive rate 4.8%" (header W-01 table and the inline table
      above the permutation function) and "5.0%": NOT REPRODUCED. Measured
      5.5-7.0% at every seed, amplitude and shape tried (TABLE 1). About 25-40%
      above nominal, consistently. This is the main reason S7D/S7E's default
      moved to 24+12.
    - "12+12 false-positive rate 5.5% / 6.2% / 4.3-6.9%": NOT REPRODUCED.
      Measured 1.1-4.1% at a fixed gap and 3.9-4.7% at a random gap - i.e.
      CONSERVATIVE and strongly gap-dependent, not slightly high. The
      configuration is below the v10.160 floor and no longer runs.
    - THE PRIOR-WEIGHT SCAN (claimed worst cell 6.1 / 6.2 / 7.2 / 9.4 / 11.4 /
      13.6% at w = 0.3 / 0.5 / 0.75 / 1.0 / 1.5 / 2.0): NOT REPRODUCED as
      stated, and the claim never recorded its amplitude, gap policy, shuffle
      count or replicate count. Re-derived - worst cell over {18+12, 24+12,
      36+12} x {smooth, sawtooth, spike}, 1500 reps/cell, gap 0-3, 300
      shuffles, nominal 5%:
        seasonal amplitude 2:  w=0 6.7  0.3 7.1  0.5 7.0  0.75 7.2  1.0 7.3
                               1.5 7.3  2.0 7.1  4.0 6.3     -> FLAT
        seasonal amplitude 10: w=0 6.7  0.3 6.9  0.5 6.9  0.75 6.7  1.0 7.0
                               1.5 8.7  2.0 10.2  4.0 13.3   -> rises above w=1
      So the scan does NOT show 0.5 is better than 0 or 0.3; what it shows is
      that w <= 1 is safe and w >= 1.5 is not, and only at large seasonal
      amplitude. CLIM_HARMONIC_PRIOR_WEIGHT stays 0.5 - it is inside the flat
      region at both amplitudes - but the JUSTIFICATION is restated to what was
      actually measured.
    - "Power 18+12 x5: 11.5 / 17.1 / 61.5" and the header's "12.4 / 16.3 /
      63.0": these are v10.156/v10.158/v10.159 comparisons. The v10.160 value
      measured here is 63.7-74.8% (TABLE 1); the v10.156 and v10.158 columns
      are historical and were NOT re-derived this session, so they are marked
      as such wherever they appear rather than repeated as current fact.
    - The v10.158-vs-v10.159 comparison columns throughout the v10.159 entry
      below (13.5%, 24.0%, 25.4%, 39.8%, 10.4%, and the 10.4-39.8% range built
      from them) describe the behaviour of code that is no longer in this file.
      They were not re-derived in v10.160. They are left in place as the
      historical record of why the v10.159 change was made, and are NOT quoted
      anywhere in v10.160 as a current property of anything.
    Where a number appears on screen it now comes from the same tables above.

  RESIDUAL RISK - WHAT A LIVE EARTH ENGINE SESSION STILL HAS TO CHECK.
    Every change in v10.160 is client-side JS, unit-tested in Node; no ee.*
    call was touched and the per-button EE call count is unchanged (STEP 4 is
    17, now derived from one expression instead of three copies). Unverified
    without a live session:
    (1) That reduceRegions() really does emit a feature for a fully-masked
        month with the band property absent. BLOCKER 1's fix depends on it, and
        the harness asserts it by construction rather than observing it. If EE
        instead OMITS those features, fssValidMonthCount() equals the feature
        count and the fix degrades to v10.159's behaviour - it does not break,
        but the DUPLICATE logic stops firing. Worth one live check.
    (2) Where OISST actually ends, and therefore whether 54mo/60mo really do
        collapse onto 48mo for the AFTER start the sidebar suggests. The
        DUPLICATE logic is correct for whatever the real span is; the "39
        months" in the worked example is still an assumption.
    (3) Whether real Sentinel-2 FAI masking resembles the clustered synthetic
        closely enough for the emission rates under S5 to transfer. The new
        24+12 default was chosen with that uncertainty in mind - it sits 10
        months clear of the floor rather than on it - but the real masking
        pattern at a real reef has not been observed here.
    (4) On-screen layout at real string lengths. BLOCKER 2's NOT TESTABLE
        headline and S7D's per-note SERIES USED block are both longer than what
        they replace.
    KNOWN RESIDUALS, MEASURED, NOT FIXED:
    - The AC1 test is calibrated above its floor but still weak: TABLE 4 shows
      a real phi 0.2 -> 0.7 change missed 5 times in 6 at 24+24 and 2 times in
      3 at 36+36. Raising the floor further would trade that for refusing
      almost every real record. The table is on screen.
    - S8's heteroscedasticity is disclosed, not fixed.
    - permutationTestAC1Fixed's ONE-SIDED null remains anti-conservative for
      the reason set out under v10.159 S2 (v10.159 measured ~10.7% at n=48
      against a nominal 5%; not re-derived here). It is now additionally gated
      by the AC1 floor, which keeps it off short records, but the null itself
      still needs a block permutation. The function remains unreferenced.
    - The variance-ratio rate at the 24+12 default is 4.5-5.9%: at or slightly
      above nominal, not below it.

v10.159 FIX 17: three blockers - one of them a CALIBRATION REGRESSION this
  series introduced in v10.158 - plus eight should-fix defects and four nits.
  Every number below was produced in a standalone Node harness this session,
  against the pure-JS functions extracted from this file at v10.158 (the
  "before" column) and at v10.159 (the "after"). Earth Engine cannot be run
  from that harness; what that leaves unverified is listed under RESIDUAL
  RISK at the end of this entry. Where a fix is partial, it says so.

  W-01 (BLOCKER, REGRESSION INTRODUCED BY v10.158) - WE REPLACED A
    CALIBRATED TEST WITH ONE RUNNING AT UP TO 8x ITS NOMINAL RATE.
    v10.158 required a COMPLETE 12/12 climatology before deseasonalizing and
    ran the permutation test on RAW values when it could not get one. The
    shipped S7D/S7E month boxes pool 30 months and S7F's pooled 22 - all
    below the 36-month floor - so all three took the raw path ALWAYS, and
    printed "p=0.021 *** likely real (p<0.05)" out of it with no warning.
    MEASURED false-positive rate of the VARIANCE-RATIO test on pure
    seasonality + white noise with NO real change (3000 reps/cell, start
    phase and BEFORE/AFTER gap both randomised, nominal 5%). These are Monte
    Carlo estimates: the standard error is about 0.4 percentage points near
    5% and about 0.9 near 25%, so an individual cell moves by up to roughly
    one point between runs. The 3-5x gaps below are far larger than that; the
    third-decimal agreement of any single cell is not claimed:
      config                 seasonal shape   v10.156  v10.158  v10.159
      18+12 (S7D/S7E default) smooth            4.8%    13.5%     4.8%
      18+12                   sawtooth          4.7%    25.4%     5.6%
      18+12                   summer spike      5.1%    14.4%     4.5%
      12+12 (S7F new default) smooth          no p-val  24.0%     6.2%
      12+12                   sawtooth        no p-val  39.8%     6.0%
      12+12                   summer spike    no p-val  10.4%     5.6%
      12+10 (S7F old default) smooth           12.8%    14.3%   no p-val
      36+12 (above the floor) smooth            4.8%     5.4%     5.0%
      36+12                   sawtooth          4.4%     4.7%     5.5%
    ("sawtooth" and "summer spike" are deliberately adversarial non-harmonic
    cycles, included so the fix could not win just by assuming its own model.)
    ROUTE TAKEN: the reviewer's option (a) first, option (b) as the backstop.
    (a) computeUsableClimatology() IMPUTES the calendar months that cannot be
        estimated from their own samples. An order-3 seasonal harmonic model
        is least-squares fitted to every valid reading at once, and each
        calendar month's climatology is its own samples blended with that
        fitted cycle at CLIM_HARMONIC_PRIOR_WEIGHT=0.5 pseudo-observations.
        v10.160 S8 CORRECTION - THE HEADLINE MECHANISM QUOTED HERE IS INERT.
        "A month with no samples is exactly the fitted value" is true of the
        returned object and irrelevant to every caller. EVERY caller in this
        file builds the climatology from a SUPERSET of the series it then
        deseasonalizes (permutationTestDeltaFixed pools BEFORE+AFTER and
        deseasonalizes each; permutationTestAC1Fixed and jsNodeStatsFixed's
        own-window branch use the same series for both; permutationTestCorrDelta
        pools per node; fssPooledStats and pooledClimFor() pool the pair), so a
        calendar month with zero own samples contains no reading to look up.
        MEASURED by instrumenting deseasonalizeSeries() and driving every entry
        point with 20000 random ragged series: 80784 deseasonalize calls, 3541
        of them against a climatology containing at least one fully-imputed
        month, 2220827 readings deseasonalized, and ZERO of those readings ever
        looked up a zero-own-sample month.
        WHAT ACTUALLY REMOVED THE HOLES is the OTHER half of the change: a
        calendar month with 1 or 2 own samples used to be REJECTED by
        computeMonthlyClimatology()'s >=3 rule, so deseasonalizeSeries() nulled
        every reading in it and the survivors had gaps. Blending those months
        with the fitted cycle keeps them, and that is what gives the series no
        holes - which is what W-02/W-03 needed. The zero-sample imputation only
        makes climatologyIsComplete() return true.
    (b) Below CLIM_MIN_TOTAL_SAMPLES (24) valid months or
        [v10.160: that constant is now 26 - see S9 in the v10.160 entry above -
         and the AC1 statistic additionally needs CSD_AC1_MIN_POOLED_MONTHS=48]
        CLIM_MIN_DISTINCT_MONTHS (9) distinct calendar months, even that is
        not usable - see the measured power table under W-01's floor note -
        the permutation tests return NO p-value, with permStatus:'unavailable'
        and a reason, instead of a raw-series one. No PERMUTATION TEST in this
        file can now emit a p-value from a series whose seasonal cycle was not
        removed. (The STEP 5 sliding-window panel still runs a Mann-Kendall
        trend test on a raw trajectory, but it shows the raw and the
        deseasonalized rows side by side and labels both, which is what the
        permutation tests were failing to do.)
        BE PRECISE ABOUT WHAT THAT FLOOR IS FOR. Re-measured against the
        shipped function with the floor lowered to 8, the false-positive rate
        at pooled totals of 12-22 months is 3.1-6.8% - about the same as
        above the floor - so the floor is NOT buying calibration. What it
        buys is power and a sane parameter count: a REAL x5 variance change
        is detected 7.3% / 19.8% / 34.3% of the time at pooled 12 / 16 / 22
        against 45.7% at 24 and 71.3% at 48, and the order-3 fit spends 7
        parameters, leaving 5 residual degrees of freedom at pooled 12.
        [v10.160 AUDIT: SUPERSEDED by TABLE 2 in the v10.160 entry above, which
         measures the same thing at 3000 reps/cell with the gap policy stated
         and gap 0 separated from gap 0-3. The qualitative claim - that the
         floor buys power, not calibration - survives. The specific cells do
         not all reproduce, and the floor itself has moved from 24 to 26.]
    BOTH CONSTANTS WERE SCANNED, NOT CHOSEN BY TASTE. Prior weight w against
    the shipped function, 1500 reps/cell, worst cell across {smooth,
    sawtooth, spike} x {12+12, 18+12, 24+12, 36+12}:
      w      0.3    0.5    0.75    1.0    1.5    2.0
      worst  6.1%   6.2%   7.2%    9.4%  11.4%  13.6%
    [v10.160 AUDIT: WITHDRAWN - does not reproduce, and the scan never stated
     its seasonal amplitude, gap policy, shuffle count or replicate count. See
     the CHANGELOG AUDIT in the v10.160 entry above, and the re-derived scan in
     the comment above fitSeasonalHarmonics(). w=0.5 is kept, on the re-derived
     evidence rather than on this table.]
    0.5 sits inside the flat part. A larger w is NOT safer: it pushes the
    climatology toward a cycle the harmonics can represent and leaves what
    they cannot in the residual, which is what inflates the rate. A linear
    TREND term in the fit was also tried on the same 12 cells and was no
    better - worst cell 6.5% with it against 6.2% without - besides being
    the wrong thing to fit, since the pooled record's apparent trend is
    partly the before/after difference the test exists to measure.
    POWER WAS CHECKED, so this is not calibration bought by suppressing
    everything. AFTER window's second-half noise SD raised x3 or x5 - a REAL
    variance change - 2000 reps, fraction detected at p<0.05:
      config    v10.156  v10.158  v10.159
      18+12 x3   10.3%    13.9%    37.9%
      18+12 x5   12.4%    16.3%    63.0%
      12+12 x5  no p-val  26.6%    48.0%
      36+12 x3   46.3%    47.0%    47.8%
      36+12 x5   72.2%    71.8%    71.3%
      48+24 x3   70.5%    68.4%    70.3%
    [v10.160 AUDIT: the v10.156 and v10.158 columns are historical and were NOT
     re-derived. The v10.159 x5 column re-measures as 63.7-74.8% at 18+12 and
     86.3-87.1% at 48+24 (TABLE 1 in the v10.160 entry above); 12+12 is below
     the v10.160 floor and no longer runs at all. Use TABLE 1.]
    Power is up sharply at the short configurations and unchanged at the long
    ones. The AC1 statistic's false-positive rate also came down (18+12 on a
    sawtooth cycle: 8.4% -> 2.8%; 36+12 sawtooth: 4.6% -> 4.6%).
    [v10.160 BLOCKER 3: THIS SENTENCE IS THE DEFECT, NOT A RESULT. A nominal-5%
     test measured at 2.8% has not been calibrated, it has been switched off.
     Re-measured across pooled totals the AC1 rate is 0.2% at 24 pooled months
     and 2.0% at 30, with 0.7-3.5% power against a real AR(1) change - dead at
     every configuration this file shipped. See TABLE 3 in the v10.160 entry
     above and CSD_AC1_MIN_POOLED_MONTHS.]
    ITEM 2 - THE FLAG WAS RETURNED AND NEVER READ. deseasonNote was produced
    on every return path and rendered in exactly ONE of six places (STEP 3
    COMPARE). STEP 4's sweet-spot table, S7D's per-node AC1 column, S7D's
    pair-coupling column, S7E and S7F all printed a bare "p=...". Every
    p-value in this file now goes through permP()/permVerdictTag(), and
    permUsable() - which tests the machine-readable `deseasonalized` flag,
    not prose - gates every significance decision, including the ones that
    feed verdicts (ac1Sig/varSig in STEP 3, sAC1Sig..cVarSig in STEP 4,
    nAc1Tested/nCorrSig and the Bonferroni hit count in S7D, _sig() in S7F).
    ITEM 3 - TWO LINES ASSERTED THE OPPOSITE OF THE TRUTH. STEP 4 and S7E
    both said the permutation p-values "remain the only figures here with a
    known false-positive rate" while those same p-values were the 10.4-39.8%
    ones in the table above. Both lines are corrected, not deleted: the claim
    now holds only where a p-value is actually reported.
    ITEM 4 - THE DEFAULT CONFIGURATION CAN NO LONGER REACH AN UNCALIBRATED
    PATH. S7D/S7E/S7F month boxes are validated 12-36 instead of 4-36 (below
    one full annual cycle per window there is nothing to estimate the
    seasonal cycle from), so the pooled record is always >= 24 months, and
    S7F's AFTER default is 12 instead of 10 - at 12+10 the SHIPPED DEFAULT
    would otherwise have produced no p-value at all.

  W-02 (BLOCKER) - THE 12/12 RULE MADE DESEASONALIZING UNREACHABLE ON GAPPY
    DATA, AND ITS PROMISED DISCLOSURE DID NOT EXIST.
    MEASURED pass rate of climatologyIsComplete(), 10000 draws per cell, by
    record length and per-month drop rate, against the v10.159 rule:
      N=36   0%     5%      10%     15%     20%
      v158  100.0   15.9     2.1     0.3     0.0
      v159  100.0  100.0   100.0    99.8    98.2
      N=48  100.0   84.4    52.2    24.5     9.3   (v158)
      N=48  100.0  100.0   100.0   100.0   100.0   (v159)
      N=60  100.0   98.5    90.1    72.5    49.5   (v158)
      N=60  100.0  100.0   100.0   100.0   100.0   (v159)
    With seasonally-CLUSTERED gaps at a 20% base rate - how cloud and ice
    actually behave - the v10.158 rule passed 0.0% of 36-month records, 0.0%
    of 48-month and 0.3% of 60-month; v10.159 passes 81.5 / 99.9 / 100.0%.
    OISST (STEP 2/3/4) is interpolated and gap-free so it passed in practice,
    but FAI (S7C/S7D/S7E/S7F, Sentinel-2 at CLOUDY_PIXEL_PERCENTAGE<20) is
    heavily masked: those four modules were on the raw path at every window
    length, always. Fixed by the same imputation as W-01.
    AND THE DISCLOSURE: jsNodeStatsFixed() wrote climatologySource,
    climatologyNote, ac1PairsUsed, ac1PairsDropped and ac1MaxGapMonths, and
    a grep of the whole file confirmed NO CALLER READ ANY OF THEM. The file
    claimed the reason was "on screen" (W-02 entry) and that the function
    "says so, loudly" (jsNodeStatsFixed header). BOTH CLAIMS WERE FALSE WHEN
    WRITTEN; both are corrected in place. nodeStatsDisclosure() and
    nodeStatsDisclosureLines() now render all five fields at every consumer -
    S7C, S7D, S7E, S7F (both paths), S12, STEP 3 COMPARE and STEP 4 - and
    de-duplicate identical messages so a 9-node network prints one line, not
    nine.

  W-03 (BLOCKER) - STEP 4's 54- AND 60-MONTH ROWS WERE NOT INDEPENDENT TESTS,
    AND THE v10.158 RE-CUT IS WHAT MADE THEM SO.
    maxWindowLen is 60 and mkMoSSTRange() builds those months unconditionally;
    months past the end of OISST come back masked and are dropped. NOTHING
    validated that the AFTER start leaves 60 months of record - the date-box
    check only tests that the box is non-empty. With the AFTER start the
    sidebar's worked example tells users to type, OISST runs about 39 months:
      AFTER record  nominal 12/24/36/48/54/60 -> actual months
       39 months     12  24  36  39  39  39
       45 months     12  24  36  45  45  45
       60 months     12  24  36  48  54  60
    rowUnderpowered tested the NOMINAL w, so at 39 months all four of
    36/48/54/60 counted as powered: nPoweredWindows=4 gave a Bonferroni bar
    of 0.0125 for TWO distinct hypotheses, and testedRows / localLeaningRows /
    ac1TestedRows / uncorrectedLocalCount triple-counted one result into
    "3 of 4 windows lean LOCAL". Under v10.156's list the powered rows
    resolved to 24/36/39 - genuinely distinct - so the re-cut made it worse.
    FIXED: the real valid month count is computed once per row and shared by
    the stats table and the permutation table; a row whose real span equals a
    SHORTER row's is marked DUPLICATE and excluded from nPoweredWindows, the
    Bonferroni divisor, uncorrectedLocalCount, testedRows, ac1TestedRows,
    bestW and bestAc1W, and is labelled as a duplicate of that shorter row
    rather than silently dropped. Verified on the transcribed block: at a
    39-month AFTER record, 4 powered rows / alpha 0.0125 becomes 2 powered
    rows / alpha 0.0250 with 54mo and 60mo marked duplicates of 48mo; at 36
    months it becomes 1 powered row / alpha 0.0500; at 60 months nothing
    changes. The table also prints the real span in brackets after the
    nominal one whenever they differ.

  S1 - THE DISPLAYED AC1 AND THE TESTED AC1 WERE DIFFERENT STATISTICS.
    jsNodeStatsFixed() passed timestamps to jsLag1AC1() (so non-adjacent
    pairs were skipped) while statAC1ForPerm() did not, so S7C-F printed a
    dAC1 from one estimator next to a p= from another. Measured on the SAME
    raw gapped series, 12% of months masked, comparing jsNodeStatsFixed's
    realAC1 against statAC1ForPerm: v10.158 disagreed on 949 of 960 draws at
    a nominal n=24 (worst |difference| 0.5447, e.g. -0.7187 displayed vs
    -0.1739 tested), 816 of 819 at n=25 and 59 of 59 at n=30. v10.159
    disagrees on 0 of 960, 0 of 819 and 0 of 59, worst difference 0.0000.
    Fixed by passing the time vector to the observed statistic - and, because
    that alone would import S2, to every shuffle as well.

  S2 - PAIR-COUNT ASYMMETRY BIASED p DOWNWARD. If the observed statistic uses
    timestamps and the null does not, the observed uses FEWER lag-1 pairs than
    each shuffle, so the null is less spread than the statistic it is judging.
    Isolated on raw white noise (true AC1=0, n=30, one-sided, 400 shuffles,
    1500 reps, nominal 5%), P(p<0.05) by per-month masking:
      masking          0%     10%    20%    30%
      shuffles get no times   5.1%   6.7%   7.3%   9.1%   (v10.158 rule)
      shuffles get the times  5.5%   5.1%   4.5%   4.8%   (v10.159 rule)
    Fixed by shuffling the VALUES and keeping the SAME times, so observed and
    null use the identical pair mask. permutationTestAC1Fixed's p-value
    denominator is also now the number of shuffles that produced a statistic
    rather than nPerm, which matters now that jsLag1AC1 can return null.
    PARTIAL, AND SAID PLAINLY: permutationTestAC1Fixed's ONE-SIDED null is
    still anti-conservative for a DIFFERENT reason - subtracting a
    climatology estimated from the same data makes readings in one calendar
    month negatively dependent, and a free shuffle can place them adjacent.
    Measured at true AC1=0, n=48, no masking: 10.7% (v10.158) and 10.8%
    (v10.159) against a nominal 5%. That needs a restricted (block or
    whole-year) permutation, which is a different design and is NOT attempted
    here. The function is currently UNREFERENCED, so this is latent - but it
    must be fixed before anything calls it. The two-sided before/after
    permutationTestDeltaFixed() is NOT affected (4.5-6.2% across every
    configuration measured above).
    v10.160 BLOCKER 3 CORRECTION: that last sentence was TRUE OF THE VARIANCE
    STATISTIC ONLY. The AC1 statistic run through the same function was
    measured at 0.2% / 2.0% / 2.1-3.2% at 24 / 30 / 36 pooled months against a
    nominal 5%, with 0.7-13.1% power - dead, not calibrated. See the v10.160
    entry at the top of this file and CSD_AC1_MIN_POOLED_MONTHS.

  S3 - STEP 4 WAS MISSING ONE OF SIX VARIANCE-ARTIFACT CHECKS. The study side
    ran three, the control side two: the control's own BEFORE-window ratio
    check, isVarRatioArtifact(cBVF, cBVarW), was absent, so v10.157's "this
    mirrors STEP 3 exactly" was not true. Reproduced with cBVF=0.0005,
    cBVar=20.0x, cVF=0.58, cVR=18.0 (so cDVar=-2.0): v10.158 raised no flag
    and scored divScored=+10.20, reaching "LOCAL CSD + AC1 CONFIRMED" off a
    control baseline STEP 3 would have excluded outright; v10.159 flags it
    and divScored becomes null. Both sides now run the identical three
    checks, so the claim is true.

  S4 - fmtDepth() PRINTED THE NO-DATA SENTINEL AS A REAL DEPTH. Its guard was
    the GEBCO physical range (v<-11500 || v>9500), and -9999 sits inside it.
    Verified: fmtDepth(-9999,0) returned the string "-9999", rendered as
    "-9999 m (GEBCO) very deep ocean" - a fabricated 10 km depth presented as
    a measurement. The old fmt() rejected v<-900 for exactly this reason.
    Explicit sentinel guards (-9999, 9999, -32768) added to BOTH fmtDepth()
    and classifyDepthLabel(), so the number and the label cannot disagree.
    After: all three return "n/a" and "" respectively; -1500 still returns
    "-1500" / " deep ocean".

  S5 - THE CHLOROPHYLL INVERSION WAS NON-MONOTONE AT EXACTLY 0, AND ONE
    READOUT CONTRADICTED IT. (i) the `cv>0` guard left s2 at the mid default
    50 when cv===0, so chl=0 scored CCS 34 while chl=0.02 scored CCS 28 - the
    header's "CCS rises monotonically with chlorophyll" was false at the one
    value where it is easiest to check. Verified across chl = 0, 0.02, 0.05,
    0.09, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 1, 1.5, 2, 3, 8 at sv=28.5,
    tv=0.02, nv=3e-5, turv=0.02, dhw=1: v10.158 gives CCS 34, 28, 28, 28, 29,
    29, 32, 32, 36, 41, 44, 44, 46, 46, 46 (NOT monotone); v10.159 gives 28,
    28, 28, 28, 29, 29, 32, 32, 36, 41, 44, 44, 46, 46, 46 (monotone).
    (ii) the sidebar label was still on the OLD direction - at 1.5 mg/m3 it
    printed "[good]" while s2 scored that same value 80/100 STRESS on the
    same screen. Relabelled to match s2's own breakpoints: <0.1 oligotrophic
    (reference reef state), <0.2 very low, <0.45 below the GBR guideline, <1
    above the 0.45 GBR guideline, <2 enriched, else bloom.

  S6 - THE CENSUS AND THE SUB-SCORES DISAGREED AT EXACTLY 0. _has(0) is TRUE,
    so a measured 0 counted as an input, but `cv>0` and `nv>0` kept s2 and s6
    at their no-data defaults. Verified: computeScore(28.5, 0, 0.02, 0, 0.02,
    1, ...) returned s2:50, s6:20, nInputs:6, dataCompleteness:100 - two
    sub-scores sitting at their "no data" values inside a score reported as
    100% complete. Reachable: TROPOMI NO2 column density can come back at or
    below zero over clean water. Guards dropped; <=0 now takes the lowest
    bin. After: s2:5, s6:10, nInputs:6, dataCompleteness:100.

  S7 - TOE-02's DISCLOSED REACH WAS OVERSTATED. See the corrected TOE-02
    entry below: "binds at k=4,6,8,12,16, stops at k=24" holds only if the k
    valid years are CONSECUTIVE. Scattered across the nominal window - the
    actual "8 of 32 salinity pixel" case - Sxx is 10-20x larger at small k
    and the gate stops binding at k=11 (even scatter) or k=14 (random-draw
    median), not 24. The honest statement is THE GATE BINDS AT k<=8, and
    that is what the header and the startup block now say.

  S8 - toeNum()'s COUNT PROBE ACCEPTED AN UNRELATED NUMERIC KEY. Its fallback
    scan took ANY numeric property, and the count call passes excludeSub=null.
    Verified on the shipped v10.158 function: toeNum({sst_count:null,
    other:7}, <candidates>, null) returned 7, with countSource:'actual' and
    no caveat - defeating the whole point of TOE-01. The fallback now requires
    the key to END IN '_count'. After: that same input returns null (and the
    verdict line then says 'nominal-fallback'); {sst_count:null,
    foo_count:12} still returns 12 and {sst_count:31} still returns 31. The
    correlation probe was re-checked and is safe - its candidate list matches
    EE's real output key and its excludeSub='value' keeps a p-value out.

  N1 - THE "13.17 -> 1.61" FIGURE WAS ATTRIBUTED TO TOE-01 ALONE. Re-derived
    in Node against this file's own tCrit95 table, for the 8-of-32 salinity
    pixel at snr=3.20: v10.156 t=13.17 vs bar 2.042; TOE-01 ALONE (n 32 -> 8
    on the derived path) t=0.58 vs bar 2.447; TOE-01 + TOE-02's MEASURED
    r=0.550 t=1.61 vs bar 2.447. Both fixes are needed to reach 1.61.
  N2 - "the 6/13/24 display strings are now derived too" WAS FALSE in
    v10.158: three STEP 4 panel labels and one startup line still held the
    literals. The window list is hoisted to CSD_SWEET_SPOT_WINDOWS at module
    scope, the onClick handler reads the same array, and all four strings are
    derived from it.
  N3 - 54 and 60 are NOT getCalibratedThresholds() keys (24/36/48 are), and
    nearest-match sends both to 48. That is the conservative direction - the
    null spread narrows as the window lengthens, so a 48mo cutoff at 54/60
    months is WIDER than the right one and fires less often - but the basis
    string did not say the request was outside the calibrated range at all.
    It does now, and STEP 4 prints which rows used an extrapolated cutoff.
  N4 - jsNodeStatsFixed()'s gap-awareness is PARTIAL and is now documented as
    partial rather than implied to be complete: only realAC1 is gap-aware;
    jsVarianceHalves() and jsSkewness() still treat kept values as evenly
    spaced. Deliberately NOT changed - the v10.153 variance-ratio cutoffs in
    CSD_NULL_CALIBRATION were measured against the current index-based
    split-half estimator, and swapping the estimator without re-running that
    13.6-year Scripps calibration would silently invalidate every varr
    threshold. ac1MaxGapMonths is returned and now shown on screen so a
    reader can see how far from evenly-spaced the series is. The n>=48 gate
    in the self-computed branch is GONE, which removes the contradiction
    between that literal and the v10.158 note claiming the rule "implies 36
    by construction": the branch is now gated on one thing only, whether
    computeUsableClimatology() can produce a climatology at all.

  RESIDUAL RISK - WHAT A LIVE EARTH ENGINE SESSION STILL HAS TO CHECK.
    Every change here is client-side JS, unit-tested in Node; no ee.* call
    was touched and the EE call count per button is unchanged. Unverified
    without a live session: (1) that OISST really does end where the W-03
    worked example assumes, and therefore that 54mo/60mo really do collapse
    onto 48mo for the suggested AFTER start - the DUPLICATE logic is correct
    for whatever the real span turns out to be, but the specific 39 is an
    assumption; (2) that the real FAI masking pattern resembles the clustered
    synthetic used above closely enough for the 81.5%/99.9% imputation pass
    rates to transfer; (3) the on-screen layout of the new SERIES USED /
    SERIES QUALITY blocks at real string lengths.
    KNOWN RESIDUAL, MEASURED: at 12+12 - the shortest configuration that can
    now produce a p-value at all - the variance-ratio rate sits a little above
    nominal, 5.3-6.2% across the three seasonal shapes over repeated 2000-3000
    rep runs, against v10.158's 10.4-39.8% in the same cells. It is at nominal
    from 18+12 upward.
    [v10.160 AUDIT: BOTH HALVES WITHDRAWN. 12+12 re-measures at 1.1-4.1% (gap 0)
     and 3.9-4.7% (gap 0-3) - BELOW nominal, not above - and 18+12 re-measures
     at 5.5-7.0%, not "at nominal". 12+12 pools 24 months, below the v10.160
     floor of 26, and no longer produces a p-value. See TABLE 1 and TABLE 2.] At 25-35% seasonally-clustered masking the imputation refuses
    most records, so few p-values are emitted at all and the rate cannot be
    estimated from the ones that are. And permutationTestAC1Fixed's one-sided
    null remains at ~10.7% for the reason set out under S2; it is dead code
    today and must not be wired up until that is fixed.

v10.158 FIX 16: two blockers created by the v10.156 BUG-08 fix, plus eight
  smaller defects, found by a global coordinate sweep of the whole file.
  Every claim below was reproduced BEFORE and re-checked AFTER in standalone
  Node harnesses run against the extracted pure-JS functions. Earth Engine
  cannot be run from that harness; the EE-side assumptions are listed under
  RESIDUAL RISK at the end of this entry.

  W-01 (BLOCKER) - THE ENFORCED MINIMUM WINDOW WAS EXACTLY WHERE THE
    PERMUTATION TEST CANNOT RUN. v10.156 BUG-08 set CSD_MIN_WINDOW_MONTHS
    = 24 and pointed users at it. But computeMonthlyClimatology() requires
    CLIM_MIN_SAMPLES_PER_MONTH = 3 samples in each calendar month (v10.151
    FIX 2), and 24 months supplies exactly 2 - so 0 of 12 calendar months
    survive, deseasonalizeSeries() nulls every point, and the tool's only
    real significance test returned pValue:null / 'statistic unavailable'.
    Measured across totals 20->60: permutationTestDeltaFixed worked at 20,
    22 and 28+, and returned NULL at 24 and 26; permutationTestAC1Fixed
    died outright at n=24 with 'insufficient data after processing (n=0)';
    permutationTestCorrDelta had the same >=24 gate. FIXED:
    (a) CSD_MIN_WINDOW_MONTHS 24 -> 36, and every place that number is
        written - input validation (the stale nMonthsIn<4 lower bound is
        gone), the refusal message, the "24-60" range text, the months
        label and placeholder, the worked example (which told users to
        type 24), the UNDERPOWERED band boundary and the two
        getCalibratedThresholds(...||24) fallbacks.
    (b) NEW CSD_DESEASON_MIN_MONTHS = 12 * CLIM_MIN_SAMPLES_PER_MONTH,
        i.e. DERIVED (=36), not typed. All three permutation tests gate on
        it, so the gate and the climatology floor cannot disagree again.
        The climatology is ALSO checked for completeness before use, so a
        ragged 40-month record with a hole fails too.
    (c) Below the floor, or on an incomplete climatology, the tests now RUN
        ON RAW (non-deseasonalized) values instead of returning null, and
        every return carries a deseasonNote saying which series was used.
        Verified: no NULL p-value at any total from 20 to 60 months.
        CAVEAT, stated plainly: a raw p-value carries the v10.101 seasonal-
        cycle inflation (the harness shows raw AC1 p=0.000 on pure seasonal
        noise below the floor). That is WHY the floor exists; the raw path
        is a disclosed fallback, not an endorsement of short windows.
    (d) FIND SWEET SPOT's window list re-cut [6,9,12,24,36,48] ->
        [12,24,36,48,54,60]. The old list left only TWO rows above the new
        floor. Now four powered rows (36,48,54,60) and two sub-floor
        diagnostic rows (12, and 24 specifically because 24 was the v10.156
        floor, so the table shows what it produced). Sub-floor rows stay
        excluded from both sweet-spot picks and every tally, unchanged.
        EE CALL BUDGET - CHECKED: it does NOT depend on the literal 6.
        multiTotal is windowLengths.length*2+1+4 and the fire loop is a
        forEach, both derived; only on-screen strings said "6"/"13"/"24
        tests". The list is kept at 6 entries, so the per-run EE cost is
        unchanged.
        v10.159 N2 CORRECTION: "and those are now derived too" WAS FALSE.
        The literals were still sitting in two STEP 4 panel labels, in the
        STEP 4 progress text and in one startup print. v10.159 hoists the
        list to CSD_SWEET_SPOT_WINDOWS at module scope and derives all four
        strings from it, so the claim is now true.
    (e) BONFERRONI DIVISOR - the previous review's finding CONFIRMED still
        wrong: alpha was 0.05/windowLengths.length (6) while the family
        actually counted is only the POWERED windows (3 before, 4 now) -
        sub-floor rows are explicitly never counted for or against. The
        divisor is now the same nPoweredWindows the panel reports as its
        denominator, computed once and used in both places.

  W-02 (BLOCKER) - LAG-1 AC1 WAS COMPUTED OVER NON-ADJACENT MONTHS.
    deseasonalizeSeries() returns v:null for any month whose calendar month
    missed the climatology floor; callers filtered those out and handed the
    survivors to jsLag1AC1() as a plain array, which paired ARRAY-adjacent
    values and called it lag-1. Measured: n=28 -> 12 of 28 months kept and
    2 of 11 "lag-1" pairs were 9 REAL MONTHS apart; n=26 -> 6 of 26 kept,
    2 of 5 pairs 11 months apart; n=30 -> 18 of 30, gaps of 7. Clean only
    at n>=36. It was returned as deseasonalized:true with no warning.
    FIXED, and deliberately NOT by relying on the 36-month floor alone,
    because jsNodeStatsFixed() and other callers can still reach it:
    (a) jsLag1AC1(resid, times) takes OPTIONAL timestamps and skips any
        pair that is not one calendar month apart. Callers whose order is
        meaningless by construction (the shuffled nulls inside the
        permutation tests) pass nothing and are unchanged.
    (b) NEW jsLag1PairAudit() reports how many pairs are fake and how big
        the worst gap is; jsNodeStatsFixed returns ac1PairsUsed /
        ac1PairsDropped / ac1MaxGapMonths and an explicit climatologyNote.
    (c) A climatology may now only be SUBTRACTED when it is complete
        (12/12 calendar months), so on the normal path no month is ever
        dropped and the gap cannot arise in the first place. Note the
        arithmetic: 12 months x 3 samples = 36 = CSD_MIN_WINDOW_MONTHS.
    A series with dropped months can no longer silently report an intact
    lag-1 AC1: it reports the gap, uses only genuine pairs, and returns
    null if fewer than 3 genuine pairs remain.

  W-03 - POOLED CLIMATOLOGY BYPASSED THE 48-MONTH FLOOR (new failure mode
    introduced by the v10.156 BUG-07 pooling fix). jsNodeStatsFixed's
    `if (sharedClimatology)` branch had NO length or completeness check at
    all; the n>=48 guard sat only on the `else if`. Reproduced: with a
    30-month pooled record (6 of 12 usable climatology months) the 6-, 9-
    and 12-month windows ALL collapsed onto the same 6 surviving points and
    returned one IDENTICAL AC1 (-0.790 in the harness, -0.646 in the
    audit's run) against -0.211/-0.468/-0.500 unpooled, labelled
    climatologySource:'shared (pooled)' with no warning. RULE PICKED, and
    it is now the same rule in BOTH branches: a climatology is usable only
    if COMPLETE - all 12 calendar months cleared the >=3-samples floor,
    which implies at least 36 months of record, the same floor as W-01.
    On failure it does not return a collapsed number: it falls back to the
    raw computation and says so in climatologySource and climatologyNote.
    After the fix the same three windows give 0.213 / 0.597 / 0.625 - three
    different windows, three different numbers.
    v10.159 CORRECTION - "with the reason on screen" WAS FALSE WHEN WRITTEN.
    jsNodeStatsFixed() wrote climatologySource, climatologyNote, ac10PairsUsed,
    ac1PairsDropped and ac1MaxGapMonths, and NOT ONE CALLER READ ANY OF THEM -
    confirmed by grep across the whole file. The reason was computed and
    discarded. v10.159 adds nodeStatsDisclosure()/nodeStatsDisclosureLines()
    and calls them at every jsNodeStatsFixed consumer (S7C, S7D, S7E, S7F,
    S12, STEP 3 COMPARE and STEP 4), so the claim is now true.
    The 12/12 RULE ITSELF is also superseded - see W-01/W-02 v10.159 below:
    refusing to deseasonalize was measured to make the permutation test's
    false-positive rate WORSE, not better.

  W-08 - CHLOROPHYLL DIRECTION WAS BACKWARDS FOR CORAL (user-decided).
    computeScore's s2 scored clear water as HIGH stress (cv<0.05 -> 85) and
    a bloom as LOW stress (cv>2.0 -> 5), and s5 added +25 for cv<0.3.
    Measured: raising chlorophyll 0.04 -> 1.5 mg/m3 LOWERED CCS by 16
    points (48 -> 32) at every site tested. That mapping is correct for S8
    computeAquaculture(), where chlorophyll is FOOD for a seaweed crop, and
    appears to have been carried across from it - but it is backwards for a
    coral reef stress score: oligotrophy is the natural reef state and
    nutrient enrichment is a documented reef stressor (macroalgal
    overgrowth, reduced calcification, higher bleaching and disease
    susceptibility). REVERSED. New s2 breakpoints and why: <0.1 -> 5 (clear
    oligotrophic reef water, the reference state); <0.2 -> 15; <0.45 -> 35;
    <1.0 -> 60; <2.0 -> 80; else 95. 0.45 mg/m3 is the GBR annual-mean
    chlorophyll water-quality guideline (De'ath & Fabricius 2010), so it is
    the first real signal rather than an invented cutoff. The s5 term flips
    with it, to cv>0.45 -> +25. After: CCS rises monotonically with
    chlorophyll, 0.04 -> 1.5 now goes 31 -> 47. THIS REVERSES A LONG-
    STANDING OUTPUT: any CCS recorded before v10.158 at a site with
    non-trivial chlorophyll is not comparable with one recorded after.
    DISCLOSED LIMIT: satellite chlorophyll is unreliable in optically
    complex nearshore water (CDOM and resuspended sediment inflate it),
    which is exactly where reefs sit - hence 15% of the composite, not
    more, with turbidity scored separately alongside it.
    computeAquaculture() is deliberately UNCHANGED - its direction is
    correct for A. taxiformis / kelp and is a different question.

  W-19 - `undefined` STILL REACHED THE HARDCODED ccs=30 THAT v10.156 BUG-05
    EXISTED TO REMOVE. The six sub-score branches tested `if(sv!==null)`,
    which `undefined` and `NaN` both pass, so they produced NaN while the
    BUG-05 input census (which uses _has(), rejecting both) reported the
    input as PRESENT. Reproduced exactly: computeScore(undefined, 0.5,
    0.03, 7e-5, 0.1, 5, ...) returned s1:NaN and ccs:30 while reporting
    nInputs:5, dataCompleteness:85 - because `if(isNaN(ccs))ccs=csat;
    if(isNaN(ccs))ccs=30;` resurrected the exact default BUG-05 removed.
    FIXED: every sub-score guard now uses the SAME _has() predicate as the
    census, so the two cannot disagree; and the isNaN(ccs) fallback no
    longer invents 30 - it routes to the insufficient-data path (ccs:null,
    insufficientData:true), which is now a single shared helper so the two
    "no score" returns cannot drift apart. Not reachable from the live call
    site today; the guards must not disagree regardless.

  W-14 - DEAD SMALL-SAMPLE CAVEAT. computeToESignal's nCaveat fired only at
    df<10, i.e. n<12, but the function returns early at `if(n<12)` - so it
    was unreachable (confirmed: n=14..34 all gave nCaveat:null). Threshold
    moved to df<30, where t_crit still differs materially from the ~1.96
    the caveat text refers to (2.042 at df=30, 2.201 at df=11); it now
    fires across n=12..31. computeToESignal() is UNREFERENCED dead code -
    the live path is calcToE() - and is KEPT, not deleted, because it is a
    correct standalone implementation; this makes it correct-if-used rather
    than leaving a branch that provably cannot execute.

  W-17 - v10.156 BUG-04 missed two on-screen version markers. All four
    self-identifying markers are now v10.158: the sidebar title, the S13
    section header (was v10.149), the footer (was 'v10.149 + GEM') and the
    per-click console header. Historical changelog lines describing what
    v10.149 actually changed are untouched, as they should be.

  PRE-EXISTING, one line each:
  W-04 - predictBleachingProbability was called with raw Math.abs(bedrock),
    no depth cap and no reef-zone gate. Depth is the model's largest
    coefficient (+0.345) on a scaler with mean 7.04 m / scale 4.21 m, so
    measured: 50 m -> 97.3%, 200 m -> 100.0%, 1000 m -> 100.0%, and a
    4000 m click sits 947 SD outside the training range and still printed
    "100.0%" in confident red - at Svalbard, the Antarctic Peninsula and
    both poles. FIXED: training range stated as mean +/- 3 SD = 0-20 m;
    20-50 m is capped to 20 m and flagged as a BOUND, not an estimate;
    beyond 50 m (the file's own shallow-reef-zone bin) the model REFUSES
    with an explanation instead of returning a number. The whole row is
    additionally gated to tropical non-EBUS reef latitudes.
  W-05 - any click deeper than 900 m rendered as "n/a m (GEBCO) very deep
    ocean" - the number suppressed by fmt()'s `v < -900` no-data sentinel
    while the label derived from that same number survived. NEW fmtDepth()
    keeps the real depth and uses the actual GEBCO physical limits
    (-11,500 / +9,500 m) as its no-data test. Used for the depth readout
    and the depth line in the console dump.
  W-06 - computeAquaculture did `var sstAnnual = sv || 25, sstPk = sv_peak
    || sv || 25;` and 0 and NaN are both FALSY. Reproduced: a measured 0
    deg C became 25 deg C and was reported "DECENT - SST and chlorophyll
    both favorable" at 100% confidence; all-NaN input produced a status
    string containing "NaN mg/m3", also at 100% confidence. Replaced with
    an explicit _aqHas() null/NaN test throughout, matching what
    computeScore already does. After: 0 deg C correctly vetoes as too cold
    (score 5) and all-NaN correctly returns INSUFFICIENT DATA at 0%
    confidence.
  W-18 (latent) - predictBleachingProbability guarded null/undefined but
    not NaN, so p:NaN printed as "NaN%" coloured GREEN (NaN > 0.6 is
    false). NaN is now rejected at the input, the output is checked for
    finiteness, and the render checks for a real number rather than merely
    the absence of an error string.
  W-20 (latent) - computeNonOverlappingTrajectory(values, 0) was an
    infinite loop (`i += 0`), which in the GEE editor locks the tab with no
    error at all. Confirmed by a timed-out child process. Any windowSize
    below 1 (or non-finite) now returns an empty trajectory.
  ES5 SWEEP - Object.assign (ES6) at getInSituBaseline() replaced with an
    explicit hasOwnProperty copy loop; this file already documents ES6
    builtins failing in the GEE sandbox (String.prototype.repeat, v10.149).
    A whole-file grep for Array.from / .includes( / .startsWith( /
    .endsWith( / .repeat( / .find( / Object.values / Object.entries /
    let / const / arrow functions / template literals found NOTHING else.
    Only ONE Object.assign existed, not the two the sweep reported.
  NULL-p WORDING - classifyToolkitConfidence said "the permutation test has
    not returned a p-value yet" even when the test HAD run and structurally
    could not produce one. It now takes a permStatus and says "NOT TESTABLE
    AT THIS WINDOW LENGTH ... waiting will not change it" in that case.

  RESIDUAL RISK / NEEDS A LIVE EARTH ENGINE SESSION:
    - Every change above is pure client-side JS and was unit-tested in
      Node. NOTHING server-side (ee.*) was touched.
    - The re-cut window list makes the single AFTER fetch 60 months instead
      of 48. Same number of EE calls, one longer ImageCollection. A user
      whose AFTER start date leaves fewer than 60 months of OISST before
      today will get a shorter series for the longest rows - the slicing
      already handles that, but it has not been observed live.
    - The 36-month floor means FIND SWEET SPOT now needs the BEFORE window
      and 60 months of AFTER data to fill all four powered rows. Whether
      real sites routinely have that has not been checked against EE.
    - Raising the floor to 36 makes some previously-accepted STEP 2 runs
      refuse. That is intended, but it changes the tool's reachable
      behaviour for existing users and has not been exercised in the UI.
    - The bleaching-probability reef gate uses the same tropical/EBUS test
      the DHW label already uses; it has not been checked against a live
      click at a borderline site (e.g. lat 29).
    - W-08's direction reversal changes every CCS at a site with
      measurable chlorophyll. Only the pure function was tested; the
      downstream Bowl Depth / omega0 / tau / p5yr readouts derive from ccs
      and will move with it. They were NOT separately re-validated.

v10.157 FIX 15: two ToE blockers and one v10.156 regression, found by an
  adversarial review. Every claim below was reproduced BEFORE and checked
  AFTER in standalone Node harnesses run against the extracted functions.
  Earth Engine cannot be run from that harness, so the EE-side assumptions
  are listed explicitly under RESIDUAL RISK at the end of this entry.

  TOE-01 (BLOCKER) - the ToE sample-size penalty used a NOMINAL record
    length, not the actual valid-point count. Every caller of calcToE()
    passed a hardcoded constant (44/27/32/7/4/32) and that constant drove
    BOTH Sxx = n(n^2-1)/12 and df = n-2. This file's own S17b block already
    documented the constant as routinely wrong: "Salinity is nominally a
    32yr record (1993-2024) in S17 above, but only N years had real, valid
    HYCOM data at this exact point". At a pixel with 8 valid annual values
    an n of 32 inflates Sxx by 65x, understates se(slope) by 8.1x and gives
    the t-test df=30 instead of df=6 - BOTH anti-conservative, so the very
    gate v10.156 added to stop short records emerging could itself pass an
    8-point record. Fixed: ee.Reducer.count() on each annual band gives the
    REAL per-pixel valid-year count; it is fetched inside the SAME
    ee.Dictionary that was already being evaluated, so it costs ZERO extra
    .evaluate() round trips (the property v10.89/v10.141 built those split
    dictionaries to protect). Sxx, SStot and df now come from that count;
    the nominal length is kept only for display and for the amplitude term
    signal=|slope*nominal_years|, which is correctly a calendar-time
    quantity. Every verdict now reads "n=8 of 32 nominal yr". GUARDS:
    count missing -> fall back to nominal AND say so on screen (the df may
    then be far too generous); count<3 -> df<1, no slope test exists,
    CANNOT emerge; count>nominal -> clamped and noted. Unit-tested: the
    8-of-32 salinity pixel at snr=3.20 went from EMERGED (t=13.17 vs bar
    2.04) to "not yet".
    v10.159 N1 CORRECTION - THE ATTRIBUTION OF THAT NUMBER WAS WRONG. The
    v10.158 text credited the landing value t=1.61 to TOE-01. TOE-01 ALONE -
    i.e. n 32 -> 8 on the DERIVED path, everything else unchanged - gives
    t=0.58 against a bar of 2.447. The 1.61 needs TOE-01 AND TOE-02's
    MEASURED r=0.550: 0.550*sqrt(6/(1-0.550^2)) = 1.61. All three figures
    re-derived in Node against this file's own tCrit95 table: 13.17 vs 2.042
    (v10.156, n=32 derived), 0.58 vs 2.447 (TOE-01 only), 1.61 vs 2.447
    (TOE-01 + TOE-02). Both fixes are needed to reach 1.61; neither alone is.
    CONSERVATISM NOTE: Sxx=n(n^2-1)/12 assumes the n valid years are
    CONSECUTIVE. Scattered across a longer window the true Sxx is larger,
    so this understates Sxx, overstates se(slope) and understates the
    DERIVED t - the conservative direction. The measured-r path (TOE-02)
    does not make that assumption at all.

  TOE-02 (BLOCKER) - the v10.156 t-gate was mathematically VACUOUS for 4 of
    the 6 ToE variables. Because SStot was ASSUMED as n*noise^2 rather than
    measured, t was a deterministic function of snr and n: with
    A = slope^2*Sxx/SStot, algebra gives A = snr^2*(n^2-1)/(12n^2) and
    t^2 = (n-2)A/(1-A), so at snr=2 exactly, t ~= sqrt((n-2)/2). Verified
    against the shipped v10.156 code: n=4 -> t=0.953 vs bar 4.303 (binds);
    n=7 -> 1.557 vs 2.571 (binds); n=12 -> 2.224 vs 2.228 (marginal);
    n=13 -> 2.335 vs 2.201, n=27 -> 3.532 vs 2.064, n=32 -> 3.870 vs 2.042,
    n=44 -> 4.581 vs 2.021 (all VACUOUS). So for SST(44), CHL(27), SAL(32)
    and DO(32) the condition snr>=2 ALREADY implied t>=t_crit and the new
    gate could never change a verdict; it bound only NO2(7) and pH(4). The
    t carried no information independent of snr because snr is itself built
    from the TOTAL SD, which already contains the signal.
    Fixed: ee.Reducer.pearsonsCorrelation() on the SAME (t, value) band
    pair already fed to linearFit gives the MEASURED per-pixel correlation,
    and the textbook exact slope t-test t = |r|*sqrt((n-2)/(1-r^2)) is used
    instead of re-deriving t from snr. Also fetched inside the SAME already-
    evaluated dictionary - no extra .evaluate(). The derived path is kept
    ONLY as a fallback when no measured r arrives, and every verdict says
    which produced it ("MEASURED r=0.550" vs "DERIVED from SNR"). |r|>=1
    (a perfect fit, infinite t) is marked unreliable and cannot emerge,
    exactly as the existing SSres<=0 branch is. EMERGED remains
    (snr>=2.0 && slopeSignificant), unchanged.
    HOW FAR THIS ACTUALLY GOES - stated plainly, because it is only a
    PARTIAL fix. The verdict now depends on a quantity the code measures
    rather than one it re-derives, and the minimum |r| the gate demands is
    0.950 at n=4, 0.755 at n=7, 0.576 at n=12, 0.382 at n=27, 0.349 at
    n=32 and 0.298 at n=44; at every one of those n there exist records
    with snr>=2 whose verdict the gate flips (checked numerically). BUT: on
    a record that is COMPLETE and gap-free (count == nominal, evenly spaced
    annual steps) and where EE's stdDev really is the population SD, the
    identity r^2 = snr^2*(n^2-1)/(12n^2) still holds EXACTLY, so r is not
    free, the measured-r t reproduces the derived t to the digit, and the
    gate is STILL VACUOUS for n=27, 32 and 44 (implied r at snr=2 is 0.577
    at every n, which clears the 0.382/0.349/0.298 bars). It stops being
    vacuous exactly where the record is INCOMPLETE - i.e. wherever TOE-01
    finds count < nominal, because then snr is measured over nominal
    calendar years while r and df come from the valid years, and the two
    decouple.
    v10.159 S7 CORRECTION - THE v10.158 REACH CLAIM WAS OVERSTATED. It said
    "binds at k=4,6,8,12,16 and stops binding at k=24". That is true only
    under the code's OWN assumption that the k valid years are CONSECUTIVE,
    which is what Sxx = k(k^2-1)/12 encodes - and the case that motivated
    both ToE blockers is the opposite: 8 valid years SCATTERED across a
    nominal 32-year window. Recomputed at nominal=32, snr=2.0 exactly, using
    the real Sxx of the actual year positions (t = |r|*sqrt(df/(1-r^2)) with
    r^2 = 4*Sxx/(k*nominal^2), against this file's own tCrit95 table):
      k     Sxx consec   Sxx scattered   ratio   t consec  t scattered  tCrit
       4         5.0           533.9     106.8      0.099       1.476   4.303
       8        42.0           823.7      19.6      0.354       2.009   2.447
      10        82.5           978.8      11.9      0.516       2.225   2.306
      12       143.0          1135.7       7.9      0.699       2.422   2.228
      16       340.0          1452.2       4.3      1.126       2.773   2.145
      24      1150.0          2089.1       1.8      2.251       3.367   2.074
    ("scattered" = k years spread evenly across the 32; with k years drawn
    at random from the 32 the median Sxx sits between the two columns.)
    So Sxx really is 10-20x larger at small k, t is correspondingly larger,
    and the gate stops binding at k=11 on the evenly-scattered pattern and
    k=14 on the random-draw median - NOT at k=24. The statement that holds
    under every pattern tested is: THE GATE BINDS AT k<=8. Above that it
    depends on how the valid years are distributed, and above k~14 it is
    vacuous again under all of them. This residual is inherent to
    defining snr from the TOTAL SD on a complete evenly-spaced series -
    there are only three numbers (slope, SD, n) and any statistic built
    from them is a function of the other two. Removing it entirely would
    mean redefining snr, which is deliberately NOT done here. So: for a
    100%-coverage SST/Chl/Salinity/DO pixel the gate still adds nothing;
    for a sparse pixel - the case that motivated both blockers - it binds.

  REGRESSION (from v10.156 BUG-03) - a CONTROL-site variance artifact
    suppressed the study's AC1 signal in STEP 4 FIND SWEET SPOT. v10.156
    collapsed FIVE checks - study AND control - into one rowVarArtifact
    flag and then used that flag to gate AC1 as well, so a near-zero
    first-half variance at the open-ocean CONTROL discarded a genuine study
    dAC1 of e.g. +0.40 from bestAc1W, localLeaningRows, testedRows, the
    ROBUSTNESS tally and the headline - and the tool then reported "NO
    WINDOW COULD SUPPORT A VERDICT" on a real signal. STEP 3 COMPARE had
    always done this correctly (studyVarArtifact and ctrlVarArtifact
    separate, only the corresponding variance term nulled, AC1 never
    touched). STEP 4 now mirrors STEP 3 exactly: rowStudyVarArtifact nulls
    only sDVarScored, rowCtrlVarArtifact nulls only cDVarScored, divScored
    requires BOTH sides intact, and neither touches sDAC1 or ac1Rose. A
    variance artifact still excludes the row from all VARIANCE-based claims
    and from the Scheffer check (which needs BOTH indicators) - that part
    was right. rowUnderpowered is UNCHANGED and still excludes the row from
    everything, AC1 included. Every downstream consumer was re-checked:
    localLeaningRows and testedRows (variance verdicts) now use varExcluded;
    a new AC1 tally and a new AC1-only headline branch use underpowered
    only; the flagged-window list distinguishes "VARIANCE ARTIFACT ...
    dAC1 STILL USABLE" from "UNDERPOWERED ... NOTHING usable". The
    permutation local-signal count never consumed rowExcluded (it tests
    w < CSD_MIN_WINDOW_MONTHS directly) and so was already correct.
    Unit-tested against the LITERAL shipped block from both versions:
    control-only artifact - AC1 kept, study variance kept, control variance
    nulled; study-only artifact - the mirror image; underpowered - still
    excludes everything in both versions; a clean powered row - bit-
    identical output in both versions.

  CHANGELOG HONESTY - the v10.156 BUG-02 entry claimed "across all 126,864
    untied orderings at n=4..10". That figure matched no enumeration; the
    sum of n! for n=4..10 is 4,037,904. Re-derived here by exact
    enumeration: over all 4,037,904 untied orderings, exactly 1,234
    verdicts change and every one LOSES significance (2 of 24 at n=4, 28 of
    720 at n=6, 1,204 of 40,320 at n=8; none at n=5, 7, 9 or 10). The 1,234
    and the directional claim were both correct - only the total was wrong.
    Corrected in the header entry and in the startup banner, with the error
    disclosed rather than silently overwritten. The other specific numbers
    in that entry were re-checked and DO hold: n=4 tau=+1 gives exactly
    0.083333 exact vs 0.041540 old; n=7 tau=+1 gives 0.000397 vs 0.001611;
    and mkExactTailP matches brute-force enumeration of all n! orderings
    for n=4..8 with max abs error 0.0e+0.

  RESIDUAL RISK - what could NOT be verified without a live Earth Engine
    session, and is therefore a GUESS until someone opens this in the
    browser:
    (a) ee.Reducer.count() band naming. The new count images are reduced
        with .reduce(ee.Reducer.count()) on a single selected band, which
        should yield "<band>_count" by the same convention that already
        makes the working stdDev bands "sst_stdDev" etc. NOT confirmed.
    (b) ee.Reducer.pearsonsCorrelation() band naming. Expected to be
        "correlation" (plus a "p-value" band), unprefixed, by the same
        convention that already makes linearFit's outputs plain "scale"
        and "offset" in this file. NOT confirmed - and v10.88 recorded
        exactly this uncertainty for the same reducer.
    MITIGATION: because (a) and (b) are guesses, NEITHER new entry does
    .select() or .get() on a guessed band name. Both pass the WHOLE
    reduceRegion dictionary through and the value is extracted CLIENT-side
    by toeNum(), which tries the expected keys, then any other numeric key
    (skipping p-value), then returns null. A wrong band name therefore
    degrades to "count unavailable -> nominal fallback, stated on screen"
    or "no measured r -> derived fallback, stated on screen", instead of
    throwing the hard server-side error that would blank all four core
    indicators. That is the same defensive pattern v10.88 adopted for
    computeSpatialEWS/extractSpatialAC1Detail.
    (c) Whether ee.Reducer.count(), stdDev, linearFit and
        pearsonsCorrelation all see the IDENTICAL set of valid pixels at a
        point. The code assumes they do. If they do not, the existing
        "reducers disagree" unreliable branch catches the derived path, and
        the measured-r path would use a slightly wrong n in (n-2).
    (d) Whether EE's stdDev is the population or sample SD. The file
        asserts population; if it is the sample form, SStot is overstated,
        the derived t is understated (conservative), and the r-implied
        algebra above shifts by under 2% at n=44.
    (e) None of the STEP 4 or ToE display strings have been rendered in a
        real ui.Panel - only their construction was exercised in Node.


v10.156 FIX 14: eight statistical bugs found by an external, unit-tested
  audit. Every fix below was reproduced BEFORE and corrected AFTER in
  standalone Node harnesses run against the extracted functions.

  BUG-01 (CRITICAL) - the regime-shift index was INVERTED.
    computeScore() fed ccs (the Coastal Cancer Score = STRESS, higher is
    worse) straight into mu, the RESILIENCE parameter of the Waddington
    double well U(q;mu)=0.25q^4-0.5*mu*q^2 whose barrier is mu^2/4 - so a
    larger mu is a DEEPER, more stable well and stress was being read as
    stability. Swept against real inputs the old code returned 81% at
    CCS=24, 74% at CCS=46 and 66% at CCS=56: monotonically backwards.
    mu is now (100-ccs)/100, the same quantity as bowl depth B. sig (the
    noise intensity) legitimately rises with stress and is UNCHANGED at
    0.04+ccs/2000. That alone still turned over above ccs~76 because the
    attempt frequency w0=sqrt(mu/me) collapses as the well flattens, so
    w0's curvature is floored at the barrier<=noise crossover mu=2*sqrt(sig),
    where Kramers escape stops applying and diffusion-limited escape
    saturates instead of vanishing. Verified: p5yr is now monotone
    non-decreasing across the FULL 0-100 CCS range (1% -> 97%), and across
    the SST 24->34 / DHW 0->20 sweep it reads 22% at CCS=24, 68% at 48 and
    79% at 55, against the old 81%/72%/67%.

  BUG-02 (HIGH) - mannKendallTest() reported impossible p-values at small n.
    It used only the normal approximation. At n=4 with tau=+1 it returned
    p=0.0416, but with 4!=24 orderings the smallest ATTAINABLE two-sided p
    at n=4 is 2/24=0.0833 - it claimed significance where significance
    cannot exist. Now computes the EXACT null distribution of S for n<=10
    with no ties, via the standard Mahonian (inversion-count) recursion,
    and adds the missing TIE CORRECTION plus the standard continuity
    correction to the normal branch it falls back to. Validated against
    brute-force enumeration of all n! orderings for n=4..8: exact match,
    max abs error 0.0e+0. n=4 tau=+1 now returns 0.0833 (was 0.0416) and
    n=7 tau=+1 returns 0.000397 (was 0.0016). (v10.157 CORRECTION: this
    entry originally said "across all 126,864 untied orderings at n=4..10".
    That figure was wrong and matched no enumeration - the sum of n! for
    n=4..10 is 4,037,904. Re-derived by exact enumeration in v10.157: over
    all 4,037,904 untied orderings at n=4..10, exactly 1,234 significance
    verdicts change and every one LOSES significance - none gains it. The
    1,234 and the directional claim were both correct; only the total was
    not. Per-n: 2 of 24 at n=4, 28 of 720 at n=6, 1,204 of 40,320 at n=8,
    and none at n=5, 7, 9 or 10.) Return shape preserved; method,
    nTieGroups added, and every on-screen verdict now says which null
    produced its p-value.

  BUG-03 (HIGH) - S13 had no near-zero-denominator guard.
    varTrendRatio = secondHalfVar / max(firstHalfVar, 1e-6), so a near-zero
    first half inflates the ratio without any genuine surge - real runs
    produced +1276%, +2876% and +34.22x at a 6-month window, and that fed a
    false "LOCAL CSD SIGNAL DETECTED" banner. S7C (v10.105) and S7D
    (v10.109) already had the guard; S13 never got it. The SAME rule is now
    ported to STEP 2, STEP 3 COMPARE, STEP 4 FIND SWEET SPOT and the S12
    sidebar: firstHalfVar < 0.001 AND |ratio or delta| > 5 -> ARTIFACT,
    warned on screen and EXCLUDED from the verdict, the divergence test,
    the toolkit tally, the sweet-spot pick and the robustness tally.
    Unit-tested: all three real cases flag, genuine surges on real variance
    (0.25 -> 2.0, an 8x rise) do not.

  BUG-04 (LOW) - version markers disagreed. The header said v10.154, the
    sidebar v10.155 and the per-click console banner still said v10.67.
    All four self-identification markers now say v10.156. Historical
    changelog entries describing what v10.154/v10.155 changed are left
    alone - they are history, not self-identification.

  BUG-05 (MEDIUM) - computeScore() returned a score from no data.
    Every sub-score falls back to a hardcoded mid default when its input is
    null (s1=50, s2=50, s3=35, s4b=30, s6=20), so with all six satellite
    inputs null it returned a confident CCS=30/100 plus a full Bowl Depth /
    omega0 / tau / regime-shift readout: "not measured" was indistinguish-
    able from "measured and benign". Now matches computeAquaculture()'s
    existing contract - ccs:null, insufficientData:true, an explicit
    human-readable dataNote - and EVERY return carries nInputs and
    dataCompleteness (the weighted share of the 15/15/15/25/20/10 composite
    backed by a real measurement), surfaced in a new sidebar row and two new
    CSV columns, with lowConfidence raised below 50%. THRESHOLD: the hard
    gate is nInputs===0, the only case where 100% of the composite is
    default and "not a score" is unarguable. scoreColors() gained a null
    branch (without it a null CCS fell through every band and painted
    itself CRITICAL red), and all 20+ downstream displays are guarded.

  BUG-06 (MEDIUM) - ToE SNR had no sample-size penalty.
    snr = |slope*recordYears| / residualSD carries no information about how
    well the slope itself is pinned down, so a 4-point record could be
    declared EMERGED as readily as a 44-point one, and S17 displayed that
    with full visual weight. Fixed in BOTH ToE paths (computeToESignal and
    the calcToE path S17 actually renders) by dividing by the STANDARD
    ERROR OF THE SLOPE, se = residualSD / sqrt(Sxx), and comparing against
    t_crit(df=n-2) rather than a flat 2.0. For the annual series S17 uses,
    Sxx = n(n^2-1)/12 and the residual follows from the fitted slope and
    the annual SD already fetched - no new Earth Engine call. EMERGED now
    requires BOTH the original amplitude criterion (snr>=2) AND a slope
    distinguishable from zero at its own df, and every verdict shows n, df,
    t and the 95% bar inline. Unit-tested at equal SNR=3.0: the 4-year
    record now reads "not yet" (t=2.18 vs bar 4.30) while 7/12/27/44-year
    records still EMERGE - the old rule called all five EMERGED.

  BUG-07 (MEDIUM) - S13 STEP 4/5 climatology was not pooled.
    COMPARE (STEP 3 / runControlCSD) has recomputed both windows against
    ONE pooled climatology since v10.151 FIX 4b; FIND SWEET SPOT never did,
    so the two panels measured different quantities - a real Bocas del Toro
    run had COMPARE at dAC1=+0.524 and FIND SWEET SPOT at +0.198 for the
    SAME site and SAME window. STEP 4 now threads the identical pooled
    recompute (computeMonthlyClimatology + jsNodeStatsFixed(series, clim))
    using the raw series the permutation test already fetches - ZERO extra
    EE calls - and falls back to the per-window figures with an explicit
    note if those series do not arrive. Verified: the STEP 4 pooled path is
    bit-identical to COMPARE's on the same inputs, and on synthetic data
    per-window vs pooled differ by up to 0.46 in dAC1, the same class of
    discrepancy the Bocas run showed. STEP 5's OVERLAPPING trajectory is
    computed server-side by computeSlidingWindowCSD(), which globally
    detrends but never deseasonalizes; pooling it needs an ee.Array
    group-by-calendar-month restructure that cannot be verified without a
    live GEE run, so instead BOTH outputs now state explicitly that they
    measure different quantities and must not be compared. STEP 5's
    independent-window test already pooled (v10.129) and is comparable.

  BUG-08 (HIGH, design) - power is very low at the window lengths the UI
    encouraged. Monte Carlo power of permutationTestDeltaFixed (BEFORE
    phi=0.2 vs AFTER phi): 24mo = 5%/19%/28% at phi=0.5/0.7/0.9; 48mo =
    23%/63%/87%. At 24 months the test has essentially no power, yet STEP 2
    accepted 4-month windows and FIND SWEET SPOT tested 6, 9 and 12-month
    windows and could name one of them the "sweet spot". The estimator is
    UNCHANGED; the guard rails are new. STEP 2 now REFUSES anything below
    24 months outright and labels 24-47 months UNDERPOWERED. FIND SWEET
    SPOT's window list changes from [6,9,12,15,18,24] to [6,9,12,24,36,48]
    - same count, so the EE call budget and the Bonferroni divisor are
    unchanged, but 15 and 18 (underpowered and redundant) are replaced by
    the first lengths that can actually detect a real change. 6, 9 and 12
    are kept for diagnostics, labelled UNDERPOWERED, and excluded from the
    sweet-spot pick, the ROBUSTNESS tally and the significance counts. The
    power table itself is printed on screen in STEP 2, STEP 3 and STEP 4 so
    the user chooses a window against the real numbers.

  NOT FIXED, deliberately: the Bowl Depth composite's 15/15/15/25/20/10
  weights and the invented effective-mass term (me = anem_N*0.2 +
  urchin_N*0.05) remain uncalibrated - BUG-01 corrects the SIGN and the
  monotonicity of the escape rate, it does not turn the index into a
  calibrated probability. The panel's own v10.143/v10.155 disclosures on
  that point still stand and are unchanged.

v10.149 NEW: S12b - Real Significance Test, converts S12 from
  heuristic to probabilistic
  Direct answer to a real question: S12 is real math (detrended AC1)
  but has zero statistical testing attached - a real number with a
  heuristic interpretation layered on top, not a genuine significance
  test. Added a real permutation/surrogate-data test: shuffles THIS
  SAME window's 24 months into random order 500 times, asking "would
  randomly-ordered data produce an AC1 this high by chance alone?" -
  a genuinely different question from S13's before/after delta test.
  This is a real, established method in nonlinear time-series
  analysis (surrogate-data testing for genuine autocorrelation).
  Built standalone (button-triggered) since S12's own computation
  runs entirely server-side in EE Arrays and never sends raw monthly
  values to the client - needed one new, small EE call, reusing
  mkMoSST() (the EXACT same Jan2023-Dec2024 collection S12 itself
  uses) and the same jsNodeStats() formula already proven throughout
  this tool, so the real AC1 shown matches S12's own number exactly.
  Verified before shipping with two direct tests: structured
  (seasonal-cycle-like) data correctly comes back highly significant
  (p=0.0000); pure random noise correctly comes back not significant
  (p=0.76).
  HONEST DISCLOSURE built directly into the tool: a significant
  result here does NOT mean AC1 has changed or a reef is destabilizing
  - real ocean temperature data has genuine seasonal structure, which
  this test will correctly flag as "significant" every time, since
  scrambling destroys that real seasonal pattern. This tests whether
  observed memory is distinguishable from noise, not whether resilience
  is being lost - that separate, more specific question is S13's job.

v10.148 FIX: S7 no longer shows coral-reef alarm language at
  cold-water kelp sites
  Direct response to a real, correct catch: a real Nuuk test showed
  S7 displaying "EXTREME - massive bloom, RED" - alarm language built
  for tropical coral-algae phase shifts - at a site S8 (and the Field
  Data section) independently identify as likely HEALTHY cold-water
  kelp habitat, not a warning. Same class of gap already fixed for
  DHW (explicitly hidden above 55 deg N/S as ecologically meaningless
  there) - S7 never got that same climate-zone awareness until now.
  FIXED: moved computeAquaculture() (which already detects cold-water
  kelp opportunities) to run BEFORE S7's algae display instead of
  after, so its kelp-detection result can contextualize S7's text.
  When kelp is detected AND the algae signal is genuinely elevated,
  S7 now shows "LIKELY COLD-WATER KELP... not a coral-algae bloom
  warning" instead of the old tropical-framed alarm text.
  Verified with real values from both directions before shipping:
  the real Nuuk test values (FAI=0.117) correctly now show the kelp-
  context message; the real Florida Keys test values (FAI=0.0796,
  kelp NOT detected there) correctly still show the normal "CRITICAL
  - dense mat" alarm, unchanged - confirms the fix only applies where
  it should, not a blanket suppression of real algae warnings.

v10.147 BUGFIX: S18 STILL broken after v10.146 - found the REAL root
  cause via a second real test
  v10.146's "empty collection" fix didn't solve it - a real test at
  One Tree Reef (36m deep, open lagoon, nowhere near a coastline)
  FALSIFIED the v10.146 hypothesis that this was about Florida Keys'
  coastal proximity - it failed identically there too. Compared
  directly against S17's OWN pH fetch, which succeeded at BOTH real
  test sites using the exact same real CAR asset, and found two real
  differences: (1) S17 uses a bare ee.Geometry.Point with NO buffer;
  S18 was using `study`, a buffered polygon - against a coarse ~25km
  global grid, a small buffered area can produce a genuinely empty
  reduceRegion result in edge cases where a bare point cleanly
  resolves to one pixel. (2) S18 used a narrow 90-day window vs
  S17's multi-year annual collection. Fixed both: switched to the
  identical bare-point geometry pattern already proven at both real
  test sites, and widened the window to 2 years (730 days).
  Honest note: this is now the tool's best real attempt at fixing
  S18, built by directly comparing against proven-working code rather
  than guessing - but it has NOT yet been re-tested live. If it still
  fails, the new diagnostic message explicitly says so rather than
  reusing the old generic "no valid pixel" text.

v10.146 BUGFIX: S18's v10.145 fix was itself broken - found via real
  testing at One Tree Reef
  A real test showed S18's new pH/pCO2/Salinity all failing with
  "error fetching real BGC assets", even though S17's pH (same real
  CAR asset) succeeded in the same test. Root cause found by
  comparing against this codebase's own established pattern: every
  other real-asset fetch (mkMoDHW, dhwProper, etc.) guards against an
  empty ImageCollection with ee.Algorithms.If(col.size().gt(0), ...,
  null) before calling .select()/.mean() on it - the v10.145 S18
  rebuild skipped that guard. If the 60-day window genuinely returned
  zero images at some location/time, that produced a real computation
  error instead of a graceful null. Fixed using the identical,
  already-proven defensive pattern; also widened the window 60->90
  days to further reduce the chance of a genuinely empty result.

v10.145 MULTI-FIX: S18 real assets, SNR-vs-Mann-Kendall diagnosis,
  S16 ECI-vs-B honesty fix, S20e wired into Bowl Depth
  1) S18 FIXED for pH + pCO2 + Salinity using the same real, confirmed
  Copernicus/HYCOM assets already used for S17/S20e - split into 3
  independent fetches so one still-dead field (O2) can never again
  blank the 3 real ones (same class of fix as v10.89/v10.141).
  pCO2 correctly converted from the real asset's native Pa to uatm
  (x9.86923), sanity-checked against realistic ocean pCO2 ranges.
  2) DO/O2: searched a THIRD time this session, still no confirmed
  working sub-collection/band name - disclosed honestly, not guessed.
  3) DIAGNOSED why SNR and S17b's Mann-Kendall sometimes disagree
  (e.g. NO2/Salinity "EMERGED" by SNR but not significant in S17b):
  the SNR formula multiplies slope by record_years with NO correction
  for how uncertain a slope estimate is at few data points - unlike a
  real p-value, it structurally over-triggers on short records. Now
  disclosed directly in S17's UI. Also found and disclosed a real,
  separate finding: Salinity's nominal 32yr record often has far
  fewer REAL valid HYCOM points at a given site (genuine data
  sparsity, not a bug) - now flagged in S17b's own output.
  4) FIXED S16's ECI-vs-B "agreement" check - investigated and found
  it was comparing a PURE physics metric (depth-only wave energy)
  against the broad ecological Fused score (SST/Chl/DHW/biology/NO2).
  These measure different things with no reason to agree - "DISAGREE"
  was implying an error that didn't exist. Relabeled honestly.
  5) S20e NOW WIRED into Bowl Depth - fulfills the v10.143 disclosed
  promise. New row shows the uncalibrated heuristic "Regime-shift
  index" directly next to S20e's real, held-out-validated P(bleaching)
  - verified with real Bocas del Toro values (78% vs 94.8%, a real
  16.8-point gap surfaced).
  NOT built this version (real scope decision, disclosed rather than
  rushed): S20e's 6-feature extension (distance-to-shore, windspeed,
  cyclone frequency) needs 3 new real global data sources not yet
  identified in GEE - a genuine, separate next step.

v10.144 NEW: S20e - Real Bleaching Probability, a genuinely fitted
  logistic regression model
  Direct build-out of the real database found and downloaded this
  session: Global Coral-Bleaching Database (van Woesik & Kratochwill
  2022) - 41,361 raw rows, cleaned to 32,716 real, complete rows
  (dropped rows with "nd" sentinel missing-data values, matching
  BCO-DMO's own documented convention). Fit a REAL logistic
  regression on 3 features that map exactly onto this tool's own
  already-computed metrics: DHW (SSTA_DHW), Turbidity, Depth.
  REAL, HONEST VALIDATION - not training-set accuracy: 70/30 train/
  test split, evaluated on a genuinely held-out 9,815-row test set
  never seen during fitting. Held-out AUC=0.620 - real, better than
  chance (0.5), but genuinely WEAK - disclosed directly in the UI,
  not hidden. Consistent with Arias-Ortiz et al. 2024's independent
  finding that DHW-style models need ~23 metrics to predict well.
  TWO REAL, INDEPENDENTLY-CONFIRMED FINDINGS from the fitted
  coefficients: (1) Depth's coefficient is POSITIVE - deeper reefs
  more likely to bleach in this real 32,716-row fit, independently
  confirming Arias-Ortiz et al.'s counter-intuitive finding already
  disclosed in S20. (2) Turbidity's coefficient is negative, matching
  a real, separate citation (Sully & Woesik 2020, "Turbid reefs
  moderate coral bleaching").
  Zero new EE calls - reuses dhwv/turv/bv, already computed for the
  main sidebar. Verified with a manual-formula-vs-sklearn sanity
  check (exact match) before porting to JS.
  Predicts BLEACHING probability, explicitly NOT collapse - see
  v10.143's terminology fix for why that distinction matters.

v10.143 FIX: terminology corrected - "collapse" language reviewed
  throughout, fixed where it overclaimed
  Direct response to a real question about a possible Bayesian/
  logistic upgrade to the Fused Score, sourced from a real, large-
  scale database (Global Coral-Bleaching Database, van Woesik &
  Kratochwill 2022 - 34,846 records, 14,405 real sites, confirmed via
  direct search; data hosted on Figshare, not yet pulled into this
  tool - would need the actual file, not just its existence, to fit a
  real model). Since any such model would honestly predict BLEACHING
  probability (what the database records), not COLLAPSE, reviewed
  every "collapse" reference in this tool for accuracy:
  FIXED (real overclaiming): "P(flip 5yr)" relabeled "Regime-shift
  index (5yr, uncalibrated)" - this heuristic was never a calibrated
  probability of anything; a disclosed caveat now points to the real
  GCBD-based upgrade path. Also fixed STEP 2's "AFTER event (post-
  collapse)" dropdown label to "AFTER event (suspect period)" - most
  AFTER windows tested this session showed NO signal at all, so
  presuming verified collapse by default was misleading.
  LEFT UNCHANGED (genuinely accurate uses, confirmed real): Levitan &
  Edmunds' real urchin population collapse citation; Byrne et al.'s
  real, field-observed coral colony collapse at One Tree Reef; S13's
  description of Scheffer et al.'s own critical-slowing-down theory,
  which is genuinely about approaching collapse/regime shifts.

v10.142 NEW: S17b - Real Mann-Kendall Significance Test
  Direct answer to a real question: since the Mann-Kendall engine
  already exists (built and fixed for STEP 5), can S17's fixed
  SNR>=2.0 heuristic be replaced with a genuine significance test?
  Built as a standalone, button-triggered tool rather than rewriting
  S17's core flow directly - that flow was already touched once this
  session (the pH fix) and is deeply async/fragile; adding a new tool
  alongside it is the same lower-risk pattern already used for
  COMPARE's statistically-valid banner (v10.139). Reuses the SAME
  real annual-value collections S17 already built (_annSSTColl etc.)
  and the SAME mannKendallTest() engine already validated for STEP 5 -
  no new architecture, genuinely real math throughout.
  Verified directly before shipping: correctly flags short records
  (n<4, e.g. a hypothetical 3-point NO2 case) as insufficient rather
  than guessing; correctly computes a real, strongly significant
  result on a longer record with an embedded trend (SST-style, 44
  points, p<0.0001).
  DISCLOSED: pH (~4yr) and NO2 (7yr) have very few real annual points
  - Kendall's tau needs real data, so these will often show
  "insufficient" here even though S17's SNR-based check above can
  still report a number for them - a real, honest limitation of
  annual-resolution testing on short records, not a bug.

v10.141 FIX: S17 pH now uses a REAL, confirmed, working Copernicus
  asset - DO still unresolved
  Direct answer to a real question: can a working pH/DO database be
  integrated? Searched and confirmed: the old dead asset ID was simply
  outdated - the real, current GEE catalog splits this Copernicus
  product into per-variable sub-collections. Confirmed directly via
  Google's own catalog pages: COPERNICUS/MARINE/
  GLOBAL_ANALYSISFORECAST_BGC_001_028/CAR, band ph_depth1 - real,
  working, current surface pH.
  IMPORTANT CAUGHT-BEFORE-SHIPPING bug: this real asset's data only
  starts 2021-10-01, not 1993 as the old 32-year assumption required.
  Verified directly: using the old 32yr constant against the real
  ~4yr record would have overstated SNR by 8x - easily turning a
  genuinely non-significant trend into a false "EMERGED" result.
  Fixed: annual pH stack rescoped to 2022-2025, record-years constant
  corrected 32->4, confidence relabeled 32yr/MARGINAL -> ~4yr/LOW.
  ALSO FIXED: pH and DO were still bundled into one evaluate() call -
  since DO's asset remains dead, that would have made DO's failure
  blank the now-working pH too (same class of bug already fixed for
  SST/Chl/Salinity/NO2 vs BGC in v10.89). Split into two fully
  independent evaluates with two independent availability flags.
  DISCLOSED: dissolved oxygen (DO) is NOT fixed this version - no
  working sub-collection/band name was confirmed for O2 after two
  direct searches - left disclosed as unavailable rather than guessed.

v10.140 NEW: real Mesoamerican Reef entry + real methodological
  caveat from a published predictive-model study
  Direct response to a real slide shown (Global Tipping Points Report
  2025, healthyreefs.org): confirmed the Healthy Reefs Initiative is
  real (70+ partner orgs, monitoring since ~2008) - report-based, not
  a live API, same access pattern as this tool's other citations.
  HONEST CORRECTION made in the process: nearly merged this regional
  data into the existing Bocas del Toro (Panama) entry - caught that
  they are genuinely DIFFERENT reef systems (Bocas sits on Panama's
  southern Caribbean coast; the Mesoamerican Reef spans Mexico/Belize/
  Guatemala/Honduras further north) - added as a real, separate 7th
  site instead, confirmed non-overlapping with Bocas del Toro
  (1042km apart). Real 2023 stats: ~40% of corals severely affected,
  regional cover fell 19%->17% (GCRMN Mesoamerican Report Card 2024).
  Also found and disclosed a real, relevant finding: Arias-Ortiz et
  al. 2024 (Communications Biology, doi:10.1038/s42003-024-07128-y)
  found a combined 23-metric model explains 75% of bleaching-severity
  variance, versus DHW alone explaining less - and, surprisingly,
  DEEPER reefs with MORE diverse coral communities showed HIGHER
  vulnerability. Added as a disclosed caveat on S20's simple DHW>=4
  threshold - not a full model port (the paper's actual fitted
  coefficients were not available from what was verified), but an
  honest flag that a more sophisticated real model exists.

v10.139 NEW: COMPARE now shows a statistically-valid PRIMARY verdict
  Direct answer to a real, fair question: since a genuine 500-shuffle
  permutation test already exists, why was the >0.01 AC1/>0.15 var
  HEURISTIC threshold still the primary displayed verdict? Real
  reason found: the permutation test fires its OWN separate EE call,
  resolving LATER than the threshold verdict - a sequencing fact, not
  a design choice. Rather than risk rewriting the large, deeply-nested
  async chain that builds the threshold verdict, added a NEW banner
  positioned ABOVE it - populated by the permutation test's own
  completion callback once real p-values are ready, so the
  statistically valid answer is now the first, most prominent thing
  shown. Same AC1-primary weighting philosophy already used
  throughout this tool (Dakos et al. 2012), now driven by p<0.05
  significance instead of the old fixed thresholds. The original
  threshold-based box remains below, explicitly relabeled as the
  secondary, faster-arriving heuristic check.
  Verified against 3 real sites already tested this session (One Tree
  Reef, Bocas del Toro, Nuuk) - all three correctly resolve to "NO
  SIGNIFICANT SIGNAL," matching what was already independently
  established as the rigorous answer at each.

v10.138 CONSOLIDATION: merged S20/S20b/S20c into ONE unified tool
  Direct response to real, fair feedback: three separate buttons for
  what's conceptually one question ("what's the combined risk, at some
  point in time?") was confusing design, not a missing capability -
  both past-date sync and live-present already worked correctly, just
  split across separate tools. Now ONE date field: leave BLANK for
  LIVE (today, real current OISST data), or type a month (YYYY-MM,
  2023-01 to 2024-12) for a HISTORICAL replay of that exact month -
  one button, one result, clearly labeled which mode ran. The
  original always-on S20 (auto-computed on every click, frozen to the
  2023-2024 peak) is retired as a separate behavior - S20 is now
  button-triggered like the rest of this family, so the click-handler
  code that used to auto-populate it on every map click was removed.

v10.137 NEW: S20d - Genus Growth-Form Lookup (Coral Trait Database)
  Direct answer to a real question: are there real, established
  global species databases available? Searched and confirmed one -
  the Coral Trait Database (coraltraits.org, Madin et al. 2016,
  Scientific Data), a real, actively-growing, peer-reviewed database
  (166,245 observations, 5,112 species as of the check). Confirmed a
  real, working direct-CSV access pattern for individual traits.
  DISCLOSED: could not directly fetch the live CSV from inside this
  GEE script - outbound access to external domains from the Code
  Editor sandbox is unverified, unlike the separate companion Python
  scripts which run outside GEE entirely. Built the honest version
  instead: GENUS_GROWTH_FORM, a small, real starter table, with each
  entry explicitly labeled either [CTD-confirmed] (individually
  verified against a live Coral Trait Database entry during this
  build) or [literature pattern] (the same well-established branching-
  vs-massive framework, Loya et al. 2001, not individually cross-
  checked for that exact genus). New S20d UI lets you type any genus
  name and see its real growth form + general vulnerability tier -
  a GENERAL, global fallback, distinct from S20's specific, real
  event-outcome mortality data at 6 named sites.

v10.136 NEW: S20c - Live Current DHW Check
  Direct answer to a real question: is live (not frozen-2023-2024)
  DHW possible? Checked directly rather than assumed - confirmed via
  Google's own Earth Engine catalog page that NOAA/CDR/OISST/V2_1
  (the exact dataset getOISSTColl() already uses everywhere in this
  tool) is real-time, updated daily, with data through essentially
  today (1-day preliminary lag, 14-day final). This was never a data
  availability limit - only the hardcoded 2023-2024 window prevented
  it. New button computes a genuine rolling recent-window DHW ending
  at today's actual date, reusing the identical monthly-max*4.33
  formula already used throughout this tool - same math, same data
  source, just pointed at "now" instead of a fixed past window.
  DISCLOSED, and important: this is an OBSERVATION of the present, not
  a forecast of the future. Genuine forecasting would need real ocean/
  climate forecast MODEL output (e.g. NOAA Coral Reef Watch's actual
  bleaching outlook) - a different kind of data source, not verified
  as available in this tool's GEE architecture, and NOT built here.

v10.135 NEW: S20b - Historical Month Check (retrospective test)
  Direct build-out of a real request: test S20's combined risk
  retrospectively, as it would have shown in early 2024. Uncovered a
  real architectural fact first: the main click flow's DHW (and most
  S1-S11 measurements) was ALREADY hardcoded to a fixed 2023-01-01 to
  2024-12-31 window this whole session, not a live "as of today"
  calculation - so the DHW shown was already the PEAK across that
  entire window, not a present-day reading. There was no "current"
  value to roll back from.
  Built the honest version instead: a real month-by-month replay
  within that same real, already-available 2023-2024 data window.
  Pick any month (YYYY-MM, 2023-01 to 2024-12), and it computes the
  REAL DHW for just that month at the last-clicked location, then runs
  it through the exact same combineSpeciesAndDHW() logic S20 uses -
  reusing getOISSTColl() and MMM_perpixel, both already proven
  elsewhere in the tool, same math as mkMoDHW(), just isolated to one
  chosen month instead of the full 24-month collection.
  DISCLOSED SCOPE: this is NOT a true arbitrary-date time machine -
  months outside 2023-2024 correctly return an explicit range error,
  not a guess, since no real data exists there in this tool's current
  architecture.

v10.134 NEW: 3 real entries added proactively (Bocas del Toro,
  Andaman/Nicobar, Brazil)
  Direct response to real feedback: proactively check sites already
  tested this session and other named regions, rather than waiting to
  be told each one individually. Searched and verified real literature
  for all three:
  1. Bocas del Toro, Panama - the exact site tested with S13/CSD this
     session. Neal et al. 2017 (Ecology and Evolution, doi:10.1002/
     ece3.2706) tracked 3 massive coral genera through the real 2005
     bleaching event for 8 years - all showed continued net tissue
     loss with no full recovery, even typically stress-tolerant forms.
     DISCLOSED: a second Bocas-specific source was found but rejected
     - its search snippet mixed content from an apparently unrelated
     Indian coral-management paper mid-sentence, a reliability red
     flag, so only the clean, verified citation was used.
  2. Andaman & Nicobar Islands, India - real 2010 mass-bleaching
     survey with clean genus-level mortality: Acropora 43%, Montipora
     22%, Porites 14%.
  3. Southwestern Atlantic Reefs, Brazil - Mussismilia harttii, Brazil's
     major endemic reef-builder, described as suffering "unprecedented
     erosion" (Braz et al. 2022, Coral Reefs) with no recovery three
     years after a major bleaching event (Corazza et al. 2024, Mar.
     Biol.). DISCLOSED: no single mortality percentage was verified
     for this exact event - reported as the real qualitative finding
     rather than an invented number.

v10.133 NEW: real third species-vulnerability entry (Maldives)
  Direct answer to a real question: can this be expanded to Maldives,
  Philippines? Searched and verified real, current literature for
  both. Maldives had a strong, clean match: "Outcomes of the fourth
  global coral bleaching (2023-2024) in the Maldives," Coral Reefs
  (Springer Nature), 2026, doi:10.1007/s00338-026-02850-x - a real
  survey of 18 reefs across central and southern atolls. Central
  atolls averaged >40% live coral cover loss, up to 57% at the
  worst-hit reef, disproportionately affecting Acropora; the
  southernmost Huvadhoo Atoll showed high resistance, with Porites-
  dominated reef flats retaining high live coral cover.
  DISCLOSED DIFFERENCE from the two existing entries: this is REGIONAL
  data (18 reefs, central/southern atolls) not single-reef data like
  One Tree Reef/Florida Keys - larger radius (150km), and the numbers
  are described as regional patterns, not one precise reef figure.
  Philippines was searched too but NOT added this version - the real
  papers found don't cleanly fit the database's format (one reports
  bleaching PREVALENCE, not mortality; another is not peer-reviewed
  and lacks a clean genus-level breakdown) - disclosed here rather
  than forcing in weaker-quality data just to have an entry.

v10.132 NEW: shared species-vulnerability database + real second entry
  (Florida Keys)
  Direct answer to a real question: can species/genus vulnerability
  data be pulled from literature and made available across the whole
  tool, not just S20? Refactored the single-site lookup into
  SPECIES_VULNERABILITY_SITES, an array any module can query via
  lookupSpeciesVulnerability(lat, lon) - genuinely shared, not
  S20-specific.
  Added a real, independently-verified SECOND entry after directly
  searching for and confirming current literature: Florida Keys/Dry
  Tortugas. Manzello et al. 2025 (Science, doi:10.1126/
  science.adx7825) confirms 97.8-100% mortality of Acropora palmata
  and A. cervicornis - described in the paper as functional
  extinction of both species at this site. Florida FWC's 2023-24
  CREMP Post-Bleaching Quick Look Report adds real genus-level data
  for Orbicella (extensive mortality, Lower Keys), Siderastrea
  siderea and Montastraea cavernosa (both ~1% decline - comparatively
  resilient).
  DISCLOSED SCOPE, stated plainly and unchanged in spirit from
  v10.131: this is explicitly NOT a comprehensive global database and
  will not claim to become one automatically - it grows one real,
  independently-verified citation at a time. Every site not
  explicitly listed still correctly returns "NOT CHECKED," never a
  guess.

v10.131 NEW: S20 - Species-Vulnerability-Weighted Combined Risk
  Direct build-out of an idea discussed: combine the tool's existing,
  real DHW measurement with real, published dominant-species
  vulnerability data, where it exists - rather than treating heat
  stress and species fragility as two separate, unconnected numbers.
  Reuses the EXACT SAME real citation already embedded in S15 (Byrne
  et al. 2025) - not a new data source, just restructured as
  programmatic genus-level data (Acropora/Goniopora HIGH vulnerability,
  95%/73% documented mortality; Porites/Pocillopora MODERATE-LOW,
  31%/23%) so it can be combined with DHW instead of only shown as
  narrative text.
  DISCLOSED SCOPE, checked directly before building: species
  composition cannot currently be measured by satellite at reliable
  accuracy - current published literature on hyperspectral coral
  species classification reports only ~56-70% accuracy at genus/
  species level (vs ~88-90% for basic coral-vs-other detection), too
  unreliable to use as ground truth. This module is explicitly NOT a
  remote-sensing capability - it's a real field-survey lookup, same
  honest pattern as S15/S5: populated only where a real published
  colony-level survey exists (currently just One Tree Reef), explicit
  "NOT CHECKED" everywhere else, never a filled-in guess.

v10.130 FIX: STEP 5's own input cap silently blocked the fix for its
  own bug
  Direct follow-up to v10.129: recommended fixing a real Panama test
  that showed "insufficient points (n=1)" on the new valid significance
  test by using 96 total months with a 12-month window (8 independent
  windows) - but STEP 5's own validation capped "Total months" at 84,
  which would have silently REJECTED that exact fix with "Total months
  must be 24-84." Raised the cap to 144 (12 years) so the fix
  recommended for this exact problem can actually be entered.
  Also added a proactive warning, checked and shown BEFORE the (slow)
  Earth Engine call fires: if the chosen total-months/window-size
  combination would leave fewer than 4 independent windows (as the
  default 50%-of-total auto-window always does, by construction - it
  only ever gives ~2), the run now says so upfront, instead of only
  finding out via "insufficient points" after the full computation
  completes, the way the real Panama test did.

v10.129 NEW: real significance test + deseasonalizing for STEP 5
  (sliding window trend test)
  Direct port of a real bug found and fixed this session in the
  companion Python deseasonalizing script: STEP 5's Mann-Kendall test
  runs on OVERLAPPING windows (step=1 month) - consecutive positions
  share nearly all their underlying data, which violates Mann-Kendall's
  independence assumption. Confirmed as a real problem, not
  theoretical: pure random noise with zero real trend gave p<0.0001
  ("significant") using overlapping windows, but p=0.109 ("not
  significant") using proper non-overlapping windows.
  Applies directly to this tool's own results: a real Bocas del Toro,
  Panama test this session reported "SIGNIFICANT CSD TREND - both AC1
  and variance rising (p<0.05)" using the (now known to be invalid)
  overlapping-window test - that specific claim needed independent
  verification, which is exactly what this new section provides.
  NEW: fetches the raw monthly series client-side (ONE new EE call,
  reusing the same collection/geometry already used for the existing
  trajectory), builds a real calendar-month climatology, deseasonalizes,
  then runs the valid Mann-Kendall test on non-overlapping windows for
  BOTH raw and deseasonalized data - shown alongside (not replacing)
  the existing overlapping-window trajectory/chart, which remains
  useful for visualizing the shape of the trend even though its own
  p-value should now be read with caution.

v10.128 FIX: deseasonalized comparison gave a factually wrong message,
  and was missing the more informative number
  Caught on the very first real test of v10.127 (One Tree Reef): the
  interpretation text said "deseasonalized AC1 is similar in size to
  the raw AC1" when the actual deseasonalized delta (-0.089) was 7.4x
  LARGER than the raw delta (-0.012), not similar - the original logic
  only ever checked for "much smaller," so any other case silently fell
  into a message that didn't match the numbers.
  Fixed the binary check into three real cases (much smaller / much
  larger / similar), and added something the same real test showed was
  actually more important than the delta comparison: the ABSOLUTE
  baseline AC1 level itself. That test showed raw AC1=0.843 vs
  deseasonalized AC1=0.434 - nearly half - suggesting a large share of
  what looks like "high thermal memory" in raw AC1 readings throughout
  this tool may reflect the seasonal cycle itself, not real resilience
  loss. Now flagged explicitly whenever the deseasonalized baseline is
  meaningfully lower than the raw one, not just when the CHANGE differs.

v10.127 NEW: deseasonalized AC1/variance comparison for STEP 3 COMPARE
  Direct fix for a caveat disclosed since v10.101 and never actually
  built: AC1/variance computed on raw monthly values can be inflated
  purely by the shared seasonal cycle (every site warmer in summer,
  colder in winter), regardless of any real dynamical change. Evidence
  this was a real problem, not theoretical: AC1 landed in a narrow
  ~0.78-0.89 band across nearly every GBR site tested this session -
  more consistent with seasonal autocorrelation dominating than genuine
  site-specific signal.
  New computeMonthlyClimatology()/deseasonalizeSeries() functions build
  each site's own calendar-month average from its raw fetched series,
  subtract it to get anomalies, then run the SAME jsNodeStats() used
  everywhere else on those anomalies instead of raw values - shown
  directly alongside a "raw" recomputation (using the identical
  client-side method, for a clean apples-to-apples comparison) so the
  seasonal cycle's actual contribution to the raw AC1 signal becomes
  visible.
  EFFICIENCY: zero new EE calls - reuses the exact same raw study
  series already fetched for the permutation test (v10.122).
  DISCLOSED LIMIT: with typical 10-24 month windows, the climatology is
  built from only 1-2 samples per calendar month - a genuinely noisy
  estimate, explicitly stated in the output, not a robust multi-decade
  seasonal average. Shown as an exploratory comparison, not a
  validated replacement for the raw AC1 already used in the verdict.

v10.126 FIX: GLOBAL SIGNAL could hide a genuinely large local
  amplification
  Caught from a real Nuuk test: study variance rose +3.40x, control
  variance rose only +0.22x - both individually crossed the 0.15x
  threshold, so the tool called it "GLOBAL SIGNAL, not a local warning"
  even though the divergence between them (+3.18x) was over 21 times
  larger than the 0.15x divergence bar this SAME tool already uses
  elsewhere (the MARGINAL LOCAL SIGNAL branch, built earlier this
  session for a Florida Keys case). That branch was only ever reachable
  when NEITHER site crossed the threshold individually - a case where
  BOTH sites genuinely rose, but by wildly different amounts, had no
  branch that could describe it honestly. The divergence value was
  already being computed and displayed in this exact verdict text, just
  never used to decide it.
  Fixed: when both sites cross the threshold AND the divergence between
  them also exceeds the divergence threshold, the verdict is now
  "GLOBAL SIGNAL, WITH STRONG LOCAL AMPLIFICATION" - explicitly stating
  that a real regional signal is present, but the study site is
  amplifying it well beyond what the open-ocean control shows. Plain
  "GLOBAL SIGNAL" is unchanged for the case where both sites rise by
  comparable amounts.

v10.125 DIAGNOSTIC: S7E hang investigation (no fix yet)
  Caught on a real run: S7E stayed stuck at "1/2 batched calls done"
  even after a page reload (ruling out an expired auth session) AND
  after shortening to 12/10-month windows (ruling out the same window-
  length slowness already fixed for S7D at 18/12mo in prior versions).
  Neither previously-suggested fix resolved it - a genuinely new,
  undiagnosed failure mode, not one of the two already-known causes.
  Rather than guess a third fix blind, ported the exact diagnostic
  pattern already proven for S7D in v10.108: both of S7E's raw-fetch
  calls now print their ACTUAL error text on failure, or their real
  feature count on success, instead of the hang staying uninformative.
  DISCLOSED LIMIT: this can only help once a call eventually returns
  something (success or a real error) - it cannot reveal anything if
  the call is a genuine infinite hang with zero response from Earth
  Engine's servers. Re-run S7E and check the Console for "S7E [" lines
  to see exactly what happens this time.

v10.124 NEW: real permutation-test p-values for S7E (algae-based
  local vs regional classification)
  Direct fix for a gap flagged in discussion: the permutation-test
  engine (v10.122/123) only ever protected the SST-based S13 findings -
  S7E's algae-based LOCAL ANOMALY / REGIONAL SIGNAL classification still
  relied entirely on the same >0.01/>0.15 fixed thresholds as before.
  Ports the exact same permutationTestDelta() engine to S7E, computing
  real p-values for both study and reference sites (AC1 and variance),
  so the same "would random reshuffling of these months produce a delta
  this big by chance?" question now applies to algae, not just
  temperature.
  EFFICIENCY: needs ZERO new Earth Engine calls - S7E already fetches
  the raw monthly FAI series for both study and reference sites to
  compute its existing statistics; this just runs the proven engine on
  that same already-fetched data, purely client-side.
  Scoped to the standalone S7E tool for this release (not yet ported to
  S7D's within-reef network or S7F's combined view).

v10.123 NEW: real, Bonferroni-corrected permutation-test p-values for
  FIND SWEET SPOT
  Direct response to a real question: if window length changes the
  result (confirmed by a real Nuuk test - a 10-month AFTER window showed
  a rising-AC1 signal, a 24-month window at the SAME coordinate showed
  nothing), how can you draw a reliable conclusion from testing several
  window lengths and picking the best-looking one? That is a classic
  multiple-comparisons trap: testing 6 windows means SOME window can
  look significant by pure chance alone, even with nothing real
  happening - the existing ROBUSTNESS note (v10.95) only ever warned
  about this in words, never corrected for it mathematically.
  Now runs the same permutationTestDelta() engine validated in STEP 3
  COMPARE (v10.122) at each of the 6 tested window lengths (study AC1,
  study variance, control AC1, control variance = 24 tests total, 300
  shuffles each), and reports an explicit Bonferroni-corrected
  significance bar (0.05 / 6 windows = p<0.0083) alongside the
  uncorrected count - so a genuinely robust signal (surviving the
  stricter bar) can be told apart from a lucky single window among six.
  EFFICIENCY: adds only 4 new EE calls, not 12 - all 6 AFTER windows
  share the same start date and differ only in length, so a shorter
  window's raw values are always a prefix of the 24-month fetch;
  fetched once per site and sliced client-side for the other 5 lengths.
  Shown as a new, clearly-separated ADDITIVE section - does not alter
  the existing threshold-based table or ROBUSTNESS note above it.

v10.122 NEW: real permutation-test p-values for STEP 3 COMPARE
  Direct response to a design discussion: most of this tool's
  LOCAL/REGIONAL/signal decisions rely on fixed thresholds (>0.01 AC1,
  >0.15 variance) copy-pasted across STEP 3, STEP 4, S7D, S7E, and S7F -
  only STEP 5's Mann-Kendall test ever answered "is this bigger than
  random noise would produce" with a real p-value. A fixed threshold
  cannot tell a genuinely surprising 0.2 AC1 jump (site where months
  normally barely move) from an unremarkable one (site where ordinary
  weather noise wobbles AC1 by 0.3 anyway) - both get treated identically.
  Built a new, reusable permutation-test engine (permutationTestDelta):
  pools the BEFORE+AFTER raw monthly values, randomly reshuffles which
  months get relabeled BEFORE/AFTER (keeping window sizes fixed),
  recomputes the same AC1/variance statistic on each of 500 shuffles,
  and reports what fraction of PURELY RANDOM reshuffles produce a delta
  at least as large as the one actually observed - a genuine, honest
  p-value. Pure client-side JS, no new EE architecture: reuses the exact
  same jsLinearDetrendResiduals/jsLag1AC1/jsVarianceHalves building
  blocks already proven this session, and the same raw-value-fetch
  pattern (extractMultiNodeSeries) already used for S7D/S7E/S7F.
  Wired into STEP 3 COMPARE first (the highest-priority, most heavily
  tested tool) as an ADDITIVE result shown alongside - not replacing -
  the existing threshold-based verdict, in its own panel, via 2 new
  independent EE calls. Deliberately scoped to COMPARE only for this
  release; S7D/S7E/S7F/STEP 4 can reuse the same engine once this is
  validated in real testing.

v10.121 FIX: COMPARE could silently show a stale result from an earlier
  run
  Caught from a real multi-site test: COMPARE showed AC1=0.794 for the
  AFTER window at a new coordinate, while STEP 2's own "AFTER stored"
  label correctly showed 0.431 for the exact same window - and 0.794
  turned out to be an exact match for a PREVIOUS site's AFTER value from
  earlier in the same session. Confirmed by direct testing (re-clicking
  COMPARE fixed it immediately, showing the correct 0.431): this was NOT
  a computation bug - COMPARE's result panel is a one-time snapshot
  taken at click time, and never auto-refreshes if STEP 2 is re-run
  afterward. The underlying math was always correct; nothing warned the
  user the on-screen result had gone stale.
  Fixed: csdCompareRanWithB/A remember which stored BEFORE/AFTER objects
  COMPARE last used. If STEP 2 is re-run afterward (creating new stored
  objects), a visible warning now appears immediately: "STALE RESULT
  BELOW... press COMPARE again to refresh." The warning clears
  automatically the next time COMPARE actually runs.

v10.120 FIX: the v10.119 fix only covered the headline box, not the
  detail panel underneath it
  Caught immediately on the very next real Nuuk test: the headline
  correctly showed the v10.119 CONFLICTS WITH AC1 message, but the
  "4-WAY CSD COMPARISON (detail)" panel below it STILL showed the raw
  unreconciled "Control variance rose but the study site stayed flat...
  This is a POSITIVE result" text - the same contradiction, just in a
  second location v10.119 did not reach.
  Root cause: this detail panel is built early and SYNCHRONOUSLY, using
  the raw variance-only verdict text, before the AC1-weighted spatial
  toolkit finishes computing (an async callback that runs later and only
  rewrites the headline box, never this detail panel). Structurally
  stuck showing pre-AC1-weighting text.
  Fixed using the same conflict-detection approach as v10.119, but with
  studyAC1Rose (already computed earlier in the code, before this panel
  is built) since fullTally is not yet available at this point. Also
  fixed the standalone console print of the same raw text right after.

v10.119 FIX: COMPARE headline could contradict its own "Regional context"
  Same class of bug already caught and fixed for FIND SWEET SPOT in
  v10.94 (two separately-computed classifications sitting in one box
  with no reconciliation), found living on in COMPARE's main verdict box
  too. Caught from a real Nuuk test: the headline correctly said
  "STRONG SIGNAL, LIKELY REGIONAL" (AC1 rose +0.234 study vs +0.242
  control - nearly identical, genuinely regional per Dakos et al.
  weighting), but the SAME box then said "Regional context: STUDY SITE
  MORE STABLE THAN CONTROL (positive result)" - the opposite claim,
  driven purely by variance (fell at study, rose at control) with no
  awareness that AC1 - the primary indicator - told a different story.
  Root cause: "Regional context: "+vTitle pasted the OLD, pre-toolkit,
  variance-only classification text directly after the NEW AC1-weighted
  combinedTitle headline, unreconciled. Fixed: the conflict is now
  detected explicitly (AC1 agrees with a rising trend AND vTitle claims
  "study more stable") and replaced with a clear statement of the
  conflict instead of silently presenting two disagreeing verdicts as
  if they were consistent.

v10.118 CRITICAL FIX: S13 COMPARE could call zero data a "positive result"
  The same bug already caught and fixed for S7E in v10.114 existed in
  S13's core COMPARE verdict too - the single most heavily used part of
  this whole tool, hit in nearly every test throughout this session.
  Caught testing a genuinely new location (Nuuk, 64.13,-51.38): the
  study site returned AC1=n/a, Var=n/a in BOTH the BEFORE and AFTER
  windows (zero valid OISST months, likely ice/coastal-pixel masking),
  yet the headline confidently declared "STUDY SITE MORE STABLE THAN
  CONTROL (positive result)" and "the study site stayed flat... This is
  a POSITIVE result for this reef." The site was never measured at all.
  Root cause: studyVarRose/studyAC1Rose both silently default to false
  when the underlying deltas are null, making "no data" indistinguishable
  from "genuinely flat" - both fell into the same reassuring branch.
  Fixed: three new explicit checks run BEFORE any of the existing GLOBAL/
  LOCAL/ANOMALOUS branches - "CANNOT ASSESS (study has no data)",
  "CANNOT ASSESS (control has no data)", "CANNOT ASSESS (neither site
  has data)" - each stating plainly that this is not a stable/flat
  result, just an unmeasurable one, with likely causes and next steps.
  Checked FIND SWEET SPOT for the same bug: already safe (its own
  per-window logic already null-checks first, before this fix existed).

v10.117 FIX: S7F combined summary was missing data density entirely
  Caught directly: asked the user to check S7F's data density output to
  validate a LOCAL ANOMALY finding, and there was nothing to show - the
  density numbers were already computed internally (needed for the
  CANNOT CLASSIFY check) but never actually pushed into the displayed
  summary. S7D and S7E each show density in their own individual
  output; S7F's combined view silently dropped it. Now shown explicitly
  for both study and reference sites, with the same INSUFFICIENT flag
  used elsewhere, so a combined-run result can be checked for data
  quality without needing to separately re-run S7D/S7E individually.

v10.116 NEW: reference-site exclusion + S7F combined orchestrator
  Two features built from a real design discussion:
  1) REFERENCE EXCLUSION: S7E could auto-select a reference reef that
     coincided with a site already tested as a STUDY location earlier
     in the same session (caught directly: the reference for one test
     landed on Low Isles, which had already shown a real signal as a
     study site). A reference is supposed to be an uninvolved baseline.
     New global s7StudySiteHistory (via recordStudySite(), called from
     S7B/S7C/S7D/S7E/S7F) tracks every coordinate tested as a study
     site this session; S7E's candidate search now skips any GEBCO
     match within 5km of a prior study site, falling through to the
     next-closest genuine candidate, with the exclusion count shown in
     the results.
  2) S7F - RUN ALL (S7D + S7E COMBINED): direct response to a request
     for full automation. Built a SCOPED version, not the literal
     "everything automatically" version - explained why in the UI
     itself: this session already hit real EE account-level concurrency
     limits once (required a manual tier upgrade) and had a single S7D
     call hang 5+ minutes with no way to cancel it (no setTimeout in
     this sandbox). An orchestrator auto-trying multiple window lengths
     across every S7 tool would multiply both risks severalfold. S7F
     instead takes ONE shared Lat/Lon + BEFORE/AFTER input (entered
     once, not twice) and runs S7D + S7E together (6 EE calls total,
     same bound as running them separately), rendering one combined
     summary table with both headline verdicts side by side. Reuses
     the exact same proven functions (extractMultiNodeSeries,
     jsNodeStats, jsPairCorrelation, groupSeriesByLabel) S7D/S7E
     already use - no new statistics code, only new orchestration.

v10.115 NEW: S7E shows the reference site's actual coordinates
  Direct response to a real question: a user asked whether a map
  landmark (Low Isles) was the auto-found reference site. Previously
  impossible to answer - S7E only ever displayed distance and depth
  ("20km away, depth=-25.0m"), never the actual lat/lon, even though
  best.lat/best.lon were already computed at candidate-selection time
  and simply never surfaced. Now shown directly in the results line, so
  the reference candidate can be identified on the map (or its depth
  checked against known landmarks - e.g. an emergent cay/island reads
  near 0m elevation, not -25m, so a candidate at meaningful negative
  depth is a submerged patch, not a visible island).

v10.114 FIX: S7E asserted "LOCAL ANOMALY DETECTED" against an untested
  reference site
  A real run showed the headline verdict "LOCAL ANOMALY DETECTED" while
  the reference site's own numbers read ΔAC1=n/a, ΔVar=n/a - the
  reference wasn't tested and found stable, it simply had no computable
  data (its AFTER period had fewer than 4 valid months, even though its
  COMBINED before+after density of 55% looked fine and didn't trip the
  existing low-confidence flag, which only checks the combined average).
  "Untested" and "tested and found stable" are different findings, but
  the decision rule collapsed both into refSignal=false, letting a real
  study-site signal get reported as a CONFIRMED local anomaly against a
  reference that never actually confirmed anything.
  Fixed: added an explicit per-period (not combined) insufficient-data
  check for BOTH study and reference. When either site has <4 valid
  months in either individual period, the verdict is now "CANNOT
  CLASSIFY" with the specific period named, instead of silently
  defaulting into one of the four normal LOCAL/REGIONAL/ANOMALOUS/NO
  SIGNAL classifications as if the comparison were valid.

v10.113 FIX: S7E was hiding a coherent regional pattern inside "NO SIGNAL"
  The same gap already caught and fixed for STEP 5 in v10.97 (a real,
  significant DECLINE was being mislabelled as "no significant trend")
  existed in S7E too, uncaught until a real run exposed it directly:
  study AC1 fell -0.675, reference fell -0.636 - nearly identical, a
  coherent pattern across two independent sites. The SIGNAL check
  correctly only flags RISING AC1/variance (per Dakos et al., a decline
  is not itself a CSD warning), so this got silently folded into a bare
  "NO SIGNAL AT EITHER SITE" with no mention of the pattern underneath.
  Fixed: when both sites show a similar, substantial AC1 decline
  (<-0.1 each, within 0.15 of each other), the verdict now says so
  explicitly - labelled a REGIONAL DECLINE, clearly distinguished from
  a CSD warning direction, rather than describing real, structured data
  as indistinguishable from noise.

v10.112 FIX: S7E now shows STUDY site data density, not just reference
  A real run at a cloudy rainforest coastline (Daintree/Cape Tribulation)
  returned NO SIGNAL AT EITHER SITE with BOTH study and reference showing
  n/a for AC1/variance. The reference site's data density was already
  shown (44%, explaining ITS n/a), but the study site's own density was
  never displayed - no way to tell whether the study site failed from
  near-zero Sentinel-2 coverage (likely, given the location) or something
  else. Now shown symmetrically: study site density (X% of BEFORE/AFTER
  months valid) alongside the reference's, with the same LOW-CONFIDENCE
  flag logic applied to both. A sparse-data "no signal" result is now
  distinguishable from a genuine "checked and found nothing" result at
  the study site too, not just the reference.

v10.111 NEW: S7E - LOCAL vs REGIONAL auto-classification
  Closes a real gap: S7D's 9 nodes all sit within 0.3-3km of each other -
  they test spatial coupling WITHIN one reef, but every node is still
  local to that same reef, so S7D alone could never actually distinguish
  a local event big enough to blanket the reef from a genuine regional
  signal. S7E adds the missing genuinely-independent reference:
  - Auto-searches outward in expanding rings (20/40/70/110/160km, 8
    bearings each = 40 candidates) for the nearest point that is
    genuinely shallow water per GEBCO (-50 to 0m, the same threshold
    used for shallowMask elsewhere in this tool) - all 40 candidates
    checked in ONE batched reduceRegions() call.
  - Runs the SAME FAI-based BEFORE/AFTER AC1/variance extraction at that
    reference site as at the study site (reusing mkMoFAIRange/
    extractMultiNodeSeries/jsNodeStats unchanged - 2 more batched calls,
    3 total for the whole tool).
  - Applies an explicit decision rule, mirroring STEP 3's already-proven
    LOCAL/GLOBAL/ANOMALOUS/NO SIGNAL classification for SST: study
    signal without reference signal = LOCAL ANOMALY; both = REGIONAL;
    reference-only = ANOMALOUS; neither = NO SIGNAL. This is the actual
    automated classification layer that was missing - not just S7D's
    within-reef synchronization count, but a real core-vs-independent-
    reference comparison with an automated verdict.
  DISCLOSED LIMIT (not solved, only partially mitigated): a GEBCO
  shallow-water match is not guaranteed to be a real reef with
  comparable ecology - could be a bare sandbar with no algae community.
  Partially checked via the reference site's own data density as a
  proxy (a sparse signal there is flagged LOW CONFIDENCE, not silently
  trusted) - this does not prove ecological comparability, only flags
  the most obvious failure mode (a site with essentially no data at all).

v10.110 FIX: v10.109's NDVI fix over-corrected - fixed buffer size
  The v10.109 fix (mean()->max(), 150m->500m buffer) correctly resolved
  Center's false "LIKELY NOT on-reef" reading, but the fixed 500m buffer
  was itself a new bug: at the default 0.5km ring radius, a real run
  showed Center/N/NE/E all reading the IDENTICAL value 0.84 - strong
  evidence their 500m buffers were overlapping so much (Center-to-ring
  spacing is only 500m at that radius; adjacent ring points are even
  closer, ~383m) that multiple "different" nodes were just picking up
  the same peak pixel, silently defeating the purpose of sampling 9
  distinct locations.
  Fixed: NDVI check buffer now SCALES with the user's chosen ring
  radius (30% of node spacing, clamped 80-250m) instead of a fixed
  500m, verified mathematically to stay clear of overlap across the
  entire valid 0.3-3km radius range. The buffer size actually used is
  now shown in the on-screen results table for transparency. The FAI
  buffers used for the actual coupling statistics remain unchanged
  (150m fixed) - only the on-reef check's buffer changed.
  Honest caveat carried forward: whether the identical-value pattern in
  the real run was purely a buffer-overlap artifact, or partly a
  genuine large uniform algae patch, could not be determined with
  certainty from that one run - re-testing with this fix will show
  whether the values differentiate.

v10.109 FIX vs v10.108: two real bugs caught from an actual completed run
  The v10.108 diagnostics worked as intended: a smaller test run (6mo/6mo
  instead of 18mo/12mo) completed normally with matching feature counts
  (BEFORE=54, AFTER=54, NDVI=9 - exactly 6x9 and 9 as expected),
  confirming the earlier multi-minute hang was a size/complexity issue
  with larger month counts in one batched request, not a hard bug. That
  completed run then surfaced two real, separate problems:
  1) NDVI ON-REEF CHECK INCONSISTENCY: Center (the exact study
     coordinate, confirmed on-reef with strongly POSITIVE NDVI in every
     prior test this session) flipped to NEGATIVE (-0.32) and got
     flagged "LIKELY NOT on-reef". Root cause: mean() over a small 150m
     buffer at 20m scale is sensitive to small-scale heterogeneity - a
     reef patch smaller than the buffer gets its average pulled negative
     by adjacent clear water. Fixed by switching to max() over a larger
     500m buffer specifically for the on-reef check, matching the
     established convention already used by the main panel and S7B (the
     FAI buffers used for the actual coupling statistics are unchanged).
  2) MISSING ARTIFACT FLAG: the same real run showed Var=+52.31x and
     +11.99x - even more extreme than the case that prompted S7C's
     v10.105 near-zero-denominator artifact warning, which was never
     ported to S7D. Now added: same threshold logic (1st-half variance
     <0.001 AND |delta|>5x triggers a warning icon on that row).

v10.108 DIAGNOSTIC: S7D real-run failure needs actual evidence, not a
  guess. A real run at the CORRECT coordinates (-23.51, 152.09, confirmed
  by the map pin) returned n/a for every field on every node -
  "COULD NOT ASSESS COUPLING" with only "1 of 3 batched calls errored"
  reported, which doesn't fully explain why ALL fields (not just
  correlation) came back empty. Rather than guess at a fix blind (the
  v10.107 batching rewrite is new, untested-live code - reduceRegions()'s
  exact output property naming for a single-band mean reducer was
  assumed, not confirmed), this adds real diagnostics:
  - New s7dDiagnose(): prints the ACTUAL error text for any of the 3
    batched calls that fails (previously only an aggregate "N errored"
    counter existed, with no detail on WHICH call or WHY).
  - On success, prints the feature count and the first feature's raw
    property names - directly confirms or refutes whether
    reduceRegions() names its output the way extractReduceRegionsValue()
    assumes, instead of leaving that assumption unverified.
  - Raw feature counts (BEFORE/AFTER/NDVI) now shown directly in the
    on-screen results table too, not just the console - distinguishes
    "the fetch itself returned nothing" (a real data/query problem) from
    "data came back but extraction is broken" (a property-naming bug),
    which need different fixes.
  No fix applied yet - this version is instrumentation only, so the next
  real run's console output will show exactly what's happening.

v10.107 PERFORMANCE: S7D rebuilt from ~43 EE calls down to 3
  Not a quantum-computing question (the actual bottleneck - raster
  compositing, API round trips, server queuing - has nothing to do with
  the kind of problems quantum algorithms accelerate) - a boring,
  effective classical fix instead: batching.
  - New extractMultiNodeSeries(): for a monthly image collection, builds
    ONE FeatureCollection covering every node's whole time series at
    once, via reduceRegions() (samples all 9 points against ONE image)
    + flatten() across all months into a single flat table. One
    .evaluate() call now returns what used to require 9 separate calls
    (one per node, each re-triggering the full monthly Sentinel-2
    compositing graph from scratch).
  - New client-side (plain JS) statistics: jsNodeStats (linear detrend,
    lag-1 AC1, variance halves, skewness) and jsPairCorrelation (Pearson
    correlation between two aligned series) - identical formulas to the
    EE versions (computeRealCSD/computeZonalSyncCSD), just computed
    locally on the small already-fetched tables instead of triggering
    new server-side graph evaluations. Same "fetch once, compute
    locally" pattern already proven safe for STEP 5's Mann-Kendall test.
  - Total EE calls: 9 nodes x 2 periods x AC1/variance + 8 pairs x 2
    periods x correlation + 9 NDVI checks (~43 calls, 2-5 min) is now
    just 2 batched series fetches (BEFORE, AFTER) + 1 batched NDVI fetch
    = 3 calls, ~20-60s expected.
  - Output format, verdict logic, and UI inputs are UNCHANGED - this is
    a pure performance rebuild, not a new feature or a different result.

v10.106 NEW: S7D - full 9-node algae coupling network (BEFORE/AFTER)
  Scales S7C's 3-node proof-of-concept to the full 9-node ring (matching
  S7B's compass geometry), with a real BEFORE/AFTER comparison instead
  of a single snapshot window.
  - ~43 Earth Engine calls in one click: 9 nodes x 2 periods x AC1/
    variance (18), 8 ring points x 2 periods x correlation-vs-Center
    (16), 9 nodes x 1 NDVI-water on-reef check (9). Deliberately scoped
    to correlation-vs-Center only (8 pairs), not the full 36-pair
    matrix, to keep runtime bounded (~2-5 min instead of much longer).
  - NDVI-water fetched once per node (current composite) specifically
    to flag likely off-reef points, per the real S7B finding that North
    showed strongly negative NDVI there. A node/pair involving a
    flagged off-reef point is marked explicitly, not silently trusted
    as a real algae-dynamics comparison.
  - Headline verdict counts how many ON-REEF ring points show rising
    correlation with Center (the Dakos et al. 2011 hyper-synchronization
    direction), separately from raw AC1-rising node count.
  - Reuses computeRealCSD/computeZonalSyncCSD/mkMoFAIRange unchanged -
    no new statistics code, only new orchestration across more nodes
    and two time periods.
  Explicitly disclosed as an empirical correlation network, not a
  mechanistic J_ij interaction matrix - same honesty framing as every
  other coupling-style indicator already in this tool.

v10.105 FIX: S7C variance ratio could look like a huge surge when it was
  actually a near-zero-denominator artifact. Caught on a real run: East
  node showed Var=15.80x, far beyond anything seen elsewhere in this
  tool (SST-based ratios never exceeded ~7x). The ratio is computed as
  secondHalfVariance / max(firstHalfVariance, 1e-6) - if the first half
  had FAI sitting near-constant (plausible for algae presence, which can
  be genuinely near-zero for long stretches), even a modest second-half
  variance inflates the ratio dramatically without any real surge having
  occurred. Now shows the raw first-half/second-half variance alongside
  the ratio, and auto-flags the case (first-half variance < 0.001 AND
  ratio > 5) with an explicit "LIKELY ARTIFACT" warning - no more trusting
  a single derived number without the components that produced it.

v10.104 NEW: S7C - algae/AC1/variance coupling proof-of-concept
  Tests whether the SAME toolkit already applied to SST in S13 (AC1,
  variance, spatial correlation) can be meaningfully computed from ALGAE
  (FAI) data instead - arguably more theoretically appropriate, since in
  Scheffer's bistable-state framework coral-vs-algae cover is the STATE
  VARIABLE that actually flips between stable states, while SST is
  closer to the external control parameter driving the system toward a
  threshold.
  - New mkMoFAIRange(): the missing monthly Sentinel-2 FAI time series
    builder - S7 previously only ever had ONE fixed 2023-2024 composite,
    never a time series, so this statistic was not computable before.
  - computeRealCSD() and computeZonalSyncCSD() needed ZERO changes - both
    already took a generic bandName parameter, so passing 'fai' instead
    of 'sst' just works.
  - Deliberately scoped to 3 nodes (Center/North/East) and one window,
    not the full 9-node ring, because Sentinel-2's cloud-masking data
    density at monthly resolution is genuinely untested until now -
    reports an explicit DATA DENSITY assessment (% of requested months
    that actually had valid data) before anything else, so a sparse-data
    result is caught and flagged rather than silently producing unreliable
    statistics.
  - If density looks good, the UI explicitly recommends scaling to the
    full 9-node ring (matching S7B's geometry) with a proper BEFORE/AFTER
    comparison next - not built yet, pending this proof-of-concept result.

v10.103 NEW: real S15 field validation for One Tree Reef
  Found a real, published paper tracking the EXACT coordinate this tool
  has been tested against: Byrne et al. 2025 (Limnol. Oceanogr. Lett.,
  doi:10.1002/lol2.10456) tracked 462 individual coral colonies at One
  Tree Reef (23.51S, 152.09E) through the 2023-24 heatwave - 66% bleached
  by Feb 2024, 80% by April, up to 52% mortality by July, with genus-
  level detail (Acropora 95% mortality/rapid collapse to rubble;
  Goniopora 73% mortality via black band disease; more resilient genera
  like Porites/Pocillopora showing partial recovery). Data publicly
  available at Sydney eScholarship (doi:10.25910/p5rq-cw63).
  getEcologicalRecoveryValidation() now takes lat/lon and checks PROXIMITY
  (haversine distance, 60km radius) to this exact site before falling
  back to the coarse region-based lookup - deliberately NOT keyed to the
  whole "Great Barrier Reef" region bucket, since that bounding box spans
  from northern reefs (Lizard Island, hit hardest in 2016) to this
  southern site (largely spared until 2024) with very different
  bleaching histories; applying this southern-GBR-specific finding
  region-wide would have misrepresented the northern reefs.
  All 3 call sites (S13 single-window test, main click S15 display,
  export log) updated to pass lat/lon through.

v10.102 NEW: S7B - MULTI-POINT ALGAE SCAN
  Direct response to a real finding: two points 1.1km apart at One Tree
  Reef showed "WATCH - mild signal" vs "EXTREME - massive bloom" for the
  same S7 macroalgae indicators - a single click can badly misrepresent
  a patchy bloom in either direction (missing a real one, or overstating
  an isolated one).
  - Samples FAI, NDCI, and NDVI-water at 8 compass points (N/NE/E/SE/S/
    SW/W/NW) plus the centre, at a user-chosen radius (0.5-10km), using
    the SAME raw images and methodology (max reducer, 100m scale) as the
    existing single-point S7 panel - so results are directly comparable.
  - Built as ONE Earth Engine call via reduceRegions() over a 9-feature
    collection, not 9 separate calls.
  - Reports a WARNING SIGN COUNT (how many of the 9 points show elevated
    FAI/NDCI/NDVI) and an explicit PATTERN classification: WIDESPREAD
    (>=70% of points elevated - likely a genuine regional bloom),
    PATCHY/ISOLATED (<=2 points - could be a real small patch OR a
    single-pixel/cloud artifact, explicitly flagged not to generalize),
    or MIXED.
  - Also reports the numeric spread (max-min) across points as a direct,
    explicit measure of the kind of spatial heterogeneity that motivated
    this feature.
  - Has its own Lat/Lon input + "Use last clicked location" button
    (reusing the lastClickLat/lastClickLon pattern from S13 STEP 1), so
    it works independently of the main click flow or S13.
  Honest limit disclosed in the UI: a 9-point compass ring is a simple
  sampling pattern, not exhaustive coverage - a bloom could sit between
  sample points. Recommends a smaller radius or the S7 map layers for
  full visual coverage.

v10.101 CAVEAT DISCLOSURE (no computation changed): a direct question
  about whether v10.100's synchronization indicator could have detected
  the GBR event earlier prompted a re-check that surfaced a real design
  concern, not caught before shipping: computeZonalSyncCSD() correlates
  RAW monthly SST, not deseasonalized anomalies. Two ocean points a short
  distance apart share a strong seasonal cycle regardless of any real
  dynamical coupling change, which likely pushes the baseline correlation
  toward a high ceiling and blunts this indicator's sensitivity.
  Supporting evidence found while checking this: across every real GBR
  test run in this tool's history, the core AC1 statistic (which has the
  same raw/linear-detrend-only limitation) consistently landed in a
  narrow ~0.78-0.89 band regardless of site or period - more consistent
  with shared seasonal month-to-month autocorrelation dominating the
  signal than with genuine site-specific critical-slowing-down dynamics.
  A proper fix (per-zone monthly climatology subtraction before
  correlating) needs live testing to verify before shipping, so it is
  NOT implemented blind here. Instead: the limitation is now disclosed
  explicitly in the STEP 3 UI, in the indicator's own display name, and
  in the source code comments, with guidance to read the DELTA rather
  than the absolute value until a tested deseasonalized version exists.
  No answer this version gives should be read as "this would have
  caught it earlier" - that claim has zero empirical support yet.

v10.100 NEW: study-control temporal synchronization indicator (STEP 3)
  Implements the "empirical interaction/covariance network" idea in a
  deliberately minimal, low-risk form: rather than inventing new zone
  geometries (reef flat / mangrove / channel, which would need habitat
  datasets not verified reliable here), this reuses the deep-water
  control site STEP 3 already auto-selects as the second node. Computes
  the Pearson correlation between the study reef's and the control
  site's month-to-month SST fluctuations, for BEFORE and AFTER
  separately - tracking whether the reef is becoming MORE locked in sync
  with the open ocean (rising correlation - losing local independence/
  buffering, the leading-indicator direction reported by Dakos et al.
  2011, Am Nat 177:E153-E166) or staying decoupled.
  - New computeZonalSyncCSD(): the actual correlation math, reusing the
    same array-based approach already proven safe elsewhere in the file.
  - Added as a new SUPPORTING (not primary) indicator in STEP 3's
    toolkit, alongside spatial variance and spatial autocorrelation -
    fetched as one more follow-up in the same non-blocking chain, so a
    failure here can't delay or break the temporal/spatial results above.
  Disclosed honestly in the code comments: this is a correlation, not a
  mechanistic interaction coefficient (no J_ij matrix); shared external
  forcing (a heatwave hitting both sites) will raise this number with or
  without any real internal dynamics - one more piece of evidence, not
  proof on its own, consistent with how every other indicator in this
  toolkit is already framed.

v10.99 CRITICAL FIX: "Cannot read property 'trim' of undefined" crash
  Real GEE Code Editor error, caught live: clicking RUN SLIDING WINDOW
  ANALYSIS with the (intentionally optional, "blank = auto") window-size
  field never typed into threw a hard crash at csdSlideWindowInput.
  getValue().trim() - GEE's ui.Textbox.getValue() returns undefined (not
  an empty string) for a field that has never been interacted with, even
  with only a placeholder set. Calling .trim() on undefined throws and
  stops the whole button handler.
  This exact pattern - X.getValue().trim() with no guard - existed in 8
  places across the file (GO TO COORDINATES, STEP 2 RUN CSD TEST, the
  ADVANCED manual control override, FIND SWEET SPOT's AFTER date, and
  all three STEP 5 fields). Any of them could crash identically if
  clicked before typing into that specific field - the STEP 5 window
  field just happened to be the one a real user hit first, since leaving
  it blank is the intended, documented way to use it.
  Fixed everywhere at once: every X.getValue().trim() is now
  (X.getValue()||'').trim(), so a never-touched field safely reads as an
  empty string and falls through to the existing "please fill this in"
  validation message instead of crashing.

v10.98 FIX: STEP 5 silently used fewer valid months than requested
  Caught on a real Lizard Island run: requested 65 months, but only 51
  were valid (37 sliding positions with a 15-month window = 51 valid
  months), with no explanation anywhere in the output for the missing 14.
  Root cause: the requested range (2022-06 + 65 months) ran into 2027-10,
  past today's real date - OISST has no satellite observations for
  future months, so they came back null and were silently dropped from
  the valid-months count. The math was correct throughout; only the
  transparency was missing.
  Fixed two ways:
  - Proactive check BEFORE firing the Earth Engine call: if the requested
    start+total runs past today's date, an explicit warning appears
    immediately (client-side, instant) explaining how many months will
    be excluded and why - no need to wait 10-40s to find out.
  - After the run: if fewer valid months were found than requested for
    ANY reason (future dates or an ordinary data gap), the results table
    now states this explicitly with the exact count, instead of leaving
    the user to reverse-engineer a position-count mismatch themselves.
  Also: FIND SWEET SPOT's ROBUSTNESS note read "only 0 of 6 windows lean
  LOCAL" for the zero case - now reads "NONE of 6" for that case.

v10.97 FIX: STEP 5 headline could claim "no significant trend" when one
  genuinely existed. Caught on a real GBR run: variance showed tau=-0.273,
  p=0.017 - a real, statistically significant DECLINING trend - but the
  headline said "NO SIGNIFICANT TREND DETECTED (p>=0.05 for both AC1 and
  variance)", which was simply false about the number sitting right below
  it. Root cause: the verdict logic only checked for RISING significant
  trends (tau>0 && p<0.05); a significantly FALLING trend (tau<0 &&
  p<0.05, a real and reportable finding) was silently lumped into the
  same bucket as a genuinely flat, noisy series with no real trend at
  all. Two very different findings, one wrong label.
  Fixed: now classifies each indicator into rising-significant /
  falling-significant / not-significant separately. A significant
  decline now gets its own accurate headline (e.g. "VARIANCE
  SIGNIFICANTLY FALLING (p<0.05) - no CSD signature... this is an
  absence-of-warning-sign result, not a 'nothing found' result"),
  distinguishing it from both a rising CSD signal and a genuinely flat
  series. "NO SIGNIFICANT TREND DETECTED" is now only used when BOTH
  indicators are actually non-significant in either direction.
  The per-window trendVerdict() line text (which already correctly showed
  direction and significance separately) was not affected - only the
  summary headline classification had this bug.

v10.96 NEW: STEP 5 - SLIDING WINDOW TREND TEST (Dakos et al. 2012 method)
  Direct implementation of the actual gold-standard method, replacing
  discrete BEFORE/AFTER chunk testing with a continuous rolling window:
  - computeSlidingWindowCSD(): detrends the whole requested series ONCE
    (linear detrend, same method as computeRealCSD - Dakos et al.'s own
    toolbox typically uses Gaussian kernel smoothing instead, which is
    more flexible but not implemented here; disclosed in the UI), then
    slides a fixed-size window forward ONE MONTH AT A TIME across the
    residuals, computing AC1 and variance at every position. Built as a
    SINGLE server-side ee.List.map() graph, so the whole scan (which can
    be dozens of window positions) costs one Earth Engine round trip.
  - mannKendallTest(): the actual Kendall's-tau / Mann-Kendall trend
    significance test (standard normal approximation for the S-statistic
    variance), run client-side on the small, already-evaluated result -
    answers "is this metric moving consistently in one direction, or
    just bouncing randomly?" with a real p-value, not an arbitrary
    numeric threshold like the rest of S13 uses.
  - New STEP 5 UI: full series start date + total months (window size
    auto-defaults to 50% of total, the Dakos et al. standard, overridable).
  - Prints the actual AC1(t) and variance(t) trajectories as line charts
    to the console (with a linear trendline) - the "smoothly climbing"
    visualization Dakos et al.'s own figures show, instead of a single
    pass/fail verdict.
  - Verdict weights AC1 as primary (consistent with v10.90/91): AC1
    significant + variance not = still a strong headline; variance
    significant alone = explicitly downgraded per Dakos et al.'s own
    finding that variance is the less robust indicator.
  This does not replace STEP 3/4 (which remain useful for the LOCAL vs
  REGIONAL control-site comparison, something the classic sliding-window
  method doesn't do on its own) - it's a genuinely different, more
  statistically rigorous complementary tool for the specific question
  "is there a consistent trend at all," addressing the multiple-
  comparisons weakness of picking a single best window out of 6.

v10.95 NEW: ROBUSTNESS check for FIND SWEET SPOT - "which window do I trust?"
  Direct response to a real usage pattern: three separate test runs at the
  same GBR coordinates, testing slightly different AFTER windows, produced
  three different verdicts (LOCAL / MARGINAL / GLOBAL). Testing 6 window
  lengths and reporting whichever looks most dramatic is a classic
  multiple-comparisons trap - with 6 tries, SOME window will cross a
  threshold by chance alone, even in pure noise. The tool had no way to
  distinguish an isolated single-window spike from a signal that holds up
  across neighbouring windows.
  - New table footer: counts how many of the 6 tested windows lean LOCAL
    (LOCAL CSD or MARGINAL LOCAL) vs how many don't. 0-1 is flagged as an
    ISOLATED result with an explicit multiple-comparisons caution and a
    recommendation to re-test using an AFTER window chosen from
    independent evidence, not from this scan.
  - The headline verdict box itself now carries a short "(N of 6 windows
    agree)" or "\u26A0 ISOLATED" flag for LOCAL CSD / MARGINAL LOCAL results,
    so the caution is visible without reading the full table.
  - Added an explicit caution note to the S13 STEP 4 intro text.
  This does not add a real statistical significance test (the tool still
  has no surrogate/null-model comparison, unlike the literature it cites)
  - it only prevents the tool from presenting an isolated, likely-chance
    result with the same confidence as a signal that is consistent across
    multiple window lengths.

v10.94 FIX: FIND SWEET SPOT headline could contradict its own breakdown
  Discovered directly from a real GBR test: the top verdict box said
  "LOCAL CSD SIGNAL DETECTED" (strong, red) for a window whose own
  detailed breakdown said "MARGINAL LOCAL" and "Scheffer NOT MET -
  neither indicator rose at the study site". Root cause: the headline
  was classified from a SEPARATE, looser check (bestDiv>THRESH alone -
  pure variance divergence magnitude, ignoring whether the study site's
  OWN variance crossed the threshold, and ignoring AC1 entirely), while
  the per-window table/breakdown used a stricter, correct classification
  (bestRow.verdict: requires studyVarRose AND !ctrlVarRose for "LOCAL
  CSD", separately labelling divergence-only cases "MARGINAL LOCAL").
  Two different rules, two different answers, same window.
  Fixed: the headline is now derived directly from bestRow.verdict (the
  SAME classification shown in the table and the sweet-spot breakdown),
  so it can no longer disagree with the detail below it. Also added an
  explicit AC1-divergence-from-control check to the headline text
  (missing before - the multi-window headline was the last remaining
  place in S13 not reflecting the AC1-primacy work from v10.90/91).

v10.93 CRITICAL FIX: Southern Hemisphere reefs showed false "no heat
  stress" because the peak/DHW calculation used a single hardcoded
  window (Jun-Oct 2023) - Northern Hemisphere summer only. Discovered
  testing a real, independently-documented site: One Tree Reef, southern
  GBR, which suffered severe bleaching in the 2023-24 Southern Hemisphere
  summer (Nov-Apr) - the tool showed DHW=0.00 deg C-wks for this exact
  site/period because it was reading Southern Hemisphere WINTER SST as
  the "peak". This one number feeds: the main DHW badge, the S4 cancer-
  score component (25% weight), the map's bleaching-risk overlay, and the
  intervention engine's thermal-stress trigger - so heat-stress detection
  was silently broken for every Southern Hemisphere reef (GBR, Ningaloo,
  S. Indian Ocean, S. Brazil, S. Africa - roughly half the world's reef
  area by latitude coverage).
  Fix: compute BOTH a Northern Hemisphere peak window (Jun-Oct) and a
  Southern Hemisphere one (Nov-Apr, same 2023-24 stress year), then
  select per-click based on the clicked latitude's sign - in
  analyzeLocation() for the sidebar/score calculation, and in
  loadLayers() for the map layers (both already had lat as a parameter).
  MMM_perpixel (the seasonal baseline) was already hemisphere-agnostic -
  only the single "current peak" snapshot had this bug. Labels/legend
  updated to show which season window applied for each click.

v10.92 NEW: FIND SWEET SPOT can now surface AC1-up/variance-down windows
  Prompted by a direct question: "where can I find AC1 rising while
  variance falls, in a real tipping ecosystem?" Two real gaps surfaced:
  1) The 6-window results table only ever showed \u0394Var per window, never
     \u0394AC1 - so it was impossible to scan across windows for this exact
     pattern without manually re-running STEP 2 for each one by hand.
  2) The "sweet spot" window pick was ENTIRELY variance-divergence-based,
     inconsistent with v10.90/91 establishing AC1 as the primary
     indicator. A window with a strong AC1 rise but falling variance
     could be ranked below - or never surfaced above - a weaker
     variance-only window.
  Fixed:
  - Table now has a Study \u0394AC1 column alongside \u0394Var.
  - Per-window verdict text now appends "+ AC1 CONFIRMED (var down)" or
    "+ AC1 CONFIRMED" so this pattern is visually flagged in the table.
  - Added a SEPARATE AC1-based ranking (bestAc1W) alongside the existing
    variance-divergence ranking (bestW); table marks both with distinct
    arrows (VAR SWEET SPOT vs AC1 SWEET SPOT) when they disagree.
  - When the two rankings pick different windows, the console now prints
    an explicit cross-check note, citing Dakos et al. 2012 Fig. 2c/4 when
    it's specifically an AC1-up/variance-down window.
  No coordinates are hardcoded or guessed anywhere - this only makes the
  pattern findable by actually running the tool against real data.

v10.91 FIX vs v10.90: AC1-up + variance-down is not weak evidence
  v10.90 correctly downweighted "variance rises but AC1 doesn't" - but
  introduced an unwanted asymmetric bug in the OTHER direction: when AC1
  WAS rising and a supporting indicator (variance/spatial) was available
  but did NOT agree (e.g. variance actively falling), the confidence
  label said "no supporting indicator corroborates it yet" - wording that
  implies pending/weak evidence. That's wrong per Dakos et al. 2012: they
  specifically document AC1 rising while variance FALLS near a genuine
  transition (their Fig. 2c "decreasing sensitivity" and Fig. 4
  "freezing" cases) - autocorrelation "remains solely dependent on the
  dominant eigenvalue" and rises "regardless of the responsiveness of
  the ecosystem" to noise. Variance disagreeing does NOT weaken AC1's
  signal; it just means variance isn't a useful witness at that site.
  - classifyToolkitConfidence(): the "AC1 rising, support available but
    disagreeing" case now explicitly cites this documented pattern
    instead of implying the signal is pending confirmation. Still scored
    MODERATE (same as AC1-alone-with-no-data) - not upgraded to HIGH,
    since variance disagreeing isn't corroboration either, just neutral.
  - STEP 3's verdict title for this case no longer says "not yet
    corroborated" - now "AC1-confirmed - the primary indicator per Dakos
    et al. 2012".
  - Added an explicit S13 intro note about this exact scenario.

v10.90 FIX vs v10.89: AC1-weighted toolkit (Dakos et al. 2012 fidelity)
  The STEP 3 toolkit tally and FIND SWEET SPOT's Scheffer check both
  previously treated every indicator (AC1, temporal variance, spatial
  variance, spatial autocorrelation) as ONE EQUAL VOTE. Dakos et al. 2012
  (Ecology 93:264-271) found this isn't justified: autocorrelation
  "appears a relatively robust indicator... regardless of the source of
  noise" across every scenario they tested, while variance "may
  sometimes decrease close to a transition" for well-documented reasons
  (their Fig. 4 - parameter-noise sensitivity changes, and a "freezing"
  effect in slow-responding systems). A variance-only rise (AC1 not
  rising) is objectively weaker evidence than an AC1-confirmed rise, and
  the code did not encode that asymmetry anywhere.
  - buildToolkitTally() now tags AC1 as the PRIMARY indicator; variance
    and both spatial indicators are SUPPORTING evidence that raises or
    lowers confidence around it, not equal votes.
  - New classifyToolkitConfidence(): HIGH only when AC1 rises AND is
    corroborated; MODERATE when AC1 rises alone; LOW-MODERATE when
    variance/spatial rise but AC1 does NOT (explicitly flagged as weaker,
    citing Dakos et al.'s finding that variance can behave unexpectedly);
    LOW when nothing rises.
  - STEP 3's final combined verdict rebuilt around this same logic
    instead of a majority-of-indicators-agree rule.
  - FIND SWEET SPOT's "Scheffer 2009 validation" line now distinguishes
    AC1-only-rising (still meaningful) from variance-only-rising (weaker,
    flagged) instead of a flat "only one rose - weaker evidence" for
    both cases equally.
  No changes to the underlying AC1/variance/skewness/spatial math itself
  (computeRealCSD, computeSpatialEWS) - this is purely a re-weighting of
  how the existing numbers are interpreted into a verdict.

v10.89 FIX vs v10.88:
  S17 (Time of Emergence) was bundling all 6 variables (SST, Chl,
  Salinity, NO2, pH, DO) into ONE ee.Dictionary and evaluating them in a
  single .evaluate() call. pH and DO both depend on the Copernicus asset
  'COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/MFC_001_028', which is currently
  returning "not found" in the GEE catalog - and because all 6 were
  bundled together, that ONE dead asset made the entire evaluate() call
  fail, blanking ALL SIX indicators to n/a - including SST/Chl/Salinity/
  NO2, which come from completely different, healthy datasets and have
  nothing to do with the BGC asset. Split into rToeCore (SST/Chl/Salinity/
  NO2) and rToeBGC (pH/DO), evaluated independently: the 4 healthy
  variables now populate normally regardless of whether the BGC asset is
  reachable, and only pH/DO show the (accurate) "dataset unavailable"
  message. Compound status now reports "X of N AVAILABLE variables",
  with N excluding pH/DO when the BGC asset is down, instead of a flat
  "n/a - error" for the whole panel.
  S18 (Biogeochemistry Snapshot) was already correctly showing all-n/a,
  since all 4 of its values legitimately come from that same one dead
  asset - that part was not a bug. Its error message and the S17/S18
  section headers now state the specific dataset/reason plainly instead
  of a bare "n/a", so this is diagnosable without reading the console.

v10.88 NEW: genuine multi-indicator "toolkit" approach in STEP 3 (COMPARE)
  Addresses standard early-warning-signal guidance (Dakos et al. 2012 and
  related literature): don't rely on a single indicator; AC1 can be more
  reliable than variance; spatial patterns can be a stronger signal than
  temporal ones; limited data makes any one indicator unreliable.
  - computeRealCSD() now also returns skewness of the detrended residuals
    (reported for context; its direction is system-dependent so it is
    NOT auto-scored as CSD-consistent either way).
  - New computeSpatialEWS(): spatial variance and a spatial-autocorrelation
    proxy (correlation between each pixel and its 3x3 neighbourhood mean)
    computed across a 15km buffer around the study site for the BEFORE and
    AFTER periods - a genuinely different indicator family from anything
    computed through time at a point.
  - New TOOLKIT SUMMARY panel (csdToolkitV) tallies how many of up to 4
    scored indicators (temporal AC1, temporal variance, spatial variance,
    spatial autocorrelation) actually agree, with a confidence label
    (HIGH/MODERATE/LOW) based on agreement level and data sufficiency -
    instead of one threshold on one indicator deciding everything.
  - The verdict box now renders in two honest stages: a PRELIMINARY
    temporal-only verdict appears immediately (no added wait for the fast
    path), then upgrades in place once the spatial indicators arrive. A
    spatial-fetch failure degrades gracefully back to the temporal-only
    verdict rather than breaking anything.
  - Final combined verdict distinguishes "MULTI-INDICATOR LOCAL CSD SIGNAL"
    (majority of available indicators agree AND it's locally-flavoured),
    "...LIKELY REGIONAL" (majority agree but matches the control site too),
    and "WEAK/SINGLE-INDICATOR SIGNAL ONLY" (toolkit does not corroborate)
    - so a single indicator can no longer produce a strong-sounding verdict
    on its own.

v10.87 FIX vs v10.86:
  S12 (main click panel) had a misleading message: "Both AC1 and
  variance rising - classic CSD pattern". This fires whenever AC1 is
  ABOVE 0.5 in a single fixed window (Jan 2023-Dec 2024, no BEFORE
  baseline at all) - that is AC1 being CURRENTLY elevated, not AC1
  RISING. This directly contradicted S13's proper Scheffer 2009
  BEFORE-vs-AFTER test, which computes a real AC1 delta against a
  stored baseline and can correctly show AC1 FALLING even while S12's
  snapshot looks "elevated". Reworded the note to say what it actually
  measures, and added an upfront caveat on the S12 header pointing to
  S13 for a validated before/after comparison. No math changed - only
  the wording, which was actively misleading.

v10.86 FIX vs v10.85:
  STEP 3 (COMPARE) messaging was misleading and had a silent hang risk:
    - The verdict box jumped straight to "Computing verdict..." on click,
      while the detail box right below it still said "Step 1/3: Finding
      control site..." - the two boxes contradicted each other. Both now
      show the same real stage (finding control site -> control site
      confirmed, running CSD -> final colour-coded verdict), and the
      verdict box only claims "verdict" once one actually exists.
    - The control-site depth lookup (GEBCO reduceRegion) previously
      ignored its own evaluate() error parameter completely - if it
      failed or stalled there was no error message at all, just an
      indefinite "Step 1/3: Finding control site...". Now wrapped with
      explicit error handling using the same friendlyEEError() translator
      as the rest of S13.
    - Clarified in STEP 3's instructions that COMPARE runs in two
      sequential stages (find/validate control site, THEN run the actual
      comparison), not as one simultaneous action.

v10.85 FIX vs v10.84:
  Sidebar Depth label contradicted the S3 map legend. The legend splits
  depth into 7 bins (land / intertidal / shallow reef / continental
  shelf / continental slope / deep ocean / very deep ocean), but the
  Depth row's text used a crude 2-way "deeper than -50m = deep ocean"
  rule - so a -60m continental-shelf point (teal on the map) was
  labelled "deep ocean" in the sidebar right next to it. New
  classifyDepthLabel() applies the exact same 7 bins as the legend, used
  by both the Depth row and the "click closer to shore" warning, so the
  text and the map colour always agree.

v10.84 CHANGES vs v10.83:
  Hardened S13 error handling end-to-end:
    - New friendlyEEError() translator recognizes common transient Earth
      Engine errors (e.g. "Unknown reference to value named ''...",
      "Failed to contact Earth Engine servers", timeouts, rate limits)
      and tells the user plainly this is a server hiccup to retry, while
      still showing the raw error underneath for debugging.
    - STEP 2 (RUN CSD TEST) now echoes the exact parsed lat/lon/dates
      before firing, so a truncated/mistyped coordinate is obvious
      immediately instead of surfacing later as a cryptic failure.
    - STEP 3 (COMPARE / runControlCSD) previously ignored the evaluate()
      error parameter entirely for the control site - a failed control
      computation just silently produced blank/n-a numbers. Now wrapped
      in try/catch with explicit error surfacing in the verdict box.
    - STEP 4 (FIND SWEET SPOT)'s catch block now uses the same friendly
      translator instead of a raw error dump.
  No changes to the underlying math anywhere.

v10.83 CRITICAL FIX vs v10.82:
  FIND SWEET SPOT was throwing "TypeError: (intermediate value)..." and
  silently dying right after "Step 3/3: Analysing results..." because it
  called the ES6 method String.prototype.repeat() to draw a separator
  line ('─'.repeat(64)) - and GEE's server-side script sandbox does not
  implement that method. Replaced with a loop-based repeatChar(ch,n)
  helper everywhere a repeated-character string was built. No other
  behavior changed.

v10.82 CHANGES vs v10.81:
  FIND SWEET SPOT (STEP 4) rebuilt:
    - The control site's "before" state is now a REAL computed value
      (one shared control-BEFORE test using the same window as STEP 2),
      not a hardcoded 1.0 placeholder. Study and control deltas are now
      apples-to-apples, same maths as COMPARE.
    - Live progress counter ("N / 13 sub-tests done") while the 13
      parallel Earth Engine calls run, so the panel never looks frozen.
    - New csdMultiSweetSpotV panel: an explicit LOCAL vs REGIONAL vs
      Scheffer-2009-validation breakdown for the winning window (AC1 and
      variance deltas at both sites, stated in plain language).
    - try/catch around the final analysis so a computation problem shows
      a red error box instead of silently stopping at "Analysing results...".
    - Verdict box now also reports whether the Scheffer validation passed.

v10.81 CHANGES vs v10.80:
  REBUILT S13 (CSD Early Warning Test) for clarity:
    - Added numbered STEP 1/2/3/4 workflow headers so the panel
      reads top-to-bottom instead of as a pile of unrelated buttons.
    - Manual control-site override moved to a clearly labelled
      "ADVANCED (OPTIONAL)" section after the normal workflow;
      defaults to blank/AUTO, with a one-line explanation of why
      it exists. Nobody has to touch it to use S13.
    - Added a short, bold, colour-coded VERDICT box (csdCompareVerdictV)
      for the COMPARE button that renders ABOVE the detailed 4-way
      numeric breakdown - mirrors the verdict-above-table pattern
      already used by FIND SWEET SPOT (csdMultiStatusV).
    - Added a "Use last clicked location" button in S13 Step 1 that
      reuses the coordinates from your last map click / GO TO
      COORDINATES, via new globals lastClickLat / lastClickLon.
    - Clearer inline status messages after each RUN (tells you what
      to do next: switch to AFTER, then press COMPARE or FIND SWEET SPOT).
  All underlying math is unchanged: computeRealCSD(), getSmartControlSite(),
  runControlCSD(), and the CSD/verdict thresholds are untouched.

v10.67 CHANGES vs v10.66:
  NEW: GEM MarineBasis Greenland stations added to MODULE A11 (S19)
    - MarineBasis Nuuk GF3, Godthåbsfjord SW Greenland (64.13N,51.38W, r=40km)
    - MarineBasis Zackenberg Young Sound NE Greenland (74.315N,20.279W, r=30km)
    - API: api.g-e-m.dk | Key: GeoMarineAnalysis | License: CC BY-SA 4.0
    - DOIs: 10.17897/KMEK-TK21 (Nuuk CTD) | 10.17897/8GPS-CE70 (Zackenberg)
    - Stats are PLACEHOLDER nulls - run gem_fetch_and_clean.py to populate
  NEW: Greenland region entries added to getRegion()
  FIX (v10.19): Hemisphere suffix parsing in GO TO COORDINATES and S13
    - parseFloat("13.5S") silently returned 13.5 (positive) - now fixed
    - parseCoordPart() helper handles N/S/E/W suffixes correctly
  FIX (v10.20): DHW blue/orange layer contradiction in loadLayers()
    - "DHW blue = no stress" was painted over ALL ocean including stressed pixels
    - Now correctly masked to only show where DHW = 0 (genuinely no stress)
  FIX (v10.20): S13 BEFORE/AFTER CSD persistent state labels
    - Added "BEFORE stored:" / "AFTER stored:" labels that update in same callback
    - Added "Clear stored BEFORE/AFTER" button to reset stale state
    - Previously no way to confirm whether storage succeeded before clicking COMPARE

v10.66 key features (unchanged from v10.66):
  v10.66: GEBCO depth layer in map (matches sidebar depth value)
  v10.65: GEBCO replaces ETOPO1 for sidebar depth readout
  v10.64: honest "insufficient data" message for aquaculture (not just "land")
  v10.63: peak SST checked in Gate 1, not just annual mean (Persian Gulf fix)
  v10.60: parallel evaluate() for S17/S18
  v10.58: disasterty GDIS property name fix (was silently always 0)
  v10.57: real cited A. taxiformis SST thresholds (Statton 2024, 17-21 optimal)
  v10.56: two-gate aquaculture architecture replacing opaque weighted average
  v10.54: OISST *0.01 conversion fix in ToE annual SST stack
  v10.42: dead soil texture asset replaced with safe masked constant
  v10.41: defensive S2 band check before normalizedDifference()
  v10.39: removed setTimeout (not available in GEE sandbox) - harmless no-op
  v10.37/38: error handling + full sidebar reset on every new click
  v10.36: S12-S19, thermal recovery, ecological validation, ECI, ToE, BGC
  v10.13: OISST replaces dead MODIS dataset

SCIENTIFIC APPROACH:
  Waddington Landscape (PNAS 2025) + Scheffer 2009 CSD +
  Ramamurthy 2024 (Bunodosoma) + Levitan 2023 (Diadema) +
  Lozano-Bilbao 2020-2024 (metals) + Peixoto 2025 (Red Sea)
============================================================

============================================================
MODULE A - DATASETS (all with ocean mask)
============================================================
v10.161 FIX 19 / S1 - THE VERSION MARKERS ARE NOW DERIVED, NOT TYPED.
A version marker has now been missed in THREE CONSECUTIVE ROUNDS (found at
v10.149 in round 3, "all five reconciled" in round 5, and the sidebar footer
still read v10.159 in the live v10.160 browser run). Patching the literals
has failed every time, so the LITERALS ARE GONE. TOOL_VERSION below is the
only place in this file where the RUNNING version is written down; every
self-identifying marker on screen - sidebar title, S13 section header, the
sidebar footer, the per-click console banner and the startup READY line -
is built by concatenating it. Bumping a version is now a one-line edit, and
a marker CANNOT fall behind because there is nothing left to fall behind.
SCOPE, stated so the next round does not "fix" the wrong thing: this covers
markers that assert "the tool you are running is version X". It deliberately
does NOT cover (a) changelog entries, which are prose about what some past
version changed and must keep their own version numbers, (b) provenance
annotations like "(v10.160, Monte Carlo)" or "S7D - ... (v10.106)", which
record WHEN a figure was measured or a module was introduced and are wrong
if they move, or (c) the file name, which is deliberately stale (see README).
```
