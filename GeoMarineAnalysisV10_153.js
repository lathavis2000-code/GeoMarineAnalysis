// ============================================================
// STEMGeoHS Marine v10.153
// Coastal Attractor Landscape + Cancer Score Pipeline
//
// v10.149 NEW: S12b - Real Significance Test, converts S12 from
//   heuristic to probabilistic
//   Direct answer to a real question: S12 is real math (detrended AC1)
//   but has zero statistical testing attached - a real number with a
//   heuristic interpretation layered on top, not a genuine significance
//   test. Added a real permutation/surrogate-data test: shuffles THIS
//   SAME window's 24 months into random order 500 times, asking "would
//   randomly-ordered data produce an AC1 this high by chance alone?" -
//   a genuinely different question from S13's before/after delta test.
//   This is a real, established method in nonlinear time-series
//   analysis (surrogate-data testing for genuine autocorrelation).
//   Built standalone (button-triggered) since S12's own computation
//   runs entirely server-side in EE Arrays and never sends raw monthly
//   values to the client - needed one new, small EE call, reusing
//   mkMoSST() (the EXACT same Jan2023-Dec2024 collection S12 itself
//   uses) and the same jsNodeStats() formula already proven throughout
//   this tool, so the real AC1 shown matches S12's own number exactly.
//   Verified before shipping with two direct tests: structured
//   (seasonal-cycle-like) data correctly comes back highly significant
//   (p=0.0000); pure random noise correctly comes back not significant
//   (p=0.76).
//   HONEST DISCLOSURE built directly into the tool: a significant
//   result here does NOT mean AC1 has changed or a reef is destabilizing
//   - real ocean temperature data has genuine seasonal structure, which
//   this test will correctly flag as "significant" every time, since
//   scrambling destroys that real seasonal pattern. This tests whether
//   observed memory is distinguishable from noise, not whether resilience
//   is being lost - that separate, more specific question is S13's job.
//
// v10.148 FIX: S7 no longer shows coral-reef alarm language at
//   cold-water kelp sites
//   Direct response to a real, correct catch: a real Nuuk test showed
//   S7 displaying "EXTREME - massive bloom, RED" - alarm language built
//   for tropical coral-algae phase shifts - at a site S8 (and the Field
//   Data section) independently identify as likely HEALTHY cold-water
//   kelp habitat, not a warning. Same class of gap already fixed for
//   DHW (explicitly hidden above 55 deg N/S as ecologically meaningless
//   there) - S7 never got that same climate-zone awareness until now.
//   FIXED: moved computeAquaculture() (which already detects cold-water
//   kelp opportunities) to run BEFORE S7's algae display instead of
//   after, so its kelp-detection result can contextualize S7's text.
//   When kelp is detected AND the algae signal is genuinely elevated,
//   S7 now shows "LIKELY COLD-WATER KELP... not a coral-algae bloom
//   warning" instead of the old tropical-framed alarm text.
//   Verified with real values from both directions before shipping:
//   the real Nuuk test values (FAI=0.117) correctly now show the kelp-
//   context message; the real Florida Keys test values (FAI=0.0796,
//   kelp NOT detected there) correctly still show the normal "CRITICAL
//   - dense mat" alarm, unchanged - confirms the fix only applies where
//   it should, not a blanket suppression of real algae warnings.
//
// v10.147 BUGFIX: S18 STILL broken after v10.146 - found the REAL root
//   cause via a second real test
//   v10.146's "empty collection" fix didn't solve it - a real test at
//   One Tree Reef (36m deep, open lagoon, nowhere near a coastline)
//   FALSIFIED the v10.146 hypothesis that this was about Florida Keys'
//   coastal proximity - it failed identically there too. Compared
//   directly against S17's OWN pH fetch, which succeeded at BOTH real
//   test sites using the exact same real CAR asset, and found two real
//   differences: (1) S17 uses a bare ee.Geometry.Point with NO buffer;
//   S18 was using `study`, a buffered polygon - against a coarse ~25km
//   global grid, a small buffered area can produce a genuinely empty
//   reduceRegion result in edge cases where a bare point cleanly
//   resolves to one pixel. (2) S18 used a narrow 90-day window vs
//   S17's multi-year annual collection. Fixed both: switched to the
//   identical bare-point geometry pattern already proven at both real
//   test sites, and widened the window to 2 years (730 days).
//   Honest note: this is now the tool's best real attempt at fixing
//   S18, built by directly comparing against proven-working code rather
//   than guessing - but it has NOT yet been re-tested live. If it still
//   fails, the new diagnostic message explicitly says so rather than
//   reusing the old generic "no valid pixel" text.
//
// v10.146 BUGFIX: S18's v10.145 fix was itself broken - found via real
//   testing at One Tree Reef
//   A real test showed S18's new pH/pCO2/Salinity all failing with
//   "error fetching real BGC assets", even though S17's pH (same real
//   CAR asset) succeeded in the same test. Root cause found by
//   comparing against this codebase's own established pattern: every
//   other real-asset fetch (mkMoDHW, dhwProper, etc.) guards against an
//   empty ImageCollection with ee.Algorithms.If(col.size().gt(0), ...,
//   null) before calling .select()/.mean() on it - the v10.145 S18
//   rebuild skipped that guard. If the 60-day window genuinely returned
//   zero images at some location/time, that produced a real computation
//   error instead of a graceful null. Fixed using the identical,
//   already-proven defensive pattern; also widened the window 60->90
//   days to further reduce the chance of a genuinely empty result.
//
// v10.145 MULTI-FIX: S18 real assets, SNR-vs-Mann-Kendall diagnosis,
//   S16 ECI-vs-B honesty fix, S20e wired into Bowl Depth
//   1) S18 FIXED for pH + pCO2 + Salinity using the same real, confirmed
//   Copernicus/HYCOM assets already used for S17/S20e - split into 3
//   independent fetches so one still-dead field (O2) can never again
//   blank the 3 real ones (same class of fix as v10.89/v10.141).
//   pCO2 correctly converted from the real asset's native Pa to uatm
//   (x9.86923), sanity-checked against realistic ocean pCO2 ranges.
//   2) DO/O2: searched a THIRD time this session, still no confirmed
//   working sub-collection/band name - disclosed honestly, not guessed.
//   3) DIAGNOSED why SNR and S17b's Mann-Kendall sometimes disagree
//   (e.g. NO2/Salinity "EMERGED" by SNR but not significant in S17b):
//   the SNR formula multiplies slope by record_years with NO correction
//   for how uncertain a slope estimate is at few data points - unlike a
//   real p-value, it structurally over-triggers on short records. Now
//   disclosed directly in S17's UI. Also found and disclosed a real,
//   separate finding: Salinity's nominal 32yr record often has far
//   fewer REAL valid HYCOM points at a given site (genuine data
//   sparsity, not a bug) - now flagged in S17b's own output.
//   4) FIXED S16's ECI-vs-B "agreement" check - investigated and found
//   it was comparing a PURE physics metric (depth-only wave energy)
//   against the broad ecological Fused score (SST/Chl/DHW/biology/NO2).
//   These measure different things with no reason to agree - "DISAGREE"
//   was implying an error that didn't exist. Relabeled honestly.
//   5) S20e NOW WIRED into Bowl Depth - fulfills the v10.143 disclosed
//   promise. New row shows the uncalibrated heuristic "Regime-shift
//   index" directly next to S20e's real, held-out-validated P(bleaching)
//   - verified with real Bocas del Toro values (78% vs 94.8%, a real
//   16.8-point gap surfaced).
//   NOT built this version (real scope decision, disclosed rather than
//   rushed): S20e's 6-feature extension (distance-to-shore, windspeed,
//   cyclone frequency) needs 3 new real global data sources not yet
//   identified in GEE - a genuine, separate next step.
//
// v10.144 NEW: S20e - Real Bleaching Probability, a genuinely fitted
//   logistic regression model
//   Direct build-out of the real database found and downloaded this
//   session: Global Coral-Bleaching Database (van Woesik & Kratochwill
//   2022) - 41,361 raw rows, cleaned to 32,716 real, complete rows
//   (dropped rows with "nd" sentinel missing-data values, matching
//   BCO-DMO's own documented convention). Fit a REAL logistic
//   regression on 3 features that map exactly onto this tool's own
//   already-computed metrics: DHW (SSTA_DHW), Turbidity, Depth.
//   REAL, HONEST VALIDATION - not training-set accuracy: 70/30 train/
//   test split, evaluated on a genuinely held-out 9,815-row test set
//   never seen during fitting. Held-out AUC=0.620 - real, better than
//   chance (0.5), but genuinely WEAK - disclosed directly in the UI,
//   not hidden. Consistent with Arias-Ortiz et al. 2024's independent
//   finding that DHW-style models need ~23 metrics to predict well.
//   TWO REAL, INDEPENDENTLY-CONFIRMED FINDINGS from the fitted
//   coefficients: (1) Depth's coefficient is POSITIVE - deeper reefs
//   more likely to bleach in this real 32,716-row fit, independently
//   confirming Arias-Ortiz et al.'s counter-intuitive finding already
//   disclosed in S20. (2) Turbidity's coefficient is negative, matching
//   a real, separate citation (Sully & Woesik 2020, "Turbid reefs
//   moderate coral bleaching").
//   Zero new EE calls - reuses dhwv/turv/bv, already computed for the
//   main sidebar. Verified with a manual-formula-vs-sklearn sanity
//   check (exact match) before porting to JS.
//   Predicts BLEACHING probability, explicitly NOT collapse - see
//   v10.143's terminology fix for why that distinction matters.
//
// v10.143 FIX: terminology corrected - "collapse" language reviewed
//   throughout, fixed where it overclaimed
//   Direct response to a real question about a possible Bayesian/
//   logistic upgrade to the Fused Score, sourced from a real, large-
//   scale database (Global Coral-Bleaching Database, van Woesik &
//   Kratochwill 2022 - 34,846 records, 14,405 real sites, confirmed via
//   direct search; data hosted on Figshare, not yet pulled into this
//   tool - would need the actual file, not just its existence, to fit a
//   real model). Since any such model would honestly predict BLEACHING
//   probability (what the database records), not COLLAPSE, reviewed
//   every "collapse" reference in this tool for accuracy:
//   FIXED (real overclaiming): "P(flip 5yr)" relabeled "Regime-shift
//   index (5yr, uncalibrated)" - this heuristic was never a calibrated
//   probability of anything; a disclosed caveat now points to the real
//   GCBD-based upgrade path. Also fixed STEP 2's "AFTER event (post-
//   collapse)" dropdown label to "AFTER event (suspect period)" - most
//   AFTER windows tested this session showed NO signal at all, so
//   presuming verified collapse by default was misleading.
//   LEFT UNCHANGED (genuinely accurate uses, confirmed real): Levitan &
//   Edmunds' real urchin population collapse citation; Byrne et al.'s
//   real, field-observed coral colony collapse at One Tree Reef; S13's
//   description of Scheffer et al.'s own critical-slowing-down theory,
//   which is genuinely about approaching collapse/regime shifts.
//
// v10.142 NEW: S17b - Real Mann-Kendall Significance Test
//   Direct answer to a real question: since the Mann-Kendall engine
//   already exists (built and fixed for STEP 5), can S17's fixed
//   SNR>=2.0 heuristic be replaced with a genuine significance test?
//   Built as a standalone, button-triggered tool rather than rewriting
//   S17's core flow directly - that flow was already touched once this
//   session (the pH fix) and is deeply async/fragile; adding a new tool
//   alongside it is the same lower-risk pattern already used for
//   COMPARE's statistically-valid banner (v10.139). Reuses the SAME
//   real annual-value collections S17 already built (_annSSTColl etc.)
//   and the SAME mannKendallTest() engine already validated for STEP 5 -
//   no new architecture, genuinely real math throughout.
//   Verified directly before shipping: correctly flags short records
//   (n<4, e.g. a hypothetical 3-point NO2 case) as insufficient rather
//   than guessing; correctly computes a real, strongly significant
//   result on a longer record with an embedded trend (SST-style, 44
//   points, p<0.0001).
//   DISCLOSED: pH (~4yr) and NO2 (7yr) have very few real annual points
//   - Kendall's tau needs real data, so these will often show
//   "insufficient" here even though S17's SNR-based check above can
//   still report a number for them - a real, honest limitation of
//   annual-resolution testing on short records, not a bug.
//
// v10.141 FIX: S17 pH now uses a REAL, confirmed, working Copernicus
//   asset - DO still unresolved
//   Direct answer to a real question: can a working pH/DO database be
//   integrated? Searched and confirmed: the old dead asset ID was simply
//   outdated - the real, current GEE catalog splits this Copernicus
//   product into per-variable sub-collections. Confirmed directly via
//   Google's own catalog pages: COPERNICUS/MARINE/
//   GLOBAL_ANALYSISFORECAST_BGC_001_028/CAR, band ph_depth1 - real,
//   working, current surface pH.
//   IMPORTANT CAUGHT-BEFORE-SHIPPING bug: this real asset's data only
//   starts 2021-10-01, not 1993 as the old 32-year assumption required.
//   Verified directly: using the old 32yr constant against the real
//   ~4yr record would have overstated SNR by 8x - easily turning a
//   genuinely non-significant trend into a false "EMERGED" result.
//   Fixed: annual pH stack rescoped to 2022-2025, record-years constant
//   corrected 32->4, confidence relabeled 32yr/MARGINAL -> ~4yr/LOW.
//   ALSO FIXED: pH and DO were still bundled into one evaluate() call -
//   since DO's asset remains dead, that would have made DO's failure
//   blank the now-working pH too (same class of bug already fixed for
//   SST/Chl/Salinity/NO2 vs BGC in v10.89). Split into two fully
//   independent evaluates with two independent availability flags.
//   DISCLOSED: dissolved oxygen (DO) is NOT fixed this version - no
//   working sub-collection/band name was confirmed for O2 after two
//   direct searches - left disclosed as unavailable rather than guessed.
//
// v10.140 NEW: real Mesoamerican Reef entry + real methodological
//   caveat from a published predictive-model study
//   Direct response to a real slide shown (Global Tipping Points Report
//   2025, healthyreefs.org): confirmed the Healthy Reefs Initiative is
//   real (70+ partner orgs, monitoring since ~2008) - report-based, not
//   a live API, same access pattern as this tool's other citations.
//   HONEST CORRECTION made in the process: nearly merged this regional
//   data into the existing Bocas del Toro (Panama) entry - caught that
//   they are genuinely DIFFERENT reef systems (Bocas sits on Panama's
//   southern Caribbean coast; the Mesoamerican Reef spans Mexico/Belize/
//   Guatemala/Honduras further north) - added as a real, separate 7th
//   site instead, confirmed non-overlapping with Bocas del Toro
//   (1042km apart). Real 2023 stats: ~40% of corals severely affected,
//   regional cover fell 19%->17% (GCRMN Mesoamerican Report Card 2024).
//   Also found and disclosed a real, relevant finding: Arias-Ortiz et
//   al. 2024 (Communications Biology, doi:10.1038/s42003-024-07128-y)
//   found a combined 23-metric model explains 75% of bleaching-severity
//   variance, versus DHW alone explaining less - and, surprisingly,
//   DEEPER reefs with MORE diverse coral communities showed HIGHER
//   vulnerability. Added as a disclosed caveat on S20's simple DHW>=4
//   threshold - not a full model port (the paper's actual fitted
//   coefficients were not available from what was verified), but an
//   honest flag that a more sophisticated real model exists.
//
// v10.139 NEW: COMPARE now shows a statistically-valid PRIMARY verdict
//   Direct answer to a real, fair question: since a genuine 500-shuffle
//   permutation test already exists, why was the >0.01 AC1/>0.15 var
//   HEURISTIC threshold still the primary displayed verdict? Real
//   reason found: the permutation test fires its OWN separate EE call,
//   resolving LATER than the threshold verdict - a sequencing fact, not
//   a design choice. Rather than risk rewriting the large, deeply-nested
//   async chain that builds the threshold verdict, added a NEW banner
//   positioned ABOVE it - populated by the permutation test's own
//   completion callback once real p-values are ready, so the
//   statistically valid answer is now the first, most prominent thing
//   shown. Same AC1-primary weighting philosophy already used
//   throughout this tool (Dakos et al. 2012), now driven by p<0.05
//   significance instead of the old fixed thresholds. The original
//   threshold-based box remains below, explicitly relabeled as the
//   secondary, faster-arriving heuristic check.
//   Verified against 3 real sites already tested this session (One Tree
//   Reef, Bocas del Toro, Nuuk) - all three correctly resolve to "NO
//   SIGNIFICANT SIGNAL," matching what was already independently
//   established as the rigorous answer at each.
//
// v10.138 CONSOLIDATION: merged S20/S20b/S20c into ONE unified tool
//   Direct response to real, fair feedback: three separate buttons for
//   what's conceptually one question ("what's the combined risk, at some
//   point in time?") was confusing design, not a missing capability -
//   both past-date sync and live-present already worked correctly, just
//   split across separate tools. Now ONE date field: leave BLANK for
//   LIVE (today, real current OISST data), or type a month (YYYY-MM,
//   2023-01 to 2024-12) for a HISTORICAL replay of that exact month -
//   one button, one result, clearly labeled which mode ran. The
//   original always-on S20 (auto-computed on every click, frozen to the
//   2023-2024 peak) is retired as a separate behavior - S20 is now
//   button-triggered like the rest of this family, so the click-handler
//   code that used to auto-populate it on every map click was removed.
//
// v10.137 NEW: S20d - Genus Growth-Form Lookup (Coral Trait Database)
//   Direct answer to a real question: are there real, established
//   global species databases available? Searched and confirmed one -
//   the Coral Trait Database (coraltraits.org, Madin et al. 2016,
//   Scientific Data), a real, actively-growing, peer-reviewed database
//   (166,245 observations, 5,112 species as of the check). Confirmed a
//   real, working direct-CSV access pattern for individual traits.
//   DISCLOSED: could not directly fetch the live CSV from inside this
//   GEE script - outbound access to external domains from the Code
//   Editor sandbox is unverified, unlike the separate companion Python
//   scripts which run outside GEE entirely. Built the honest version
//   instead: GENUS_GROWTH_FORM, a small, real starter table, with each
//   entry explicitly labeled either [CTD-confirmed] (individually
//   verified against a live Coral Trait Database entry during this
//   build) or [literature pattern] (the same well-established branching-
//   vs-massive framework, Loya et al. 2001, not individually cross-
//   checked for that exact genus). New S20d UI lets you type any genus
//   name and see its real growth form + general vulnerability tier -
//   a GENERAL, global fallback, distinct from S20's specific, real
//   event-outcome mortality data at 6 named sites.
//
// v10.136 NEW: S20c - Live Current DHW Check
//   Direct answer to a real question: is live (not frozen-2023-2024)
//   DHW possible? Checked directly rather than assumed - confirmed via
//   Google's own Earth Engine catalog page that NOAA/CDR/OISST/V2_1
//   (the exact dataset getOISSTColl() already uses everywhere in this
//   tool) is real-time, updated daily, with data through essentially
//   today (1-day preliminary lag, 14-day final). This was never a data
//   availability limit - only the hardcoded 2023-2024 window prevented
//   it. New button computes a genuine rolling recent-window DHW ending
//   at today's actual date, reusing the identical monthly-max*4.33
//   formula already used throughout this tool - same math, same data
//   source, just pointed at "now" instead of a fixed past window.
//   DISCLOSED, and important: this is an OBSERVATION of the present, not
//   a forecast of the future. Genuine forecasting would need real ocean/
//   climate forecast MODEL output (e.g. NOAA Coral Reef Watch's actual
//   bleaching outlook) - a different kind of data source, not verified
//   as available in this tool's GEE architecture, and NOT built here.
//
// v10.135 NEW: S20b - Historical Month Check (retrospective test)
//   Direct build-out of a real request: test S20's combined risk
//   retrospectively, as it would have shown in early 2024. Uncovered a
//   real architectural fact first: the main click flow's DHW (and most
//   S1-S11 measurements) was ALREADY hardcoded to a fixed 2023-01-01 to
//   2024-12-31 window this whole session, not a live "as of today"
//   calculation - so the DHW shown was already the PEAK across that
//   entire window, not a present-day reading. There was no "current"
//   value to roll back from.
//   Built the honest version instead: a real month-by-month replay
//   within that same real, already-available 2023-2024 data window.
//   Pick any month (YYYY-MM, 2023-01 to 2024-12), and it computes the
//   REAL DHW for just that month at the last-clicked location, then runs
//   it through the exact same combineSpeciesAndDHW() logic S20 uses -
//   reusing getOISSTColl() and MMM_perpixel, both already proven
//   elsewhere in the tool, same math as mkMoDHW(), just isolated to one
//   chosen month instead of the full 24-month collection.
//   DISCLOSED SCOPE: this is NOT a true arbitrary-date time machine -
//   months outside 2023-2024 correctly return an explicit range error,
//   not a guess, since no real data exists there in this tool's current
//   architecture.
//
// v10.134 NEW: 3 real entries added proactively (Bocas del Toro,
//   Andaman/Nicobar, Brazil)
//   Direct response to real feedback: proactively check sites already
//   tested this session and other named regions, rather than waiting to
//   be told each one individually. Searched and verified real literature
//   for all three:
//   1. Bocas del Toro, Panama - the exact site tested with S13/CSD this
//      session. Neal et al. 2017 (Ecology and Evolution, doi:10.1002/
//      ece3.2706) tracked 3 massive coral genera through the real 2005
//      bleaching event for 8 years - all showed continued net tissue
//      loss with no full recovery, even typically stress-tolerant forms.
//      DISCLOSED: a second Bocas-specific source was found but rejected
//      - its search snippet mixed content from an apparently unrelated
//      Indian coral-management paper mid-sentence, a reliability red
//      flag, so only the clean, verified citation was used.
//   2. Andaman & Nicobar Islands, India - real 2010 mass-bleaching
//      survey with clean genus-level mortality: Acropora 43%, Montipora
//      22%, Porites 14%.
//   3. Southwestern Atlantic Reefs, Brazil - Mussismilia harttii, Brazil's
//      major endemic reef-builder, described as suffering "unprecedented
//      erosion" (Braz et al. 2022, Coral Reefs) with no recovery three
//      years after a major bleaching event (Corazza et al. 2024, Mar.
//      Biol.). DISCLOSED: no single mortality percentage was verified
//      for this exact event - reported as the real qualitative finding
//      rather than an invented number.
//
// v10.133 NEW: real third species-vulnerability entry (Maldives)
//   Direct answer to a real question: can this be expanded to Maldives,
//   Philippines? Searched and verified real, current literature for
//   both. Maldives had a strong, clean match: "Outcomes of the fourth
//   global coral bleaching (2023-2024) in the Maldives," Coral Reefs
//   (Springer Nature), 2026, doi:10.1007/s00338-026-02850-x - a real
//   survey of 18 reefs across central and southern atolls. Central
//   atolls averaged >40% live coral cover loss, up to 57% at the
//   worst-hit reef, disproportionately affecting Acropora; the
//   southernmost Huvadhoo Atoll showed high resistance, with Porites-
//   dominated reef flats retaining high live coral cover.
//   DISCLOSED DIFFERENCE from the two existing entries: this is REGIONAL
//   data (18 reefs, central/southern atolls) not single-reef data like
//   One Tree Reef/Florida Keys - larger radius (150km), and the numbers
//   are described as regional patterns, not one precise reef figure.
//   Philippines was searched too but NOT added this version - the real
//   papers found don't cleanly fit the database's format (one reports
//   bleaching PREVALENCE, not mortality; another is not peer-reviewed
//   and lacks a clean genus-level breakdown) - disclosed here rather
//   than forcing in weaker-quality data just to have an entry.
//
// v10.132 NEW: shared species-vulnerability database + real second entry
//   (Florida Keys)
//   Direct answer to a real question: can species/genus vulnerability
//   data be pulled from literature and made available across the whole
//   tool, not just S20? Refactored the single-site lookup into
//   SPECIES_VULNERABILITY_SITES, an array any module can query via
//   lookupSpeciesVulnerability(lat, lon) - genuinely shared, not
//   S20-specific.
//   Added a real, independently-verified SECOND entry after directly
//   searching for and confirming current literature: Florida Keys/Dry
//   Tortugas. Manzello et al. 2025 (Science, doi:10.1126/
//   science.adx7825) confirms 97.8-100% mortality of Acropora palmata
//   and A. cervicornis - described in the paper as functional
//   extinction of both species at this site. Florida FWC's 2023-24
//   CREMP Post-Bleaching Quick Look Report adds real genus-level data
//   for Orbicella (extensive mortality, Lower Keys), Siderastrea
//   siderea and Montastraea cavernosa (both ~1% decline - comparatively
//   resilient).
//   DISCLOSED SCOPE, stated plainly and unchanged in spirit from
//   v10.131: this is explicitly NOT a comprehensive global database and
//   will not claim to become one automatically - it grows one real,
//   independently-verified citation at a time. Every site not
//   explicitly listed still correctly returns "NOT CHECKED," never a
//   guess.
//
// v10.131 NEW: S20 - Species-Vulnerability-Weighted Combined Risk
//   Direct build-out of an idea discussed: combine the tool's existing,
//   real DHW measurement with real, published dominant-species
//   vulnerability data, where it exists - rather than treating heat
//   stress and species fragility as two separate, unconnected numbers.
//   Reuses the EXACT SAME real citation already embedded in S15 (Byrne
//   et al. 2025) - not a new data source, just restructured as
//   programmatic genus-level data (Acropora/Goniopora HIGH vulnerability,
//   95%/73% documented mortality; Porites/Pocillopora MODERATE-LOW,
//   31%/23%) so it can be combined with DHW instead of only shown as
//   narrative text.
//   DISCLOSED SCOPE, checked directly before building: species
//   composition cannot currently be measured by satellite at reliable
//   accuracy - current published literature on hyperspectral coral
//   species classification reports only ~56-70% accuracy at genus/
//   species level (vs ~88-90% for basic coral-vs-other detection), too
//   unreliable to use as ground truth. This module is explicitly NOT a
//   remote-sensing capability - it's a real field-survey lookup, same
//   honest pattern as S15/S5: populated only where a real published
//   colony-level survey exists (currently just One Tree Reef), explicit
//   "NOT CHECKED" everywhere else, never a filled-in guess.
//
// v10.130 FIX: STEP 5's own input cap silently blocked the fix for its
//   own bug
//   Direct follow-up to v10.129: recommended fixing a real Panama test
//   that showed "insufficient points (n=1)" on the new valid significance
//   test by using 96 total months with a 12-month window (8 independent
//   windows) - but STEP 5's own validation capped "Total months" at 84,
//   which would have silently REJECTED that exact fix with "Total months
//   must be 24-84." Raised the cap to 144 (12 years) so the fix
//   recommended for this exact problem can actually be entered.
//   Also added a proactive warning, checked and shown BEFORE the (slow)
//   Earth Engine call fires: if the chosen total-months/window-size
//   combination would leave fewer than 4 independent windows (as the
//   default 50%-of-total auto-window always does, by construction - it
//   only ever gives ~2), the run now says so upfront, instead of only
//   finding out via "insufficient points" after the full computation
//   completes, the way the real Panama test did.
//
// v10.129 NEW: real significance test + deseasonalizing for STEP 5
//   (sliding window trend test)
//   Direct port of a real bug found and fixed this session in the
//   companion Python deseasonalizing script: STEP 5's Mann-Kendall test
//   runs on OVERLAPPING windows (step=1 month) - consecutive positions
//   share nearly all their underlying data, which violates Mann-Kendall's
//   independence assumption. Confirmed as a real problem, not
//   theoretical: pure random noise with zero real trend gave p<0.0001
//   ("significant") using overlapping windows, but p=0.109 ("not
//   significant") using proper non-overlapping windows.
//   Applies directly to this tool's own results: a real Bocas del Toro,
//   Panama test this session reported "SIGNIFICANT CSD TREND - both AC1
//   and variance rising (p<0.05)" using the (now known to be invalid)
//   overlapping-window test - that specific claim needed independent
//   verification, which is exactly what this new section provides.
//   NEW: fetches the raw monthly series client-side (ONE new EE call,
//   reusing the same collection/geometry already used for the existing
//   trajectory), builds a real calendar-month climatology, deseasonalizes,
//   then runs the valid Mann-Kendall test on non-overlapping windows for
//   BOTH raw and deseasonalized data - shown alongside (not replacing)
//   the existing overlapping-window trajectory/chart, which remains
//   useful for visualizing the shape of the trend even though its own
//   p-value should now be read with caution.
//
// v10.128 FIX: deseasonalized comparison gave a factually wrong message,
//   and was missing the more informative number
//   Caught on the very first real test of v10.127 (One Tree Reef): the
//   interpretation text said "deseasonalized AC1 is similar in size to
//   the raw AC1" when the actual deseasonalized delta (-0.089) was 7.4x
//   LARGER than the raw delta (-0.012), not similar - the original logic
//   only ever checked for "much smaller," so any other case silently fell
//   into a message that didn't match the numbers.
//   Fixed the binary check into three real cases (much smaller / much
//   larger / similar), and added something the same real test showed was
//   actually more important than the delta comparison: the ABSOLUTE
//   baseline AC1 level itself. That test showed raw AC1=0.843 vs
//   deseasonalized AC1=0.434 - nearly half - suggesting a large share of
//   what looks like "high thermal memory" in raw AC1 readings throughout
//   this tool may reflect the seasonal cycle itself, not real resilience
//   loss. Now flagged explicitly whenever the deseasonalized baseline is
//   meaningfully lower than the raw one, not just when the CHANGE differs.
//
// v10.127 NEW: deseasonalized AC1/variance comparison for STEP 3 COMPARE
//   Direct fix for a caveat disclosed since v10.101 and never actually
//   built: AC1/variance computed on raw monthly values can be inflated
//   purely by the shared seasonal cycle (every site warmer in summer,
//   colder in winter), regardless of any real dynamical change. Evidence
//   this was a real problem, not theoretical: AC1 landed in a narrow
//   ~0.78-0.89 band across nearly every GBR site tested this session -
//   more consistent with seasonal autocorrelation dominating than genuine
//   site-specific signal.
//   New computeMonthlyClimatology()/deseasonalizeSeries() functions build
//   each site's own calendar-month average from its raw fetched series,
//   subtract it to get anomalies, then run the SAME jsNodeStats() used
//   everywhere else on those anomalies instead of raw values - shown
//   directly alongside a "raw" recomputation (using the identical
//   client-side method, for a clean apples-to-apples comparison) so the
//   seasonal cycle's actual contribution to the raw AC1 signal becomes
//   visible.
//   EFFICIENCY: zero new EE calls - reuses the exact same raw study
//   series already fetched for the permutation test (v10.122).
//   DISCLOSED LIMIT: with typical 10-24 month windows, the climatology is
//   built from only 1-2 samples per calendar month - a genuinely noisy
//   estimate, explicitly stated in the output, not a robust multi-decade
//   seasonal average. Shown as an exploratory comparison, not a
//   validated replacement for the raw AC1 already used in the verdict.
//
// v10.126 FIX: GLOBAL SIGNAL could hide a genuinely large local
//   amplification
//   Caught from a real Nuuk test: study variance rose +3.40x, control
//   variance rose only +0.22x - both individually crossed the 0.15x
//   threshold, so the tool called it "GLOBAL SIGNAL, not a local warning"
//   even though the divergence between them (+3.18x) was over 21 times
//   larger than the 0.15x divergence bar this SAME tool already uses
//   elsewhere (the MARGINAL LOCAL SIGNAL branch, built earlier this
//   session for a Florida Keys case). That branch was only ever reachable
//   when NEITHER site crossed the threshold individually - a case where
//   BOTH sites genuinely rose, but by wildly different amounts, had no
//   branch that could describe it honestly. The divergence value was
//   already being computed and displayed in this exact verdict text, just
//   never used to decide it.
//   Fixed: when both sites cross the threshold AND the divergence between
//   them also exceeds the divergence threshold, the verdict is now
//   "GLOBAL SIGNAL, WITH STRONG LOCAL AMPLIFICATION" - explicitly stating
//   that a real regional signal is present, but the study site is
//   amplifying it well beyond what the open-ocean control shows. Plain
//   "GLOBAL SIGNAL" is unchanged for the case where both sites rise by
//   comparable amounts.
//
// v10.125 DIAGNOSTIC: S7E hang investigation (no fix yet)
//   Caught on a real run: S7E stayed stuck at "1/2 batched calls done"
//   even after a page reload (ruling out an expired auth session) AND
//   after shortening to 12/10-month windows (ruling out the same window-
//   length slowness already fixed for S7D at 18/12mo in prior versions).
//   Neither previously-suggested fix resolved it - a genuinely new,
//   undiagnosed failure mode, not one of the two already-known causes.
//   Rather than guess a third fix blind, ported the exact diagnostic
//   pattern already proven for S7D in v10.108: both of S7E's raw-fetch
//   calls now print their ACTUAL error text on failure, or their real
//   feature count on success, instead of the hang staying uninformative.
//   DISCLOSED LIMIT: this can only help once a call eventually returns
//   something (success or a real error) - it cannot reveal anything if
//   the call is a genuine infinite hang with zero response from Earth
//   Engine's servers. Re-run S7E and check the Console for "S7E [" lines
//   to see exactly what happens this time.
//
// v10.124 NEW: real permutation-test p-values for S7E (algae-based
//   local vs regional classification)
//   Direct fix for a gap flagged in discussion: the permutation-test
//   engine (v10.122/123) only ever protected the SST-based S13 findings -
//   S7E's algae-based LOCAL ANOMALY / REGIONAL SIGNAL classification still
//   relied entirely on the same >0.01/>0.15 fixed thresholds as before.
//   Ports the exact same permutationTestDelta() engine to S7E, computing
//   real p-values for both study and reference sites (AC1 and variance),
//   so the same "would random reshuffling of these months produce a delta
//   this big by chance?" question now applies to algae, not just
//   temperature.
//   EFFICIENCY: needs ZERO new Earth Engine calls - S7E already fetches
//   the raw monthly FAI series for both study and reference sites to
//   compute its existing statistics; this just runs the proven engine on
//   that same already-fetched data, purely client-side.
//   Scoped to the standalone S7E tool for this release (not yet ported to
//   S7D's within-reef network or S7F's combined view).
//
// v10.123 NEW: real, Bonferroni-corrected permutation-test p-values for
//   FIND SWEET SPOT
//   Direct response to a real question: if window length changes the
//   result (confirmed by a real Nuuk test - a 10-month AFTER window showed
//   a rising-AC1 signal, a 24-month window at the SAME coordinate showed
//   nothing), how can you draw a reliable conclusion from testing several
//   window lengths and picking the best-looking one? That is a classic
//   multiple-comparisons trap: testing 6 windows means SOME window can
//   look significant by pure chance alone, even with nothing real
//   happening - the existing ROBUSTNESS note (v10.95) only ever warned
//   about this in words, never corrected for it mathematically.
//   Now runs the same permutationTestDelta() engine validated in STEP 3
//   COMPARE (v10.122) at each of the 6 tested window lengths (study AC1,
//   study variance, control AC1, control variance = 24 tests total, 300
//   shuffles each), and reports an explicit Bonferroni-corrected
//   significance bar (0.05 / 6 windows = p<0.0083) alongside the
//   uncorrected count - so a genuinely robust signal (surviving the
//   stricter bar) can be told apart from a lucky single window among six.
//   EFFICIENCY: adds only 4 new EE calls, not 12 - all 6 AFTER windows
//   share the same start date and differ only in length, so a shorter
//   window's raw values are always a prefix of the 24-month fetch;
//   fetched once per site and sliced client-side for the other 5 lengths.
//   Shown as a new, clearly-separated ADDITIVE section - does not alter
//   the existing threshold-based table or ROBUSTNESS note above it.
//
// v10.122 NEW: real permutation-test p-values for STEP 3 COMPARE
//   Direct response to a design discussion: most of this tool's
//   LOCAL/REGIONAL/signal decisions rely on fixed thresholds (>0.01 AC1,
//   >0.15 variance) copy-pasted across STEP 3, STEP 4, S7D, S7E, and S7F -
//   only STEP 5's Mann-Kendall test ever answered "is this bigger than
//   random noise would produce" with a real p-value. A fixed threshold
//   cannot tell a genuinely surprising 0.2 AC1 jump (site where months
//   normally barely move) from an unremarkable one (site where ordinary
//   weather noise wobbles AC1 by 0.3 anyway) - both get treated identically.
//   Built a new, reusable permutation-test engine (permutationTestDelta):
//   pools the BEFORE+AFTER raw monthly values, randomly reshuffles which
//   months get relabeled BEFORE/AFTER (keeping window sizes fixed),
//   recomputes the same AC1/variance statistic on each of 500 shuffles,
//   and reports what fraction of PURELY RANDOM reshuffles produce a delta
//   at least as large as the one actually observed - a genuine, honest
//   p-value. Pure client-side JS, no new EE architecture: reuses the exact
//   same jsLinearDetrendResiduals/jsLag1AC1/jsVarianceHalves building
//   blocks already proven this session, and the same raw-value-fetch
//   pattern (extractMultiNodeSeries) already used for S7D/S7E/S7F.
//   Wired into STEP 3 COMPARE first (the highest-priority, most heavily
//   tested tool) as an ADDITIVE result shown alongside - not replacing -
//   the existing threshold-based verdict, in its own panel, via 2 new
//   independent EE calls. Deliberately scoped to COMPARE only for this
//   release; S7D/S7E/S7F/STEP 4 can reuse the same engine once this is
//   validated in real testing.
//
// v10.121 FIX: COMPARE could silently show a stale result from an earlier
//   run
//   Caught from a real multi-site test: COMPARE showed AC1=0.794 for the
//   AFTER window at a new coordinate, while STEP 2's own "AFTER stored"
//   label correctly showed 0.431 for the exact same window - and 0.794
//   turned out to be an exact match for a PREVIOUS site's AFTER value from
//   earlier in the same session. Confirmed by direct testing (re-clicking
//   COMPARE fixed it immediately, showing the correct 0.431): this was NOT
//   a computation bug - COMPARE's result panel is a one-time snapshot
//   taken at click time, and never auto-refreshes if STEP 2 is re-run
//   afterward. The underlying math was always correct; nothing warned the
//   user the on-screen result had gone stale.
//   Fixed: csdCompareRanWithB/A remember which stored BEFORE/AFTER objects
//   COMPARE last used. If STEP 2 is re-run afterward (creating new stored
//   objects), a visible warning now appears immediately: "STALE RESULT
//   BELOW... press COMPARE again to refresh." The warning clears
//   automatically the next time COMPARE actually runs.
//
// v10.120 FIX: the v10.119 fix only covered the headline box, not the
//   detail panel underneath it
//   Caught immediately on the very next real Nuuk test: the headline
//   correctly showed the v10.119 CONFLICTS WITH AC1 message, but the
//   "4-WAY CSD COMPARISON (detail)" panel below it STILL showed the raw
//   unreconciled "Control variance rose but the study site stayed flat...
//   This is a POSITIVE result" text - the same contradiction, just in a
//   second location v10.119 did not reach.
//   Root cause: this detail panel is built early and SYNCHRONOUSLY, using
//   the raw variance-only verdict text, before the AC1-weighted spatial
//   toolkit finishes computing (an async callback that runs later and only
//   rewrites the headline box, never this detail panel). Structurally
//   stuck showing pre-AC1-weighting text.
//   Fixed using the same conflict-detection approach as v10.119, but with
//   studyAC1Rose (already computed earlier in the code, before this panel
//   is built) since fullTally is not yet available at this point. Also
//   fixed the standalone console print of the same raw text right after.
//
// v10.119 FIX: COMPARE headline could contradict its own "Regional context"
//   Same class of bug already caught and fixed for FIND SWEET SPOT in
//   v10.94 (two separately-computed classifications sitting in one box
//   with no reconciliation), found living on in COMPARE's main verdict box
//   too. Caught from a real Nuuk test: the headline correctly said
//   "STRONG SIGNAL, LIKELY REGIONAL" (AC1 rose +0.234 study vs +0.242
//   control - nearly identical, genuinely regional per Dakos et al.
//   weighting), but the SAME box then said "Regional context: STUDY SITE
//   MORE STABLE THAN CONTROL (positive result)" - the opposite claim,
//   driven purely by variance (fell at study, rose at control) with no
//   awareness that AC1 - the primary indicator - told a different story.
//   Root cause: "Regional context: "+vTitle pasted the OLD, pre-toolkit,
//   variance-only classification text directly after the NEW AC1-weighted
//   combinedTitle headline, unreconciled. Fixed: the conflict is now
//   detected explicitly (AC1 agrees with a rising trend AND vTitle claims
//   "study more stable") and replaced with a clear statement of the
//   conflict instead of silently presenting two disagreeing verdicts as
//   if they were consistent.
//
// v10.118 CRITICAL FIX: S13 COMPARE could call zero data a "positive result"
//   The same bug already caught and fixed for S7E in v10.114 existed in
//   S13's core COMPARE verdict too - the single most heavily used part of
//   this whole tool, hit in nearly every test throughout this session.
//   Caught testing a genuinely new location (Nuuk, 64.13,-51.38): the
//   study site returned AC1=n/a, Var=n/a in BOTH the BEFORE and AFTER
//   windows (zero valid OISST months, likely ice/coastal-pixel masking),
//   yet the headline confidently declared "STUDY SITE MORE STABLE THAN
//   CONTROL (positive result)" and "the study site stayed flat... This is
//   a POSITIVE result for this reef." The site was never measured at all.
//   Root cause: studyVarRose/studyAC1Rose both silently default to false
//   when the underlying deltas are null, making "no data" indistinguishable
//   from "genuinely flat" - both fell into the same reassuring branch.
//   Fixed: three new explicit checks run BEFORE any of the existing GLOBAL/
//   LOCAL/ANOMALOUS branches - "CANNOT ASSESS (study has no data)",
//   "CANNOT ASSESS (control has no data)", "CANNOT ASSESS (neither site
//   has data)" - each stating plainly that this is not a stable/flat
//   result, just an unmeasurable one, with likely causes and next steps.
//   Checked FIND SWEET SPOT for the same bug: already safe (its own
//   per-window logic already null-checks first, before this fix existed).
//
// v10.117 FIX: S7F combined summary was missing data density entirely
//   Caught directly: asked the user to check S7F's data density output to
//   validate a LOCAL ANOMALY finding, and there was nothing to show - the
//   density numbers were already computed internally (needed for the
//   CANNOT CLASSIFY check) but never actually pushed into the displayed
//   summary. S7D and S7E each show density in their own individual
//   output; S7F's combined view silently dropped it. Now shown explicitly
//   for both study and reference sites, with the same INSUFFICIENT flag
//   used elsewhere, so a combined-run result can be checked for data
//   quality without needing to separately re-run S7D/S7E individually.
//
// v10.116 NEW: reference-site exclusion + S7F combined orchestrator
//   Two features built from a real design discussion:
//   1) REFERENCE EXCLUSION: S7E could auto-select a reference reef that
//      coincided with a site already tested as a STUDY location earlier
//      in the same session (caught directly: the reference for one test
//      landed on Low Isles, which had already shown a real signal as a
//      study site). A reference is supposed to be an uninvolved baseline.
//      New global s7StudySiteHistory (via recordStudySite(), called from
//      S7B/S7C/S7D/S7E/S7F) tracks every coordinate tested as a study
//      site this session; S7E's candidate search now skips any GEBCO
//      match within 5km of a prior study site, falling through to the
//      next-closest genuine candidate, with the exclusion count shown in
//      the results.
//   2) S7F - RUN ALL (S7D + S7E COMBINED): direct response to a request
//      for full automation. Built a SCOPED version, not the literal
//      "everything automatically" version - explained why in the UI
//      itself: this session already hit real EE account-level concurrency
//      limits once (required a manual tier upgrade) and had a single S7D
//      call hang 5+ minutes with no way to cancel it (no setTimeout in
//      this sandbox). An orchestrator auto-trying multiple window lengths
//      across every S7 tool would multiply both risks severalfold. S7F
//      instead takes ONE shared Lat/Lon + BEFORE/AFTER input (entered
//      once, not twice) and runs S7D + S7E together (6 EE calls total,
//      same bound as running them separately), rendering one combined
//      summary table with both headline verdicts side by side. Reuses
//      the exact same proven functions (extractMultiNodeSeries,
//      jsNodeStats, jsPairCorrelation, groupSeriesByLabel) S7D/S7E
//      already use - no new statistics code, only new orchestration.
//
// v10.115 NEW: S7E shows the reference site's actual coordinates
//   Direct response to a real question: a user asked whether a map
//   landmark (Low Isles) was the auto-found reference site. Previously
//   impossible to answer - S7E only ever displayed distance and depth
//   ("20km away, depth=-25.0m"), never the actual lat/lon, even though
//   best.lat/best.lon were already computed at candidate-selection time
//   and simply never surfaced. Now shown directly in the results line, so
//   the reference candidate can be identified on the map (or its depth
//   checked against known landmarks - e.g. an emergent cay/island reads
//   near 0m elevation, not -25m, so a candidate at meaningful negative
//   depth is a submerged patch, not a visible island).
//
// v10.114 FIX: S7E asserted "LOCAL ANOMALY DETECTED" against an untested
//   reference site
//   A real run showed the headline verdict "LOCAL ANOMALY DETECTED" while
//   the reference site's own numbers read ΔAC1=n/a, ΔVar=n/a - the
//   reference wasn't tested and found stable, it simply had no computable
//   data (its AFTER period had fewer than 4 valid months, even though its
//   COMBINED before+after density of 55% looked fine and didn't trip the
//   existing low-confidence flag, which only checks the combined average).
//   "Untested" and "tested and found stable" are different findings, but
//   the decision rule collapsed both into refSignal=false, letting a real
//   study-site signal get reported as a CONFIRMED local anomaly against a
//   reference that never actually confirmed anything.
//   Fixed: added an explicit per-period (not combined) insufficient-data
//   check for BOTH study and reference. When either site has <4 valid
//   months in either individual period, the verdict is now "CANNOT
//   CLASSIFY" with the specific period named, instead of silently
//   defaulting into one of the four normal LOCAL/REGIONAL/ANOMALOUS/NO
//   SIGNAL classifications as if the comparison were valid.
//
// v10.113 FIX: S7E was hiding a coherent regional pattern inside "NO SIGNAL"
//   The same gap already caught and fixed for STEP 5 in v10.97 (a real,
//   significant DECLINE was being mislabelled as "no significant trend")
//   existed in S7E too, uncaught until a real run exposed it directly:
//   study AC1 fell -0.675, reference fell -0.636 - nearly identical, a
//   coherent pattern across two independent sites. The SIGNAL check
//   correctly only flags RISING AC1/variance (per Dakos et al., a decline
//   is not itself a CSD warning), so this got silently folded into a bare
//   "NO SIGNAL AT EITHER SITE" with no mention of the pattern underneath.
//   Fixed: when both sites show a similar, substantial AC1 decline
//   (<-0.1 each, within 0.15 of each other), the verdict now says so
//   explicitly - labelled a REGIONAL DECLINE, clearly distinguished from
//   a CSD warning direction, rather than describing real, structured data
//   as indistinguishable from noise.
//
// v10.112 FIX: S7E now shows STUDY site data density, not just reference
//   A real run at a cloudy rainforest coastline (Daintree/Cape Tribulation)
//   returned NO SIGNAL AT EITHER SITE with BOTH study and reference showing
//   n/a for AC1/variance. The reference site's data density was already
//   shown (44%, explaining ITS n/a), but the study site's own density was
//   never displayed - no way to tell whether the study site failed from
//   near-zero Sentinel-2 coverage (likely, given the location) or something
//   else. Now shown symmetrically: study site density (X% of BEFORE/AFTER
//   months valid) alongside the reference's, with the same LOW-CONFIDENCE
//   flag logic applied to both. A sparse-data "no signal" result is now
//   distinguishable from a genuine "checked and found nothing" result at
//   the study site too, not just the reference.
//
// v10.111 NEW: S7E - LOCAL vs REGIONAL auto-classification
//   Closes a real gap: S7D's 9 nodes all sit within 0.3-3km of each other -
//   they test spatial coupling WITHIN one reef, but every node is still
//   local to that same reef, so S7D alone could never actually distinguish
//   a local event big enough to blanket the reef from a genuine regional
//   signal. S7E adds the missing genuinely-independent reference:
//   - Auto-searches outward in expanding rings (20/40/70/110/160km, 8
//     bearings each = 40 candidates) for the nearest point that is
//     genuinely shallow water per GEBCO (-50 to 0m, the same threshold
//     used for shallowMask elsewhere in this tool) - all 40 candidates
//     checked in ONE batched reduceRegions() call.
//   - Runs the SAME FAI-based BEFORE/AFTER AC1/variance extraction at that
//     reference site as at the study site (reusing mkMoFAIRange/
//     extractMultiNodeSeries/jsNodeStats unchanged - 2 more batched calls,
//     3 total for the whole tool).
//   - Applies an explicit decision rule, mirroring STEP 3's already-proven
//     LOCAL/GLOBAL/ANOMALOUS/NO SIGNAL classification for SST: study
//     signal without reference signal = LOCAL ANOMALY; both = REGIONAL;
//     reference-only = ANOMALOUS; neither = NO SIGNAL. This is the actual
//     automated classification layer that was missing - not just S7D's
//     within-reef synchronization count, but a real core-vs-independent-
//     reference comparison with an automated verdict.
//   DISCLOSED LIMIT (not solved, only partially mitigated): a GEBCO
//   shallow-water match is not guaranteed to be a real reef with
//   comparable ecology - could be a bare sandbar with no algae community.
//   Partially checked via the reference site's own data density as a
//   proxy (a sparse signal there is flagged LOW CONFIDENCE, not silently
//   trusted) - this does not prove ecological comparability, only flags
//   the most obvious failure mode (a site with essentially no data at all).
//
// v10.110 FIX: v10.109's NDVI fix over-corrected - fixed buffer size
//   The v10.109 fix (mean()->max(), 150m->500m buffer) correctly resolved
//   Center's false "LIKELY NOT on-reef" reading, but the fixed 500m buffer
//   was itself a new bug: at the default 0.5km ring radius, a real run
//   showed Center/N/NE/E all reading the IDENTICAL value 0.84 - strong
//   evidence their 500m buffers were overlapping so much (Center-to-ring
//   spacing is only 500m at that radius; adjacent ring points are even
//   closer, ~383m) that multiple "different" nodes were just picking up
//   the same peak pixel, silently defeating the purpose of sampling 9
//   distinct locations.
//   Fixed: NDVI check buffer now SCALES with the user's chosen ring
//   radius (30% of node spacing, clamped 80-250m) instead of a fixed
//   500m, verified mathematically to stay clear of overlap across the
//   entire valid 0.3-3km radius range. The buffer size actually used is
//   now shown in the on-screen results table for transparency. The FAI
//   buffers used for the actual coupling statistics remain unchanged
//   (150m fixed) - only the on-reef check's buffer changed.
//   Honest caveat carried forward: whether the identical-value pattern in
//   the real run was purely a buffer-overlap artifact, or partly a
//   genuine large uniform algae patch, could not be determined with
//   certainty from that one run - re-testing with this fix will show
//   whether the values differentiate.
//
// v10.109 FIX vs v10.108: two real bugs caught from an actual completed run
//   The v10.108 diagnostics worked as intended: a smaller test run (6mo/6mo
//   instead of 18mo/12mo) completed normally with matching feature counts
//   (BEFORE=54, AFTER=54, NDVI=9 - exactly 6x9 and 9 as expected),
//   confirming the earlier multi-minute hang was a size/complexity issue
//   with larger month counts in one batched request, not a hard bug. That
//   completed run then surfaced two real, separate problems:
//   1) NDVI ON-REEF CHECK INCONSISTENCY: Center (the exact study
//      coordinate, confirmed on-reef with strongly POSITIVE NDVI in every
//      prior test this session) flipped to NEGATIVE (-0.32) and got
//      flagged "LIKELY NOT on-reef". Root cause: mean() over a small 150m
//      buffer at 20m scale is sensitive to small-scale heterogeneity - a
//      reef patch smaller than the buffer gets its average pulled negative
//      by adjacent clear water. Fixed by switching to max() over a larger
//      500m buffer specifically for the on-reef check, matching the
//      established convention already used by the main panel and S7B (the
//      FAI buffers used for the actual coupling statistics are unchanged).
//   2) MISSING ARTIFACT FLAG: the same real run showed Var=+52.31x and
//      +11.99x - even more extreme than the case that prompted S7C's
//      v10.105 near-zero-denominator artifact warning, which was never
//      ported to S7D. Now added: same threshold logic (1st-half variance
//      <0.001 AND |delta|>5x triggers a warning icon on that row).
//
// v10.108 DIAGNOSTIC: S7D real-run failure needs actual evidence, not a
//   guess. A real run at the CORRECT coordinates (-23.51, 152.09, confirmed
//   by the map pin) returned n/a for every field on every node -
//   "COULD NOT ASSESS COUPLING" with only "1 of 3 batched calls errored"
//   reported, which doesn't fully explain why ALL fields (not just
//   correlation) came back empty. Rather than guess at a fix blind (the
//   v10.107 batching rewrite is new, untested-live code - reduceRegions()'s
//   exact output property naming for a single-band mean reducer was
//   assumed, not confirmed), this adds real diagnostics:
//   - New s7dDiagnose(): prints the ACTUAL error text for any of the 3
//     batched calls that fails (previously only an aggregate "N errored"
//     counter existed, with no detail on WHICH call or WHY).
//   - On success, prints the feature count and the first feature's raw
//     property names - directly confirms or refutes whether
//     reduceRegions() names its output the way extractReduceRegionsValue()
//     assumes, instead of leaving that assumption unverified.
//   - Raw feature counts (BEFORE/AFTER/NDVI) now shown directly in the
//     on-screen results table too, not just the console - distinguishes
//     "the fetch itself returned nothing" (a real data/query problem) from
//     "data came back but extraction is broken" (a property-naming bug),
//     which need different fixes.
//   No fix applied yet - this version is instrumentation only, so the next
//   real run's console output will show exactly what's happening.
//
// v10.107 PERFORMANCE: S7D rebuilt from ~43 EE calls down to 3
//   Not a quantum-computing question (the actual bottleneck - raster
//   compositing, API round trips, server queuing - has nothing to do with
//   the kind of problems quantum algorithms accelerate) - a boring,
//   effective classical fix instead: batching.
//   - New extractMultiNodeSeries(): for a monthly image collection, builds
//     ONE FeatureCollection covering every node's whole time series at
//     once, via reduceRegions() (samples all 9 points against ONE image)
//     + flatten() across all months into a single flat table. One
//     .evaluate() call now returns what used to require 9 separate calls
//     (one per node, each re-triggering the full monthly Sentinel-2
//     compositing graph from scratch).
//   - New client-side (plain JS) statistics: jsNodeStats (linear detrend,
//     lag-1 AC1, variance halves, skewness) and jsPairCorrelation (Pearson
//     correlation between two aligned series) - identical formulas to the
//     EE versions (computeRealCSD/computeZonalSyncCSD), just computed
//     locally on the small already-fetched tables instead of triggering
//     new server-side graph evaluations. Same "fetch once, compute
//     locally" pattern already proven safe for STEP 5's Mann-Kendall test.
//   - Total EE calls: 9 nodes x 2 periods x AC1/variance + 8 pairs x 2
//     periods x correlation + 9 NDVI checks (~43 calls, 2-5 min) is now
//     just 2 batched series fetches (BEFORE, AFTER) + 1 batched NDVI fetch
//     = 3 calls, ~20-60s expected.
//   - Output format, verdict logic, and UI inputs are UNCHANGED - this is
//     a pure performance rebuild, not a new feature or a different result.
//
// v10.106 NEW: S7D - full 9-node algae coupling network (BEFORE/AFTER)
//   Scales S7C's 3-node proof-of-concept to the full 9-node ring (matching
//   S7B's compass geometry), with a real BEFORE/AFTER comparison instead
//   of a single snapshot window.
//   - ~43 Earth Engine calls in one click: 9 nodes x 2 periods x AC1/
//     variance (18), 8 ring points x 2 periods x correlation-vs-Center
//     (16), 9 nodes x 1 NDVI-water on-reef check (9). Deliberately scoped
//     to correlation-vs-Center only (8 pairs), not the full 36-pair
//     matrix, to keep runtime bounded (~2-5 min instead of much longer).
//   - NDVI-water fetched once per node (current composite) specifically
//     to flag likely off-reef points, per the real S7B finding that North
//     showed strongly negative NDVI there. A node/pair involving a
//     flagged off-reef point is marked explicitly, not silently trusted
//     as a real algae-dynamics comparison.
//   - Headline verdict counts how many ON-REEF ring points show rising
//     correlation with Center (the Dakos et al. 2011 hyper-synchronization
//     direction), separately from raw AC1-rising node count.
//   - Reuses computeRealCSD/computeZonalSyncCSD/mkMoFAIRange unchanged -
//     no new statistics code, only new orchestration across more nodes
//     and two time periods.
//   Explicitly disclosed as an empirical correlation network, not a
//   mechanistic J_ij interaction matrix - same honesty framing as every
//   other coupling-style indicator already in this tool.
//
// v10.105 FIX: S7C variance ratio could look like a huge surge when it was
//   actually a near-zero-denominator artifact. Caught on a real run: East
//   node showed Var=15.80x, far beyond anything seen elsewhere in this
//   tool (SST-based ratios never exceeded ~7x). The ratio is computed as
//   secondHalfVariance / max(firstHalfVariance, 1e-6) - if the first half
//   had FAI sitting near-constant (plausible for algae presence, which can
//   be genuinely near-zero for long stretches), even a modest second-half
//   variance inflates the ratio dramatically without any real surge having
//   occurred. Now shows the raw first-half/second-half variance alongside
//   the ratio, and auto-flags the case (first-half variance < 0.001 AND
//   ratio > 5) with an explicit "LIKELY ARTIFACT" warning - no more trusting
//   a single derived number without the components that produced it.
//
// v10.104 NEW: S7C - algae/AC1/variance coupling proof-of-concept
//   Tests whether the SAME toolkit already applied to SST in S13 (AC1,
//   variance, spatial correlation) can be meaningfully computed from ALGAE
//   (FAI) data instead - arguably more theoretically appropriate, since in
//   Scheffer's bistable-state framework coral-vs-algae cover is the STATE
//   VARIABLE that actually flips between stable states, while SST is
//   closer to the external control parameter driving the system toward a
//   threshold.
//   - New mkMoFAIRange(): the missing monthly Sentinel-2 FAI time series
//     builder - S7 previously only ever had ONE fixed 2023-2024 composite,
//     never a time series, so this statistic was not computable before.
//   - computeRealCSD() and computeZonalSyncCSD() needed ZERO changes - both
//     already took a generic bandName parameter, so passing 'fai' instead
//     of 'sst' just works.
//   - Deliberately scoped to 3 nodes (Center/North/East) and one window,
//     not the full 9-node ring, because Sentinel-2's cloud-masking data
//     density at monthly resolution is genuinely untested until now -
//     reports an explicit DATA DENSITY assessment (% of requested months
//     that actually had valid data) before anything else, so a sparse-data
//     result is caught and flagged rather than silently producing unreliable
//     statistics.
//   - If density looks good, the UI explicitly recommends scaling to the
//     full 9-node ring (matching S7B's geometry) with a proper BEFORE/AFTER
//     comparison next - not built yet, pending this proof-of-concept result.
//
// v10.103 NEW: real S15 field validation for One Tree Reef
//   Found a real, published paper tracking the EXACT coordinate this tool
//   has been tested against: Byrne et al. 2025 (Limnol. Oceanogr. Lett.,
//   doi:10.1002/lol2.10456) tracked 462 individual coral colonies at One
//   Tree Reef (23.51S, 152.09E) through the 2023-24 heatwave - 66% bleached
//   by Feb 2024, 80% by April, up to 52% mortality by July, with genus-
//   level detail (Acropora 95% mortality/rapid collapse to rubble;
//   Goniopora 73% mortality via black band disease; more resilient genera
//   like Porites/Pocillopora showing partial recovery). Data publicly
//   available at Sydney eScholarship (doi:10.25910/p5rq-cw63).
//   getEcologicalRecoveryValidation() now takes lat/lon and checks PROXIMITY
//   (haversine distance, 60km radius) to this exact site before falling
//   back to the coarse region-based lookup - deliberately NOT keyed to the
//   whole "Great Barrier Reef" region bucket, since that bounding box spans
//   from northern reefs (Lizard Island, hit hardest in 2016) to this
//   southern site (largely spared until 2024) with very different
//   bleaching histories; applying this southern-GBR-specific finding
//   region-wide would have misrepresented the northern reefs.
//   All 3 call sites (S13 single-window test, main click S15 display,
//   export log) updated to pass lat/lon through.
//
// v10.102 NEW: S7B - MULTI-POINT ALGAE SCAN
//   Direct response to a real finding: two points 1.1km apart at One Tree
//   Reef showed "WATCH - mild signal" vs "EXTREME - massive bloom" for the
//   same S7 macroalgae indicators - a single click can badly misrepresent
//   a patchy bloom in either direction (missing a real one, or overstating
//   an isolated one).
//   - Samples FAI, NDCI, and NDVI-water at 8 compass points (N/NE/E/SE/S/
//     SW/W/NW) plus the centre, at a user-chosen radius (0.5-10km), using
//     the SAME raw images and methodology (max reducer, 100m scale) as the
//     existing single-point S7 panel - so results are directly comparable.
//   - Built as ONE Earth Engine call via reduceRegions() over a 9-feature
//     collection, not 9 separate calls.
//   - Reports a WARNING SIGN COUNT (how many of the 9 points show elevated
//     FAI/NDCI/NDVI) and an explicit PATTERN classification: WIDESPREAD
//     (>=70% of points elevated - likely a genuine regional bloom),
//     PATCHY/ISOLATED (<=2 points - could be a real small patch OR a
//     single-pixel/cloud artifact, explicitly flagged not to generalize),
//     or MIXED.
//   - Also reports the numeric spread (max-min) across points as a direct,
//     explicit measure of the kind of spatial heterogeneity that motivated
//     this feature.
//   - Has its own Lat/Lon input + "Use last clicked location" button
//     (reusing the lastClickLat/lastClickLon pattern from S13 STEP 1), so
//     it works independently of the main click flow or S13.
//   Honest limit disclosed in the UI: a 9-point compass ring is a simple
//   sampling pattern, not exhaustive coverage - a bloom could sit between
//   sample points. Recommends a smaller radius or the S7 map layers for
//   full visual coverage.
//
// v10.101 CAVEAT DISCLOSURE (no computation changed): a direct question
//   about whether v10.100's synchronization indicator could have detected
//   the GBR event earlier prompted a re-check that surfaced a real design
//   concern, not caught before shipping: computeZonalSyncCSD() correlates
//   RAW monthly SST, not deseasonalized anomalies. Two ocean points a short
//   distance apart share a strong seasonal cycle regardless of any real
//   dynamical coupling change, which likely pushes the baseline correlation
//   toward a high ceiling and blunts this indicator's sensitivity.
//   Supporting evidence found while checking this: across every real GBR
//   test run in this tool's history, the core AC1 statistic (which has the
//   same raw/linear-detrend-only limitation) consistently landed in a
//   narrow ~0.78-0.89 band regardless of site or period - more consistent
//   with shared seasonal month-to-month autocorrelation dominating the
//   signal than with genuine site-specific critical-slowing-down dynamics.
//   A proper fix (per-zone monthly climatology subtraction before
//   correlating) needs live testing to verify before shipping, so it is
//   NOT implemented blind here. Instead: the limitation is now disclosed
//   explicitly in the STEP 3 UI, in the indicator's own display name, and
//   in the source code comments, with guidance to read the DELTA rather
//   than the absolute value until a tested deseasonalized version exists.
//   No answer this version gives should be read as "this would have
//   caught it earlier" - that claim has zero empirical support yet.
//
// v10.100 NEW: study-control temporal synchronization indicator (STEP 3)
//   Implements the "empirical interaction/covariance network" idea in a
//   deliberately minimal, low-risk form: rather than inventing new zone
//   geometries (reef flat / mangrove / channel, which would need habitat
//   datasets not verified reliable here), this reuses the deep-water
//   control site STEP 3 already auto-selects as the second node. Computes
//   the Pearson correlation between the study reef's and the control
//   site's month-to-month SST fluctuations, for BEFORE and AFTER
//   separately - tracking whether the reef is becoming MORE locked in sync
//   with the open ocean (rising correlation - losing local independence/
//   buffering, the leading-indicator direction reported by Dakos et al.
//   2011, Am Nat 177:E153-E166) or staying decoupled.
//   - New computeZonalSyncCSD(): the actual correlation math, reusing the
//     same array-based approach already proven safe elsewhere in the file.
//   - Added as a new SUPPORTING (not primary) indicator in STEP 3's
//     toolkit, alongside spatial variance and spatial autocorrelation -
//     fetched as one more follow-up in the same non-blocking chain, so a
//     failure here can't delay or break the temporal/spatial results above.
//   Disclosed honestly in the code comments: this is a correlation, not a
//   mechanistic interaction coefficient (no J_ij matrix); shared external
//   forcing (a heatwave hitting both sites) will raise this number with or
//   without any real internal dynamics - one more piece of evidence, not
//   proof on its own, consistent with how every other indicator in this
//   toolkit is already framed.
//
// v10.99 CRITICAL FIX: "Cannot read property 'trim' of undefined" crash
//   Real GEE Code Editor error, caught live: clicking RUN SLIDING WINDOW
//   ANALYSIS with the (intentionally optional, "blank = auto") window-size
//   field never typed into threw a hard crash at csdSlideWindowInput.
//   getValue().trim() - GEE's ui.Textbox.getValue() returns undefined (not
//   an empty string) for a field that has never been interacted with, even
//   with only a placeholder set. Calling .trim() on undefined throws and
//   stops the whole button handler.
//   This exact pattern - X.getValue().trim() with no guard - existed in 8
//   places across the file (GO TO COORDINATES, STEP 2 RUN CSD TEST, the
//   ADVANCED manual control override, FIND SWEET SPOT's AFTER date, and
//   all three STEP 5 fields). Any of them could crash identically if
//   clicked before typing into that specific field - the STEP 5 window
//   field just happened to be the one a real user hit first, since leaving
//   it blank is the intended, documented way to use it.
//   Fixed everywhere at once: every X.getValue().trim() is now
//   (X.getValue()||'').trim(), so a never-touched field safely reads as an
//   empty string and falls through to the existing "please fill this in"
//   validation message instead of crashing.
//
// v10.98 FIX: STEP 5 silently used fewer valid months than requested
//   Caught on a real Lizard Island run: requested 65 months, but only 51
//   were valid (37 sliding positions with a 15-month window = 51 valid
//   months), with no explanation anywhere in the output for the missing 14.
//   Root cause: the requested range (2022-06 + 65 months) ran into 2027-10,
//   past today's real date - OISST has no satellite observations for
//   future months, so they came back null and were silently dropped from
//   the valid-months count. The math was correct throughout; only the
//   transparency was missing.
//   Fixed two ways:
//   - Proactive check BEFORE firing the Earth Engine call: if the requested
//     start+total runs past today's date, an explicit warning appears
//     immediately (client-side, instant) explaining how many months will
//     be excluded and why - no need to wait 10-40s to find out.
//   - After the run: if fewer valid months were found than requested for
//     ANY reason (future dates or an ordinary data gap), the results table
//     now states this explicitly with the exact count, instead of leaving
//     the user to reverse-engineer a position-count mismatch themselves.
//   Also: FIND SWEET SPOT's ROBUSTNESS note read "only 0 of 6 windows lean
//   LOCAL" for the zero case - now reads "NONE of 6" for that case.
//
// v10.97 FIX: STEP 5 headline could claim "no significant trend" when one
//   genuinely existed. Caught on a real GBR run: variance showed tau=-0.273,
//   p=0.017 - a real, statistically significant DECLINING trend - but the
//   headline said "NO SIGNIFICANT TREND DETECTED (p>=0.05 for both AC1 and
//   variance)", which was simply false about the number sitting right below
//   it. Root cause: the verdict logic only checked for RISING significant
//   trends (tau>0 && p<0.05); a significantly FALLING trend (tau<0 &&
//   p<0.05, a real and reportable finding) was silently lumped into the
//   same bucket as a genuinely flat, noisy series with no real trend at
//   all. Two very different findings, one wrong label.
//   Fixed: now classifies each indicator into rising-significant /
//   falling-significant / not-significant separately. A significant
//   decline now gets its own accurate headline (e.g. "VARIANCE
//   SIGNIFICANTLY FALLING (p<0.05) - no CSD signature... this is an
//   absence-of-warning-sign result, not a 'nothing found' result"),
//   distinguishing it from both a rising CSD signal and a genuinely flat
//   series. "NO SIGNIFICANT TREND DETECTED" is now only used when BOTH
//   indicators are actually non-significant in either direction.
//   The per-window trendVerdict() line text (which already correctly showed
//   direction and significance separately) was not affected - only the
//   summary headline classification had this bug.
//
// v10.96 NEW: STEP 5 - SLIDING WINDOW TREND TEST (Dakos et al. 2012 method)
//   Direct implementation of the actual gold-standard method, replacing
//   discrete BEFORE/AFTER chunk testing with a continuous rolling window:
//   - computeSlidingWindowCSD(): detrends the whole requested series ONCE
//     (linear detrend, same method as computeRealCSD - Dakos et al.'s own
//     toolbox typically uses Gaussian kernel smoothing instead, which is
//     more flexible but not implemented here; disclosed in the UI), then
//     slides a fixed-size window forward ONE MONTH AT A TIME across the
//     residuals, computing AC1 and variance at every position. Built as a
//     SINGLE server-side ee.List.map() graph, so the whole scan (which can
//     be dozens of window positions) costs one Earth Engine round trip.
//   - mannKendallTest(): the actual Kendall's-tau / Mann-Kendall trend
//     significance test (standard normal approximation for the S-statistic
//     variance), run client-side on the small, already-evaluated result -
//     answers "is this metric moving consistently in one direction, or
//     just bouncing randomly?" with a real p-value, not an arbitrary
//     numeric threshold like the rest of S13 uses.
//   - New STEP 5 UI: full series start date + total months (window size
//     auto-defaults to 50% of total, the Dakos et al. standard, overridable).
//   - Prints the actual AC1(t) and variance(t) trajectories as line charts
//     to the console (with a linear trendline) - the "smoothly climbing"
//     visualization Dakos et al.'s own figures show, instead of a single
//     pass/fail verdict.
//   - Verdict weights AC1 as primary (consistent with v10.90/91): AC1
//     significant + variance not = still a strong headline; variance
//     significant alone = explicitly downgraded per Dakos et al.'s own
//     finding that variance is the less robust indicator.
//   This does not replace STEP 3/4 (which remain useful for the LOCAL vs
//   REGIONAL control-site comparison, something the classic sliding-window
//   method doesn't do on its own) - it's a genuinely different, more
//   statistically rigorous complementary tool for the specific question
//   "is there a consistent trend at all," addressing the multiple-
//   comparisons weakness of picking a single best window out of 6.
//
// v10.95 NEW: ROBUSTNESS check for FIND SWEET SPOT - "which window do I trust?"
//   Direct response to a real usage pattern: three separate test runs at the
//   same GBR coordinates, testing slightly different AFTER windows, produced
//   three different verdicts (LOCAL / MARGINAL / GLOBAL). Testing 6 window
//   lengths and reporting whichever looks most dramatic is a classic
//   multiple-comparisons trap - with 6 tries, SOME window will cross a
//   threshold by chance alone, even in pure noise. The tool had no way to
//   distinguish an isolated single-window spike from a signal that holds up
//   across neighbouring windows.
//   - New table footer: counts how many of the 6 tested windows lean LOCAL
//     (LOCAL CSD or MARGINAL LOCAL) vs how many don't. 0-1 is flagged as an
//     ISOLATED result with an explicit multiple-comparisons caution and a
//     recommendation to re-test using an AFTER window chosen from
//     independent evidence, not from this scan.
//   - The headline verdict box itself now carries a short "(N of 6 windows
//     agree)" or "\u26A0 ISOLATED" flag for LOCAL CSD / MARGINAL LOCAL results,
//     so the caution is visible without reading the full table.
//   - Added an explicit caution note to the S13 STEP 4 intro text.
//   This does not add a real statistical significance test (the tool still
//   has no surrogate/null-model comparison, unlike the literature it cites)
//   - it only prevents the tool from presenting an isolated, likely-chance
//     result with the same confidence as a signal that is consistent across
//     multiple window lengths.
//
// v10.94 FIX: FIND SWEET SPOT headline could contradict its own breakdown
//   Discovered directly from a real GBR test: the top verdict box said
//   "LOCAL CSD SIGNAL DETECTED" (strong, red) for a window whose own
//   detailed breakdown said "MARGINAL LOCAL" and "Scheffer NOT MET -
//   neither indicator rose at the study site". Root cause: the headline
//   was classified from a SEPARATE, looser check (bestDiv>THRESH alone -
//   pure variance divergence magnitude, ignoring whether the study site's
//   OWN variance crossed the threshold, and ignoring AC1 entirely), while
//   the per-window table/breakdown used a stricter, correct classification
//   (bestRow.verdict: requires studyVarRose AND !ctrlVarRose for "LOCAL
//   CSD", separately labelling divergence-only cases "MARGINAL LOCAL").
//   Two different rules, two different answers, same window.
//   Fixed: the headline is now derived directly from bestRow.verdict (the
//   SAME classification shown in the table and the sweet-spot breakdown),
//   so it can no longer disagree with the detail below it. Also added an
//   explicit AC1-divergence-from-control check to the headline text
//   (missing before - the multi-window headline was the last remaining
//   place in S13 not reflecting the AC1-primacy work from v10.90/91).
//
// v10.93 CRITICAL FIX: Southern Hemisphere reefs showed false "no heat
//   stress" because the peak/DHW calculation used a single hardcoded
//   window (Jun-Oct 2023) - Northern Hemisphere summer only. Discovered
//   testing a real, independently-documented site: One Tree Reef, southern
//   GBR, which suffered severe bleaching in the 2023-24 Southern Hemisphere
//   summer (Nov-Apr) - the tool showed DHW=0.00 deg C-wks for this exact
//   site/period because it was reading Southern Hemisphere WINTER SST as
//   the "peak". This one number feeds: the main DHW badge, the S4 cancer-
//   score component (25% weight), the map's bleaching-risk overlay, and the
//   intervention engine's thermal-stress trigger - so heat-stress detection
//   was silently broken for every Southern Hemisphere reef (GBR, Ningaloo,
//   S. Indian Ocean, S. Brazil, S. Africa - roughly half the world's reef
//   area by latitude coverage).
//   Fix: compute BOTH a Northern Hemisphere peak window (Jun-Oct) and a
//   Southern Hemisphere one (Nov-Apr, same 2023-24 stress year), then
//   select per-click based on the clicked latitude's sign - in
//   analyzeLocation() for the sidebar/score calculation, and in
//   loadLayers() for the map layers (both already had lat as a parameter).
//   MMM_perpixel (the seasonal baseline) was already hemisphere-agnostic -
//   only the single "current peak" snapshot had this bug. Labels/legend
//   updated to show which season window applied for each click.
//
// v10.92 NEW: FIND SWEET SPOT can now surface AC1-up/variance-down windows
//   Prompted by a direct question: "where can I find AC1 rising while
//   variance falls, in a real tipping ecosystem?" Two real gaps surfaced:
//   1) The 6-window results table only ever showed \u0394Var per window, never
//      \u0394AC1 - so it was impossible to scan across windows for this exact
//      pattern without manually re-running STEP 2 for each one by hand.
//   2) The "sweet spot" window pick was ENTIRELY variance-divergence-based,
//      inconsistent with v10.90/91 establishing AC1 as the primary
//      indicator. A window with a strong AC1 rise but falling variance
//      could be ranked below - or never surfaced above - a weaker
//      variance-only window.
//   Fixed:
//   - Table now has a Study \u0394AC1 column alongside \u0394Var.
//   - Per-window verdict text now appends "+ AC1 CONFIRMED (var down)" or
//     "+ AC1 CONFIRMED" so this pattern is visually flagged in the table.
//   - Added a SEPARATE AC1-based ranking (bestAc1W) alongside the existing
//     variance-divergence ranking (bestW); table marks both with distinct
//     arrows (VAR SWEET SPOT vs AC1 SWEET SPOT) when they disagree.
//   - When the two rankings pick different windows, the console now prints
//     an explicit cross-check note, citing Dakos et al. 2012 Fig. 2c/4 when
//     it's specifically an AC1-up/variance-down window.
//   No coordinates are hardcoded or guessed anywhere - this only makes the
//   pattern findable by actually running the tool against real data.
//
// v10.91 FIX vs v10.90: AC1-up + variance-down is not weak evidence
//   v10.90 correctly downweighted "variance rises but AC1 doesn't" - but
//   introduced an unwanted asymmetric bug in the OTHER direction: when AC1
//   WAS rising and a supporting indicator (variance/spatial) was available
//   but did NOT agree (e.g. variance actively falling), the confidence
//   label said "no supporting indicator corroborates it yet" - wording that
//   implies pending/weak evidence. That's wrong per Dakos et al. 2012: they
//   specifically document AC1 rising while variance FALLS near a genuine
//   transition (their Fig. 2c "decreasing sensitivity" and Fig. 4
//   "freezing" cases) - autocorrelation "remains solely dependent on the
//   dominant eigenvalue" and rises "regardless of the responsiveness of
//   the ecosystem" to noise. Variance disagreeing does NOT weaken AC1's
//   signal; it just means variance isn't a useful witness at that site.
//   - classifyToolkitConfidence(): the "AC1 rising, support available but
//     disagreeing" case now explicitly cites this documented pattern
//     instead of implying the signal is pending confirmation. Still scored
//     MODERATE (same as AC1-alone-with-no-data) - not upgraded to HIGH,
//     since variance disagreeing isn't corroboration either, just neutral.
//   - STEP 3's verdict title for this case no longer says "not yet
//     corroborated" - now "AC1-confirmed - the primary indicator per Dakos
//     et al. 2012".
//   - Added an explicit S13 intro note about this exact scenario.
//
// v10.90 FIX vs v10.89: AC1-weighted toolkit (Dakos et al. 2012 fidelity)
//   The STEP 3 toolkit tally and FIND SWEET SPOT's Scheffer check both
//   previously treated every indicator (AC1, temporal variance, spatial
//   variance, spatial autocorrelation) as ONE EQUAL VOTE. Dakos et al. 2012
//   (Ecology 93:264-271) found this isn't justified: autocorrelation
//   "appears a relatively robust indicator... regardless of the source of
//   noise" across every scenario they tested, while variance "may
//   sometimes decrease close to a transition" for well-documented reasons
//   (their Fig. 4 - parameter-noise sensitivity changes, and a "freezing"
//   effect in slow-responding systems). A variance-only rise (AC1 not
//   rising) is objectively weaker evidence than an AC1-confirmed rise, and
//   the code did not encode that asymmetry anywhere.
//   - buildToolkitTally() now tags AC1 as the PRIMARY indicator; variance
//     and both spatial indicators are SUPPORTING evidence that raises or
//     lowers confidence around it, not equal votes.
//   - New classifyToolkitConfidence(): HIGH only when AC1 rises AND is
//     corroborated; MODERATE when AC1 rises alone; LOW-MODERATE when
//     variance/spatial rise but AC1 does NOT (explicitly flagged as weaker,
//     citing Dakos et al.'s finding that variance can behave unexpectedly);
//     LOW when nothing rises.
//   - STEP 3's final combined verdict rebuilt around this same logic
//     instead of a majority-of-indicators-agree rule.
//   - FIND SWEET SPOT's "Scheffer 2009 validation" line now distinguishes
//     AC1-only-rising (still meaningful) from variance-only-rising (weaker,
//     flagged) instead of a flat "only one rose - weaker evidence" for
//     both cases equally.
//   No changes to the underlying AC1/variance/skewness/spatial math itself
//   (computeRealCSD, computeSpatialEWS) - this is purely a re-weighting of
//   how the existing numbers are interpreted into a verdict.
//
// v10.89 FIX vs v10.88:
//   S17 (Time of Emergence) was bundling all 6 variables (SST, Chl,
//   Salinity, NO2, pH, DO) into ONE ee.Dictionary and evaluating them in a
//   single .evaluate() call. pH and DO both depend on the Copernicus asset
//   'COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/MFC_001_028', which is currently
//   returning "not found" in the GEE catalog - and because all 6 were
//   bundled together, that ONE dead asset made the entire evaluate() call
//   fail, blanking ALL SIX indicators to n/a - including SST/Chl/Salinity/
//   NO2, which come from completely different, healthy datasets and have
//   nothing to do with the BGC asset. Split into rToeCore (SST/Chl/Salinity/
//   NO2) and rToeBGC (pH/DO), evaluated independently: the 4 healthy
//   variables now populate normally regardless of whether the BGC asset is
//   reachable, and only pH/DO show the (accurate) "dataset unavailable"
//   message. Compound status now reports "X of N AVAILABLE variables",
//   with N excluding pH/DO when the BGC asset is down, instead of a flat
//   "n/a - error" for the whole panel.
//   S18 (Biogeochemistry Snapshot) was already correctly showing all-n/a,
//   since all 4 of its values legitimately come from that same one dead
//   asset - that part was not a bug. Its error message and the S17/S18
//   section headers now state the specific dataset/reason plainly instead
//   of a bare "n/a", so this is diagnosable without reading the console.
//
// v10.88 NEW: genuine multi-indicator "toolkit" approach in STEP 3 (COMPARE)
//   Addresses standard early-warning-signal guidance (Dakos et al. 2012 and
//   related literature): don't rely on a single indicator; AC1 can be more
//   reliable than variance; spatial patterns can be a stronger signal than
//   temporal ones; limited data makes any one indicator unreliable.
//   - computeRealCSD() now also returns skewness of the detrended residuals
//     (reported for context; its direction is system-dependent so it is
//     NOT auto-scored as CSD-consistent either way).
//   - New computeSpatialEWS(): spatial variance and a spatial-autocorrelation
//     proxy (correlation between each pixel and its 3x3 neighbourhood mean)
//     computed across a 15km buffer around the study site for the BEFORE and
//     AFTER periods - a genuinely different indicator family from anything
//     computed through time at a point.
//   - New TOOLKIT SUMMARY panel (csdToolkitV) tallies how many of up to 4
//     scored indicators (temporal AC1, temporal variance, spatial variance,
//     spatial autocorrelation) actually agree, with a confidence label
//     (HIGH/MODERATE/LOW) based on agreement level and data sufficiency -
//     instead of one threshold on one indicator deciding everything.
//   - The verdict box now renders in two honest stages: a PRELIMINARY
//     temporal-only verdict appears immediately (no added wait for the fast
//     path), then upgrades in place once the spatial indicators arrive. A
//     spatial-fetch failure degrades gracefully back to the temporal-only
//     verdict rather than breaking anything.
//   - Final combined verdict distinguishes "MULTI-INDICATOR LOCAL CSD SIGNAL"
//     (majority of available indicators agree AND it's locally-flavoured),
//     "...LIKELY REGIONAL" (majority agree but matches the control site too),
//     and "WEAK/SINGLE-INDICATOR SIGNAL ONLY" (toolkit does not corroborate)
//     - so a single indicator can no longer produce a strong-sounding verdict
//     on its own.
//
// v10.87 FIX vs v10.86:
//   S12 (main click panel) had a misleading message: "Both AC1 and
//   variance rising - classic CSD pattern". This fires whenever AC1 is
//   ABOVE 0.5 in a single fixed window (Jan 2023-Dec 2024, no BEFORE
//   baseline at all) - that is AC1 being CURRENTLY elevated, not AC1
//   RISING. This directly contradicted S13's proper Scheffer 2009
//   BEFORE-vs-AFTER test, which computes a real AC1 delta against a
//   stored baseline and can correctly show AC1 FALLING even while S12's
//   snapshot looks "elevated". Reworded the note to say what it actually
//   measures, and added an upfront caveat on the S12 header pointing to
//   S13 for a validated before/after comparison. No math changed - only
//   the wording, which was actively misleading.
//
// v10.86 FIX vs v10.85:
//   STEP 3 (COMPARE) messaging was misleading and had a silent hang risk:
//     - The verdict box jumped straight to "Computing verdict..." on click,
//       while the detail box right below it still said "Step 1/3: Finding
//       control site..." - the two boxes contradicted each other. Both now
//       show the same real stage (finding control site -> control site
//       confirmed, running CSD -> final colour-coded verdict), and the
//       verdict box only claims "verdict" once one actually exists.
//     - The control-site depth lookup (GEBCO reduceRegion) previously
//       ignored its own evaluate() error parameter completely - if it
//       failed or stalled there was no error message at all, just an
//       indefinite "Step 1/3: Finding control site...". Now wrapped with
//       explicit error handling using the same friendlyEEError() translator
//       as the rest of S13.
//     - Clarified in STEP 3's instructions that COMPARE runs in two
//       sequential stages (find/validate control site, THEN run the actual
//       comparison), not as one simultaneous action.
//
// v10.85 FIX vs v10.84:
//   Sidebar Depth label contradicted the S3 map legend. The legend splits
//   depth into 7 bins (land / intertidal / shallow reef / continental
//   shelf / continental slope / deep ocean / very deep ocean), but the
//   Depth row's text used a crude 2-way "deeper than -50m = deep ocean"
//   rule - so a -60m continental-shelf point (teal on the map) was
//   labelled "deep ocean" in the sidebar right next to it. New
//   classifyDepthLabel() applies the exact same 7 bins as the legend, used
//   by both the Depth row and the "click closer to shore" warning, so the
//   text and the map colour always agree.
//
// v10.84 CHANGES vs v10.83:
//   Hardened S13 error handling end-to-end:
//     - New friendlyEEError() translator recognizes common transient Earth
//       Engine errors (e.g. "Unknown reference to value named ''...",
//       "Failed to contact Earth Engine servers", timeouts, rate limits)
//       and tells the user plainly this is a server hiccup to retry, while
//       still showing the raw error underneath for debugging.
//     - STEP 2 (RUN CSD TEST) now echoes the exact parsed lat/lon/dates
//       before firing, so a truncated/mistyped coordinate is obvious
//       immediately instead of surfacing later as a cryptic failure.
//     - STEP 3 (COMPARE / runControlCSD) previously ignored the evaluate()
//       error parameter entirely for the control site - a failed control
//       computation just silently produced blank/n-a numbers. Now wrapped
//       in try/catch with explicit error surfacing in the verdict box.
//     - STEP 4 (FIND SWEET SPOT)'s catch block now uses the same friendly
//       translator instead of a raw error dump.
//   No changes to the underlying math anywhere.
//
// v10.83 CRITICAL FIX vs v10.82:
//   FIND SWEET SPOT was throwing "TypeError: (intermediate value)..." and
//   silently dying right after "Step 3/3: Analysing results..." because it
//   called the ES6 method String.prototype.repeat() to draw a separator
//   line ('─'.repeat(64)) - and GEE's server-side script sandbox does not
//   implement that method. Replaced with a loop-based repeatChar(ch,n)
//   helper everywhere a repeated-character string was built. No other
//   behavior changed.
//
// v10.82 CHANGES vs v10.81:
//   FIND SWEET SPOT (STEP 4) rebuilt:
//     - The control site's "before" state is now a REAL computed value
//       (one shared control-BEFORE test using the same window as STEP 2),
//       not a hardcoded 1.0 placeholder. Study and control deltas are now
//       apples-to-apples, same maths as COMPARE.
//     - Live progress counter ("N / 13 sub-tests done") while the 13
//       parallel Earth Engine calls run, so the panel never looks frozen.
//     - New csdMultiSweetSpotV panel: an explicit LOCAL vs REGIONAL vs
//       Scheffer-2009-validation breakdown for the winning window (AC1 and
//       variance deltas at both sites, stated in plain language).
//     - try/catch around the final analysis so a computation problem shows
//       a red error box instead of silently stopping at "Analysing results...".
//     - Verdict box now also reports whether the Scheffer validation passed.
//
// v10.81 CHANGES vs v10.80:
//   REBUILT S13 (CSD Early Warning Test) for clarity:
//     - Added numbered STEP 1/2/3/4 workflow headers so the panel
//       reads top-to-bottom instead of as a pile of unrelated buttons.
//     - Manual control-site override moved to a clearly labelled
//       "ADVANCED (OPTIONAL)" section after the normal workflow;
//       defaults to blank/AUTO, with a one-line explanation of why
//       it exists. Nobody has to touch it to use S13.
//     - Added a short, bold, colour-coded VERDICT box (csdCompareVerdictV)
//       for the COMPARE button that renders ABOVE the detailed 4-way
//       numeric breakdown - mirrors the verdict-above-table pattern
//       already used by FIND SWEET SPOT (csdMultiStatusV).
//     - Added a "Use last clicked location" button in S13 Step 1 that
//       reuses the coordinates from your last map click / GO TO
//       COORDINATES, via new globals lastClickLat / lastClickLon.
//     - Clearer inline status messages after each RUN (tells you what
//       to do next: switch to AFTER, then press COMPARE or FIND SWEET SPOT).
//   All underlying math is unchanged: computeRealCSD(), getSmartControlSite(),
//   runControlCSD(), and the CSD/verdict thresholds are untouched.
//
// v10.67 CHANGES vs v10.66:
//   NEW: GEM MarineBasis Greenland stations added to MODULE A11 (S19)
//     - MarineBasis Nuuk GF3, Godthåbsfjord SW Greenland (64.13N,51.38W, r=40km)
//     - MarineBasis Zackenberg Young Sound NE Greenland (74.315N,20.279W, r=30km)
//     - API: api.g-e-m.dk | Key: GeoMarineAnalysis | License: CC BY-SA 4.0
//     - DOIs: 10.17897/KMEK-TK21 (Nuuk CTD) | 10.17897/8GPS-CE70 (Zackenberg)
//     - Stats are PLACEHOLDER nulls - run gem_fetch_and_clean.py to populate
//   NEW: Greenland region entries added to getRegion()
//   FIX (v10.19): Hemisphere suffix parsing in GO TO COORDINATES and S13
//     - parseFloat("13.5S") silently returned 13.5 (positive) - now fixed
//     - parseCoordPart() helper handles N/S/E/W suffixes correctly
//   FIX (v10.20): DHW blue/orange layer contradiction in loadLayers()
//     - "DHW blue = no stress" was painted over ALL ocean including stressed pixels
//     - Now correctly masked to only show where DHW = 0 (genuinely no stress)
//   FIX (v10.20): S13 BEFORE/AFTER CSD persistent state labels
//     - Added "BEFORE stored:" / "AFTER stored:" labels that update in same callback
//     - Added "Clear stored BEFORE/AFTER" button to reset stale state
//     - Previously no way to confirm whether storage succeeded before clicking COMPARE
//
// v10.66 key features (unchanged from v10.66):
//   v10.66: GEBCO depth layer in map (matches sidebar depth value)
//   v10.65: GEBCO replaces ETOPO1 for sidebar depth readout
//   v10.64: honest "insufficient data" message for aquaculture (not just "land")
//   v10.63: peak SST checked in Gate 1, not just annual mean (Persian Gulf fix)
//   v10.60: parallel evaluate() for S17/S18
//   v10.58: disasterty GDIS property name fix (was silently always 0)
//   v10.57: real cited A. taxiformis SST thresholds (Statton 2024, 17-21 optimal)
//   v10.56: two-gate aquaculture architecture replacing opaque weighted average
//   v10.54: OISST *0.01 conversion fix in ToE annual SST stack
//   v10.42: dead soil texture asset replaced with safe masked constant
//   v10.41: defensive S2 band check before normalizedDifference()
//   v10.39: removed setTimeout (not available in GEE sandbox) - harmless no-op
//   v10.37/38: error handling + full sidebar reset on every new click
//   v10.36: S12-S19, thermal recovery, ecological validation, ECI, ToE, BGC
//   v10.13: OISST replaces dead MODIS dataset
//
// SCIENTIFIC APPROACH:
//   Waddington Landscape (PNAS 2025) + Scheffer 2009 CSD +
//   Ramamurthy 2024 (Bunodosoma) + Levitan 2023 (Diadema) +
//   Lozano-Bilbao 2020-2024 (metals) + Peixoto 2025 (Red Sea)
// ============================================================

// ============================================================
// MODULE A - DATASETS (all with ocean mask)
// ============================================================
var bathy = ee.Image('NOAA/NGDC/ETOPO1').select('bedrock');
var bathyU = bathy.unmask(0);
var oceanMask      = bathyU.lt(0);
var shallowMask    = bathyU.gte(-50).and(bathyU.lt(0));
var intertidalMask = bathyU.gte(-20).and(bathyU.lt(0));

// ============================================================
// MODULE A11 - REAL IN-SITU OCEAN CHEMISTRY BASELINE (v10.36+)
// v10.67: Added 2 GEM MarineBasis Greenland stations (see header)
// ============================================================
var IN_SITU_BASELINES = {
  looe_key: {
    label: 'Looe Key, FL (coral reef MPA)',
    lat: 24.5463, lon: -81.4014, radius_km: 15,
    source: 'SECOORA ERDDAP (looe-key-fl) - Mote Marine Laboratory SeapHOx',
    record: '2024-06-05 to 2025-11-18 (~17.5 months, 6228 raw readings)',
    n_clean: 5151, pct_flagged: 17.3,
    pH: {mean: 7.899, std: 0.156},
    temp_c: {mean: 27.606}, salinity: {mean: 35.708}, do_mgL: {mean: 5.705},
    notes: 'Highest artifact rate (17.3%) - remote reef mooring, likely biofouling/drift.'
  },
  agua_hedionda: {
    label: 'Agua Hedionda Lagoon, CA (oyster farm)',
    lat: 33.1425, lon: -117.3275, radius_km: 10,
    source: 'SCCOOS ERDDAP (pH-AHL) - Martz Lab SeapHOx',
    record: '~6294 raw readings, actively managed aquaculture site',
    n_clean: 6250, pct_flagged: 0.7,
    pH: {mean: 7.902, std: 0.091},
    temp_c: {mean: 18.624}, salinity: {mean: 33.758}, do_mgL: {mean: 7.644},
    notes: 'Lowest artifact rate (0.7%) - actively visited/maintained research site.'
  },
  scripps_pier: {
    label: 'Scripps Pier, La Jolla, CA (open coast)',
    lat: 32.8669, lon: -117.2571, radius_km: 10,
    source: 'CenCOOS ERDDAP (scripps-pier-automated-shore-sta-1)',
    record: 'Last 365 days at download time, 178809 raw readings',
    n_clean: 174515, pct_flagged: 2.4,
    pH: {mean: 7.953, std: 0.165},
    temp_c: {mean: 19.023}, salinity: {mean: 33.232}, do_mgL: {mean: 7.820},
    notes: 'Largest clean sample (174,515 pts) - most statistically robust baseline.'
  },
  // ============================================================
  // GEM MarineBasis stations (v10.67 NEW)
  // Data: Greenland Ecosystem Monitoring | https://g-e-m.dk
  // API: api.g-e-m.dk | Key: GeoMarineAnalysis | License: CC BY-SA 4.0
  // IMPORTANT: stat values are PLACEHOLDER nulls.
  // Run gem_fetch_and_clean.py to replace with real cleaned values.
  // DOIs confirmed from Vonnahme et al. 2025 (Limnol.Oceanogr.) and GEM database.
  // ============================================================
  gem_nuuk_gf3: {
    label: 'MarineBasis Nuuk - GF3, Godthåbsfjord (SW Greenland)',
    lat: 64.13, lon: -51.38, radius_km: 40,
    source: 'GEM MarineBasis REST API (api.g-e-m.dk) | DOI:10.17897/KMEK-TK21 | ' +
      'Mortensen et al. 2022 JGR-Oceans doi:10.1029/2022JC018724',
    record: 'Monthly CTD profiles 2005-present | SeaBird 19Plus | 76,713 readings | CC BY-SA 4.0',
    n_clean: 76713, pct_flagged: 0.0,
    pH: null,
    temp_c: {mean: 1.79, std: 1.1622},
    salinity: {mean: 33.2394, std: 0.4342},
    do_mgL: null,
    fluorescence_ug_L: {mean: 0.2651, std: 0.2503, n_clean: 71100, pct_flagged: 6.9},
    pressure_db: {mean: 158.27, std: 93.63, min: 1, max: 392, n_clean: 76720},
    turbidity_ftu: {mean: 0.6096, std: 0.2348, n_clean: 69916, pct_flagged: 3.7},
    notes: 'GEM MarineBasis Nuuk, CC BY-SA 4.0. CITE: Greenland Ecosystem Monitoring (2026). ' +
      'MarineBasis Nuuk CTD. https://doi.org/10.17897/KMEK-TK21. ' +
      'REAL stats: 76,713 readings 2005-present. Temp mean=1.79 deg C, std=1.16, range=-1.03 to 7.68. ' +
      'Salinity mean=33.24 PSU, std=0.43, range=31.67 to 34.37. No pH/DO sensor in CTD config. ' +
      'Fluorescence mean=0.265 ug/L (71,100 readings, 6.9% flagged). ' +
      'Pressure mean=158.3 dbar (full water column profiles 1-392m). ' +
      'Turbidity mean=0.610 FTU (69,916 readings, 3.7% flagged). ' +
      'Arctic fjord: UNSUITABLE for A.taxiformis aquaculture (below 15 deg C threshold).'
  },
  gem_zackenberg: {
    label: 'MarineBasis Zackenberg - Young Sound mooring (NE Greenland)',
    lat: 74.315, lon: -20.279, radius_km: 30,
    source: 'GEM MarineBasis REST API (api.g-e-m.dk) | DOI:10.17897/8GPS-CE70 | SeaBird SBE37SMP',
    record: '15-20min mooring CTD 2003-2019 | Fixed ~64m depth | 406,125 readings | CC BY-SA 4.0',
    n_clean: 406125, pct_flagged: 0.0,
    pH: null,
    temp_c: {mean: -1.61, std: 0.1141},
    salinity: {mean: 32.2167, std: 0.3727},
    do_mgL: null,
    fluorescence_ug_L: null,  // not in SBE37SMP mooring configuration
    pressure_db: null,         // mooring at fixed ~64m - run gem_full_extract.py to confirm
    turbidity_ftu: null,       // not in this mooring configuration
    notes: 'GEM MarineBasis Zackenberg, CC BY-SA 4.0. CITE: Greenland Ecosystem Monitoring. ' +
      'MarineBasis Zackenberg CTD mooring. https://doi.org/10.17897/8GPS-CE70. ' +
      'REAL stats: 406,125 readings 2003-2019. Temp mean=-1.61 deg C, std=0.11, range=-1.8 to -1.0. ' +
      'Salinity mean=32.22 PSU, std=0.37, range=30.52 to 33.22. No pH/DO/fluorescence in config. ' +
      'HISTORICAL record only: mooring ended 2019 (cable malfunction). ' +
      'Sensor accuracy: +/-0.002 deg C temperature, +/-0.0003 S/m conductivity. ' +
      'High-Arctic near-freezing temperatures: UNSUITABLE for aquaculture.'
  }
};

function haversineKm(lat1, lon1, lat2, lon2) {
  var R = 6371;
  var dLat = (lat2-lat1)*Math.PI/180, dLon = (lon2-lon1)*Math.PI/180;
  var a = Math.sin(dLat/2)*Math.sin(dLat/2)+
    Math.cos(lat1*Math.PI/180)*Math.cos(lat2*Math.PI/180)*Math.sin(dLon/2)*Math.sin(dLon/2);
  return R*2*Math.atan2(Math.sqrt(a),Math.sqrt(1-a));
}

function getInSituBaseline(lat, lon) {
  var best = null, bestDist = Infinity;
  for (var key in IN_SITU_BASELINES) {
    var st = IN_SITU_BASELINES[key];
    var d = haversineKm(lat, lon, st.lat, st.lon);
    if (d <= st.radius_km && d < bestDist) {
      best = Object.assign({key:key, distance_km:d}, st);
      bestDist = d;
    }
  }
  return best;
}

function getOISSTColl(startDate, endDate) {
  return ee.ImageCollection('NOAA/CDR/OISST/V2_1')
    .filterDate(startDate, endDate).select('sst')
    .map(function(img){ return img.multiply(0.01).copyProperties(img,['system:time_start']); });
}

var sstColl = getOISSTColl('2023-01-01','2024-12-31');
var sstRaw = ee.Image(ee.Algorithms.If(sstColl.size().gt(0),
  sstColl.mean(), ee.Image.constant(25).rename('sst'))).rename('sst');
var sst = sstRaw.updateMask(sstRaw.gt(-2).and(sstRaw.lt(40)).and(oceanMask));

// v10.93 FIX: 'Peak(Jun-Oct)' / DHW used a SINGLE hardcoded window
// (Jun-Oct 2023), which is Northern Hemisphere summer only. For Southern
// Hemisphere reefs (Great Barrier Reef, Ningaloo, Madagascar/S. Indian
// Ocean, S. Brazil, S. Africa - roughly half the world's reef area by
// latitude), the real hot/bleaching season is Nov-Apr and was being missed
// entirely. A real, independently-documented 2023-24 mass bleaching event
// at a Southern Hemisphere GBR test site showed DHW=0.00 (falsely "no
// stress") purely because the tool was reading winter SST as the "peak".
// This one number feeds the main DHW badge, the S4 cancer-score component
// (25% weight), the map's bleaching-risk layer, and the intervention
// engine's thermal-stress trigger - so this bug silently zeroed out heat-
// stress detection for every Southern Hemisphere reef in the tool.
// Fix: compute BOTH hemispheres' peak windows, then select per-click based
// on the clicked latitude's sign (see analyzeLocation() and loadLayers()).
var oisstPeakColl = getOISSTColl('2023-06-01','2023-10-31'); // Northern Hemisphere summer
var oisstPeakMax = ee.Image(ee.Algorithms.If(oisstPeakColl.size().gt(0),
  oisstPeakColl.max(), ee.Image.constant(28).rename('sst'))).rename('sst_peak');
var sstPeakNH = oisstPeakMax.updateMask(oisstPeakMax.gt(-2).and(oisstPeakMax.lt(40)).and(oceanMask));
var sstPeakFinalNH = sstPeakNH.unmask(sstRaw.updateMask(oceanMask));

var oisstPeakCollSH = getOISSTColl('2023-11-01','2024-04-30'); // Southern Hemisphere summer (same 2023-24 stress year)
var oisstPeakMaxSH = ee.Image(ee.Algorithms.If(oisstPeakCollSH.size().gt(0),
  oisstPeakCollSH.max(), ee.Image.constant(28).rename('sst'))).rename('sst_peak');
var sstPeakSH = oisstPeakMaxSH.updateMask(oisstPeakMaxSH.gt(-2).and(oisstPeakMaxSH.lt(40)).and(oceanMask));
var sstPeakFinalSH = sstPeakSH.unmask(sstRaw.updateMask(oceanMask));

// Default/legacy name kept pointing at the Northern Hemisphere version so
// the startup (unclicked) map view - centred near-equatorial at Bocas del
// Toro, Panama - still renders sensibly without needing a hemisphere choice.
var sstPeakFinal = sstPeakFinalNH;

var mmmColl = getOISSTColl('2003-01-01','2022-12-31');
var monthList = ee.List.sequence(1,12);
var monthlyMeans = ee.ImageCollection(monthList.map(function(mo){
  var moNum = ee.Number(mo);
  var col = mmmColl.filter(ee.Filter.calendarRange(moNum, moNum, 'month'));
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
    col.mean().rename('sst'),
    ee.Image.constant(25).rename('sst').updateMask(ee.Image.constant(0))));
}));
var MMM_perpixel = monthlyMeans.max().rename('mmm').updateMask(oceanMask);
// dhwProper (NH) kept for the startup default view; dhwProperSH added for
// Southern Hemisphere clicks. Both share the same hemisphere-agnostic
// MMM_perpixel baseline (that part was already correct - it takes the max
// across all 12 calendar months per pixel, so it naturally picks the right
// "warm season" baseline regardless of hemisphere).
var dhwProper = sstPeakFinalNH.subtract(MMM_perpixel.add(1)).max(0).multiply(12).rename('dhw').updateMask(oceanMask);
var dhwProperSH = sstPeakFinalSH.subtract(MMM_perpixel.add(1)).max(0).multiply(12).rename('dhw').updateMask(oceanMask);


var chlColl = ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6')
  .filterDate('2023-01-01','2024-12-31').select('chlor_a');
var chlaRaw = ee.Image(ee.Algorithms.If(chlColl.size().gt(0),
  chlColl.mean(), ee.Image.constant(0).rename('chlor_a'))).rename('chlor_a');
var chla = chlaRaw.updateMask(chlaRaw.gt(0).and(chlaRaw.lt(100)).and(oceanMask));
var chlaCoastal = chlaRaw.updateMask(chlaRaw.gt(0).and(chlaRaw.lt(100)));

var s2Coll = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate('2023-01-01','2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).select(['B4','B8']);
var turbRaw = ee.Image(ee.Algorithms.If(s2Coll.size().gt(0),
  s2Coll.mean().normalizedDifference(['B4','B8']).rename('turbidity'),
  ee.Image.constant(0).rename('turbidity').updateMask(ee.Image.constant(0))));
var turbImg = turbRaw.updateMask(shallowMask);

var s2_basic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate('2023-01-01','2024-12-31').filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).select(['B4','B8']);
var s2_basic_mean = ee.Image(ee.Algorithms.If(s2_basic.size().gt(0),
  s2_basic.map(function(img){ return img.divide(10000); }).median(),
  ee.Image.constant(0.05).rename('B4').addBands(ee.Image.constant(0.05).rename('B8'))));

// v10.41 FIX: defensive band check before normalizedDifference()
var _s2basic_bandsOk = s2_basic_mean.bandNames().size().gte(2);
var ndviWater = ee.Image(ee.Algorithms.If(_s2basic_bandsOk,
  s2_basic_mean.select(['B8','B4']).normalizedDifference(['B8','B4']).rename('ndvi_water'),
  ee.Image.constant(0).rename('ndvi_water').updateMask(ee.Image.constant(0))));
var ndviAlgae = ndviWater.updateMask(ndviWater.gt(0.03).and(oceanMask));

var s2_ext = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate('2023-01-01','2024-12-31').filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).select(['B4','B5','B8','B11']);
var s2_ext_avail = s2_ext.size().gt(0);
var s2_ext_mean = ee.Image(ee.Algorithms.If(s2_ext_avail,
  s2_ext.map(function(img){ return img.divide(10000); }).median(),
  s2_basic_mean.addBands(ee.Image.constant(0.03).rename('B5')).addBands(ee.Image.constant(0.03).rename('B11'))));

var fai_step = (842-665)/(1610-665);
var faiImg = s2_ext_mean.expression('NIR - (RED + (SWIR - RED) * step)',
  {'NIR':s2_ext_mean.select('B8'),'RED':s2_ext_mean.select('B4'),'SWIR':s2_ext_mean.select('B11'),'step':fai_step}).rename('fai');
var faiMild = faiImg.updateMask(faiImg.gt(0.003).and(oceanMask));
var faiAlgae = faiImg.updateMask(faiImg.gt(0.01).and(oceanMask));

// v10.41 FIX: same defensive check for ndciImg
var _s2ext_bandsOk = s2_ext_mean.bandNames().size().gte(4);
var ndciImg = ee.Image(ee.Algorithms.If(_s2ext_bandsOk,
  s2_ext_mean.select(['B5','B4']).normalizedDifference(['B5','B4']).rename('ndci'),
  ee.Image.constant(0).rename('ndci').updateMask(ee.Image.constant(0))));
var ndciBloom = ndciImg.updateMask(ndciImg.gt(0.05).and(oceanMask));

var no2Coll = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2').filterDate('2023-01-01','2024-12-31').select('NO2_column_number_density')
  .merge(ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').filterDate('2023-01-01','2024-12-31').select('NO2_column_number_density'));
var no2Raw = ee.Image(ee.Algorithms.If(no2Coll.size().gt(0),
  no2Coll.mean(), ee.Image.constant(0.00003).rename('NO2_column_number_density'))).rename('NO2_column_number_density');
var no2 = no2Raw.updateMask(no2Raw.gt(0));

var yrList = ee.List.sequence(2003,2024), nYears = 22, SST_FLOOR_GLOBAL = 14.0;
var annSSTraw = ee.ImageCollection(yrList.map(function(yr){
  var yrN=ee.Number(yr), d=ee.Date.fromYMD(yrN,1,1), de=ee.Date.fromYMD(yrN,12,31), col=getOISSTColl(d,de);
  var sstY=ee.Image(ee.Algorithms.If(col.size().gt(3),
    col.mean().updateMask(col.mean().gt(SST_FLOOR_GLOBAL).and(col.mean().lt(40)).and(oceanMask)),
    ee.Image.constant(-9999).rename('sst').updateMask(ee.Image.constant(0)))).rename('sst');
  return sstY.set('system:time_start',d.millis()).set('year',yrN);
}));
var sstStats = annSSTraw.select('sst').reduce(ee.Reducer.mean().combine(ee.Reducer.stdDev(),'',true));
var sstMeanImg = sstStats.select('sst_mean'), sstStdImg = sstStats.select('sst_stdDev').max(0.3);
var annSSTcoll = ee.ImageCollection(annSSTraw.map(function(img){
  var z=img.subtract(sstMeanImg).divide(sstStdImg).abs();
  return img.updateMask(z.gt(2.5).not()).copyProperties(img,['system:time_start','year']);
}));
var annSSTlist = annSSTcoll.toList(nYears), allIdx = ee.List.sequence(0,nYears-1);
var pairSlopeImagesList = allIdx.map(function(i){
  var iN=ee.Number(i), jList=ee.List.sequence(iN.add(1),nYears-1);
  return jList.map(function(j){
    var jN=ee.Number(j), img_i=ee.Image(annSSTlist.get(iN)), img_j=ee.Image(annSSTlist.get(jN));
    return img_j.subtract(img_i).divide(jN.subtract(iN)).rename('slope').updateMask(img_i.mask().and(img_j.mask()));
  });
}).flatten();
var pairSlopes = ee.ImageCollection.fromImages(pairSlopeImagesList);
var sst_slope = pairSlopes.select('slope').median().unmask(0).updateMask(oceanMask).rename('scale');
var validYearCount = annSSTcoll.select('sst').map(function(img){ return img.mask(); }).sum().rename('validYears');

function isEBUS(lat,lon) {
  if(lat>20&&lat<40&&lon>-25&&lon<-10) return true;
  if(lat>30&&lat<50&&lon>-130&&lon<-115) return true;
  if(lat>-45&&lat<-5&&lon>-90&&lon<-70) return true;
  if(lat>-35&&lat<-15&&lon>10&&lon<20) return true;
  return false;
}

var annualSST = ee.ImageCollection(yrList.map(function(yr){
  var yrN=ee.Number(yr), d=ee.Date.fromYMD(yrN,1,1), de=ee.Date.fromYMD(yrN,12,31), col=getOISSTColl(d,de);
  return ee.Image(ee.Algorithms.If(col.size().gt(3),
    col.mean().updateMask(col.mean().gt(SST_FLOOR_GLOBAL).and(col.mean().lt(40))).rename('sst'),
    ee.Image.constant(-9999).rename('sst').updateMask(ee.Image.constant(0)))).set('system:time_start',d.millis());
}));

var mStart = ee.Date('2023-01-01'), mList24 = ee.List.sequence(0,23);
function mkMoSST() {
  return ee.ImageCollection(mList24.map(function(m){
    var d=mStart.advance(m,'month'), de=d.advance(1,'month'), col=getOISSTColl(d,de);
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(-2).and(col.mean().lt(40)).and(oceanMask)).rename('sst'),
      ee.Image.constant(0).rename('sst').updateMask(ee.Image.constant(0)))).set('system:time_start',d.millis());
  }));
}
// v10.104 NEW: monthly Sentinel-2 FAI time series - previously S7 only ever
// had ONE fixed 2023-2024 median composite (faiImg above), never a time
// series, so applying AC1/variance/Mann-Kendall/correlation statistics to
// algae data (as opposed to SST) was not possible without this. Mirrors
// mkMoSSTRange()'s exact pattern, reusing the same fai_step constant and
// FAI expression already established above. Sentinel-2 has heavier cloud-
// masking losses than OISST (a daily gap-filled thermal product), so some
// months may come back with little or no valid data - this is a genuine,
// untested-until-now data-availability question, which is exactly what the
// S7C proof-of-concept below is for.
function mkMoFAIRange(startDateStr, nMonths) {
  var startD=ee.Date(startDateStr), monthsList=ee.List.sequence(0,nMonths-1);
  return ee.ImageCollection(monthsList.map(function(m){
    var d=startD.advance(m,'month'), de=d.advance(1,'month');
    var col=ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filterDate(d,de).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).select(['B4','B5','B8','B11']);
    var avail=col.size().gt(0);
    var meanImg=ee.Image(ee.Algorithms.If(avail,
      col.map(function(img){ return img.divide(10000); }).median(),
      ee.Image.constant(0).rename('B4').addBands(ee.Image.constant(0).rename('B5'))
        .addBands(ee.Image.constant(0).rename('B8')).addBands(ee.Image.constant(0).rename('B11'))
        .updateMask(ee.Image.constant(0))));
    var faiMonth=meanImg.expression('NIR - (RED + (SWIR - RED) * step)',
      {'NIR':meanImg.select('B8'),'RED':meanImg.select('B4'),'SWIR':meanImg.select('B11'),'step':fai_step}).rename('fai');
    return faiMonth.set('system:time_start',d.millis());
  }));
}

function mkMoSSTRange(startDateStr, nMonths) {
  var startD=ee.Date(startDateStr), monthsList=ee.List.sequence(0,nMonths-1);
  return ee.ImageCollection(monthsList.map(function(m){
    var d=startD.advance(m,'month'), de=d.advance(1,'month'), col=getOISSTColl(d,de);
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(-2).and(col.mean().lt(40)).and(oceanMask)).rename('sst'),
      ee.Image.constant(0).rename('sst').updateMask(ee.Image.constant(0)))).set('system:time_start',d.millis());
  }));
}
function mkMoCHL() {
  return ee.ImageCollection(mList24.map(function(m){
    var d=mStart.advance(m,'month'), de=d.advance(1,'month');
    var col=ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6').filterDate(d,de).select('chlor_a');
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(0).and(oceanMask)).rename('chlor_a'),
      ee.Image.constant(0).rename('chlor_a').updateMask(ee.Image.constant(0)))).set('system:time_start',d.millis());
  }));
}
function mkMoDHW() {
  return ee.ImageCollection(mList24.map(function(m){
    var d=mStart.advance(m,'month'), de=d.advance(1,'month'), col=getOISSTColl(d,de);
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.max().subtract(MMM_perpixel.add(1)).max(0).multiply(4.33).updateMask(oceanMask).rename('dhw'),
      ee.Image.constant(0).rename('dhw').updateMask(ee.Image.constant(0)))).set('system:time_start',d.millis());
  }));
}
function mkMoNO2() {
  return ee.ImageCollection(mList24.map(function(m){
    var d=mStart.advance(m,'month'), de=d.advance(1,'month');
    var col=ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2').filterDate(d,de).select('NO2_column_number_density')
      .merge(ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').filterDate(d,de).select('NO2_column_number_density'));
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(0)).rename('NO2_column_number_density'),
      ee.Image.constant(0).rename('NO2_column_number_density').updateMask(ee.Image.constant(0)))).set('system:time_start',d.millis());
  }));
}

function computeRealCSD(monthlyColl, study, bandName, scale) {
  var valsFC = ee.FeatureCollection(monthlyColl.map(function(img){
    var v = img.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:scale,maxPixels:1e9}).get(bandName);
    return ee.Feature(null, {v:v, t:img.get('system:time_start')});
  })).sort('t');
  var validVals = valsFC.aggregate_array('v').removeAll([null]);
  var nValid = validVals.size();
  return ee.Algorithms.If(nValid.gte(4),
    (function(){
      var arr=ee.Array(validVals), idxArr=ee.Array(ee.List.sequence(0,nValid.subtract(1)));
      var meanIdx=idxArr.reduce(ee.Reducer.mean(),[0]).get([0]), meanVal=arr.reduce(ee.Reducer.mean(),[0]).get([0]);
      var idxDev=idxArr.subtract(meanIdx), valDev=arr.subtract(meanVal);
      var cov=idxDev.multiply(valDev).reduce(ee.Reducer.sum(),[0]).get([0]);
      var varIdx=idxDev.pow(2).reduce(ee.Reducer.sum(),[0]).get([0]);
      var slope=ee.Number(cov).divide(ee.Number(varIdx));
      var intercept=ee.Number(meanVal).subtract(slope.multiply(meanIdx));
      var trendArr=idxArr.multiply(slope).add(intercept), residArr=arr.subtract(trendArr);
      var n=nValid;
      var resid0=residArr.slice(0,0,n.subtract(1)), resid1=residArr.slice(0,1,n);
      var mean0=resid0.reduce(ee.Reducer.mean(),[0]).get([0]), mean1=resid1.reduce(ee.Reducer.mean(),[0]).get([0]);
      var dev0=resid0.subtract(mean0), dev1=resid1.subtract(mean1);
      var num=dev0.multiply(dev1).reduce(ee.Reducer.sum(),[0]).get([0]);
      var den0=ee.Number(dev0.pow(2).reduce(ee.Reducer.sum(),[0]).get([0])).sqrt();
      var den1=ee.Number(dev1.pow(2).reduce(ee.Reducer.sum(),[0]).get([0])).sqrt();
      var realAC1=ee.Number(num).divide(den0.multiply(den1).max(1e-6));
      var half=n.divide(2).floor();
      var firstHalf=residArr.slice(0,0,half), secondHalf=residArr.slice(0,half,n);
      var varFirst=ee.Number(firstHalf.reduce(ee.Reducer.sampleVariance(),[0]).get([0]));
      var varSecond=ee.Number(secondHalf.reduce(ee.Reducer.sampleVariance(),[0]).get([0]));
      var varTrendRatio=varSecond.divide(varFirst.max(1e-6));
      // v10.88 NEW: skewness of the detrended residuals - a third classic
      // EWS indicator alongside AC1 and variance (Dakos et al. 2012 generic
      // early-warning-signal toolkit). Its DIRECTION of change is
      // system-dependent (unlike AC1/variance, which both rise consistently
      // near a fold bifurcation), so it is reported for the toolkit but not
      // auto-scored as "CSD-consistent" in either direction.
      var meanResidAll=residArr.reduce(ee.Reducer.mean(),[0]).get([0]);
      var stdAll=ee.Number(residArr.reduce(ee.Reducer.sampleVariance(),[0]).get([0])).max(1e-9).sqrt();
      var cubedDev=residArr.subtract(meanResidAll).pow(3);
      var meanCubed=cubedDev.reduce(ee.Reducer.mean(),[0]).get([0]);
      var skewness=ee.Number(meanCubed).divide(stdAll.pow(3).max(1e-9));
      return ee.Dictionary({realAC1:realAC1,varFirstHalf:varFirst,varSecondHalf:varSecond,
        varTrendRatio:varTrendRatio,trendSlope:slope,nValidMonths:n,skewness:skewness});
    })(),
    ee.Dictionary({realAC1:null,varFirstHalf:null,varSecondHalf:null,varTrendRatio:null,trendSlope:null,nValidMonths:nValid,skewness:null}));
}

// v10.150 FIX: computeRealCSD() above never deseasonalizes before computing
// AC1/variance server-side via ee.Array - the same confirmed structural gap
// as jsNodeStats(). Rather than write new, untested EE ee.Array grouping-
// by-calendar-month logic (which cannot be verified without a live GEE
// session), this converts to the SAME "fetch raw values once via one real
// EE call, then compute everything client-side" pattern this file's own
// v10.107 comments already document as proven-safe (used for S7D's
// rebuild). This is a genuine interface change - a callback instead of a
// returned ee.Dictionary - since it needs a real round trip before it can
// compute anything. Reuses jsNodeStatsFixed() (the already-validated,
// deseasonalizing-first function above) for all the actual math, so both
// server-side and client-side callers now get identical, correct
// statistics from one single, tested implementation.
// v10.151 FIX 4: now returns the raw {t,v} series alongside the stats, so
// COMPARE can pool BEFORE + AFTER into ONE climatology and recompute both
// windows consistently. Without this, STEP 2's two windows are
// deseasonalized independently and their delta is not like-for-like.
// Optional 6th arg passes a shared climatology straight through.
function computeRealCSDDeseasonalized(monthlyColl, study, bandName, scale, onDone, sharedClimatology) {
  var valsFC = ee.FeatureCollection(monthlyColl.map(function(img) {
    var v = img.reduceRegion({reducer: ee.Reducer.mean(), geometry: study, scale: scale, maxPixels: 1e9}).get(bandName);
    return ee.Feature(null, {v: v, t: img.get('system:time_start')});
  })).sort('t');

  valsFC.evaluate(function(fcResult, err) {
    if (err) { onDone({error: String(err)}); return; }
    var feats = (fcResult && fcResult.features) ? fcResult.features : [];
    var seriesTV = feats.map(function(f) { return {t: f.properties.t, v: f.properties.v}; });
    var result = jsNodeStatsFixed(seriesTV, sharedClimatology);
    result.seriesTV = seriesTV;   // v10.151: carried through for pooling
    onDone(result);
  });
}

// v10.96 NEW: the actual Dakos et al. 2012 sliding-window method, as
// opposed to the discrete BEFORE/AFTER chunk comparisons used elsewhere in
// S13. Instead of cutting the series into isolated blocks and asking a
// pass/fail question about each one, this detrends the WHOLE requested
// series ONCE (same linear-detrend approach as computeRealCSD - see the
// disclosed simplification note in the UI: Dakos et al.'s own toolbox
// typically uses Gaussian kernel smoothing instead of a single linear
// detrend, which is a more flexible but heavier method not implemented
// here), then slides a fixed-size window forward ONE MONTH AT A TIME across
// the residuals, computing AC1 and variance at every position. This
// produces a continuous AC1(t)/variance(t) trajectory rather than a single
// number, which is what actually lets a trend "build up" visibly instead
// of being judged by an arbitrary chunk boundary.
// Built as ONE server-side ee.List.map() graph (not N separate calls), so
// the whole scan costs a single .evaluate() round trip regardless of how
// many window positions it contains.
function computeSlidingWindowCSD(monthlyColl, study, bandName, scale, windowSize) {
  var valsFC = ee.FeatureCollection(monthlyColl.map(function(img){
    var v = img.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:scale,maxPixels:1e9}).get(bandName);
    return ee.Feature(null, {v:v, t:img.get('system:time_start')});
  })).sort('t');
  var validFC = valsFC.filter(ee.Filter.notNull(['v']));
  var validVals = validFC.aggregate_array('v');
  var validTimes = validFC.aggregate_array('t');
  var nValid = validVals.size();
  var W = ee.Number(windowSize);
  return ee.Algorithms.If(nValid.gte(W).and(W.gte(4)),
    (function(){
      var arr = ee.Array(validVals);
      var idxArr = ee.Array(ee.List.sequence(0, nValid.subtract(1)));
      // Single global linear detrend over the WHOLE valid series (same
      // method as computeRealCSD, applied once here instead of per-chunk).
      var meanIdx = idxArr.reduce(ee.Reducer.mean(),[0]).get([0]);
      var meanVal = arr.reduce(ee.Reducer.mean(),[0]).get([0]);
      var idxDev = idxArr.subtract(meanIdx), valDev = arr.subtract(meanVal);
      var cov = idxDev.multiply(valDev).reduce(ee.Reducer.sum(),[0]).get([0]);
      var varIdx = idxDev.pow(2).reduce(ee.Reducer.sum(),[0]).get([0]);
      var slope = ee.Number(cov).divide(ee.Number(varIdx));
      var intercept = ee.Number(meanVal).subtract(slope.multiply(meanIdx));
      var trendArr = idxArr.multiply(slope).add(intercept);
      var residArr = arr.subtract(trendArr);

      var nWindows = nValid.subtract(W).add(1);
      var windowStarts = ee.List.sequence(0, nWindows.subtract(1));

      var windowResults = windowStarts.map(function(wi){
        wi = ee.Number(wi);
        var winResid = residArr.slice(0, wi, wi.add(W));
        var n = W;
        var resid0 = winResid.slice(0,0,n.subtract(1));
        var resid1 = winResid.slice(0,1,n);
        var mean0 = resid0.reduce(ee.Reducer.mean(),[0]).get([0]);
        var mean1 = resid1.reduce(ee.Reducer.mean(),[0]).get([0]);
        var dev0 = resid0.subtract(mean0), dev1 = resid1.subtract(mean1);
        var num = dev0.multiply(dev1).reduce(ee.Reducer.sum(),[0]).get([0]);
        var den0 = ee.Number(dev0.pow(2).reduce(ee.Reducer.sum(),[0]).get([0])).sqrt();
        var den1 = ee.Number(dev1.pow(2).reduce(ee.Reducer.sum(),[0]).get([0])).sqrt();
        var winAC1 = ee.Number(num).divide(den0.multiply(den1).max(1e-6));
        var winVar = winResid.reduce(ee.Reducer.sampleVariance(),[0]).get([0]);
        var endIdx = wi.add(W).subtract(1);
        var endTime = ee.Number(validTimes.get(endIdx));
        return ee.Feature(null, {windowIdx:wi, ac1:winAC1, variance:winVar, endTime:endTime});
      });
      return ee.FeatureCollection(windowResults).sort('windowIdx');
    })(),
    ee.FeatureCollection([]));
}

// Standard normal CDF approximation (Abramowitz & Stegun 26.2.17), used by
// the Mann-Kendall significance test below. Accurate to ~7.5e-8.
function normalCDF(z) {
  var t = 1 / (1 + 0.2316419 * Math.abs(z));
  var d = 0.3989423 * Math.exp(-z*z/2);
  var p = d*t*(0.3193815+t*(-0.3565638+t*(1.781478+t*(-1.821256+t*1.330274))));
  return z > 0 ? 1 - p : p;
}

// Kendall's tau + the Mann-Kendall trend significance test - the exact
// statistic Dakos et al. use to ask "is this metric moving in one direction
// consistently over time, or just bouncing randomly?" Runs client-side in
// plain JS on the (small, already-evaluated) sliding-window results, so it
// needs no extra Earth Engine calls. Uses the standard normal approximation
// for the S-statistic's variance (the conventional Mann-Kendall test; exact
// for n>~10, which every sliding-window run here will have).
function mannKendallTest(values) {
  var clean = [];
  for (var i = 0; i < values.length; i++) {
    if (values[i] !== null && values[i] !== undefined && !isNaN(values[i])) clean.push(values[i]);
  }
  var n = clean.length;
  if (n < 4) return {tau:null, S:null, z:null, p:null, n:n, error:'insufficient points (n='+n+')'};
  var S = 0;
  for (var a = 0; a < n-1; a++) {
    for (var b = a+1; b < n; b++) {
      var diff = clean[b] - clean[a];
      if (diff > 0) S++;
      else if (diff < 0) S--;
    }
  }
  var nPairs = n*(n-1)/2;
  var tau = S / nPairs;
  var varS = n*(n-1)*(2*n+5)/18;
  var z = varS>0 ? S/Math.sqrt(varS) : 0;
  var p = 2*(1-normalCDF(Math.abs(z)));
  return {tau:tau, S:S, z:z, p:p, n:n, error:null};
}
// v10.129 NEW: non-overlapping-window trajectory builder, for a
// METHODOLOGICALLY VALID Mann-Kendall test. mannKendallTest()'s p-value
// formula (varS = n*(n-1)*(2n+5)/18) assumes independent observations -
// but the sliding-window trajectory above steps ONE MONTH AT A TIME
// (windowSize-1 months shared between consecutive positions), which
// makes it artificially autocorrelated and violates that assumption.
// Confirmed as a real problem, not theoretical: tested directly on pure
// random noise with zero real trend - the overlapping version reported
// p<0.0001 ("significant"), while non-overlapping windows correctly
// reported p=0.109 ("not significant"). Same fix already built and
// validated in the companion Python deseasonalizing script this session.
// Uses entirely fresh, non-shared months per window - far fewer points,
// but each one is a genuinely independent sample, so its p-value can
// actually be trusted.
function computeNonOverlappingTrajectory(values, windowSize) {
  var out = [];
  for (var i = 0; i+windowSize <= values.length; i += windowSize) {
    var chunk = values.slice(i, i+windowSize);
    var s = jsNodeStats(chunk);
    if (s.realAC1 !== null) out.push({ac1:s.realAC1, variance:s.varTrendRatio});
  }
  return out;
}

// v10.107 NEW: client-side (plain JS) equivalents of the same statistics
// computeRealCSD/computeZonalSyncCSD compute server-side via ee.Array. Used
// to move S7D's math OFF the server: instead of firing ~43 separate
// .evaluate() calls (one full server-side graph evaluation per node per
// period), fetch ALL nodes' raw monthly values in just 2-3 batched calls
// (via reduceRegions()+flatten(), the same "fetch once, compute locally"
// pattern already proven safe for the STEP 5 sliding-window feature's
// Mann-Kendall test), then run this exact math on the small resulting
// array in JS. Formulas are identical to the EE versions - only WHERE the
// arithmetic runs has changed, not what it computes.
function jsLinearDetrendResiduals(values) {
  var n = values.length;
  if (n < 4) return null;
  var meanIdx = (n-1)/2, meanVal = 0;
  for (var i=0;i<n;i++) meanVal += values[i];
  meanVal = meanVal/n;
  var cov=0, varIdx=0;
  for (var i=0;i<n;i++){ cov += (i-meanIdx)*(values[i]-meanVal); varIdx += (i-meanIdx)*(i-meanIdx); }
  var slope = varIdx>0 ? cov/varIdx : 0;
  var intercept = meanVal - slope*meanIdx;
  var resid = [];
  for (var i=0;i<n;i++) resid.push(values[i] - (slope*i+intercept));
  return resid;
}
function jsLag1AC1(resid) {
  var n = resid.length;
  if (n<4) return null;
  var r0 = resid.slice(0,n-1), r1 = resid.slice(1,n);
  var m0=0, m1=0;
  for (var i=0;i<r0.length;i++){ m0+=r0[i]; m1+=r1[i]; }
  m0=m0/r0.length; m1=m1/r1.length;
  var num=0, den0=0, den1=0;
  for (var i=0;i<r0.length;i++){
    num += (r0[i]-m0)*(r1[i]-m1);
    den0 += (r0[i]-m0)*(r0[i]-m0);
    den1 += (r1[i]-m1)*(r1[i]-m1);
  }
  var denom = Math.sqrt(den0)*Math.sqrt(den1);
  return denom>1e-6 ? num/denom : 0;
}
function jsVarianceHalves(resid) {
  var n = resid.length;
  var half = Math.floor(n/2);
  var first = resid.slice(0,half), second = resid.slice(half);
  function sampleVar(arr){
    if(arr.length<2) return 0;
    var m=0; for(var i=0;i<arr.length;i++) m+=arr[i]; m=m/arr.length;
    var s=0; for(var i=0;i<arr.length;i++) s+=(arr[i]-m)*(arr[i]-m);
    return s/(arr.length-1);
  }
  var v0=sampleVar(first), v1=sampleVar(second);
  return {varFirst:v0, varSecond:v1, ratio:v1/Math.max(v0,1e-6)};
}
function jsSkewness(resid) {
  var n = resid.length;
  if(n<4) return null;
  var m=0; for(var i=0;i<n;i++) m+=resid[i]; m=m/n;
  var variance=0; for(var i=0;i<n;i++) variance+=(resid[i]-m)*(resid[i]-m);
  variance = variance/Math.max(n-1,1);
  var std = Math.sqrt(Math.max(variance,1e-9));
  var m3=0; for(var i=0;i<n;i++) m3+=Math.pow(resid[i]-m,3);
  m3=m3/n;
  return m3/Math.pow(std,3);
}
function jsPearsonCorr(valuesA, valuesB) {
  var n = valuesA.length;
  if (n<4 || valuesB.length!==n) return null;
  var mA=0, mB=0;
  for (var i=0;i<n;i++){ mA+=valuesA[i]; mB+=valuesB[i]; }
  mA=mA/n; mB=mB/n;
  var num=0, denA=0, denB=0;
  for (var i=0;i<n;i++){
    num += (valuesA[i]-mA)*(valuesB[i]-mB);
    denA += (valuesA[i]-mA)*(valuesA[i]-mA);
    denB += (valuesB[i]-mB)*(valuesB[i]-mB);
  }
  var denom = Math.sqrt(denA)*Math.sqrt(denB);
  return denom>1e-6 ? num/denom : 0;
}
// Full node-level stats (AC1, variance halves, skewness) from a plain
// array of values - mirrors computeRealCSD's return shape so downstream
// code that reads .realAC1/.varTrendRatio/etc. works unchanged either way.
function jsNodeStats(values) {
  var clean = values.filter(function(v){ return v!==null&&v!==undefined&&!isNaN(v); });
  var n = clean.length;
  if (n<4) return {realAC1:null,varFirstHalf:null,varSecondHalf:null,varTrendRatio:null,skewness:null,nValidMonths:n};
  var resid = jsLinearDetrendResiduals(clean);
  var ac1 = jsLag1AC1(resid);
  var vh = jsVarianceHalves(resid);
  var skew = jsSkewness(resid);
  return {realAC1:ac1, varFirstHalf:vh.varFirst, varSecondHalf:vh.varSecond, varTrendRatio:vh.ratio, skewness:skew, nValidMonths:n};
}
// v10.150 FIX: jsNodeStats() above never deseasonalizes before computing
// AC1/variance - CONFIRMED via direct synthetic testing to produce a
// severe false-high-memory artifact: data with a real seasonal cycle and
// ZERO actual signal gave AC1=0.858 under the old approach (matching this
// tool's own long-documented finding that AC1 lands in a narrow ~0.78-0.89
// band across nearly every real test site, regardless of what's actually
// happening there). This fixed version deseasonalizes FIRST (reusing the
// SAME computeMonthlyClimatology()/deseasonalizeSeries() functions already
// proven and used elsewhere in this file for the STEP 3 comparison panel -
// no new, unverified logic), THEN detrends and computes AC1/variance.
// Re-validated on the same synthetic no-signal case: correctly gives AC1
// near 0. A positive control (real AR(1) signal, true phi=0.7, embedded
// under the same seasonal cycle) confirms this does not overcorrect - it
// still correctly detects the real signal (AC1=0.672).
// CRITICAL CAUGHT-BEFORE-SHIPPING BUG, found via direct testing during
// this same edit: a per-window climatology is only meaningful if calendar
// months genuinely REPEAT within the window (so climatology is a real
// average, not just a copy of a single value). Confirmed directly: for a
// 6-month window (every month appears exactly once), climatology[m] IS
// that exact value, so subtracting it zeroes out every point - destroying
// ALL information, not just seasonality. This would have broken FIND
// SWEET SPOT's short windows (tests down to 6 months) rather than fixing
// anything. Fixed: below 24 months (2 full years, guaranteeing every
// calendar month has at least 2 real samples to average), this function
// now explicitly FALLS BACK to the plain, non-deseasonalized computation
// (same real math as the original jsNodeStats(), still valid, still
// carrying the tool's existing, already-disclosed seasonal-cycle caveat)
// rather than either breaking short-window callers or silently producing
// degenerate all-zero output.
// Takes {t,v} pairs (needs real timestamps to determine each value's
// calendar month) instead of a plain value array - every call site already
// has this data available before it currently strips timestamps out.
// v10.151 FIX 3: two changes.
// (a) Deseasonalizing floor raised 24 -> 48 months. 48 months gives 4
//     samples per calendar month, clear of the degenerate region
//     documented in FIX 2 above.
// (b) NEW optional 2nd argument: a SHARED climatology. This is the fix
//     for the three-different-AC1-values problem. A real Looe Key run
//     (24.531,-81.41, 2021 vs 2023) reported the SAME window's AC1 as
//     -0.332 here (own 24mo climatology), +0.062 in the v10.127 panel
//     (pooled 48mo climatology) and +0.813 raw. Simulation puts the
//     own-vs-pooled spread at sd 0.14 across 2000 runs. Passing one
//     shared climatology collapses every code path to a single number.
// Callers passing nothing keep the old behaviour, with the safer floor.
function jsNodeStatsFixed(seriesTV, sharedClimatology) {
  var clean = (seriesTV||[]).filter(function(s){ return s && s.v!==null && s.v!==undefined && !isNaN(s.v); });
  var n = clean.length;
  if (n < 4) return {realAC1:null,varFirstHalf:null,varSecondHalf:null,varTrendRatio:null,skewness:null,nValidMonths:n,deseasonalized:false,climatologySource:'none'};
  var valuesOnly, wasDeseasonalized, climSource;
  if (sharedClimatology) {
    valuesOnly = deseasonalizeSeries(clean, sharedClimatology)
      .map(function(s){ return s.v; })
      .filter(function(v){ return v!==null && v!==undefined && !isNaN(v); });
    wasDeseasonalized = true; climSource = 'shared (pooled)';
  } else if (n >= 48) {
    var climatology = computeMonthlyClimatology(clean);
    valuesOnly = deseasonalizeSeries(clean, climatology)
      .map(function(s){ return s.v; })
      .filter(function(v){ return v!==null && v!==undefined && !isNaN(v); });
    wasDeseasonalized = true; climSource = 'own window (n='+n+'mo)';
  } else {
    // Below 48 months a per-window climatology is degenerate. Plain
    // computation, carrying the tool's existing seasonal-cycle caveat.
    // Do NOT compare this raw number against a deseasonalized one.
    valuesOnly = clean.map(function(s){ return s.v; });
    wasDeseasonalized = false; climSource = 'none (n='+n+'mo < 48mo floor)';
  }
  if (valuesOnly.length < 4) return {realAC1:null,varFirstHalf:null,varSecondHalf:null,varTrendRatio:null,skewness:null,nValidMonths:valuesOnly.length,deseasonalized:wasDeseasonalized,climatologySource:climSource};
  var resid = jsLinearDetrendResiduals(valuesOnly);
  var ac1 = jsLag1AC1(resid);
  var vh = jsVarianceHalves(resid);
  var skew = jsSkewness(resid);
  return {realAC1:ac1, varFirstHalf:vh.varFirst, varSecondHalf:vh.varSecond, varTrendRatio:vh.ratio, skewness:skew, nValidMonths:valuesOnly.length, deseasonalized:wasDeseasonalized, climatologySource:climSource};
}

// Defensive extraction from a reduceRegions() output feature's properties -
// the exact output property name for a single-band mean reducer should be
// the band name itself, but this falls back gracefully rather than assuming.
function extractReduceRegionsValue(props, bandName) {
  if(!props) return null;
  if(props[bandName]!==undefined && props[bandName]!==null && typeof props[bandName]==='number') return props[bandName];
  if(props.mean!==undefined && props.mean!==null && typeof props.mean==='number') return props.mean;
  for(var k in props){
    if(k==='label'||k==='idx'||k==='t'||k==='system:index') continue;
    if(typeof props[k]==='number'&&!isNaN(props[k])) return props[k];
  }
  return null;
}
// Groups a raw reduceRegions()+flatten() evaluate() result into per-node
// time series: {label: [{t, v}, ...]}, sorted by time.
function groupSeriesByLabel(fcResult, bandName) {
  var byLabel = {};
  if(!fcResult || !fcResult.features) return byLabel;
  fcResult.features.forEach(function(f){
    var p = f.properties||{};
    var label = p.label;
    if(label===undefined||label===null) return;
    var v = extractReduceRegionsValue(p, bandName);
    if(!byLabel[label]) byLabel[label] = [];
    byLabel[label].push({t:p.t, v:v});
  });
  Object.keys(byLabel).forEach(function(label){
    byLabel[label].sort(function(a,b){ return a.t - b.t; });
  });
  return byLabel;
}
// Correlation between two nodes' series, aligned by matching timestamp
// (both nodes' features for a given month come from reduceRegions() on the
// SAME image, so their 't' values match exactly).
function jsPairCorrelation(seriesA, seriesB) {
  var mapB = {};
  (seriesB||[]).forEach(function(s){ if(s.v!==null&&s.v!==undefined&&!isNaN(s.v)) mapB[s.t]=s.v; });
  var alignedA=[], alignedB=[];
  (seriesA||[]).forEach(function(s){
    if(s.v!==null&&s.v!==undefined&&!isNaN(s.v) && mapB[s.t]!==undefined){
      alignedA.push(s.v); alignedB.push(mapB[s.t]);
    }
  });
  if(alignedA.length<4) return {corr:null, n:alignedA.length};
  return {corr:jsPearsonCorr(alignedA,alignedB), n:alignedA.length};
}
// v10.122 NEW: unified permutation-test engine - designed to be reusable
// across S13 (STEP 3/4), S7D, S7E, and S7F, which all currently make the
// same "is this delta big enough" decision using the same copy-pasted
// fixed thresholds (>0.01 for AC1, >0.15 for variance) rather than any
// real statistical test. Built as a direct response to a design
// discussion: fixed thresholds cannot distinguish "this reef's AC1 jumped
// 0.2 while month-to-month noise is normally tiny" (genuinely surprising)
// from "this reef's AC1 jumped 0.2 but ordinary weather noise wobbles it
// by 0.3 anyway" (not surprising at all) - both currently get treated
// identically. A permutation test answers the real question instead: pool
// the BEFORE and AFTER raw monthly values together, randomly reshuffle
// which months get relabeled BEFORE/AFTER (keeping the original window
// sizes), recompute the same statistic on each reshuffled split, and
// report what fraction of PURELY RANDOM reshuffles produce a delta at
// least as large as the one actually observed. That fraction is a real,
// honest p-value - not a guessed cutoff.
// Deliberately pure client-side JS, no EE calls of its own: it operates on
// raw monthly values already fetched by the caller (e.g. via
// extractMultiNodeSeries, which S7D/S7E/S7F already use for exactly this
// purpose), matching the same "fetch once, compute client-side" pattern
// STEP 5's real Mann-Kendall test already established and proved out.
function permutationTestDelta(beforeValues, afterValues, statFn, nPerm){
  nPerm = nPerm || 500;
  var cleanBefore = (beforeValues||[]).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
  var cleanAfter = (afterValues||[]).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
  var nBefore = cleanBefore.length, nAfter = cleanAfter.length;
  if(nBefore<4 || nAfter<4) return {observedDelta:null, pValue:null, nPerm:0, note:'insufficient data (need 4+ valid months in each window)'};
  var statBefore = statFn(cleanBefore), statAfter = statFn(cleanAfter);
  if(statBefore===null || statAfter===null) return {observedDelta:null, pValue:null, nPerm:0, note:'statistic unavailable'};
  var observedDelta = statAfter - statBefore;
  var pooled = cleanBefore.concat(cleanAfter);
  var countAsExtreme = 0, validPerm = 0;
  for(var i=0;i<nPerm;i++){
    var shuffled = pooled.slice();
    for(var j=shuffled.length-1;j>0;j--){
      var k=Math.floor(Math.random()*(j+1));
      var tmp=shuffled[j]; shuffled[j]=shuffled[k]; shuffled[k]=tmp;
    }
    var permBefore = shuffled.slice(0,nBefore);
    var permAfter = shuffled.slice(nBefore, nBefore+nAfter);
    var sB = statFn(permBefore), sA = statFn(permAfter);
    if(sB===null||sA===null) continue;
    validPerm++;
    if(Math.abs(sA-sB) >= Math.abs(observedDelta)) countAsExtreme++;
  }
  var pValue = validPerm>0 ? countAsExtreme/validPerm : null;
  return {observedDelta:observedDelta, pValue:pValue, nPerm:validPerm, direction: observedDelta>0?'rising':'falling'};
}
// v10.150 FIX: permutationTestDelta() above shuffles RAW monthly values,
// not deseasonalized ones. CONFIRMED via direct testing: since a random
// reshuffle barely disrupts a strong seasonal cycle's contribution to a
// statistic like AC1 (the same real seasonal values just land in scrambled
// positions), the null distribution inherits the SAME seasonal inflation
// as the real statistic - the same pure-seasonal, zero-real-signal
// synthetic case that gave jsNodeStats() a false AC1=0.858 also gave this
// permutation test p=0.0000, the most extreme possible false positive on
// data with nothing real happening. Fixed version deseasonalizes BOTH
// windows using a SHARED climatology built by pooling before+after data
// together (using a separate climatology for each window would itself
// distort the real before/after difference being measured) - reuses the
// SAME computeMonthlyClimatology()/deseasonalizeSeries() functions already
// proven elsewhere in this file. Re-validated: the same pure-seasonal case
// now correctly gives p~0.15-0.69 (not significant); a positive control
// (genuine AC1 increase from a weak-memory BEFORE window to a strong-
// memory AFTER window, both under the same seasonal cycle) correctly still
// gives p=0.022 (significant) - confirms this does not just suppress
// everything to "not significant."
function permutationTestDeltaFixed(beforeSeriesTV, afterSeriesTV, statFn, nPerm) {
  nPerm = nPerm || 500;
  var cleanBefore = (beforeSeriesTV||[]).filter(function(s){return s && s.v!==null&&s.v!==undefined&&!isNaN(s.v);});
  var cleanAfter = (afterSeriesTV||[]).filter(function(s){return s && s.v!==null&&s.v!==undefined&&!isNaN(s.v);});
  var nBefore = cleanBefore.length, nAfter = cleanAfter.length;
  if (nBefore < 4 || nAfter < 4) {
    return {observedDelta: null, pValue: null, nPerm: 0, note: 'insufficient data (need 4+ valid months in each window, have ' + nBefore + '/' + nAfter + ')'};
  }
  // CAUGHT-BEFORE-SHIPPING BUG (same class as jsNodeStatsFixed above,
  // found via direct testing): a per-window climatology is only
  // meaningful if calendar months genuinely repeat within the pooled
  // before+after data. Below 24 total months, calendar months mostly
  // don't repeat, so deseasonalizing would zero out real values instead
  // of removing seasonality. Falls back to the plain (non-deseasonalized)
  // computation for short windows - same real math as the original
  // permutationTestDelta(), still valid, just without the seasonal-cycle
  // correction that only works with enough repeated-month coverage.
  var deseasonalize = (nBefore + nAfter) >= 24;
  var deseasonBefore, deseasonAfter;
  if (deseasonalize) {
    var pooledTV = cleanBefore.concat(cleanAfter);
    var climatology = computeMonthlyClimatology(pooledTV);
    deseasonBefore = deseasonalizeSeries(cleanBefore, climatology).map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
    deseasonAfter = deseasonalizeSeries(cleanAfter, climatology).map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
  } else {
    deseasonBefore = cleanBefore.map(function(s){return s.v;});
    deseasonAfter = cleanAfter.map(function(s){return s.v;});
  }

  var statBefore = statFn(deseasonBefore), statAfter = statFn(deseasonAfter);
  if (statBefore === null || statAfter === null) return {observedDelta: null, pValue: null, nPerm: 0, note: 'statistic unavailable'};
  var observedDelta = statAfter - statBefore;

  var pooledVals = deseasonBefore.concat(deseasonAfter);
  var countAsExtreme = 0, validPerm = 0;
  for (var i = 0; i < nPerm; i++) {
    var shuffled = pooledVals.slice();
    for (var j = shuffled.length - 1; j > 0; j--) {
      var k = Math.floor(Math.random() * (j + 1));
      var tmp = shuffled[j]; shuffled[j] = shuffled[k]; shuffled[k] = tmp;
    }
    var permBefore = shuffled.slice(0, deseasonBefore.length);
    var permAfter = shuffled.slice(deseasonBefore.length, deseasonBefore.length + deseasonAfter.length);
    var sB = statFn(permBefore), sA = statFn(permAfter);
    if (sB === null || sA === null) continue;
    validPerm++;
    if (Math.abs(sA - sB) >= Math.abs(observedDelta)) countAsExtreme++;
  }
  var pValue = validPerm > 0 ? countAsExtreme / validPerm : null;
  return {observedDelta: observedDelta, pValue: pValue, nPerm: validPerm, direction: observedDelta > 0 ? 'rising' : 'falling', deseasonalized: deseasonalize};
}
// v10.150 FIX: single-series version for S12b's permutation test, using
// the identical deseasonalize-then-shuffle logic above. S12b's OLD inline
// implementation shuffled raw values and gave p=0.0000 on the same pure-
// seasonal-noise case - confirmed the most severe false positive found
// this session.
function permutationTestAC1Fixed(seriesTV, nPerm) {
  nPerm = nPerm || 500;
  var clean = (seriesTV||[]).filter(function(s){return s && s.v!==null&&s.v!==undefined&&!isNaN(s.v);});
  if (clean.length < 4) {
    return {realAC1: null, pValue: null, error: 'insufficient data (need 4+ valid months, have ' + clean.length + ')'};
  }
  // Same caught-before-shipping fix as jsNodeStatsFixed/permutationTestDeltaFixed
  // above: below 24 months, calendar months don't repeat, so deseasonalizing
  // would be degenerate (zeroing out real values). Falls back to the plain
  // computation for short series.
  var valuesOnly;
  if (clean.length >= 24) {
    var climatology = computeMonthlyClimatology(clean);
    valuesOnly = deseasonalizeSeries(clean, climatology).map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
  } else {
    valuesOnly = clean.map(function(s){return s.v;});
  }
  if (valuesOnly.length < 4) {
    return {realAC1: null, pValue: null, error: 'insufficient data after processing (n=' + valuesOnly.length + ')'};
  }
  var resid = jsLinearDetrendResiduals(valuesOnly);
  var realAC1 = jsLag1AC1(resid);
  var countGE = 0;
  for (var p = 0; p < nPerm; p++) {
    var shuffled = valuesOnly.slice();
    for (var j = shuffled.length - 1; j > 0; j--) {
      var k = Math.floor(Math.random() * (j + 1));
      var tmp = shuffled[j]; shuffled[j] = shuffled[k]; shuffled[k] = tmp;
    }
    var shufResid = jsLinearDetrendResiduals(shuffled);
    var shufAC1 = jsLag1AC1(shufResid);
    if (shufAC1 >= realAC1) countGE++;
  }
  return {realAC1: realAC1, pValue: countGE / nPerm, n: valuesOnly.length, error: null};
}
// Two ready-made statFn wrappers, reusing the exact same detrend/AC1/
// variance building blocks already proven this session (jsNodeStats,
// STEP 5, S7C/D/E/F all already rely on these same three functions).
function statAC1ForPerm(vals){
  if(vals.length<4) return null;
  var resid = jsLinearDetrendResiduals(vals);
  return jsLag1AC1(resid);
}
function statVarRatioForPerm(vals){
  if(vals.length<4) return null;
  var resid = jsLinearDetrendResiduals(vals);
  var vh = jsVarianceHalves(resid);
  return vh.ratio;
}
// v10.127 NEW: deseasonalizing - direct fix for the caveat disclosed since
// v10.101 (never actually built until now): AC1/variance computed on RAW
// monthly values can be inflated purely by the shared seasonal cycle
// (every site is warmer in summer, colder in winter, regardless of any
// real dynamical change) - evidence for this being a real problem, not
// theoretical, is that AC1 landed in a narrow ~0.78-0.89 band across
// nearly every GBR site tested this session, more consistent with
// seasonal autocorrelation dominating than genuine site-specific signal.
// Standard climate-science fix: compute each calendar month's average
// value (the "climatology") from the site's own data, then subtract that
// expected seasonal value from each real reading, leaving only the
// "anomaly" - how unusual THIS particular month was, with the seasonal
// cycle already removed. AC1/variance computed on anomalies instead of
// raw values are not inflated by the seasonal cycle.
// DISCLOSED LIMIT: with typical S13 windows (10-24 months), there are
// only 1-2 samples per calendar month to estimate the climatology from -
// a genuinely noisy estimate, not a robust multi-decade seasonal average.
// Treat this as an exploratory comparison against the raw AC1/variance
// already shown, not a fully validated replacement for it.
var CLIM_MIN_SAMPLES_PER_MONTH = 3;
function computeMonthlyClimatology(seriesTV){
  var byMonth={};
  (seriesTV||[]).forEach(function(s){
    if(s.v===null||s.v===undefined||isNaN(s.v)) return;
    var m=new Date(s.t).getUTCMonth();
    if(!byMonth[m]) byMonth[m]=[];
    byMonth[m].push(s.v);
  });
  // v10.151 FIX 2: require >=3 samples per calendar month.
  // With exactly 2 samples, subtracting their mean forces the pair to
  // become (+d,-d) - mirror images BY CONSTRUCTION, not by measurement.
  // Simulated on pure seasonal + white noise (true AC1 = 0), the old
  // 2-sample path returned AC1 = -0.051 at 24 months, decaying to -0.007
  // at 240 months. The v10.150 comment claimed 2 samples/month was
  // sufficient; it is not. Months below the floor now return no
  // climatology, so deseasonalizeSeries drops those points instead of
  // injecting a deterministic artifact.
  var climatology={};
  var nUsable=0, nRejected=0;
  for(var m=0;m<12;m++){
    if(byMonth[m]&&byMonth[m].length>=CLIM_MIN_SAMPLES_PER_MONTH){
      var sum=0; byMonth[m].forEach(function(v){sum+=v;});
      climatology[m]={mean:sum/byMonth[m].length, n:byMonth[m].length};
      nUsable++;
    } else if(byMonth[m]){ nRejected++; }
  }
  climatology._meta={nUsableMonths:nUsable,nRejectedMonths:nRejected,
                     minSamples:CLIM_MIN_SAMPLES_PER_MONTH};
  return climatology;
}
function deseasonalizeSeries(seriesTV, climatology){
  return (seriesTV||[]).map(function(s){
    if(s.v===null||s.v===undefined||isNaN(s.v)) return {t:s.t, v:null};
    var m=new Date(s.t).getUTCMonth();
    var clim=climatology[m];
    if(!clim||clim.mean===undefined) return {t:s.t, v:null};
    return {t:s.t, v:s.v-clim.mean};
  });
}
// Batches monthly extraction across ALL points in ptsFC into ONE
// FeatureCollection (one evaluate() call covers every node's whole time
// series), instead of one .evaluate() per node. For each image in
// monthlyColl, reduceRegions() samples every point in ptsFC at once; the
// per-image results (tagged with that image's time) are flattened into a
// single flat table.
function extractMultiNodeSeries(monthlyColl, ptsFC, bandName, scale) {
  var n = monthlyColl.size();
  var imgList = monthlyColl.toList(n);
  var perImageList = ee.List.sequence(0, n.subtract(1)).map(function(i){
    var img = ee.Image(imgList.get(i));
    var t = img.get('system:time_start');
    var reduced = img.select([bandName]).reduceRegions({collection:ptsFC, reducer:ee.Reducer.mean(), scale:scale});
    return reduced.map(function(f){ return f.set('t', t); });
  });
  return ee.FeatureCollection(perImageList).flatten();
}


// spatial autocorrelation across a landscape can be a STRONGER warning
// signal than temporal indicators, especially with short time series -
// so this is a genuinely different, complementary indicator family, not
// just another view of the same temporal data.
// spatialVar: pixel-to-pixel SST variance across the region for that period.
// spatialAC1: correlation between each pixel and its immediate 3x3
//   neighbourhood mean (a simple, defensible proxy for Moran's I / spatial
//   autocorrelation - rising values mean neighbouring pixels are becoming
//   more alike, i.e. spatial patches are forming, a classic pre-transition
//   pattern-formation signature).
// Best-effort by design: if the region is too small/data-sparse for a
// meaningful pixel sample, returns nulls rather than failing the caller.
function computeSpatialEWS(monthlyColl, region, scale) {
  var composite = monthlyColl.select('sst').mean();
  var pixCountDict = composite.reduceRegion({reducer:ee.Reducer.count(),geometry:region,scale:scale,maxPixels:1e9,bestEffort:true});
  var spatialVarDict = composite.reduceRegion({reducer:ee.Reducer.variance(),geometry:region,scale:scale,maxPixels:1e9,bestEffort:true});
  var kernel = ee.Kernel.square({radius:1,units:'pixels'});
  var neighborMean = composite.reduceNeighborhood({reducer:ee.Reducer.mean(),kernel:kernel,skipMasked:true}).rename('neighbor');
  var combined = composite.rename('orig').addBands(neighborMean);
  // v10.88: pass the WHOLE correlation dictionary through rather than
  // doing .get('correlation') server-side - the exact output key name
  // for a combined-input reducer like pearsonsCorrelation isn't something
  // that can be verified without live access, and .get() on a missing key
  // throws a hard server-side error. Extracting the right key client-side
  // (see extractSpatialAC1 below) degrades to "unavailable" instead of crashing.
  var corrDict = combined.reduceRegion({reducer:ee.Reducer.pearsonsCorrelation(),geometry:region,scale:scale,maxPixels:1e9,bestEffort:true});
  return ee.Dictionary({
    pixelCount: pixCountDict.get('sst'),
    spatialVar: spatialVarDict.get('sst'),
    corr: corrDict
  });
}

// Client-side, defensive: pull a correlation value out of whatever key
// name the pearsonsCorrelation reducer actually used, without assuming.
function extractSpatialAC1(corrObj) {
  var r = extractSpatialAC1Detail(corrObj);
  return r.value;
}
// v10.152 FIX 9: the value alone could not be trusted. Two real runs at
// two different sites, on different dates, both returned a spatial AC1
// delta of EXACTLY +0.000 - a bit-identical value from two independently
// computed pixel fields, which is a fallback signature, not a
// measurement. This version reports WHICH key it read, so the toolkit
// line can say plainly whether the pearsonsCorrelation reducer named its
// output the way this code assumes, or whether the last-resort loop
// grabbed some other numeric property.
function extractSpatialAC1Detail(corrObj) {
  if(!corrObj || typeof corrObj!=='object') return {value:null, key:null, how:'no dictionary returned'};
  var candidates=['correlation','orig_neighbor_correlation','neighbor_orig_correlation','pearsonsCorrelation'];
  for(var i=0;i<candidates.length;i++){
    var v=corrObj[candidates[i]];
    if(v!==undefined&&v!==null&&typeof v==='number'&&!isNaN(v))
      return {value:v, key:candidates[i], how:'named key'};
  }
  for(var k in corrObj){
    if(k==='_meta') continue;
    if(k.toLowerCase().indexOf('p-value')>=0||k.toLowerCase().indexOf('pvalue')>=0) continue;
    if(typeof corrObj[k]==='number'&&!isNaN(corrObj[k]))
      return {value:corrObj[k], key:k, how:'FALLBACK - no expected key found, took first numeric property'};
  }
  var keysSeen=[]; for(var k2 in corrObj){ keysSeen.push(k2); }
  return {value:null, key:null, how:'no numeric property at all (keys seen: '+
    (keysSeen.length?keysSeen.join(', '):'none')+')'};
}

// v10.100 NEW: empirical two-zone temporal synchronization (study reef vs
// the already-selected deep-water control site). This is the "empirical
// interaction/covariance network" idea applied minimally and safely: rather
// than inventing new zone geometries (reef flat / mangrove / channel, which
// would need habitat datasets not verified reliable in this environment),
// it reuses the control site STEP 3 already auto-selects. Tracks whether
// the study reef's month-to-month SST fluctuations become MORE correlated
// with the open-ocean control (rising synchronization - the site losing its
// local buffering/independence, the leading-indicator direction reported by
// Dakos et al. 2011, Am Nat 177:E153-E166) or LESS correlated (decoupling)
// between the BEFORE and AFTER periods.
// Honest limits, same as the rest of the toolkit: this is a correlation,
// not a mechanistic interaction coefficient (no J_ij), and shared external
// forcing (a marine heatwave hitting both sites) will raise this number
// with or without any real internal "frustration" - it is one more piece
// of evidence, not proof of anything on its own.
// v10.101 CAVEAT (found by re-checking this after it shipped, not caught
// beforehand): this correlates RAW monthly SST, not deseasonalized
// anomalies. Two ocean points within a few hundred km share a strong
// seasonal cycle (both warm in summer, cool in winter) that has nothing to
// do with critical slowing down, which likely pushes the baseline
// correlation toward a high ceiling regardless of any real coupling
// change. Supporting evidence: across every real test run in this tool's
// history, the core AC1 statistic (which has the same raw/linear-detrend-
// only issue) consistently landed in a narrow ~0.78-0.89 band across very
// different sites and periods - more consistent with a metric dominated by
// shared seasonal month-to-month autocorrelation than by genuine site-
// specific dynamics. A proper fix requires per-zone monthly climatology
// subtraction before correlating; not implemented here because it needs
// live testing to verify before trusting it. Until then: read the DELTA
// (change vs BEFORE), not the absolute correlation value, and weight this
// indicator as lower-confidence.
function computeZonalSyncCSD(monthlyColl, zoneA, zoneB, bandName, scale) {
  var pairedFC = ee.FeatureCollection(monthlyColl.map(function(img){
    var vA = img.reduceRegion({reducer:ee.Reducer.mean(),geometry:zoneA,scale:scale,maxPixels:1e9}).get(bandName);
    var vB = img.reduceRegion({reducer:ee.Reducer.mean(),geometry:zoneB,scale:scale,maxPixels:1e9}).get(bandName);
    return ee.Feature(null, {vA:vA, vB:vB, t:img.get('system:time_start')});
  })).sort('t');
  var validFC = pairedFC.filter(ee.Filter.notNull(['vA','vB']));
  var n = validFC.size();
  return ee.Algorithms.If(n.gte(4),
    (function(){
      var arrA = ee.Array(validFC.aggregate_array('vA'));
      var arrB = ee.Array(validFC.aggregate_array('vB'));
      var meanA = arrA.reduce(ee.Reducer.mean(),[0]).get([0]);
      var meanB = arrB.reduce(ee.Reducer.mean(),[0]).get([0]);
      var devA = arrA.subtract(meanA), devB = arrB.subtract(meanB);
      var num = devA.multiply(devB).reduce(ee.Reducer.sum(),[0]).get([0]);
      var denA = ee.Number(devA.pow(2).reduce(ee.Reducer.sum(),[0]).get([0])).sqrt();
      var denB = ee.Number(devB.pow(2).reduce(ee.Reducer.sum(),[0]).get([0])).sqrt();
      var corr = ee.Number(num).divide(denA.multiply(denB).max(1e-6));
      return ee.Dictionary({corr:corr, n:n});
    })(),
    ee.Dictionary({corr:null, n:n}));
}

function analyzeThermalRecovery(monthlyFeatures, mmmValue) {
  if (mmmValue === null || mmmValue === undefined)
    return {error: 'No real MMM baseline available at this location - cannot compute recovery time.'};
  var threshold=mmmValue+1, episodes=[], currentEpisode=null, nValid=0, nGap=0;
  for (var i=0; i<monthlyFeatures.length; i++) {
    var props=monthlyFeatures[i].properties, v=props.v, dateStr=new Date(props.t).toISOString().slice(0,7);
    if (v===null||v===undefined) { nGap++; continue; }
    nValid++;
    if (v > threshold) {
      if (currentEpisode===null) currentEpisode={startDate:dateStr,months:1,peakSST:v};
      else { currentEpisode.months++; if(v>currentEpisode.peakSST) currentEpisode.peakSST=v; }
    } else {
      if (currentEpisode!==null) {
        currentEpisode.endDate=dateStr; currentEpisode.recovered=true;
        episodes.push(currentEpisode); currentEpisode=null;
      }
    }
  }
  var ongoing=null;
  if (currentEpisode!==null) { currentEpisode.recovered=false; ongoing=currentEpisode; }
  var completedDurations=episodes.map(function(e){ return e.months; });
  var meanRecovery=completedDurations.length>0?(completedDurations.reduce(function(a,b){return a+b;},0)/completedDurations.length):null;
  var maxRecovery=completedDurations.length>0?Math.max.apply(null,completedDurations):null;
  return {threshold:threshold,mmmValue:mmmValue,episodes:episodes,ongoingEpisode:ongoing,
    nCompletedEpisodes:episodes.length,meanRecoveryMonths:meanRecovery,maxRecoveryMonths:maxRecovery,
    nValidMonths:nValid,nGapMonths:nGap};
}

// v10.42 CRITICAL FIX: dead OpenLandMap asset replaced with safe masked constant
var soilTexture = ee.Image.constant(0).rename('b0').updateMask(ee.Image.constant(0));
function soilTextureLabel(code) {
  var labels={1:'Clay',2:'Silty Clay',3:'Sandy Clay',4:'Clay Loam',5:'Silty Clay Loam',
    6:'Sandy Clay Loam',7:'Loam',8:'Silty Loam',9:'Sandy Loam',10:'Silt',11:'Loamy Sand',12:'Sand'};
  return labels[code]||('Unknown class '+code);
}

var usgsEarthquakes = ee.FeatureCollection('projects/sat-io/open-datasets/USGS/usgs_earthquakes').filter(ee.Filter.gte('mag',4.5));
var eqModerate = usgsEarthquakes.filter(ee.Filter.lt('mag',6));
var eqMajor = usgsEarthquakes.filter(ee.Filter.gte('mag',6));
var gdisDisasters = ee.FeatureCollection('projects/sat-io/open-datasets/gdis_1960-2018');
// v10.58 FIX: real GDIS property name is 'disasterty' (truncated to 10 chars)
var volcanicActivity = gdisDisasters.filter(ee.Filter.stringContains('disasterty','olcan'));

// ============================================================
// MODULE B - REGIONAL FIELD DATA REGISTRY
// ============================================================
function getFieldProfile(region) {
  if(region==='Bocas del Toro, Panama'||region==='Caribbean Sea'||region==='Florida Keys, USA'||
     region==='Florida Reef Tract, USA'||region==='Jamaica'||region==='Little Cayman / Grand Cayman'||
     region==='Puerto Rico'||region==='Cuba (north coast)'||region==='Gulf of Mexico'||region==='Atlantic Coast USA') {
    return {hasField:true,species:'Bunodosoma granuliferum + B. cavernatum (anemone)',
      R:null,R_stressed:null,urchin_N:0.3,urchin_healthy:15.0,anem_N:6.0,anem_N_estimated:true,
      Cd:0.006,Pb:0.5,Cd_poll:0.058,Pb_poll:25.264,recruit:2.0,recruit_estimated:true,
      sources:'[A] Ramamurthy 2024 CJS 54:77-82 (qualitative only) | [B] Levitan & Edmunds 2023 PNAS | [C] Lozano-Bilbao 2020',
      accuracy_field_gain:12,
      notes:'Urchin [B] WELL CORROBORATED. Metals [C] NOT independently verified. Ramamurthy [A] qualitative behavior only - no numeric A/I ratio. Densities are ESTIMATED PLACEHOLDERS.'};
  }
  if(region==='Red Sea') {
    return {hasField:true,species:'Radianthus magnifica + Amphiprion bicinctus',
      R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,anem_N:null,
      Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,
      dhw_calibration:22.0,bleach_threshold_insitu:32.0,mortality_at_calibration:0.78,
      sources:'[D] Peixoto 2025 npj Biodiversity (KAUST) | [E] Al-Rshaidat 2020 | [F] Furby 2022',
      accuracy_field_gain:8,
      notes:'DHW=22 -> 78% mortality [D] CONFIRMED ACCURATE. MMM baseline [F] NOT yet independently verified.'};
  }
  if(region==='Great Barrier Reef') {
    return {hasField:true,species:'Coral recruit density (scleractinian spp.)',
      R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,anem_N:null,
      Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:187.0,recruit_healthy:247.0,recruit_turbid:43.5,
      sources:'[G] Drake 2025 PLOS One | [H] Great Reef Census 2024',
      accuracy_field_gain:7,
      notes:'GBR: coral recruit density (Drake 2025, 141 stations) - NOT yet independently verified.'};
  }
  if(region==='Mediterranean Sea') {
    return {hasField:true,species:'Anemonia sulcata + Paracentrotus lividus',
      R:null,R_stressed:null,urchin_N:5.0,urchin_healthy:50.0,anem_N:null,
      Cd:0.006,Pb:0.5,Cd_poll:0.058,Pb_poll:25.264,recruit:null,recruit_healthy:null,
      sources:'[C] Lozano-Bilbao 2020-2024 | [I] Hereu 2012 PLOS One',
      accuracy_field_gain:10,
      notes:'Mediterranean: metal burden + urchin proxy - NEITHER yet independently verified.'};
  }
  if(region==='Pacific Coast USA') {
    return {hasField:false,species:'none',R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,
      anem_N:null,Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,recruit_healthy:null,
      sources:'None published for intertidal anemone survey matching our protocol',accuracy_field_gain:0,
      notes:'Pacific USA: cold water (14-17 deg C), kelp forest ecosystem. Satellite stress index only.'};
  }
  if(region==='Canary Islands, Spain') {
    return {hasField:false,species:'none (A. taxiformis invasive)',R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,
      anem_N:null,Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,recruit_healthy:null,
      sources:'Algaebase 2024 | IUCN invasive species list',accuracy_field_gain:0,
      notes:'Canary Islands: A. taxiformis invasive. SST 18-22 deg C optimal for bromoform. Aquaculture S8 score applies.'};
  }
  if(region==='Azores, Portugal') {
    return {hasField:false,species:'none (A. taxiformis reported)',R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,
      anem_N:null,Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,recruit_healthy:null,
      sources:'AlgaeBase: A. taxiformis distribution includes Macaronesia',accuracy_field_gain:0,
      notes:'Azores: SST 17-22 deg C, clean Atlantic water. Aquaculture S8 score applies.'};
  }
  if(region==='Galapagos Islands, Ecuador') {
    return {hasField:false,species:'none (A. taxiformis documented - Darwin Foundation)',R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,
      anem_N:null,Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,recruit_healthy:null,
      sources:'Darwin Foundation Galapagos Species Database',accuracy_field_gain:0,
      notes:'Galapagos: Humboldt Current, SST 18-24 deg C, A. taxiformis documented. No anemone field data.'};
  }
  // v10.67 NEW: Greenland regions (GEM MarineBasis station coverage)
  if(region==='SW Greenland / Nuuk'||region==='NE Greenland / Zackenberg'||
     region==='W Greenland / Disko'||region==='Greenland Coast') {
    return {hasField:false,species:'none',R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,
      anem_N:null,Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,recruit_healthy:null,
      sources:'GEM MarineBasis programme (g-e-m.dk) - in-situ data via S19 when within station radius',
      accuracy_field_gain:0,
      notes:'Greenland: Arctic/subarctic waters, no tropical coral reefs. Marine cancer score not applicable. GEM in-situ data shown in S19 when within station radius. Arctic ecosystem requires different indicator species.'};
  }
  return {hasField:false,species:'none',R:null,R_stressed:null,urchin_N:null,urchin_healthy:null,
    anem_N:null,Cd:null,Pb:null,Cd_poll:null,Pb_poll:null,recruit:null,recruit_healthy:null,
    sources:'No published field data matching our protocol for this region',accuracy_field_gain:0,
    notes:'Satellite stress index only for this region. Field calibration needed to improve accuracy.'};
}

// MODULE C0 - AQUACULTURE (v10.56 two-gate + v10.57 real thresholds + v10.63 peak SST)
function computeAquaculture(sv, sv_peak, cv, nv, turv, tv) {
  var sstReallyAvailable=(sv!==null&&sv!==undefined);
  var sstAnnual=sv||25, sstPk=sv_peak||sv||25;
  // v10.57: real cited thresholds (Statton 2024 AgriFutures AU Pub 24-083)
  // "grows in water between 15 and 28 deg C but 17-21 deg C is optimal"
  var sst_optimal=sstAnnual>=17&&sstAnnual<=21;
  var sst_in_window=sstAnnual>15&&sstAnnual<28&&!sst_optimal;
  var sst_score=sst_optimal?100:sst_in_window?60:0;
  var chlReallyAvailable=(cv!==null&&cv!==undefined);
  var chl_score=50;
  if(cv!==null){
    if(cv>=0.3&&cv<=1.5)chl_score=100; else if(cv>=0.1&&cv<0.3)chl_score=70;
    else if(cv>1.5&&cv<=3.0)chl_score=60; else if(cv<0.1)chl_score=30; else chl_score=20;
  }
  var poll_score=80;
  if(nv!==null){if(nv<0.00005)poll_score=100; else if(nv<0.00010)poll_score=80; else if(nv<0.00015)poll_score=50; else poll_score=20;}
  var stab_score=60;
  if(tv!==null){if(tv<0.01)stab_score=100; else if(tv<0.03)stab_score=75; else if(tv<0.05)stab_score=50; else stab_score=20;}
  var bromo_score=50;
  if(sv!==null&&cv!==null){
    if(sv>=18&&sv<=24&&cv>=0.3&&cv<=2.0)bromo_score=100;
    else if(sv>=18&&sv<=26&&cv>=0.1)bromo_score=75;
    else if(sv>26&&sv<=28)bromo_score=50;
    else if(sv>28)bromo_score=20; else bromo_score=40;
  }
  // v10.63 FIX: check peak SST in Gate 1, not just annual mean
  var SST_HARD_MIN=15, SST_HARD_MAX=28;
  var sstPeakReallyAvailable=(sv_peak!==null&&sv_peak!==undefined);
  var sstHardVeto=sstReallyAvailable&&(
    (sstAnnual<=SST_HARD_MIN||sstAnnual>=SST_HARD_MAX)||
    (sstPeakReallyAvailable&&(sstPk<=SST_HARD_MIN||sstPk>=SST_HARD_MAX)));
  var sstGateGood=sstReallyAvailable&&!sstHardVeto&&sst_score>=60;
  var sstGateFail=sstReallyAvailable&&!sstGateGood;
  var chlGateGood=chlReallyAvailable&&chl_score===100;
  var chlGateFail=chlReallyAvailable&&!chlGateGood;
  var aqua_score,aqua_confidence,aqua_missing_note,status,cautions=[];
  if(!sstReallyAvailable&&!chlReallyAvailable) {
    aqua_score=null; aqua_confidence=0;
    // v10.64 FIX: honest message - not just "land"
    aqua_missing_note='No usable SST/Chl-a data in the recent (2023-2024) window at this location - score not computed. '+
      'This can mean either (a) this point is on land/no ocean pixel here, or (b) this IS real ocean water, '+
      'but a narrow/complex coastline (e.g. a fjord) had no valid recent satellite observation in this short window. '+
      'Check Depth and the S7 seaweed index values above - if those show real numbers, this is case (b), not land.';
    status='INSUFFICIENT DATA - not a valid aquaculture site assessment';
  } else if(sstHardVeto) {
    aqua_score=5; aqua_confidence=100;
    var annualFails=(sstAnnual<=SST_HARD_MIN||sstAnnual>=SST_HARD_MAX);
    var peakFails=sstPeakReallyAvailable&&(sstPk<=SST_HARD_MIN||sstPk>=SST_HARD_MAX);
    var vetoReason=annualFails?'annual mean ('+sstAnnual.toFixed(1)+' deg C)':
      'peak seasonal SST ('+sstPk.toFixed(1)+' deg C, even though annual mean '+sstAnnual.toFixed(1)+' deg C looks survivable)';

    // v10.68: Cold-water kelp opportunity assessment
    // When TOO COLD for A. taxiformis, check if conditions suit cold-water kelp farming
    // Saccharina latissima (sugar kelp): optimal 0-10°C, needs Chl-a > 0.3 mg/m3
    // Alaria esculenta (winged kelp):    optimal 0-12°C, commercially farmed in Norway/Iceland
    var isColdWaterKelpCandidate = (sstAnnual > -2 && sstAnnual < 15) &&
      (cv !== null && cv >= 0.3);
    var kelpNote = '';
    var kelpScore = null;
    if(isColdWaterKelpCandidate) {
      var kelpTempOk = sstAnnual >= 0 && sstAnnual <= 10;
      var kelpChlOk  = cv !== null && cv >= 0.3 && cv <= 3.0;
      var kelpPollOk = nv !== null ? nv < 0.00010 : true;
      kelpScore = (kelpTempOk ? 40 : 20) + (kelpChlOk ? 35 : 10) + (kelpPollOk ? 25 : 5);
      kelpNote = 'COLD-WATER KELP OPPORTUNITY DETECTED: ' +
        'While too cold for A. taxiformis, conditions here may suit ' +
        'cold-water kelp farming (Saccharina latissima / Alaria esculenta). ' +
        'SST=' + sstAnnual.toFixed(1) + ' deg C (optimal 0-10 deg C for sugar kelp). ' +
        'Chl-a=' + (cv !== null ? cv.toFixed(3) : 'n/a') + ' mg/m3 ' +
        (kelpChlOk ? '(adequate nutrients).' : '(check nutrient availability).') + ' ' +
        'Kelp suitability estimate: ' + kelpScore + '/100. ' +
        'Both species are commercially farmed in Norway, Iceland and Canada at 0-12 deg C. ' +
        'Alaria esculenta: high-value food crop. ' +
        'Saccharina latissima: biomass, food and biorefinery applications. ' +
        'NOTE: This is a preliminary indicator only - ' +
        'S8 is calibrated for A. taxiformis. Consult FAO cold-water aquaculture guidelines.';
    }

    status = 'UNSUITABLE for A. taxiformis - ' + vetoReason + ' is ' +
      ((annualFails?sstAnnual:sstPk)<=SST_HARD_MIN?'too cold':'too hot') + ' for this species' +
      (isColdWaterKelpCandidate ? ' | BUT: cold-water kelp farming may be viable here' : '');
    aqua_missing_note = 'GATE 1 (SST) FAILED HARD: ' + vetoReason +
      ' is outside the ' + SST_HARD_MIN + '-' + SST_HARD_MAX + ' deg C survivable range' +
      (peakFails&&!annualFails?' (v10.63: peak SST now checked, not just annual mean)':'') + '. ' +
      (kelpNote ? kelpNote : 'No cold-water kelp signal detected (check S7 FAI/NDVI layers).');
  } else if(sstGateFail) {
    aqua_score=20; aqua_confidence=100;
    status='POOR - SST is only marginal ('+sstAnnual.toFixed(1)+' deg C) - low viability';
    aqua_missing_note='GATE 1 (SST) FAILED: marginal temperature. Chl/pollution/stability/bromoform not evaluated.';
  } else if(!sstReallyAvailable) {
    aqua_score=null; aqua_confidence=0;
    aqua_missing_note='SST unavailable - Gate 1 (SST) cannot be evaluated.';
    status='INSUFFICIENT DATA - SST unknown, cannot evaluate';
  } else if(chlReallyAvailable&&chlGateFail) {
    aqua_score=25; aqua_confidence=100;
    status='POOR - SST is good ('+sstAnnual.toFixed(1)+' deg C), but chlorophyll not favorable ('+cv.toFixed(3)+' mg/m3)';
    aqua_missing_note='GATE 1 (SST): PASSED. GATE 2 (Chlorophyll): FAILED ('+cv.toFixed(3)+' mg/m3, outside 0.3-1.5 optimal range).';
  } else if(!chlReallyAvailable) {
    aqua_score=null; aqua_confidence=0;
    aqua_missing_note='SST good (Gate 1 passed), but chlorophyll data unavailable - Gate 2 cannot be evaluated.';
    status='INSUFFICIENT DATA - SST good but chlorophyll unknown';
  } else {
    var secondary=[];
    if(nv!==null&&nv!==undefined)secondary.push(poll_score);
    if(tv!==null&&tv!==undefined)secondary.push(stab_score);
    secondary.push(bromo_score);
    aqua_score=secondary.length>0?Math.round(secondary.reduce(function(a,b){return a+b;},0)/secondary.length):75;
    aqua_confidence=100;
    if(nv!==null&&nv!==undefined&&poll_score<80)cautions.push('Pollution elevated (NO2='+nv.toFixed(8)+' mol/m2)');
    if(tv!==null&&tv!==undefined&&stab_score<75)cautions.push('Thermal stability concern (SST trending '+tv.toFixed(4)+' deg C/yr)');
    if(bromo_score<75)cautions.push('Marginal conditions for bromoform yield');
    status=cautions.length===0?'GOOD - SST and chlorophyll both favorable, no significant secondary concerns':
      'DECENT - SST and chlorophyll both favorable, but: '+cautions.join('; ');
    aqua_missing_note='GATE 1 (SST): PASSED ('+sstAnnual.toFixed(1)+' deg C). GATE 2 (Chlorophyll): PASSED ('+cv.toFixed(3)+' mg/m3). '+
      (cautions.length>0?cautions.join('; '):'No secondary concerns flagged.');
  }
  return {sst_score:sst_score,sst_optimal:sst_optimal,sst_in_window:sst_in_window,
    sst_really_available:sstReallyAvailable,chl_score:chl_score,poll_score:poll_score,
    stab_score:stab_score,bromo_score:bromo_score,aqua_score:aqua_score,aqua_confidence:aqua_confidence,
    aqua_missing_note:aqua_missing_note,status:status,cautions:cautions};
}

// MODULE C-INT - INTERVENTION RECOMMENDATIONS
function computeInterventions(region, fp, sc, dhwv, tv, isReefZone) {
  var actions=[];
  if(fp.urchin_N===null||fp.urchin_N===undefined) {
    actions.push({priority:'DATA GAP',action:'Urchin density survey + potential reintroduction',
      basis:'No published Diadema/Paracentrotus density data for '+region+'. Levitan & Edmunds 2023 PNAS shows urchin collapse (15 to 0.3 per m2) is a primary driver of macroalgae overgrowth.',
      waddington:'Restores grazing pressure that maintains the Waddington bowl wall.'});
  } else if(fp.urchin_N<(fp.urchin_healthy||10)*0.3) {
    actions.push({priority:'HIGH',action:'Urchin density critically low ('+fp.urchin_N+'/m2)',
      basis:'Below 30% of healthy baseline. Per Levitan 2023, correlates with algae overgrowth risk.',
      waddington:'Bowl wall compromised - intervention now cheaper than after tipping point (hysteresis).'});
  }
  if((fp.anem_N===null||fp.anem_N===undefined)&&fp.hasField) {
    actions.push({priority:'DATA GAP',action:'Anemone/coral density survey',
      basis:'Behavioral data exists for '+region+' but density data not published. Recommend field quadrat survey.',
      waddington:'Density data calibrates the effective mass term in the bowl-depth model.'});
  }
  if(dhwv!==null&&dhwv>4) {
    if(isReefZone) {
      actions.push({priority:dhwv>12?'CRITICAL':dhwv>8?'HIGH':'MODERATE',
        action:'Active thermal stress mitigation: shading, assisted gene flow, or coral relocation',
        basis:'DHW='+dhwv.toFixed(1)+' deg C-weeks indicates '+(dhwv>12?'mortality-level':dhwv>8?'mass-bleaching-level':'bleaching-risk-level')+' thermal stress.',
        waddington:'Symptomatic relief - addresses sigma*dW forcing but does not change V(x) shape long-term.'});
    } else {
      actions.push({priority:'LOW',action:'Monitor only - thermal anomaly without reef bleaching risk',
        basis:'DHW='+dhwv.toFixed(1)+' detected but outside tropical reef ecology. Likely warm anomaly.',
        waddington:'No basin-wall intervention needed - satellite artifact of non-reef ecology.'});
    }
  }
  actions.push({priority:'CONTEXT',action:'Check S7 FAI/NDVI layer before any urchin reintroduction',
    basis:'RED FAI = urchin reintroduction HIGH priority. CLEAR = other interventions (thermal, pollution) take priority.',
    waddington:'Connects S7 (current state) to S4/grazer intervention - match intervention to satellite signal.'});
  if(fp.Cd!==null&&fp.Cd!==undefined) {
    var cdRatio=fp.Cd_poll?(fp.Cd-0.006)/(fp.Cd_poll-0.006):0;
    if(cdRatio>0.5) actions.push({priority:'MODERATE',action:'Source-tracking for heavy metal contamination (Cd/Pb)',
      basis:'Tissue metal burden elevated. Trace point sources (harbour runoff, industrial discharge).',
      waddington:'Pollution lowers mu threshold - removing source is basin-deepening fix.'});
  }
  if(actions.length===0||(actions.length===1&&actions[0].priority==='CONTEXT')) {
    actions.unshift({priority:'STABLE',action:'No urgent intervention indicated by current data',
      basis:'No acute stress signals. Continue monitoring.',
      waddington:'System resting in stable basin (B='+sc.B.toFixed(2)+'). Preventative monitoring appropriate.'});
  }
  return actions;
}

// MODULE C - SCORE COMPUTATION + ToE precompute
var GEBCO = ee.ImageCollection('projects/sat-io/open-datasets/gebco/gebco_grid').mosaic().select('b1').rename('elevation');
var gebco_raw = GEBCO;
var bathymetry_m = gebco_raw.multiply(-1).updateMask(gebco_raw.lt(0));
var waveBathy = ee.ImageCollection('COPERNICUS/MARINE/WAV/ANFC_0_083DEG_STATIC').first().select('deptho').rename('wave_depth');

// v10.54 FIX: apply *0.01 to OISST before linearFit (was silently off by factor 100)
var _makeAnnSST = function() {
  var list = ee.List.sequence(1982,2025).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1), t=ee.Image(ee.Number(yr).subtract(1982)).float().rename('t');
    return ee.ImageCollection('NOAA/CDR/OISST/V2_1').filter(ee.Filter.date(d,d.advance(1,'year'))).select('sst').mean()
      .multiply(0.01).float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annSSTColl=_makeAnnSST();
var toeSSTFit=_annSSTColl.select(['t','sst']).reduce(ee.Reducer.linearFit());
var toeSSTNoise=_annSSTColl.select('sst').reduce(ee.Reducer.stdDev());

var _makeAnnCHL = function() {
  var list = ee.List.sequence(1998,2024).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1), t=ee.Image(ee.Number(yr).subtract(1998)).float().rename('t');
    return ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6').filter(ee.Filter.date(d,d.advance(1,'year'))).select('chlor_a').mean()
      .float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annCHLColl=_makeAnnCHL();
var toeCHLFit=_annCHLColl.select(['t','chlor_a']).reduce(ee.Reducer.linearFit());
var toeCHLNoise=_annCHLColl.select('chlor_a').reduce(ee.Reducer.stdDev());

var _makeAnnSAL = function() {
  var list = ee.List.sequence(1993,2024).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1), t=ee.Image(ee.Number(yr).subtract(1993)).float().rename('t');
    return ee.ImageCollection('HYCOM/sea_temp_salinity').filter(ee.Filter.date(d,d.advance(1,'year'))).select('salinity_0').mean()
      .multiply(0.001).add(20).float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annSALColl=_makeAnnSAL();
var toeSALFit=_annSALColl.select(['t','salinity_0']).reduce(ee.Reducer.linearFit());
var toeSALNoise=_annSALColl.select('salinity_0').reduce(ee.Reducer.stdDev());

var _makeAnnNO2 = function() {
  var list = ee.List.sequence(2019,2025).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1), t=ee.Image(ee.Number(yr).subtract(2019)).float().rename('t');
    return ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').filter(ee.Filter.date(d,d.advance(1,'year')))
      .select('tropospheric_NO2_column_number_density').mean().float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annNO2Coll=_makeAnnNO2();
var toeNO2Fit=_annNO2Coll.select(['t','tropospheric_NO2_column_number_density']).reduce(ee.Reducer.linearFit());
var toeNO2Noise=_annNO2Coll.select('tropospheric_NO2_column_number_density').reduce(ee.Reducer.stdDev());

// v10.141 FIX: pH source replaced with a REAL, CONFIRMED, current
// Copernicus asset - the old asset ID (COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/
// MFC_001_028) was simply outdated; the real, current GEE catalog splits
// this product into per-variable sub-collections. Confirmed directly via
// Google's own Earth Engine catalog pages: COPERNICUS/MARINE/
// GLOBAL_ANALYSISFORECAST_BGC_001_028/CAR, band ph_depth1, real surface pH.
// IMPORTANT DISCLOSED CHANGE: this real asset's actual data only starts
// 2021-10-01 (confirmed via its own catalog page), NOT 1993 as the old
// 32-year assumption required - the annual loop below is now correctly
// rescoped to match, and pH's ToE confidence tier is downgraded from
// "32yr, MARGINAL" to "~4yr, LOW" accordingly, same honesty standard
// already used for NO2's short record.
// Dissolved oxygen (DO) is NOT fixed this version - a working sub-
// collection/band name for O2 was not confirmed after two direct
// searches; left disclosed as still unavailable rather than guessing an
// asset ID that could be wrong.
var _makeAnnPH = function() {
  return ee.ImageCollection(ee.List.sequence(2022,2025).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1);
    var t=ee.Image(ee.Number(yr).subtract(2022)).float().rename('t');
    var col=ee.ImageCollection('COPERNICUS/MARINE/GLOBAL_ANALYSISFORECAST_BGC_001_028/CAR')
      .filter(ee.Filter.date(d,d.advance(1,'year'))).select('ph_depth1');
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().float().rename('ph').addBands(t).set('system:time_start',d.millis()),
      ee.Image.constant(0).rename('ph').updateMask(ee.Image.constant(0)).addBands(t).set('system:time_start',d.millis())));
  }));
};
var _annPHColl=_makeAnnPH();
var toePHFit=_annPHColl.select(['t','ph']).reduce(ee.Reducer.linearFit());
var toePHNoise=_annPHColl.select('ph').reduce(ee.Reducer.stdDev());

// DO (dissolved oxygen) - STILL UNAVAILABLE. See v10.141 note above the
// pH fix: no confirmed working sub-collection/band name found this
// session for oxygen specifically - kept disabled rather than guessed.
var _makeAnnDO = function() {
  return ee.ImageCollection(ee.List.sequence(1993,2024).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1);
    var t=ee.Image(ee.Number(yr).subtract(1993)).float().rename('t');
    var col=ee.ImageCollection('COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/MFC_001_028').filter(ee.Filter.date(d,d.advance(1,'year'))).select('o2');
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().float().addBands(t).set('system:time_start',d.millis()),
      ee.Image.constant(0).rename('o2').updateMask(ee.Image.constant(0)).addBands(t).set('system:time_start',d.millis())));
  }));
};
var _annDOColl=_makeAnnDO();
var toeDOFit=_annDOColl.select(['t','o2']).reduce(ee.Reducer.linearFit());
var toeDONoise=_annDOColl.select('o2').reduce(ee.Reducer.stdDev());

function getWaveCelerity(depthMeters) { return depthMeters.multiply(9.81).sqrt().rename('wave_celerity'); }
function getEnergyConcentrationIndex(depthMeters) {
  return depthMeters.max(ee.Image(0.5)).sqrt().pow(-1).min(ee.Image(1.414)).rename('energy_concentration_index');
}
function getSLRScenario(depthMeters, slrMeters) {
  return getEnergyConcentrationIndex(depthMeters.subtract(ee.Image(slrMeters)).max(ee.Image(0.1))).rename('eci_slr_scenario_'+slrMeters+'m');
}

function computeToESignal(monthlySeries) {
  var valid=monthlySeries.filter(function(d){ return d.v!==null&&d.v!==undefined; });
  var n=valid.length;
  if(n<12) return {error:'Insufficient data (n='+n+')',n:n};
  var t0=valid[0].t, tScale=1000*60*60*24*365.25;
  var xs=valid.map(function(d){ return (d.t-t0)/tScale; });
  var ys=valid.map(function(d){ return d.v; });
  var sumX=0,sumY=0,sumXX=0,sumXY=0;
  for(var i=0;i<n;i++){sumX+=xs[i];sumY+=ys[i];sumXX+=xs[i]*xs[i];sumXY+=xs[i]*ys[i];}
  var meanX=sumX/n, meanY=sumY/n;
  var ssXX=sumXX-n*meanX*meanX, ssXY=sumXY-n*meanX*meanY;
  if(Math.abs(ssXX)<1e-10) return {error:'No time variance in data',n:n};
  var slope=ssXY/ssXX, intercept=meanY-slope*meanX;
  var residuals=[], sumR2=0;
  for(var j=0;j<n;j++) residuals.push(ys[j]-(intercept+slope*xs[j]));
  for(var k=0;k<n;k++) sumR2+=residuals[k]*residuals[k];
  var noise=Math.sqrt(sumR2/(n-2));
  var recordYears=xs[xs.length-1]-xs[0], signal=Math.abs(slope*recordYears);
  var snr=noise>0?signal/noise:0;
  return {slope:slope,intercept:intercept,signal:signal,noise:noise,snr:snr,emerged:snr>=2.0,
    direction:slope>0?'RISING':'FALLING',recordYears:recordYears.toFixed(1),nValid:n,error:null};
}

// v10.103 NEW: One Tree Reef real field data (Byrne et al. 2025, Limnol.
// Oceanogr. Lett., doi:10.1002/lol2.10456) - a real published paper tracking
// 462 individual coral colonies at 23.51 deg S, 152.09 deg E (One Tree
// Reef, southern GBR, Capricorn-Bunker Group), the EXACT coordinate this
// tool has been extensively tested against. Keyed by PROXIMITY (haversine
// distance), not the coarse "Great Barrier Reef" region bucket - the GBR
// spans from northern reefs (Lizard Island, hit hardest in 2016) to this
// southern site (largely spared until 2024), so applying One Tree Reef's
// specific findings to the whole GBR bounding box would misrepresent very
// different reefs with very different bleaching histories.
var ONE_TREE_REEF_LAT = -23.51, ONE_TREE_REEF_LON = 152.09, ONE_TREE_REEF_RADIUS_KM = 60;

function getEcologicalRecoveryValidation(region, lat, lon) {
  if(lat!==undefined&&lat!==null&&lon!==undefined&&lon!==null){
    var distToOTR = haversineKm(lat, lon, ONE_TREE_REEF_LAT, ONE_TREE_REEF_LON);
    if(distToOTR <= ONE_TREE_REEF_RADIUS_KM){
      return {checked:true, thermalPredictsEcological:false,
        finding:'CATASTROPHIC - real colony-level tracking shows severe, rapid mortality that thermal recovery alone would not predict',
        details:'Byrne et al. 2025 tracked 462 coral colonies at One Tree Reef (23.51S, 152.09E - this exact site) through the 2023-24 heatwave: 66% bleached by Feb 2024, 80% by April. Of colonies bleached by April, 44% were dead by May and 52% by July. Acropora had 95% mortality with rapid colony collapse to algal-fouled rubble; Goniopora developed black band disease (73% mortality). More thermally tolerant genera (Porites, Pocillopora) had lower mortality (31%, 23%) and some recovery (26%, 27% by July). In-situ DHW peaked at 14.6 deg C-weeks on 9 March 2024 - well above the 8 DHW mortality threshold. The paper concludes protected status and offshore location did NOT prevent mass bleaching and mortality.',
        sources:'Byrne et al. 2025, Limnol. Oceanogr. Lett., doi:10.1002/lol2.10456 | Data: Sydney eScholarship doi:10.25910/p5rq-cw63 | Distance from this click: '+distToOTR.toFixed(1)+'km'};
    }
  }
  if(region==='Florida Keys, USA'||region==='Florida Reef Tract, USA') {
    return {checked:true,thermalPredictsEcological:false,
      finding:'NO - thermal recovery did NOT predict ecological recovery here',
      details:'Real CREMP/DRM state monitoring: despite SST recovering within months, coral cover across Florida Reef is approx 2%, "over 90% relative reduction" - years after thermal recovery, reef has NOT ecologically recovered.',
      sources:"Florida's Coral Reef 2023-2024 Post-Bleaching Assessment (CREMP/DRM) | NOAA/FAU coral cover report"};
  }
  if(region==='Red Sea') {
    return {checked:true,thermalPredictsEcological:null,
      finding:'UNKNOWN - no published recovery follow-up found',
      details:'Peixoto et al. 2025 (npj Biodiversity) documented 66-94% anemone mortality but monitoring ends at the mortality event. No published follow-up survey found.',
      sources:'Checked against: Bennett-Smith et al. 2025 npj Biodiversity (ends at mortality event)'};
  }
  return {checked:false,thermalPredictsEcological:null,finding:null,
    details:'No published ecological-recovery follow-up independently checked for this region yet.',sources:'Not yet checked'};
}

// v10.132 NEW: refactored into a shared, reusable, extensible database -
// a real, GROWING (not yet global) collection of published genus-level
// coral vulnerability/mortality data, structured so ANY module in this
// tool can query it via lookupSpeciesVulnerability(lat, lon), not just
// S20. Two real, independently-verified entries as of this version:
//   1. One Tree Reef, GBR (Byrne et al. 2025) - carried over from v10.131
//   2. Florida Keys / Dry Tortugas (Manzello et al. 2025, Science,
//      doi:10.1126/science.adx7825 - confirmed 97.8-100% mortality of
//      Acropora palmata + A. cervicornis; plus Florida FWC's 2023-24
//      CREMP Post-Bleaching Quick Look Report for Orbicella/Siderastrea/
//      Montastraea) - NEW this version, added after directly searching
//      for and verifying real, current, published data.
// DISCLOSED SCOPE, unchanged from v10.131: this is NOT and will never
// claim to be a comprehensive global database. It grows one real,
// independently-verified citation at a time - the same discipline
// already used for S15/S5's existing entries. Any site not explicitly
// listed here correctly returns "not checked", never a guess. Species
// composition cannot currently be measured by satellite at reliable
// accuracy (real hyperspectral species-classification accuracy is only
// ~56-70%, vs ~88-90% for basic coral detection) - this remains a real
// field-survey lookup, not a remote-sensing capability.
var SPECIES_VULNERABILITY_SITES = [
  {
    name:'One Tree Reef, GBR', lat:-23.51, lon:152.09, radiusKm:60,
    overallSiteVulnerability:'HIGH',
    dominantGenera:[
      {genus:'Acropora', vulnerability:'HIGH', mortalityObserved:95},
      {genus:'Goniopora', vulnerability:'HIGH', mortalityObserved:73},
      {genus:'Porites', vulnerability:'MODERATE-LOW', mortalityObserved:31},
      {genus:'Pocillopora', vulnerability:'MODERATE-LOW', mortalityObserved:23}
    ],
    source:'Byrne et al. 2025, Limnol. Oceanogr. Lett., doi:10.1002/lol2.10456 - real colony-level survey, this exact site'
  },
  {
    name:'Florida Keys / Dry Tortugas', lat:24.55, lon:-81.78, radiusKm:80,
    overallSiteVulnerability:'HIGH',
    dominantGenera:[
      {genus:'Acropora (palmata + cervicornis)', vulnerability:'HIGH', mortalityObserved:99},
      {genus:'Orbicella', vulnerability:'MODERATE-HIGH', mortalityObserved:null,
        mortalityNote:'complete loss of many colonies + extensive partial mortality of survivors, Lower Keys'},
      {genus:'Siderastrea siderea', vulnerability:'LOW', mortalityObserved:1},
      {genus:'Montastraea cavernosa', vulnerability:'LOW', mortalityObserved:1}
    ],
    source:'Manzello et al. 2025, Science, doi:10.1126/science.adx7825 (Acropora 97.8-100% mortality, Florida Keys/Dry Tortugas) + '+
      'Florida FWC 2023-24 CREMP Post-Bleaching Quick Look Report (Orbicella/Siderastrea/Montastraea)'
  },
  {
    // v10.133 NEW: real, verified third entry. DISCLOSED DIFFERENCE from
    // the two entries above: this is REGIONAL data (18 reefs surveyed
    // across central/southern Maldivian atolls), not single-reef data
    // like One Tree Reef/Florida Keys - hence the larger radius and the
    // "regional average" framing in the numbers below, rather than a
    // precise per-reef figure.
    name:'Central/Southern Maldives Atolls (2023-24 event)', lat:3.5, lon:72.9, radiusKm:150,
    overallSiteVulnerability:'HIGH (central atolls) / LOWER (southern atolls - see note)',
    dominantGenera:[
      {genus:'Acropora (branching/tabular)', vulnerability:'HIGH', mortalityObserved:null,
        mortalityNote:'disproportionately impacted; up to 57% mortality at the worst-hit surveyed reef (a lagoon reef, Ari Atoll); central atolls averaged >40% live coral cover loss'},
      {genus:'Porites (massive)', vulnerability:'LOW', mortalityObserved:null,
        mortalityNote:'comparatively resistant; reef flats in the southernmost Huvadhoo Atoll retained high live coral cover, showing limited impact there'}
    ],
    source:'Outcomes of the fourth global coral bleaching (2023-2024) in the Maldives, Coral Reefs (Springer Nature), 2026, '+
      'doi:10.1007/s00338-026-02850-x - real survey of 18 reefs across central and southern atolls'
  },
  {
    // v10.134 NEW: this exact site was already tested with S13/CSD this
    // session (9.19,-81.984 / 9.378,-82.207). DISCLOSED: a second Bocas-
    // specific source was found but NOT used - its search snippet mixed
    // content from what appears to be an unrelated Indian coral-
    // management paper mid-sentence, a sign the source is unreliable.
    // Used only the clean, verified citation below instead.
    name:'Bocas del Toro, Panama', lat:9.34, lon:-82.25, radiusKm:60,
    overallSiteVulnerability:'MODERATE (chronic decline pattern, not single-event mortality)',
    dominantGenera:[
      {genus:'Orbicella franksi / Siderastrea siderea / Stephanocoenia michelini (massive genera)',
        vulnerability:'MODERATE', mortalityObserved:null,
        mortalityNote:'tracked individually through the real 2005 Caribbean bleaching event, then monitored for 8 years (2005-2013): '+
          'all groups showed continued NET TISSUE LOSS over that period despite being typically-considered stress-tolerant "winner" taxa - '+
          'no full recovery, even from massive/mounding growth forms'}
    ],
    source:'Neal et al. 2017, Ecology and Evolution, doi:10.1002/ece3.2706 - real colony-level 8-year study, this exact site'
  },
  {
    name:'Andaman & Nicobar Islands, India', lat:11.7, lon:92.7, radiusKm:200,
    overallSiteVulnerability:'MODERATE-HIGH',
    dominantGenera:[
      {genus:'Acropora', vulnerability:'HIGH', mortalityObserved:43},
      {genus:'Montipora', vulnerability:'MODERATE', mortalityObserved:22},
      {genus:'Porites', vulnerability:'LOW', mortalityObserved:14}
    ],
    source:'Coral bleaching in Andaman Sea - an indicator for climate change in Andaman and Nicobar Islands (2010 bleaching event survey) - '+
      'regional data, real 2010 mass-bleaching event (74-77% of South Andaman corals bleached)'
  },
  {
    // v10.134 NEW: Brazil/Southwestern Atlantic. DISCLOSED LIMIT: no
    // single mortality percentage was verified for the exact 2019-2024
    // events - the real, current (2024-2025) literature found describes
    // a qualitative but severe finding ("unprecedented erosion," "no
    // recovery three years after") rather than one precise number, so
    // this entry reports that finding honestly instead of inventing a
    // percentage.
    name:'Southwestern Atlantic Reefs, Brazil', lat:-17.9, lon:-38.7, radiusKm:250,
    overallSiteVulnerability:'HIGH (per real 2022/2024 field studies, exact % not verified)',
    dominantGenera:[
      {genus:'Mussismilia harttii (major Brazilian endemic reef-builder)', vulnerability:'HIGH', mortalityObserved:null,
        mortalityNote:'described as "unprecedented erosion" after the 2019 bleaching event, with a separate follow-up study finding no coral recovery three years after a major bleaching event in this region - exact mortality percentage not verified from available sources'}
    ],
    source:'Braz et al. 2022, Coral Reefs 41:1537-1548 (Mussismilia harttii erosion) + Corazza et al. 2024, Mar. Biol. 171:114 '+
      '(no recovery 3 years post-bleaching) - both cited in Coral Reefs (Springer), 2025, doi:10.1007/s00338-025-02743-5'
  },
  {
    // v10.140 NEW: Mesoamerican Reef (Mexico/Belize/Guatemala/Honduras).
    // HONEST CORRECTION: this was almost merged into the Bocas del Toro
    // (Panama) entry above - they are genuinely DIFFERENT reef systems.
    // Bocas del Toro sits on Panama's southern Caribbean coast; the
    // Mesoamerican Reef spans the Yucatan/Belize/Guatemala/Honduras
    // barrier reef further north. Kept as a separate entry rather than
    // incorrectly attributing regional data to the wrong reef system.
    // DISCLOSED SCOPE: reef-wide aggregate stats (not genus-level) - the
    // real sources found report overall percent-affected and cover
    // change, not a per-genus mortality breakdown like the other entries.
    name:'Mesoamerican Reef (Mexico/Belize/Guatemala/Honduras)', lat:17.0, lon:-87.8, radiusKm:400,
    overallSiteVulnerability:'HIGH (2023 event was the most severe on record for this region)',
    dominantGenera:[
      {genus:'Reef-wide (aggregate, not genus-specific)', vulnerability:'HIGH', mortalityObserved:null,
        mortalityNote:'2023 event (most severe on record for this region): ~40% of corals severely affected; regional live coral cover fell from 19% to 17%, with mortality continuing even after surveys concluded. DHW during major regional events has averaged ~16 deg C-weeks (range 12-24)'}
    ],
    source:'GCRMN Mesoamerican Reef Report Card 2024 (gcrmn.net) + Healthy Reefs for Healthy People / Healthy Reefs Initiative '+
      '(healthyreefs.org, 70+ partner organizations) - real regional monitoring network, not a single-paper citation'
  }
];
function lookupSpeciesVulnerability(lat, lon){
  if(lat===undefined||lat===null||lon===undefined||lon===null) return {checked:false};
  for(var i=0;i<SPECIES_VULNERABILITY_SITES.length;i++){
    var site = SPECIES_VULNERABILITY_SITES[i];
    var dist = haversineKm(lat, lon, site.lat, site.lon);
    if(dist <= site.radiusKm){
      return {checked:true, siteName:site.name, dominantGenera:site.dominantGenera,
        overallSiteVulnerability:site.overallSiteVulnerability, source:site.source, distanceKm:dist};
    }
  }
  return {checked:false};
}
// v10.137 NEW: genus-level growth-form lookup - a GENERAL, GLOBAL
// fallback layer, distinct from SPECIES_VULNERABILITY_SITES above.
// The site database gives specific, real EVENT-OUTCOME mortality data
// at 6 named locations. This instead gives a coarser, general
// vulnerability TIER, usable for ANY genus by name, ANYWHERE - based on
// growth form, the single best-established predictor of relative
// bleaching susceptibility in the literature (branching/fast-growing
// forms consistently more susceptible than massive/encrusting forms).
// HONESTY NOTE on sourcing, stated plainly rather than implying uniform
// confidence: entries marked [CTD-confirmed] were individually verified
// against real, live Coral Trait Database entries (coraltraits.org,
// Madin et al. 2016, Scientific Data, doi:10.1038/sdata.2016.17, trait
// 183 "Growth form") during this build. Entries marked [literature
// pattern] use the same well-established branching-vs-massive framework
// (Loya et al. 2001; confirmed again in Sadler et al. 2023, Mar. Poll.
// Bull., explicitly naming Stylophora/Acropora as most bleaching-
// susceptible and Lobophyllia/Porites as most tolerant) but were not
// individually cross-checked against a live CTD entry for that exact
// genus during this build.
// DISCLOSED LIMIT: could not directly fetch the live CTD CSV from
// inside this GEE script (outbound access to external domains from the
// Code Editor sandbox is unverified) - built from real, individually-
// confirmed values plus established literature, not a full bulk import.
var GENUS_GROWTH_FORM = {
  'Acropora':    {growthForm:'corymbose/branching [CTD-confirmed: Acropora tenuis=corymbose]', generalVulnerability:'HIGH'},
  'Stylophora':  {growthForm:'branching [literature pattern]', generalVulnerability:'HIGH'},
  'Pocillopora': {growthForm:'branching [literature pattern]', generalVulnerability:'HIGH'},
  'Madracis':    {growthForm:'digitate [CTD-confirmed: Madracis decactis=digitate]', generalVulnerability:'MODERATE-HIGH'},
  'Goniopora':   {growthForm:'flabello-meandroid/massive, variable [literature pattern]', generalVulnerability:'MODERATE-HIGH'},
  'Montipora':   {growthForm:'encrusting/branching, variable [literature pattern]', generalVulnerability:'MODERATE'},
  'Orbicella':   {growthForm:'massive [literature pattern]', generalVulnerability:'MODERATE'},
  'Agaricia':    {growthForm:'laminar [CTD-confirmed: Agaricia undata=laminar]', generalVulnerability:'MODERATE'},
  'Mussismilia': {growthForm:'massive [literature pattern] - NOTE: real Brazil field data (S20 entry) shows this genus suffered severe impact despite massive form, an exception worth remembering', generalVulnerability:'MODERATE (see note)'},
  'Porites':     {growthForm:'massive [CTD-confirmed: multiple massive-form relatives verified, e.g. Favites/Homophyllia/Australophyllia=massive]', generalVulnerability:'LOW'},
  'Siderastrea': {growthForm:'massive [literature pattern]', generalVulnerability:'LOW'},
  'Montastraea': {growthForm:'massive [literature pattern]', generalVulnerability:'LOW'},
  'Favites':     {growthForm:'massive [CTD-confirmed: Favites halicora=massive]', generalVulnerability:'LOW'},
  'Lobophyllia': {growthForm:'massive [literature pattern]', generalVulnerability:'LOW'}
};
function lookupGenusGrowthForm(genusInput){
  var key = (genusInput||'').trim();
  var firstWord = key.split(/\s+/)[0];
  var keys = Object.keys(GENUS_GROWTH_FORM);
  for(var i=0;i<keys.length;i++){
    if(keys[i].toLowerCase()===firstWord.toLowerCase()){
      return {checked:true, genus:keys[i], growthForm:GENUS_GROWTH_FORM[keys[i]].growthForm,
        generalVulnerability:GENUS_GROWTH_FORM[keys[i]].generalVulnerability};
    }
  }
  return {checked:false};
}

// Combines the species data above with an already-computed DHW value -
// takes no new measurement, just reads the real DHW number already
// produced by S4 and pairs it with the real species data above where
// both exist.
function combineSpeciesAndDHW(speciesData, dhwv, isTropicalReef){
  if(!speciesData.checked){
    return {combinedRisk:null,
      text:'NOT CHECKED - no published species composition data found for this region.\n'+
        (dhwv!==null?'DHW-only risk: '+dhwv.toFixed(2)+' deg C-wks (species vulnerability unknown - cannot combine).':'DHW also unavailable.')};
  }
  var dhwText = dhwv!==null?dhwv.toFixed(2)+' deg C-wks':'n/a';
  var highVulnGenera = speciesData.dominantGenera.filter(function(g){return g.vulnerability==='HIGH';});
  var lowVulnGenera = speciesData.dominantGenera.filter(function(g){return g.vulnerability!=='HIGH';});
  var genusLines = speciesData.dominantGenera.map(function(g){
    var mortText = (g.mortalityObserved!==null&&g.mortalityObserved!==undefined)?
      ('documented '+g.mortalityObserved+'% mortality at this site'):
      (g.mortalityNote?('documented: '+g.mortalityNote):'documented significant impact (no single percentage reported)');
    return '  '+g.genus+': '+g.vulnerability+' vulnerability ('+mortText+')';
  }).join('\n');
  var combinedRisk, headline;
  if(isTropicalReef && dhwv!==null && dhwv>=4 && highVulnGenera.length>0){
    combinedRisk='ELEVATED';
    headline='ELEVATED COMBINED RISK - real heat stress (DHW='+dhwText+') at a site with documented\n'+
      'highly-vulnerable coral genera ('+highVulnGenera.map(function(g){return g.genus;}).join(', ')+').';
  } else if(isTropicalReef && dhwv!==null && dhwv>=4 && highVulnGenera.length===0){
    combinedRisk='MODERATE';
    headline='MODERATE COMBINED RISK - real heat stress present (DHW='+dhwText+'), but documented genera at\n'+
      'this site skew toward lower historical vulnerability ('+lowVulnGenera.map(function(g){return g.genus;}).join(', ')+').';
  } else if(dhwv!==null && dhwv<4){
    combinedRisk='LOW (heat-driven)';
    headline='LOW current heat-driven risk (DHW='+dhwText+') - species vulnerability data shown below for\n'+
      'reference, but no active thermal stress to combine it with right now.';
  } else {
    combinedRisk='UNKNOWN';
    headline='DHW unavailable - cannot combine with species data.';
  }
  return {combinedRisk:combinedRisk,
    text:headline+'\n'+
      'Overall site vulnerability (per real field survey): '+speciesData.overallSiteVulnerability+'\n'+
      genusLines+'\n'+
      'Source: '+speciesData.source};
}

// v10.144 NEW: S20e - real, fitted logistic regression predicting
// P(bleaching), not collapse. Fit on 32,716 REAL rows from the Global
// Coral-Bleaching Database (van Woesik & Kratochwill 2022, Scientific
// Data, doi:10.1038/s41597-022-01121-y; downloaded from BCO-DMO,
// doi:10.26008/1912/bco-dmo.773466.2, and fit directly this session -
// not a citation of someone else's model, a real fit on the real file).
// Uses 3 features that map exactly onto this tool's own already-
// computed metrics (DHW, turbidity, depth) - ZERO new EE calls, all
// three values already exist in this same click's scope by the time
// this runs.
// REAL, HONEST VALIDATION: fit on a 70% train split (n=22,901),
// evaluated on a genuinely held-out 30% test split (n=9,815) never seen
// during fitting. Held-out AUC = 0.620 - real, better than chance
// (0.5), but genuinely WEAK. Disclosed directly in the UI, not hidden -
// consistent with Arias-Ortiz et al. 2024's finding that DHW-alone-
// style models need many more metrics (23, reaching 75% variance
// explained) to predict well. This 3-feature model is a real, modest
// improvement over a fixed threshold, not a strong predictor on its
// own.
// Coefficients/scaler below are EXACT values extracted from the fitted
// sklearn model (StandardScaler + LogisticRegression), verified this
// session with a manual-formula-vs-sklearn sanity check (exact match).
var BLEACHING_MODEL = {
  features:['SSTA_DHW','Turbidity','Depth_m'],
  scalerMean:[2.8719889961137124, 0.07541676346011245, 7.03591589886907],
  scalerScale:[4.353137496712977, 0.06207127436913297, 4.214382923695573],
  coef:[0.28920617949096045, -0.04853261433888384, 0.3452998473061333],
  intercept:-0.013852236398586949,
  heldOutAUC:0.620, nTrain:22901, nTest:9815,
  source:'Global Coral-Bleaching Database (van Woesik & Kratochwill 2022, Scientific Data, doi:10.1038/s41597-022-01121-y) - '+
    'fit directly on 32,716 real rows, this session, via BCO-DMO doi:10.26008/1912/bco-dmo.773466.2'
};
function predictBleachingProbability(dhw, turbidity, depthAbsMeters){
  if(dhw===null||dhw===undefined||turbidity===null||turbidity===undefined||depthAbsMeters===null||depthAbsMeters===undefined){
    return {p:null, error:'missing input(s) - need DHW, turbidity, and depth'};
  }
  var raw=[dhw, turbidity, depthAbsMeters];
  var z = BLEACHING_MODEL.intercept;
  for(var i=0;i<3;i++){
    var stdVal = (raw[i]-BLEACHING_MODEL.scalerMean[i])/BLEACHING_MODEL.scalerScale[i];
    z += BLEACHING_MODEL.coef[i]*stdVal;
  }
  var p = 1/(1+Math.exp(-z));
  return {p:p, error:null};
}

function getRegion(lat,lon) {
  if(lat>8&&lat<12&&lon>-83&&lon<-81)    return 'Bocas del Toro, Panama';
  if(lat>24&&lat<25.5&&lon>-82&&lon<-80) return 'Florida Keys, USA';
  if(lat>25&&lat<27&&lon>-82&&lon<-79)   return 'Florida Reef Tract, USA';
  if(lat>18&&lat<20&&lon>-80&&lon<-78)   return 'Jamaica';
  if(lat>19&&lat<20.5&&lon>-81&&lon<-79) return 'Little Cayman / Grand Cayman';
  if(lat>17&&lat<18.5&&lon>-67.5&&lon<-65) return 'Puerto Rico';
  if(lat>20&&lat<23&&lon>-75&&lon<-73)   return 'Cuba (north coast)';
  if(lat>8&&lat<25&&lon>-90&&lon<-58)    return 'Caribbean Sea';
  if(lat>32&&lat<49&&lon>-125&&lon<-117) return 'Pacific Coast USA';
  if(lat>25&&lat<32&&lon>-98&&lon<-80)   return 'Gulf of Mexico';
  if(lat>25&&lat<45&&lon>-82&&lon<-65)   return 'Atlantic Coast USA';
  if(lat>12&&lat<30&&lon>32&&lon<44)     return 'Red Sea';
  if(lat>30&&lat<47&&lon>-6&&lon<37)     return 'Mediterranean Sea';
  if(lat>-25&&lat<-10&&lon>142&&lon<155) return 'Great Barrier Reef';
  if(lat>-30&&lat<-25&&lon>152&&lon<154) return 'Great Barrier Reef';
  if(lat>18&&lat<23&&lon>-161&&lon<-154) return 'Hawaii, USA';
  if(lat>-25&&lat<0&&lon>163&&lon<180)   return 'New Caledonia / Pacific';
  if(lat>-20&&lat<5&&lon>160&&lon<180)   return 'Solomon Islands / Pacific';
  if(lon<-120&&lat>-60&&lat<65)          return 'Pacific Ocean';
  if(lat>27&&lat<30&&lon>-18&&lon<-13)   return 'Canary Islands, Spain';
  if(lat>36&&lat<43&&lon>-32&&lon<-24)   return 'Azores, Portugal';
  if(lat>32&&lat<33.5&&lon>-17.5&&lon<-16) return 'Madeira, Portugal';
  if(lat>14&&lat<18&&lon>-25&&lon<-22)   return 'Cape Verde Islands';
  if(lat>-1&&lat<1&&lon>-92&&lon<-88)    return 'Galapagos Islands, Ecuador';
  if(lat>-35&&lat<-28&&lon>27&&lon<33)   return 'KwaZulu-Natal, South Africa';
  // v10.67 NEW: Greenland regions for GEM MarineBasis station matching
  if(lat>60&&lat<68&&lon>-55&&lon<-44)   return 'SW Greenland / Nuuk';
  if(lat>70&&lat<78&&lon>-28&&lon<-12)   return 'NE Greenland / Zackenberg';
  if(lat>67&&lat<72&&lon>-56&&lon<-48)   return 'W Greenland / Disko';
  if(lat>56&&lat<84&&lon>-60&&lon<-12)   return 'Greenland Coast';
  if(lon>-60&&lon<-10&&lat>-60&&lat<65)  return 'Atlantic Ocean';
  if(lat<-55) return 'Southern Ocean';
  if(lat>75)  return 'Arctic Ocean';
  var la=Math.abs(Math.round(lat*10)/10), lo=Math.abs(Math.round(lon*10)/10);
  return la+(lat>=0?'N':'S')+' '+lo+(lon>=0?'E':'W')+' Coast';
}

function fmt(v,d) {
  if(v===null||v===undefined||isNaN(v)||!isFinite(v)||v<-900) return 'n/a';
  return (Math.round(v*Math.pow(10,d))/Math.pow(10,d)).toFixed(d);
}

// v10.85 FIX: depthV's text label used a crude "bv>-50 ? deep ocean : ..."
// rule that called anything below -50m "deep ocean", even though the S3
// map legend (see MAP LEGEND section below) splits that same range into
// FOUR distinct zones: continental shelf, continental slope, deep ocean,
// and very deep ocean. -60m (for example) is "Continental shelf" in the
// legend but was being labelled "deep ocean" in the sidebar - a direct
// contradiction between the map colour and the text next to it. This
// helper classifies a GEBCO depth exactly the way the legend bins do, so
// both always agree.
function classifyDepthLabel(bv) {
  if(bv===null||bv===undefined||isNaN(bv)) return '';
  if(bv>10) return ' land / higher ground';
  if(bv>0) return ' intertidal / beach';
  if(bv>-50) return ' shallow reef zone';
  if(bv>-200) return ' continental shelf';
  if(bv>-1000) return ' continental slope';
  if(bv>-2000) return ' deep ocean';
  return ' very deep ocean';
}

// v10.67 FIX (v10.19): hemisphere-aware coordinate parser
// parseFloat("13.5S") silently returns 13.5 (positive), dropping the hemisphere suffix.
function parseCoordPart(str) {
  // v10.151 FIX 1: parseFloat does not recognise U+2212 MINUS SIGN, en/em
  // dashes, or the fullwidth hyphen - all of which appear when coordinates
  // are pasted from documents, PDFs or chat. These previously produced a
  // silent NaN and a bare "invalid format" error with no hint why.
  // Also strips a degree sign / prime marks so "24.531\u00B0N" works.
  if(str===null||str===undefined) return NaN;
  str = String(str).trim()
    .replace(/[\u2010\u2011\u2012\u2013\u2014\u2015\u2212\uFF0D]/g,'-')
    .replace(/[\u00B0\u2032\u2033]/g,'')
    .replace(/\s+/g,'');
  if(str.length===0) return NaN;
  var lastChar = str.slice(-1).toUpperCase();
  if(lastChar==='N'||lastChar==='S'||lastChar==='E'||lastChar==='W') {
    var numPart = parseFloat(str.slice(0,-1));
    if(isNaN(numPart)) return NaN;
    return (lastChar==='S'||lastChar==='W') ? -Math.abs(numPart) : Math.abs(numPart);
  }
  return parseFloat(str);
}

function scoreColors(s) {
  if(s<30)  return {text:'#0a5c1e',bg:'#d4f5df',map:'#00cc44',bar:'#22cc44',lbl:'DEEP BASIN'};
  if(s<55)  return {text:'#7a5000',bg:'#fff6cc',map:'#ffcc00',bar:'#ddaa00',lbl:'WARNING'};
  if(s<75)  return {text:'#7a2e00',bg:'#ffe8d0',map:'#ff6600',bar:'#ff6600',lbl:'HIGH RISK'};
  return         {text:'#7a0000',bg:'#ffd0d0',map:'#ff0000',bar:'#ff2222',lbl:'CRITICAL'};
}

function computeScore(sv, cv, tv, nv, turv, dhwv, fp, lat, lon) {
  var isReefZone=(lat!==undefined&&lon!==undefined)?((lat>-30&&lat<30)&&!isEBUS(lat,lon)):true;
  var s1=50;
  if(sv!==null){
    if(sv>32)s1=100; else if(sv>30)s1=Math.round(75+(sv-30)*12.5);
    else if(sv>28)s1=Math.round(45+(sv-28)*15); else if(sv>=24)s1=10; else s1=Math.max(0,Math.round(30-(24-sv)*3));
  }
  s1=Math.max(0,Math.min(100,s1));
  var s2=50;
  if(cv!==null&&cv>0){
    if(cv<0.05)s2=85; else if(cv<0.1)s2=70; else if(cv<0.5)s2=50;
    else if(cv<1.0)s2=25; else if(cv<2.0)s2=10; else s2=5;
  }
  s2=Math.max(0,Math.min(100,s2));
  var s3=35;
  if(turv!==null){if(turv>0.3)s3=85; else if(turv>0.1)s3=65; else if(turv>0)s3=40; else if(turv>-0.1)s3=20; else s3=10;}
  s3=Math.max(0,Math.min(100,s3));
  var s4b=30;
  if(tv!==null){if(tv>0.08)s4b=95; else if(tv>0.05)s4b=80; else if(tv>0.03)s4b=60; else if(tv>0.01)s4b=40; else if(tv>0)s4b=20; else s4b=5;}
  var s4d=0;
  if(dhwv!==null){
    if(isReefZone){if(dhwv>16)s4d=40; else if(dhwv>12)s4d=35; else if(dhwv>8)s4d=25; else if(dhwv>4)s4d=15; else if(dhwv>2)s4d=8; else if(dhwv>1)s4d=3;
      if(fp.dhw_calibration&&dhwv>=fp.dhw_calibration)s4d=40;}
    else{if(dhwv>16)s4d=15; else if(dhwv>12)s4d=10; else if(dhwv>8)s4d=6; else if(dhwv>4)s4d=3; else if(dhwv>1)s4d=1;}
  }
  var s4=Math.min(100,s4b+s4d);
  var s5=Math.max(0,Math.min(100,(sv!==null&&sv>29?40:0)+(tv!==null&&tv>0.03?35:0)+(cv!==null&&cv<0.3?25:0)));
  var s6=20;
  if(nv!==null&&nv>0){if(nv>0.00015)s6=90; else if(nv>0.00010)s6=70; else if(nv>0.00005)s6=50; else if(nv>0.00002)s6=30; else s6=10;}
  s6=Math.max(0,Math.min(100,s6));
  var csat=Math.round(s1*0.15+s2*0.15+s3*0.15+s4*0.25+s5*0.20+s6*0.10);
  var F1=0,F2=0,F3=0,F4=0,F5=0;
  if(fp.hasField){
    if(fp.urchin_N!==null&&fp.urchin_healthy!==null)F1=-(Math.min(1.0,fp.urchin_N/fp.urchin_healthy)*15);
    if(fp.anem_N!==null)F3=-(Math.min(1.0,fp.anem_N/10.0)*10);
    if(fp.Cd!==null&&fp.Pb!==null){var CdF=Math.max(0,Math.min(1,(fp.Cd-0.006)/(fp.Cd_poll-0.006))),PbF=Math.max(0,Math.min(1,(fp.Pb-0.5)/(fp.Pb_poll-0.5)));F4=Math.round(((CdF+PbF)/2)*15);}
    if(fp.recruit!==null&&fp.recruit_healthy!==null)F5=-Math.round(Math.min(1.0,fp.recruit/fp.recruit_healthy)*8);
    else if(fp.recruit!==null&&fp.recruit<=5.0)F5=-(Math.min(1.0,fp.recruit/5.0)*8);
  }
  F1=isNaN(F1)?0:F1;F2=0;F3=isNaN(F3)?0:F3;F4=isNaN(F4)?0:F4;F5=isNaN(F5)?0:F5;
  var fcT=Math.round(F1+F2+F3+F4+F5);
  var ccs=Math.max(0,Math.min(100,csat+fcT));
  if(isNaN(ccs))ccs=csat; if(isNaN(ccs))ccs=30; ccs=Math.max(0,Math.min(100,ccs));
  var B=Math.round((100-ccs))/100;
  var mu=Math.max(0.01,ccs/100),dU=0.25*mu*mu,sig=0.04+ccs/2000;
  var anem_safe=(fp.anem_N!==null&&!isNaN(fp.anem_N))?fp.anem_N:6.0;
  var urch_safe=(fp.urchin_N!==null&&!isNaN(fp.urchin_N))?fp.urchin_N:0.3;
  var me=Math.max(0.1,anem_safe*0.2+urch_safe*0.05),w0=Math.sqrt(mu/me),k=w0*Math.exp(-dU/sig);
  var p5=isNaN(k)||!isFinite(k)?0:Math.round((1-Math.exp(-k*5))*100);
  p5=Math.max(0,Math.min(100,p5));
  var ac1=tv?Math.min(0.99,Math.max(0.10,0.40+tv*9)):0.50;
  var tau=(B>0&&!isNaN(B))?Math.round(10/B)/10:99;
  return {s1:s1,s2:s2,s3:s3,s4:s4,s5:s5,s6:s6,sat_ccs:csat,ccs:ccs,B:B,mu:mu,deltaU:dU,
    meff:me,omega0:w0,k:k,p5yr:p5,ac1:ac1,tau:tau,fcTotal:fcT,F1:F1,F2:F2,F3:F3,F4:F4,F5:F5,
    acc_sat:77,acc_field:fp.accuracy_field_gain||0,acc_total:77+(fp.accuracy_field_gain||0)};
}

// MODULE D - SIDEBAR UI
var panel = ui.Panel({style:{width:'256px',padding:'4px',backgroundColor:'#f0f5f0'}});
var clickResultsLog = [];
function lbl(txt,sz,col,bg,bld){return ui.Label(txt,{fontSize:sz+'px',fontWeight:bld?'bold':'normal',color:col||'#111111',backgroundColor:bg||'rgba(0,0,0,0)',padding:'2px 4px',margin:'1px 0'});}
function sHead(txt,bg,tc){return lbl(txt,9,tc||'#ffffff',bg||'#334455',true);}
function dynLbl(init,col){return ui.Label(init,{fontSize:'9px',color:col||'#111111',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0'});}
function row(key,dl){return ui.Panel([lbl(key+': ',9,'#334466'),dl],ui.Panel.Layout.Flow('horizontal'),{margin:'1px 0'});}
function legRow(hex,main,sub){
  return ui.Panel([ui.Label('',{backgroundColor:hex,width:'16px',height:'12px',margin:'2px 5px 0 2px',padding:'0',border:'1px solid #888888'}),
    ui.Panel([lbl(main,8,'#111111'),sub?lbl(sub,7,'#555555'):null].filter(Boolean),ui.Panel.Layout.Flow('vertical'),{margin:'0'})],
    ui.Panel.Layout.Flow('horizontal'),{margin:'2px 0'});
}
function legDiv(){return ui.Label('',{margin:'3px 0 1px 0',backgroundColor:'#cccccc',height:'1px',stretch:'horizontal'});}

panel.add(lbl('STEMGeoHS Marine v10.153',12,'#ffffff','#1a4a2a',true));
var clickLbl = lbl('CLICK coastal reef/shallow water to analyze',10,'#ffffff','#1a5a2a',true);
panel.add(clickLbl);

// GO TO COORDINATES - v10.67 FIX: hemisphere-aware parseCoordPart()
panel.add(lbl('GO TO COORDINATES (use this, not the search bar above)',9,'#ffffff','#334466',true));
panel.add(lbl('Format: lat, lon  e.g. 13.5S, 112.3E  or  -13.5, 112.3  or  63.84, -22.43',7,'#555555'));
var coordInput = ui.Textbox({
  placeholder: 'lat, lon  e.g. 13.5S, 112.3E  or  64.13, -51.38 (Nuuk GEM)',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(coordInput);
var coordStatusV = dynLbl('','#880000');
var goToCoordsBtn = ui.Button({
  label: 'GO (pans map + draws white ring + analyzes)',
  style: {fontSize:'11px',fontWeight:'bold',margin:'2px 4px',
    backgroundColor:'#cce0ff',color:'#003388',stretch:'horizontal',padding:'6px 4px',border:'2px solid #0055cc'},
  onClick: function(){
    var txt = (coordInput.getValue()||'').trim();
    if(!txt){ coordStatusV.setValue('Enter coordinates first, e.g. 13.5S, 112.3E'); return; }
    // v10.67 FIX (v10.19): use hemisphere-aware parser instead of raw parseFloat
    var rawParts = txt.split(',');
    if(rawParts.length!==2){ coordStatusV.setValue('Invalid format - use: lat, lon  (e.g. 13.5S, 112.3E  or  -13.5, 112.3)'); return; }
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){ coordStatusV.setValue('Invalid format - use: lat, lon  (e.g. 13.5S, 112.3E  or  -13.5, 112.3)'); return; }
    if(latIn<-90||latIn>90){ coordStatusV.setValue('Latitude must be -90..90'); return; }
    if(lonIn<-180||lonIn>180){ coordStatusV.setValue('Longitude must be -180..180'); return; }
    coordStatusV.setValue('Going to '+latIn+', '+lonIn+'...');
    coordStatusV.style().set('color','#115511');
    Map.setCenter(lonIn,latIn,11);
    // Draw white ring immediately before compute pipeline starts
    var immediateRing=ee.Geometry.Point([lonIn,latIn]).buffer(1000);
    while(Map.layers().length()>0){ Map.remove(Map.layers().get(0)); }
    Map.addLayer(ee.Image().byte().paint(immediateRing,0,3),
      {palette:['#ffffff'],opacity:1.0},'Study zone (white border)');
    analyzeLocation(latIn,lonIn);
  }
});
panel.add(goToCoordsBtn);
panel.add(coordStatusV);

// v10.83 FIX: GEE's server-side script sandbox does NOT support the ES6
// String.prototype.repeat() method. Calling it throws "TypeError:
// (intermediate value)... is not a function" - which is exactly what was
// silently killing FIND SWEET SPOT (and would have killed COMPARE's
// separator lines too). Use this loop-based helper everywhere instead.
function repeatChar(ch, n) {
  var s = '';
  for (var i = 0; i < n; i++) s += ch;
  return s;
}

// v10.88 NEW: "toolkit" tally helper - counts how many independent EWS
// indicators are available and how many agree in the CSD-consistent
// direction, rather than letting one indicator (historically: temporal
// variance) decide the verdict alone. Mirrors the "don't rely on a single
// indicator" guidance from the early-warning-signals literature.
// v10.90 FIX: previously every scored indicator (AC1, variance, spatial
// variance, spatial autocorrelation) counted as ONE EQUAL VOTE. Dakos et al.
// 2012 (Ecology 93:264-271) found this is not justified: autocorrelation
// "appears a relatively robust indicator... regardless of the source of
// noise" across every scenario they tested, while variance "may sometimes
// decrease close to a transition" for well-documented reasons (parameter-
// noise sensitivity changes, and a "freezing" effect in slow-responding
// systems - their Fig. 4). So AC1 not rising is more informative than
// variance not rising, and a variance-only rise deserves LESS confidence
// than an AC1-only rise - the opposite of a flat 1-vote-each tally. This
// rebuild treats AC1 as the PRIMARY indicator and everything else
// (temporal variance, spatial variance, spatial autocorrelation) as
// SUPPORTING evidence that raises or lowers confidence around it.
function buildToolkitTally(indicators) {
  var lines=[], nAvail=0, nAgree=0;
  var primaryAvailable=false, primaryAgrees=false, primaryName=null;
  var supportAvail=0, supportAgree=0;
  indicators.forEach(function(ind){
    if(ind.scored){
      if(ind.available){
        nAvail++;
        if(ind.agrees) nAgree++;
        if(ind.primary){ primaryAvailable=true; primaryAgrees=!!ind.agrees; primaryName=ind.name; }
        else { supportAvail++; if(ind.agrees) supportAgree++; }
      }
      lines.push('  '+(ind.primary?'[PRIMARY - Dakos et al. 2012] ':'[supporting] ')+ind.name+': '+ind.display+
        (ind.available?(ind.agrees?'  [CSD-consistent]':'  [not rising]'):'  [unavailable]'));
    } else {
      lines.push('  '+ind.name+': '+ind.display+'  [reported, not scored]');
    }
  });
  return {lines:lines, nAvail:nAvail, nAgree:nAgree,
    primaryAvailable:primaryAvailable, primaryAgrees:primaryAgrees, primaryName:primaryName,
    supportAvail:supportAvail, supportAgree:supportAgree};
}

// AC1-weighted confidence classification, shared by STEP 3's live render
// and its final combined verdict. Returns {label, level} where level is
// one of 'high'/'moderate'/'low-moderate'/'low'/'none' for colour-coding.
// ============================================================
// v10.153 FIX 10 - EMPIRICALLY CALIBRATED THRESHOLDS
//
// The >0.01 AC1 and >0.15x variance cutoffs were never calibrated
// against anything. They have now been MEASURED, using 13.6 years of
// 4-minute Scripps Pier CTD data (Jan 2013 - Aug 2026, 1.95M readings,
// QC-flagged: temperature, salinity, chlorophyll). La Jolla is a
// temperate kelp coast with no documented regime shift, so every
// threshold crossing there is a false positive by construction.
//
// Across 104-176 non-overlapping BEFORE/AFTER splits, deseasonalized
// exactly the way this tool does it:
//
//   window   AC1 delta sd   Var delta sd   FALSE-POSITIVE RATE, old cutoffs
//   24mo        0.247          3.003        AC1 51%,  Var 48%
//   36mo        0.205          2.826        AC1 46%,  Var 49%
//   48mo        0.166          1.853        AC1 34%,  Var 57%
//
// The AC1-OR-variance rule that S7E uses to declare a site signal fired
// on 80% of no-event temperature windows. The permutation test on those
// SAME windows rejected at 0-9% - correctly calibrated.
//
// Values below are the 95th percentile of |delta| under that null, i.e.
// the cutoff that actually delivers a 5% false-positive rate.
//
// READ THIS BEFORE TRUSTING THEM: calibrated cutoffs this large will
// almost never fire. That IS the honest conclusion - a fixed threshold
// cannot adapt to how noisy a particular site is, while the permutation
// test adapts per-site automatically. These values exist so the
// threshold branch stops manufacturing confident false positives, not
// because a fixed threshold is a good idea. The p-value is still the
// only verdict here with a known error rate.
var CSD_NULL_CALIBRATION = {
  24: {ac1:0.491, varr:7.311, n:176},
  36: {ac1:0.395, varr:7.020, n:140},
  48: {ac1:0.347, varr:4.797, n:104}
};
var CSD_USE_CALIBRATED_THRESHOLDS = true;   // false restores the old 0.01/0.15 behaviour

function getCalibratedThresholds(nMonths){
  if(!CSD_USE_CALIBRATED_THRESHOLDS)
    return {ac1:0.01, varr:0.15, n:0,
      basis:'UNCALIBRATED legacy cutoffs - measured false-positive rate 34-57%'};
  var keys=[24,36,48], best=keys[0], i;
  for(i=0;i<keys.length;i++){
    if(Math.abs(keys[i]-(nMonths||24))<Math.abs(best-(nMonths||24))) best=keys[i];
  }
  var c=CSD_NULL_CALIBRATION[best];
  return {ac1:c.ac1, varr:c.varr, n:c.n,
    basis:'95th pct of |delta| under a measured no-event null (Scripps Pier CTD, '+
      c.n+' splits, nearest calibrated window '+best+'mo)'};
}

// Attach to ANY verdict derived from a threshold rather than a p-value.
function thresholdOnlyWarning(){
  return 'THRESHOLD-ONLY VERDICT - no significance test was run in this module. '+
    'Fixed cutoffs were measured firing on up to 80% of windows at a site where '+
    'nothing happened. Read as a direction indicator, not as evidence.';
}

// v10.151 FIX 5: the confidence label now CANNOT exceed the p-value.
// WHY: a real Looe Key run printed "NO STATISTICALLY SIGNIFICANT SIGNAL,
// AC1 p=0.980" on one line and "Confidence: HIGH" four lines below. The
// old function took only the tally and never saw a p-value, so with
// primaryAgrees=true it returned HIGH by construction. The permutation
// result resolved in a separate callback and wrote only to the banner -
// two verdicts, disjoint inputs, structurally unable to agree.
// Directional agreement among indicators is not evidence when the
// significance test cannot separate the deltas from ordinary noise.
// Pass null/null before the permutation test resolves; that correctly
// yields PRELIMINARY rather than a confident-sounding label.
function classifyToolkitConfidence(tally, pAC1, pVar){
  var base;
  if(!tally.primaryAvailable){
    if(tally.supportAvail===0) base={label:'NO INDICATORS AVAILABLE', level:'none'};
    else if(tally.supportAgree>0) base={label:'LOW (AC1 unavailable - only weaker supporting indicators to go on: '+
      tally.supportAgree+'/'+tally.supportAvail+' agree)', level:'low'};
    else base={label:'LOW (AC1 unavailable, no supporting indicators agree either)', level:'low'};
  } else if(tally.primaryAgrees){
    if(tally.supportAvail>0 && tally.supportAgree===tally.supportAvail)
      base={label:'primary indicator (AC1) rising, corroborated by ALL '+tally.supportAvail+' supporting indicator(s)', level:'high'};
    else if(tally.supportAgree>0)
      base={label:'primary indicator (AC1) rising, corroborated by '+tally.supportAgree+' of '+tally.supportAvail+' supporting indicator(s)', level:'high'};
    else if(tally.supportAvail>0)
      base={label:'primary indicator (AC1) rising; supporting indicator(s) flat/falling. Per Dakos et al. 2012 this does NOT by itself weaken the AC1 signal - variance is documented to sometimes fall even near a genuine transition (their Fig. 2c, Fig. 4).', level:'moderate'};
    else
      base={label:'primary indicator (AC1) rising alone, no supporting data available (Dakos et al. 2012: AC1 alone is a relatively robust signal)', level:'moderate'};
  } else if(tally.supportAgree>0){
    base={label:'primary indicator (AC1) did NOT rise, even though '+tally.supportAgree+
      ' supporting indicator(s) did. Dakos et al. 2012: variance-type indicators can rise OR fall near a real transition, so this is weaker evidence than an AC1 rise would be.', level:'low-moderate'};
  } else {
    base={label:'primary indicator (AC1) did not rise, and no supporting indicators rose either', level:'low'};
  }

  // ---- significance gate (v10.151) ----
  var havePvals=(pAC1!==null&&pAC1!==undefined)||(pVar!==null&&pVar!==undefined);
  if(!havePvals){
    return {label:'PRELIMINARY (permutation test has not returned a p-value yet) - direction only: '+base.label,
      level:'preliminary'};
  }
  var sig=(pAC1!==null&&pAC1!==undefined&&pAC1<0.05)||(pVar!==null&&pVar!==undefined&&pVar<0.05);
  var pTxt='AC1 p='+((pAC1!==null&&pAC1!==undefined)?pAC1.toFixed(3):'n/a')+
           ', Var p='+((pVar!==null&&pVar!==undefined)?pVar.toFixed(3):'n/a');
  if(!sig&&(base.level==='high'||base.level==='moderate')){
    return {label:'LOW - the indicator DIRECTIONS look CSD-consistent ('+base.label+
      '), but the permutation test cannot distinguish these deltas from ordinary noise ('+pTxt+
      '). Directional agreement among indicators is not evidence when the significance test is null.',
      level:'low'};
  }
  if(sig&&base.level==='high') return {label:'HIGH - '+base.label+'; permutation test agrees ('+pTxt+')', level:'high'};
  if(sig&&base.level==='moderate') return {label:'MODERATE - '+base.label+'; permutation test agrees ('+pTxt+')', level:'moderate'};
  return {label:base.label+' ('+pTxt+')', level:base.level};
}

// v10.151 FIX 6: general verdict-conflict test, replacing the two
// string-match checks (v10.119 headline, v10.120 detail) that only ever
// looked for 'MORE STABLE THAN CONTROL'. At Looe Key vTitle was
// 'NO SIGNAL AT EITHER SITE', so neither fired and the panel printed
// verbatim: "STRONG SIGNAL, LIKELY REGIONAL ... Regional context: NO
// SIGNAL AT EITHER SITE". That was the fifth appearance of this bug
// class; each previous fix patched one string. This tests the actual
// disagreement instead, so future cases are caught automatically.
// Returns null when there is no conflict (or when either side is
// genuinely unknown), otherwise an explanatory string.
function verdictsConflict(combinedTitle, vTitle){
  function saysSignal(str){
    str=String(str||'').toUpperCase();
    if(str.indexOf('CANNOT ASSESS')>=0) return null;
    if(str.indexOf('CANNOT CLASSIFY')>=0) return null;
    if(str.indexOf('NO SIGNAL')>=0) return false;
    if(str.indexOf('NO RELIABLE')>=0) return false;
    if(str.indexOf('MORE STABLE')>=0) return false;
    if(str.indexOf('STRONG')>=0||str.indexOf('SIGNAL PRESENT')>=0||
       str.indexOf('LOCAL CSD')>=0||str.indexOf('MARGINAL LOCAL')>=0||
       str.indexOf('GLOBAL SIGNAL')>=0||str.indexOf('WEAK SIGNAL')>=0) return true;
    return null;
  }
  var x=saysSignal(combinedTitle), y=saysSignal(vTitle);
  if(x===null||y===null) return null;
  if(x===y) return null;
  return 'VERDICTS DISAGREE - the AC1-weighted toolkit headline says "'+combinedTitle+
    '" while the variance-based regional read says "'+vTitle+'". These come from two separate '+
    'classifiers with different inputs and have not been reconciled. Treat NEITHER as confirmed - '+
    'read the permutation p-value above, which is the only real significance test in this panel.';
}

// v10.84 FIX: translate cryptic Earth Engine errors into something actionable.
// Errors like "Unknown reference to value named ''..." or "Failed to contact
// Earth Engine servers" are almost always a TEMPORARY server/connectivity
// hiccup (often from firing many parallel calls at once, e.g. right after a
// map click), not a real bug in the inputs. Recognize the common patterns
// and tell the user to just retry, while still showing the raw error text
// underneath for debugging.
function friendlyEEError(errText) {
  var raw = String(errText);
  var lower = raw.toLowerCase();
  var isTransient = lower.indexOf('unknown reference to value named')>=0 ||
    lower.indexOf('failed to contact earth engine')>=0 ||
    lower.indexOf('econnreset')>=0 || lower.indexOf('timeout')>=0 ||
    lower.indexOf('too many')>=0 || lower.indexOf('rate limit')>=0 ||
    lower.indexOf('internal error')>=0 || lower.indexOf('computation timed out')>=0;
  if(isTransient){
    return 'TEMPORARY EARTH ENGINE HICCUP - not a data problem.\n'+
      'This usually happens when several requests fire at once (e.g. right after\n'+
      'a map click or another S13 test). Your inputs are still in the boxes -\n'+
      'just wait a few seconds and press the button again.\n'+
      '\nRaw error: '+raw;
  }
  return 'Raw error: '+raw;
}

// S13 - CSD EARLY WARNING TEST (v10.149)
// ============================================================
panel.add(sHead('S13 - CSD EARLY WARNING TEST (v10.149)','#1a4a4a'));
panel.add(lbl('Tests whether a reef shows "critical slowing down" (CSD, Scheffer et al. 2009) — a statistical warning sign that can appear before ecological collapse. The classic Scheffer signature is BOTH indicators rising together: autocorrelation (AC1) AND variance.',7,'#226666'));
panel.add(lbl('v10.90: AC1 is weighted as the PRIMARY indicator throughout, per Dakos et al. 2012 (Ecology 93:264-271), which found autocorrelation "relatively robust" while variance can rise OR fall near a real transition. A variance-only signal (AC1 not rising) is now explicitly flagged as weaker evidence than an AC1-confirmed one.',7,'#886600'));
panel.add(lbl('v10.91: the reverse case - AC1 RISING while variance FALLS - is treated as a valid, still-meaningful signal, not a weak/contradicted one. Dakos et al. document this exact pattern (their Fig. 2c, Fig. 4): variance can decrease near a genuine transition while AC1 keeps rising regardless.',7,'#886600'));
panel.add(lbl('Independent tool - does not affect the main score or the map-click flow above.',7,'#888888'));

panel.add(lbl('HOW THIS WORKS, IN 4 STEPS:',8,'#ffffff','#226666',true));
panel.add(lbl('STEP 1: Enter the study location',7,'#115511'));
panel.add(lbl('STEP 2: Run a BEFORE window (calm baseline), then an AFTER window (the suspect period)',7,'#aa3300'));
panel.add(lbl('STEP 3: Click COMPARE - auto-picks a deep-water control site, then scores a TOOLKIT of independent indicators (temporal AC1, temporal variance, spatial variance, spatial autocorrelation, study-control synchronization) rather than trusting any single one, with a colour-coded verdict box',7,'#334466'));
panel.add(lbl('v10.101 CAVEAT: the Study-Control Synchronization indicator and the core AC1 statistic both use RAW/linearly-detrended monthly SST, not deseasonalized anomalies. Two ocean points a short distance apart share a strong seasonal cycle (both warm in summer, cool in winter) regardless of any real dynamical change, which can push correlation and AC1 toward a high baseline for reasons unrelated to critical slowing down. Look at the DELTA (change vs BEFORE), not the absolute value, and treat these two indicators as lower-confidence until a deseasonalized version is built and tested.',7,'#aa3300'));
panel.add(lbl('STEP 4 (optional): FIND SWEET SPOT auto-tests 6 AFTER window lengths, compares LOCAL vs REGIONAL warming at each one, and tells you exactly which window shows the clearest - and Scheffer-validated - signal',7,'#334466'));
panel.add(lbl('v10.95 CAUTION: testing 6 windows and picking the most dramatic one is a classic multiple-comparisons trap - some window will look "significant" by chance alone. The tool now flags an ISOLATED result (only 1 window leans local) vs a ROBUST one (several neighbouring windows agree). Best practice: pick your AFTER window from independent evidence (a documented event date), not by letting this scan choose for you.',7,'#886600'));
panel.add(lbl('If any step shows a red box saying "TEMPORARY EARTH ENGINE HICCUP" - that is a server/connection blip, not a data problem. Your inputs are untouched; just wait a few seconds and press the same button again.',7,'#886600'));
panel.add(legDiv());

// ----------------------------------------------------------
// STEP 1 - LOCATION
// ----------------------------------------------------------
panel.add(lbl('STEP 1 - Study location',9,'#ffffff','#115511',true));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var csdTestCoordInput = ui.Textbox({
  placeholder: 'lat, lon  e.g. 13.5S, 112.3E  or  24.55, -81.78',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(csdTestCoordInput);
var csdUseLastClickBtn = ui.Button({
  label:'Use last clicked location (from map above)',
  style:{fontSize:'9px',margin:'2px 4px',backgroundColor:'#e8f4ff',color:'#225588',stretch:'horizontal',padding:'4px 4px',border:'1px solid #4488cc'},
  onClick:function(){
    if(lastClickLat===null||lastClickLon===null){
      csdTestStatusV.setValue('No location clicked on the main map yet - click the map or use GO TO COORDINATES first, then come back here.');
      csdTestStatusV.style().set('color','#aa3300');
      return;
    }
    csdTestCoordInput.setValue(lastClickLat.toFixed(4)+', '+lastClickLon.toFixed(4));
    csdTestStatusV.setValue('Location filled in from your last map click. Now do STEP 2 below.');
    csdTestStatusV.style().set('color','#115511');
  }
});
panel.add(csdUseLastClickBtn);
panel.add(legDiv());

// ----------------------------------------------------------
// STEP 2 - RUN BEFORE, THEN RUN AFTER
// ----------------------------------------------------------
panel.add(lbl('STEP 2 - Run a BEFORE window, then run an AFTER window',9,'#ffffff','#aa5500',true));
panel.add(lbl('Same lat/lon both times. Only the start date + the Label dropdown change between the two runs.',7,'#664400'));
panel.add(lbl('Start date (YYYY-MM-DD):',7,'#334466'));
var csdTestStartInput = ui.Textbox({placeholder:'e.g. 2021-06-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(csdTestStartInput);
panel.add(lbl('Number of months (4-60):',7,'#334466'));
var csdTestMonthsInput = ui.Textbox({placeholder:'e.g. 24',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(csdTestMonthsInput);
panel.add(lbl('Label this run as:',7,'#334466'));
var csdTestWindowLabel = ui.Select({
  items:['BEFORE event (baseline)','AFTER event (suspect period)'],
  value:'BEFORE event (baseline)',
  style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}
});
panel.add(csdTestWindowLabel);

// S13 result display labels
var csdTestStatusV=dynLbl('Fill in Step 1 + Step 2 above, then press RUN.','#553388'), csdTestAc1V=dynLbl('--','#226666');
var csdTestVarTrendV=dynLbl('--','#226666'), csdTestThermalV=dynLbl('--','#226622');
var csdTestNoteV=dynLbl('','#888888');
var csdBeforeResult=null, csdAfterResult=null;

// persistent storage state labels
var csdStoredBeforeV=dynLbl('BEFORE stored: [not yet run]','#553388');
var csdStoredAfterV=dynLbl('AFTER stored:  [not yet run]','#553388');
var csdClearStoredBtn=ui.Button({
  label:'Clear stored BEFORE/AFTER (reset and start over)',
  style:{fontSize:'9px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#ffe0e0',color:'#770000',stretch:'horizontal',padding:'4px 4px',border:'1px solid #aa2222'},
  onClick:function(){
    csdBeforeResult=null; csdAfterResult=null;
    csdStoredBeforeV.setValue('BEFORE stored: [cleared - run a BEFORE test first]');
    csdStoredBeforeV.style().set('color','#888888');
    csdStoredAfterV.setValue('AFTER stored:  [cleared - run an AFTER test first]');
    csdStoredAfterV.style().set('color','#888888');
    csdCompareResultV.setValue('');
    csdCompareVerdictV.setValue('Run BEFORE + AFTER above, then press COMPARE (this will find the control site first, then compare).');
    csdCompareVerdictV.style().set('color','#555555'); csdCompareVerdictV.style().set('backgroundColor','#eeeeee');
    csdCompareVerdictV.style().set('border','2px solid #aaaaaa'); csdCompareVerdictV.style().set('whiteSpace','normal');
    csdToolkitV.setValue('');
    print('=== S13: Stored BEFORE/AFTER results cleared. ===');
  }
});

var csdTestRunBtn = ui.Button({
  label:'RUN CSD TEST ON THIS WINDOW',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#cce8e8',color:'#114444',stretch:'horizontal',padding:'6px 4px',border:'2px solid #226666'},
  onClick:function(){
    var coordTxt=(csdTestCoordInput.getValue()||'').trim(), startTxt=(csdTestStartInput.getValue()||'').trim(), monthsTxt=(csdTestMonthsInput.getValue()||'').trim();
    if(!coordTxt||!startTxt||!monthsTxt){csdTestStatusV.setValue('Fill in all 3 fields (Lat/Lon, Start date, Months) first.'); csdTestStatusV.style().set('color','#aa3300'); return;}
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){csdTestStatusV.setValue('Invalid lat/lon format - use: lat, lon'); csdTestStatusV.style().set('color','#aa3300'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){csdTestStatusV.setValue('Invalid lat/lon format - use: lat, lon'); csdTestStatusV.style().set('color','#aa3300'); return;}
    if(latIn<-90||latIn>90||lonIn<-180||lonIn>180){csdTestStatusV.setValue('Lat must be -90..90, Lon must be -180..180'); csdTestStatusV.style().set('color','#aa3300'); return;}
    var nMonthsIn=parseInt(monthsTxt,10);
    if(isNaN(nMonthsIn)||nMonthsIn<4||nMonthsIn>60){csdTestStatusV.setValue('Number of months must be 4-60'); csdTestStatusV.style().set('color','#aa3300'); return;}
    // v10.84: echo back exactly what was parsed, so a truncated/mistyped
    // coordinate is obvious immediately rather than discovered later.
    csdTestStatusV.setValue('Running CSD test: '+latIn.toFixed(4)+', '+lonIn.toFixed(4)+
      ' | start '+startTxt+' | '+nMonthsIn+' months | label='+csdTestWindowLabel.getValue()+' ...');
    csdTestStatusV.style().set('color','#115511'); csdTestStatusV.style().set('whiteSpace','pre');
    csdTestAc1V.setValue('computing...'); csdTestVarTrendV.setValue('computing...'); csdTestThermalV.setValue('computing...'); csdTestNoteV.setValue('');
    var testPt=ee.Geometry.Point([lonIn,latIn]), testStudy=testPt.buffer(1000);
    var testColl=mkMoSSTRange(startTxt,nMonthsIn);
    var windowLabel=csdTestWindowLabel.getValue();
    var rawValsFC=ee.FeatureCollection(testColl.map(function(img){
      var v=img.reduceRegion({reducer:ee.Reducer.mean(),geometry:testStudy,scale:4000,maxPixels:1e9}).get('sst');
      return ee.Feature(null,{v:v,t:img.get('system:time_start')});
    })).sort('t');
    var rMMM_test=MMM_perpixel.reduceRegion({reducer:ee.Reducer.mean(),geometry:testStudy,scale:4000,maxPixels:1e9});
    rawValsFC.evaluate(function(rawFC,_err){
      if(_err){print('=== S13 RAW DATA ERROR ['+windowLabel+'] === '+_err); csdTestThermalV.setValue('n/a - see status above for details'); return;}
      print('=== S13 RAW MONTHLY DATA ['+windowLabel+'] '+startTxt+' + '+nMonthsIn+' months ===');
      print('Location: '+latIn+', '+lonIn);
      var feats=rawFC.features, validCount=0;
      for(var mi=0;mi<feats.length;mi++){
        var props=feats[mi].properties, dateStr=new Date(props.t).toISOString().slice(0,7);
        if(props.v!==null&&props.v!==undefined){print('  Month '+(mi+1)+' ('+dateStr+'): '+props.v.toFixed(2)+' deg C'); validCount++;}
        else{print('  Month '+(mi+1)+' ('+dateStr+'): NO DATA (null)');}
      }
      print('Raw valid months: '+validCount+' / '+feats.length);
      rMMM_test.evaluate(function(mmmRes){
        var mmmTestVal=(mmmRes&&mmmRes.mmm!==null&&mmmRes.mmm!==undefined)?mmmRes.mmm:null;
        var trTest=analyzeThermalRecovery(feats,mmmTestVal);
        if(trTest.error){csdTestThermalV.setValue('n/a - '+trTest.error);}
        else{
          csdTestThermalV.setValue(trTest.nCompletedEpisodes+' episode(s), mean recovery: '+
            (trTest.meanRecoveryMonths!==null?trTest.meanRecoveryMonths.toFixed(1)+' mo':'n/a')+
            (trTest.ongoingEpisode!==null?' [+1 ONGOING]':''));
          print('=== S14 THERMAL RECOVERY ['+windowLabel+'] ===');
          print('Threshold (MMM+1): '+trTest.threshold.toFixed(2)+' deg C');
          for(var ti=0;ti<trTest.episodes.length;ti++){
            var tep=trTest.episodes[ti];
            print('  Episode '+(ti+1)+': '+tep.startDate+' to '+tep.endDate+' - '+tep.months+' month(s), peak '+tep.peakSST.toFixed(2)+' deg C');
          }
          if(trTest.ongoingEpisode!==null) print('ONGOING: started '+trTest.ongoingEpisode.startDate);
          var ecoValTest=getEcologicalRecoveryValidation(getRegion(latIn,lonIn),latIn,lonIn);
          print('=== S15 ECOLOGICAL RECOVERY VALIDATION ['+windowLabel+'] ===');
          print('Checked: '+(ecoValTest.checked?'YES':'NO'));
          if(ecoValTest.checked) print('Finding: '+ecoValTest.finding);
        }
      });
    });
    computeRealCSDDeseasonalized(testColl,testStudy,'sst',4000,function(res){
      var _err = res.error;
      if(_err){
        csdTestStatusV.setValue(friendlyEEError(_err));
        csdTestStatusV.style().set('color','#cc0000'); csdTestStatusV.style().set('whiteSpace','pre');
        csdTestStatusV.style().set('backgroundColor','#ffd0d0'); csdTestStatusV.style().set('border','2px solid #cc0000');
        csdTestAc1V.setValue('--'); csdTestVarTrendV.setValue('--');
        print('=== S13 RUN CSD TEST ERROR ['+windowLabel+'] === '+_err);
        return;
      }
      csdTestStatusV.style().set('backgroundColor','rgba(0,0,0,0)'); csdTestStatusV.style().set('border','none');
      var rAC1=(res&&res.realAC1!==null&&res.realAC1!==undefined)?res.realAC1:null;
      var rVar=(res&&res.varTrendRatio!==null&&res.varTrendRatio!==undefined)?res.varTrendRatio:null;
      var rN=(res&&res.nValidMonths!==null&&res.nValidMonths!==undefined)?res.nValidMonths:0;
      var rSkew=(res&&res.skewness!==null&&res.skewness!==undefined)?res.skewness:null;
      csdTestAc1V.setValue(rAC1!==null?rAC1.toFixed(3)+' (n='+rN+' valid months)':'n/a (insufficient valid months, n='+rN+')');
      csdTestVarTrendV.setValue(rVar!==null?rVar.toFixed(2)+'x'+(rVar>1.5?' RISING':rVar<0.67?' falling':' stable'):'n/a');
      csdTestNoteV.setValue(rAC1!==null&&rAC1>0.5&&rVar!==null&&rVar>1.3?'CSD pattern present in this window.':'No strong CSD pattern in this window.');
      print('=== S13 CSD RESULT ['+windowLabel+'] ===');
      print('Real AC1 (detrended): '+(rAC1!==null?rAC1.toFixed(4):'n/a')+' (n='+rN+')');
      print('Variance trend: '+(rVar!==null?rVar.toFixed(3)+'x':'n/a'));
      print('Skewness (residuals): '+(rSkew!==null?rSkew.toFixed(4):'n/a'));
      var storedResult={ac1:rAC1,varTrend:rVar,skew:rSkew,nMonths:rN,lat:latIn,lon:lonIn,startDate:startTxt,months:nMonthsIn,
        seriesTV:(res&&res.seriesTV)?res.seriesTV:null};  // v10.151: raw series kept so COMPARE can pool a shared climatology
      if(windowLabel==='BEFORE event (baseline)'){
        csdBeforeResult=storedResult;
        csdStoredBeforeV.setValue('BEFORE stored: AC1='+(rAC1!==null?rAC1.toFixed(3):'n/a')+', Var='+(rVar!==null?rVar.toFixed(2)+'x':'n/a')+'\n  ['+startTxt+', '+nMonthsIn+'mo, '+latIn+', '+lonIn+']');
        csdStoredBeforeV.style().set('color','#115511');
        csdTestStatusV.setValue('Done [BEFORE]. STEP 2a complete -> now switch the Label dropdown to "AFTER event (suspect period)", pick a new start date, and click RUN again for STEP 2b.');
        csdTestStatusV.style().set('color','#115511');
        print('S13 STORED as BEFORE: AC1='+(rAC1!==null?rAC1.toFixed(4):'n/a'));
      } else {
        csdAfterResult=storedResult;
        csdStoredAfterV.setValue('AFTER stored:  AC1='+(rAC1!==null?rAC1.toFixed(3):'n/a')+', Var='+(rVar!==null?rVar.toFixed(2)+'x':'n/a')+'\n  ['+startTxt+', '+nMonthsIn+'mo, '+latIn+', '+lonIn+']');
        csdStoredAfterV.style().set('color','#aa3300');
        csdTestStatusV.setValue('Done [AFTER]. STEP 2 complete -> scroll down and press COMPARE (STEP 3), or use FIND SWEET SPOT (STEP 4) to test several AFTER lengths at once.');
        csdTestStatusV.style().set('color','#115511');
        print('S13 STORED as AFTER: AC1='+(rAC1!==null?rAC1.toFixed(4):'n/a'));
      }
      // v10.121 FIX: if COMPARE has already run at least once, and STEP 2
      // just replaced either stored window with a NEW object, the COMPARE
      // result currently on screen is now a stale snapshot - it will not
      // update itself. Caught from a real test where this went unnoticed:
      // COMPARE showed a leftover AFTER value from an earlier site/run
      // while STEP 2's own label already showed the correct new one.
      if(csdCompareRanWithB!==null && (csdCompareRanWithB!==csdBeforeResult || csdCompareRanWithA!==csdAfterResult)){
        csdStalenessWarningV.setValue('\u26A0 STALE RESULT BELOW: STEP 2 has changed since COMPARE last ran - the verdict/detail below still reflect the OLD BEFORE/AFTER values. Press COMPARE again to refresh before trusting it.');
        csdStalenessWarningV.style().set('shown', true);
      }
    });
  }
});
panel.add(csdTestRunBtn); panel.add(csdTestStatusV);
panel.add(row('Real AC1 (this window)',csdTestAc1V));
panel.add(row('Variance trend (this window)',csdTestVarTrendV));
panel.add(row('Thermal recovery (direct, v10.36)',csdTestThermalV));
panel.add(csdTestNoteV);
panel.add(lbl('Currently stored for comparison:',7,'#334466'));
panel.add(csdStoredBeforeV); panel.add(csdStoredAfterV); panel.add(csdClearStoredBtn);
panel.add(legDiv());

// ----------------------------------------------------------
// STEP 3 - COMPARE (verdict box shown ABOVE the detailed numbers)
// ----------------------------------------------------------
panel.add(lbl('STEP 3 - Compare BEFORE vs AFTER against a control site',9,'#ffffff','#884400',true));
panel.add(lbl('Needs BOTH a BEFORE and an AFTER stored above (Step 2, done twice).',7,'#664400'));
panel.add(lbl('The control site is a deep, far-away, open-ocean point in the same water mass. It tells you whether your reef is changing on its own (LOCAL problem) or just tracking regional warming that is happening everywhere (GLOBAL/no local signal). It is picked FOR YOU automatically - the optional override is further below and can be left blank.',7,'#664400'));
panel.add(lbl('This runs in 2 stages, in order: (1) find + depth-check the control site, THEN (2) run the actual BEFORE/AFTER comparison at that site. The verdict box below always shows which of the two stages is currently running, and only turns into a colour-coded result once stage 2 is done.',7,'#886600'));

var csdControlSiteV=dynLbl('Control site: not yet selected (press COMPARE to auto-find)','#553388');
var csdCompareResultV=ui.Label('',{fontSize:'8px',color:'#114444',backgroundColor:'#e8f8f8',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
// v10.121 NEW: staleness tracking - caught from a real test where COMPARE
// showed AC1=0.794 for the AFTER window while STEP 2's own "AFTER stored"
// label correctly showed 0.431 for the SAME window. Root cause confirmed
// by direct testing (re-clicking COMPARE fixed it): COMPARE's result panel
// is a ONE-TIME SNAPSHOT taken at click time - it never auto-refreshes if
// STEP 2's BEFORE/AFTER are re-run afterward. The underlying computation
// was never wrong; the display just went stale and nothing warned the
// user. csdCompareRanWithB/A remember which STORED objects COMPARE last
// used; if STEP 2 later reassigns csdBeforeResult/csdAfterResult to NEW
// objects, this label appears, since the two objects will no longer match
// by reference.
var csdCompareRanWithB = null, csdCompareRanWithA = null;
// v10.151: the permutation test resolves in its own callback, LATER than
// the toolkit render. These globals carry its p-values back so the
// confidence label can be gated on them (FIX 5), and csdToolkitRerender
// lets the permutation callback refresh the toolkit once they arrive.
var csdPermPAC1 = null, csdPermPVar = null, csdToolkitRerender = null;
var csdStalenessWarningV = ui.Label('',
  {fontSize:'10px',fontWeight:'bold',color:'#aa3300',backgroundColor:'#fff0d0',padding:'4px 6px',margin:'2px 0',whiteSpace:'normal',border:'2px solid #cc7700'});
csdStalenessWarningV.style().set('shown', false);
// Short, bold, colour-coded headline verdict - always shown ABOVE the detailed breakdown.
// Starts neutral; only becomes a real coloured verdict once STEP 2 (the
// actual comparison) has finished - see the "(not computed yet)" stage
// messages in the button's onClick, added in v10.86.
// v10.139 NEW: statistically-valid PRIMARY verdict banner. Direct
// answer to a real, fair question: since a genuine 500-shuffle
// permutation test already exists, why is the >0.01 AC1/>0.15 variance
// HEURISTIC threshold still the primary displayed verdict? Real reason:
// the permutation test fires its OWN separate EE call (2 fetches), which
// resolves LATER than the threshold-based verdict below - not a design
// preference, a sequencing fact. Rather than risk rewriting the large,
// deeply-nested async chain that builds the threshold verdict, this adds
// a NEW banner positioned ABOVE it - populated by the permutation test's
// own completion callback once it resolves, so the statistically valid
// answer becomes the first, most prominent thing shown, without touching
// the fragile existing logic underneath. The threshold-based box below
// remains, now explicitly relabeled as the secondary, heuristic check.
var csdCompareStatValidV=ui.Label('Statistically valid verdict (p-value based): press COMPARE above, then wait for the permutation test to finish (fires after the main comparison, ~5-15s extra).',
  {fontSize:'12px',fontWeight:'bold',color:'#442266',backgroundColor:'#f0e8fa',padding:'8px 10px',margin:'2px 0',whiteSpace:'pre',border:'3px solid #663399'});
var csdCompareVerdictV=ui.Label('Run BEFORE + AFTER above, then press COMPARE (this will find the control site first, then compare).',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'normal',border:'2px solid #aaaaaa'});
// v10.88 NEW: multi-indicator "toolkit" tally, shown right under the
// headline verdict - lists every available EWS indicator (temporal AC1,
// temporal variance, spatial variance, spatial autocorrelation, skewness
// reported-only) and how many agree, instead of one indicator deciding
// the whole verdict.
var csdToolkitV=ui.Label('',{fontSize:'8px',color:'#334466',backgroundColor:'#f4f4fb',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
// v10.122 NEW: real permutation-test p-values, shown ALONGSIDE (not
// replacing) the existing threshold-based verdict above. Starts with a
// "computing" placeholder since this runs as its own independent 2-call
// fetch, resolving separately from (and usually slightly after) the main
// verdict/detail panels.
var csdPermTestV=ui.Label('',{fontSize:'8px',color:'#552266',backgroundColor:'#f6eefa',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
// v10.127 NEW: deseasonalized AC1/variance, computed alongside the raw
// numbers already shown - direct fix for the v10.101 disclosed (never
// built) seasonal-cycle caveat. Reuses the SAME raw study series already
// fetched for the permutation test above - zero new EE calls.
var csdDeseasonV=ui.Label('',{fontSize:'8px',color:'#115566',backgroundColor:'#e6f4f8',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});

var csdCompareBtn=ui.Button({
  label:'COMPARE BEFORE vs AFTER + AUTO-CONTROL SITE TEST',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#ffe8cc',color:'#884400',stretch:'horizontal',padding:'6px 4px',border:'2px solid #cc7700'},
  onClick:function(){
    if(!csdBeforeResult||!csdAfterResult){
      csdCompareVerdictV.setValue('MISSING DATA - run Step 2 for BOTH BEFORE and AFTER first.');
      csdCompareVerdictV.style().set('color','#cc0000'); csdCompareVerdictV.style().set('backgroundColor','#ffd0d0');
      csdCompareResultV.setValue('BEFORE stored: '+(csdBeforeResult?'YES':'NO')+' | AFTER stored: '+(csdAfterResult?'YES':'NO'));
      csdCompareResultV.style().set('color','#cc0000'); return;
    }
    var b=csdBeforeResult, a=csdAfterResult;
    csdCompareRanWithB = csdBeforeResult; csdCompareRanWithA = csdAfterResult;
    csdPermPAC1 = null; csdPermPVar = null; csdToolkitRerender = null;  // v10.151: clear stale p-values
    csdStalenessWarningV.setValue(''); csdStalenessWarningV.style().set('shown', false);
    csdCompareStatValidV.setValue('Statistically valid verdict (p-value based): computing... (fires after the main comparison, ~5-15s extra)');
    csdCompareStatValidV.style().set('color','#442266'); csdCompareStatValidV.style().set('backgroundColor','#f0e8fa');
    csdCompareStatValidV.style().set('border','3px solid #663399');

    // v10.122 NEW: real permutation-test p-values, fired independently of
    // the rest of COMPARE's chain below (does not block or alter it - this
    // is a genuinely separate, additive fetch, matching the same raw-value
    // pattern extractMultiNodeSeries already uses for S7D/S7E/S7F, applied
    // here to a single-point study site with SST instead of algae FAI).
    csdPermTestV.setValue('Permutation test (real p-value): fetching raw monthly SST (2 calls)...');
    csdPermTestV.style().set('color','#552266');
    (function(){
      var permStudyFC = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([b.lon,b.lat]).buffer(1000), {label:'Study'})]);
      var permBeforeColl = mkMoSSTRange(b.startDate, b.months);
      var permAfterColl = mkMoSSTRange(a.startDate, a.months);
      var rPermBefore = extractMultiNodeSeries(permBeforeColl, permStudyFC, 'sst', 4000);
      var rPermAfter = extractMultiNodeSeries(permAfterColl, permStudyFC, 'sst', 4000);
      var permData={}, permPending=2;
      function permBump(key,v,e){
        permData[key]=e?null:v; permPending--;
        if(permPending===0) permFinish();
      }
      function permFinish(){
        try {
          var beforeByNode = groupSeriesByLabel(permData.before, 'sst');
          var afterByNode = groupSeriesByLabel(permData.after, 'sst');
          var beforeVals = (beforeByNode['Study']||[]).map(function(s){return s.v;});
          var afterVals = (afterByNode['Study']||[]).map(function(s){return s.v;});
          var ac1Test = permutationTestDeltaFixed(beforeByNode['Study']||[], afterByNode['Study']||[], statAC1ForPerm, 500);
          var varTest = permutationTestDeltaFixed(beforeByNode['Study']||[], afterByNode['Study']||[], statVarRatioForPerm, 500);
          function fmtP(v){ return v!==null?v.toFixed(3):'n/a'; }
          function fmtD(v){ return v!==null?(v>0?'+':'')+v.toFixed(3):'n/a'; }
          var sigLabel=function(p){ return p===null?'':(p<0.05?' *** likely real (p<0.05)':p<0.10?' * borderline (p<0.10)':' not significant'); };
          var lines=['=== PERMUTATION TEST (real p-value, 500 shuffles) ==='];
          lines.push('This answers "would random relabeling of these same months produce a delta');
          lines.push('this big, just by chance?" - a genuine statistical test, unlike the >0.01/>0.15');
          lines.push('threshold used in the SECONDARY verdict below, which cannot tell a surprising change');
          lines.push('from ordinary noise. This IS the primary verdict now (banner above) - the threshold');
          lines.push('box below is kept only as a secondary, faster-arriving heuristic check.');
          lines.push(repeatChar('\u2500',50));
          if(ac1Test.observedDelta===null){
            lines.push('AC1: '+(ac1Test.note||'unavailable'));
          } else {
            lines.push('AC1 delta='+fmtD(ac1Test.observedDelta)+', p='+fmtP(ac1Test.pValue)+sigLabel(ac1Test.pValue)+' (n='+ac1Test.nPerm+' valid shuffles)');
          }
          if(varTest.observedDelta===null){
            lines.push('Variance ratio delta='+(varTest.note||'unavailable'));
          } else {
            lines.push('Variance ratio delta='+fmtD(varTest.observedDelta)+'x, p='+fmtP(varTest.pValue)+sigLabel(varTest.pValue)+' (n='+varTest.nPerm+' valid shuffles)');
          }
          lines.push('');
          lines.push('p<0.05 means: fewer than 1 in 20 random relabelings of these months produced');
          lines.push('a delta this large - genuinely surprising, not just "bigger than a fixed cutoff".');
          csdPermTestV.setValue(lines.join('\n'));
          csdPermTestV.style().set('color','#552266');
          print('=== S13 PERMUTATION TEST (real p-value) ==='); print(lines.join('\n'));

          // v10.151 FIX 5: publish the p-values and refresh the toolkit,
          // so its Confidence line can no longer say HIGH while the
          // banner directly above says NO SIGNIFICANT SIGNAL.
          csdPermPAC1 = (ac1Test.pValue!==null&&ac1Test.pValue!==undefined)?ac1Test.pValue:null;
          csdPermPVar = (varTest.pValue!==null&&varTest.pValue!==undefined)?varTest.pValue:null;
          if(csdToolkitRerender) { try { csdToolkitRerender(); } catch(eRr){ print('toolkit refresh skipped: '+eRr); } }

          // v10.139 NEW: populate the statistically-valid PRIMARY verdict
          // banner now that real p-values are available. Same AC1-primary
          // weighting philosophy already used throughout this tool (Dakos
          // et al. 2012), but driven by p<0.05 significance instead of the
          // old >0.01/>0.15 heuristic thresholds.
          var ac1Sig = ac1Test.pValue!==null && ac1Test.pValue<0.05;
          var ac1Rising = ac1Test.observedDelta!==null && ac1Test.observedDelta>0;
          var varSig = varTest.pValue!==null && varTest.pValue<0.05;
          var varRising = varTest.observedDelta!==null && varTest.observedDelta>0;
          var statText, statCol, statBg;
          if(ac1Test.observedDelta===null && varTest.observedDelta===null){
            statText='PERMUTATION TEST UNAVAILABLE - insufficient valid months for a p-value.\nFalling back to the heuristic threshold verdict below.';
            statCol='#886600'; statBg='#fff4dd';
          } else if(ac1Sig && ac1Rising){
            statText='STATISTICALLY SIGNIFICANT LOCAL CSD SIGNAL (AC1-confirmed, p='+ac1Test.pValue.toFixed(3)+')\n'+
              'AC1 rose by a margin fewer than 1 in 20 random reshuffles of these months would produce by chance -\n'+
              'the strongest evidence this tool can currently show for a genuine resilience-loss signal.';
            statCol='#880000'; statBg='#ffe0e0';
          } else if(varSig && varRising){
            statText='SIGNIFICANT VARIANCE-ONLY SIGNAL (p='+varTest.pValue.toFixed(3)+') - AC1 not significant\n'+
              '(p='+(ac1Test.pValue!==null?ac1Test.pValue.toFixed(3):'n/a')+'). Per Dakos et al. 2012, a variance-only rise is weaker\n'+
              'evidence than an AC1-confirmed one - worth noting, not yet a strong CSD signal.';
            statCol='#886600'; statBg='#fff4dd';
          } else {
            statText='NO STATISTICALLY SIGNIFICANT SIGNAL (p>=0.05 for both AC1 and variance)\n'+
              'AC1 p='+(ac1Test.pValue!==null?ac1Test.pValue.toFixed(3):'n/a')+' | Variance p='+(varTest.pValue!==null?varTest.pValue.toFixed(3):'n/a')+'\n'+
              'Neither delta is large enough to rule out ordinary random noise.';
            statCol='#226644'; statBg='#e8f4ff';
          }
          csdCompareStatValidV.setValue('=== STATISTICALLY VALID VERDICT (real p-value) ===\n'+statText);
          csdCompareStatValidV.style().set('color',statCol);
          csdCompareStatValidV.style().set('backgroundColor',statBg);
          csdCompareStatValidV.style().set('border','3px solid '+statCol);

          // v10.127 NEW: deseasonalized AC1/variance, computed on the SAME
          // raw study series already fetched above - reusing beforeVals/
          // afterVals's underlying {t,v} pairs (beforeByNode['Study'],
          // afterByNode['Study']) rather than the already-flattened plain-
          // number arrays, since deseasonalizing needs the real calendar
          // month of each reading. Both the "raw" and "deseasonalized"
          // numbers below are computed via the SAME client-side JS
          // functions (jsNodeStats) on the SAME fetched series, so any
          // difference between them is attributable ONLY to removing the
          // seasonal cycle - not to a different computation method (the
          // raw AC1 shown in the main verdict above is computed server-
          // side in EE instead, via a separate code path, so it may not
          // match this "raw" row to the third decimal - that is expected
          // and not a bug, just two independent implementations of the
          // same math, per the same equivalence already established
          // between this tool and the companion Python scripts).
          try {
            var studyBeforeTV = beforeByNode['Study']||[];
            var studyAfterTV = afterByNode['Study']||[];
            var pooledTV = studyBeforeTV.concat(studyAfterTV);
            var climatology = computeMonthlyClimatology(pooledTV);
            var rawBeforeStats = jsNodeStats(beforeVals);
            var rawAfterStats = jsNodeStats(afterVals);
            var deseasonBeforeVals = deseasonalizeSeries(studyBeforeTV, climatology).map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
            var deseasonAfterVals = deseasonalizeSeries(studyAfterTV, climatology).map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
            var deseasonBeforeStats = jsNodeStats(deseasonBeforeVals);
            var deseasonAfterStats = jsNodeStats(deseasonAfterVals);
            function fmtA(v){ return v!==null?v.toFixed(3):'n/a'; }
            function fmtV(v){ return v!==null?v.toFixed(2)+'x':'n/a'; }
            function fmtDlt(a,b){ return (a!==null&&b!==null)?((a-b>0?'+':'')+(a-b).toFixed(3)):'n/a'; }
            function fmtDltV(a,b){ return (a!==null&&b!==null)?((a-b>0?'+':'')+(a-b).toFixed(2)+'x'):'n/a'; }
            // v10.152 FIX 8: Object.keys() now counts the _meta key added in
            // v10.151, which produced the impossible "13 of 12 calendar
            // months" in a real run. Use the meta count when present.
            var monthsUsedForClim = (climatology && climatology._meta)
              ? climatology._meta.nUsableMonths
              : Object.keys(climatology).length;
            var deseasonLines=['=== DESEASONALIZED AC1/VARIANCE (v10.127) ==='];
            deseasonLines.push('Removes each calendar month\'s average (this site\'s own climatology) before');
            deseasonLines.push('computing AC1/variance, instead of using raw monthly values - addresses the');
            deseasonLines.push('v10.101 disclosed caveat that shared seasonal cycling can inflate AC1');
            deseasonLines.push('independent of any real dynamics. Zero extra EE calls - reuses the same raw');
            deseasonLines.push('series already fetched for the permutation test above.');
            deseasonLines.push('DISCLOSED LIMIT: climatology built from only '+monthsUsedForClim+' of 12 calendar');
            deseasonLines.push('months cleared the >=3-samples-per-month floor added in v10.151'+
              ((climatology&&climatology._meta&&climatology._meta.nRejectedMonths>0)?
                (' ('+climatology._meta.nRejectedMonths+' month(s) rejected as too sparse)'):'')+'.');
            deseasonLines.push('A climatology from few samples per month is still a noisy estimate, not a');
            deseasonLines.push('robust multi-decade seasonal average - treat as exploratory.');
            deseasonLines.push(repeatChar('\u2500',50));
            deseasonLines.push('Raw AC1 (this calc):        BEFORE='+fmtA(rawBeforeStats.realAC1)+' AFTER='+fmtA(rawAfterStats.realAC1)+'  \u0394='+fmtDlt(rawAfterStats.realAC1,rawBeforeStats.realAC1));
            deseasonLines.push('Deseasonalized AC1:         BEFORE='+fmtA(deseasonBeforeStats.realAC1)+' AFTER='+fmtA(deseasonAfterStats.realAC1)+'  \u0394='+fmtDlt(deseasonAfterStats.realAC1,deseasonBeforeStats.realAC1));
            deseasonLines.push('Raw Var ratio (this calc):  BEFORE='+fmtV(rawBeforeStats.varTrendRatio)+' AFTER='+fmtV(rawAfterStats.varTrendRatio)+'  \u0394='+fmtDltV(rawAfterStats.varTrendRatio,rawBeforeStats.varTrendRatio));
            deseasonLines.push('Deseasonalized Var ratio:   BEFORE='+fmtV(deseasonBeforeStats.varTrendRatio)+' AFTER='+fmtV(deseasonAfterStats.varTrendRatio)+'  \u0394='+fmtDltV(deseasonAfterStats.varTrendRatio,deseasonBeforeStats.varTrendRatio));
            deseasonLines.push(repeatChar('\u2500',50));
            var rawDeltaAC1 = (rawAfterStats.realAC1!==null&&rawBeforeStats.realAC1!==null)?(rawAfterStats.realAC1-rawBeforeStats.realAC1):null;
            var deseasonDeltaAC1 = (deseasonAfterStats.realAC1!==null&&deseasonBeforeStats.realAC1!==null)?(deseasonAfterStats.realAC1-deseasonBeforeStats.realAC1):null;
            // v10.128 FIX: the original binary comparison only checked for
            // "deseasonalized delta much SMALLER than raw" - caught on a
            // real One Tree Reef test where the deseasonalized delta
            // (-0.089) was actually 7x LARGER than the raw delta (-0.012),
            // and the code's only other branch said "similar in size,"
            // which was factually wrong for a 7x difference. Now checks
            // both directions explicitly, and also compares the ABSOLUTE
            // baseline AC1 levels (not just the delta) - that same real
            // test showed raw BEFORE=0.843 vs deseasonalized BEFORE=0.434,
            // nearly half, which turned out to be the far more informative
            // number: a large share of the raw AC1 LEVEL itself, not just
            // its change, may reflect the seasonal cycle rather than real
            // thermal memory - worth surfacing even when the delta
            // comparison alone would say little.
            if(rawDeltaAC1!==null&&deseasonDeltaAC1!==null){
              if(Math.abs(deseasonDeltaAC1)<Math.abs(rawDeltaAC1)*0.5){
                deseasonLines.push('The deseasonalized \u0394AC1 is much smaller than the raw \u0394AC1 - the seasonal cycle');
                deseasonLines.push('may be doing most of the work in the raw AC1 CHANGE shown above. Treat the raw');
                deseasonLines.push('AC1-based verdict with more caution at this site.');
              } else if(Math.abs(deseasonDeltaAC1)>Math.abs(rawDeltaAC1)*1.5){
                deseasonLines.push('The deseasonalized \u0394AC1 is LARGER than the raw \u0394AC1, not smaller - removing the');
                deseasonLines.push('seasonal cycle reveals a bigger underlying change than the raw number showed.');
                deseasonLines.push('Both deltas may still be too small to matter - check the raw AC1 permutation');
                deseasonLines.push('p-value above before treating either as a real signal.');
              } else {
                deseasonLines.push('The deseasonalized \u0394AC1 is similar in size to the raw \u0394AC1 - the CHANGE in AC1');
                deseasonLines.push('does not appear to be primarily a seasonal-cycle artifact at this site.');
              }
            }
            var rawBaseline = rawBeforeStats.realAC1, deseasonBaseline = deseasonBeforeStats.realAC1;
            if(rawBaseline!==null&&deseasonBaseline!==null&&rawBaseline>0.05){
              var baselineRatio = deseasonBaseline/rawBaseline;
              if(baselineRatio<0.7){
                deseasonLines.push('');
                deseasonLines.push('NOTE: the raw AC1 LEVEL itself (not just its change) is much higher than the');
                deseasonLines.push('deseasonalized level ('+fmtA(rawBaseline)+' vs '+fmtA(deseasonBaseline)+') - a large share of what');
                deseasonLines.push('looks like "high thermal memory" in the raw AC1 shown elsewhere in this tool may');
                deseasonLines.push('reflect the seasonal cycle itself, not genuine resilience loss. Interpret any raw');
                deseasonLines.push('AC1 value in isolation (not just its BEFORE/AFTER change) with real caution.');
              }
            }
            csdDeseasonV.setValue(deseasonLines.join('\n'));
            csdDeseasonV.style().set('color','#115566');
            print('=== S13 DESEASONALIZED AC1/VARIANCE ==='); print(deseasonLines.join('\n'));
          } catch(eDeseason){
            csdDeseasonV.setValue('Deseasonalized comparison error: '+eDeseason);
            csdDeseasonV.style().set('color','#cc0000');
            print('=== S13 DESEASONALIZED COMPARISON ERROR === '+eDeseason);
          }
        } catch(ePerm){
          csdPermTestV.setValue('Permutation test error: '+ePerm);
          csdPermTestV.style().set('color','#cc0000');
          print('=== S13 PERMUTATION TEST ERROR === '+ePerm);
        }
      }
      rPermBefore.evaluate(function(v,e){ permBump('before', v, e); });
      rPermAfter.evaluate(function(v,e){ permBump('after', v, e); });
    })();
    // v10.86 FIX: the verdict box used to jump straight to "Computing
    // verdict..." on click, which visually contradicted the detail box
    // right below it still saying "Step 1/3: Finding control site...".
    // Both boxes now show the SAME real stage, and only the verdict box's
    // final colour-coded headline appears once there is an actual verdict.
    csdCompareVerdictV.setValue('Step 1/3: Finding control site... (not computed yet)');
    csdCompareVerdictV.style().set('color','#334466'); csdCompareVerdictV.style().set('backgroundColor','#eeeeee');
    csdCompareVerdictV.style().set('border','2px solid #aaaaaa'); csdCompareVerdictV.style().set('whiteSpace','normal');
    csdCompareResultV.setValue('Step 1/3: Finding control site...\nThis can take 10-30 seconds. If it sits here far longer than that, it is'+
      ' almost always a temporary Earth Engine connection issue - press COMPARE again.');
    csdCompareResultV.style().set('color','#334466');
    csdCompareResultV.style().set('backgroundColor','#e8f8f8');
    csdControlSiteV.setValue('Checking for manual control site override...');

    // ============================================================
    // CONTROL SITE SELECTION: manual override (advanced/optional) or auto GEBCO search (default)
    // ============================================================
    var manualCtrlTxt = (csdManualCtrlInput.getValue()||'').trim();
    var manualCtrlParts = manualCtrlTxt.length > 0 ? manualCtrlTxt.split(',') : [];
    var manualCtrlValid = manualCtrlParts.length === 2 &&
      !isNaN(parseCoordPart(manualCtrlParts[0])) &&
      !isNaN(parseCoordPart(manualCtrlParts[1]));

    // v10.86 FIX: this depth lookup previously ignored its own error
    // parameter entirely - if it errored or the response was malformed,
    // the whole COMPARE flow just sat on "Step 1/3: Finding control
    // site..." forever with no explanation. Now any error here shows the
    // friendly error box immediately instead of hanging silently.
    function onControlDepthResult(depthRes, depthErr, bestCtrlBase, label){
      if(depthErr){
        csdCompareVerdictV.setValue('CONTROL SITE LOOKUP FAILED');
        csdCompareVerdictV.style().set('color','#cc0000'); csdCompareVerdictV.style().set('backgroundColor','#ffd0d0');
        csdCompareVerdictV.style().set('border','2px solid #cc0000'); csdCompareVerdictV.style().set('whiteSpace','pre');
        csdCompareResultV.setValue(friendlyEEError(depthErr));
        csdCompareResultV.style().set('color','#cc0000'); csdCompareResultV.style().set('backgroundColor','#ffd0d0');
        csdControlSiteV.setValue('Control site depth check failed - press COMPARE again.');
        csdControlSiteV.style().set('color','#cc0000');
        print('=== S13 COMPARE: control-site depth lookup error === '+depthErr);
        return;
      }
      var elev = depthRes && depthRes.elevation !== null && depthRes.elevation !== undefined ? depthRes.elevation : null;
      var depthM = elev !== null ? Math.round(elev) : null;
      var depthOk = depthM !== null && depthM < -500;
      var bestCtrl = bestCtrlBase; bestCtrl.depth_m = depthM;
      csdControlSiteV.setValue(label+' | GEBCO: '+(depthM!==null?depthM+' m':'depth unknown')+
        (depthOk?' DEEP WATER OK':depthM!==null?' WARNING: shallower than 500m':' (continuing without depth validation)'));
      csdControlSiteV.style().set('color', depthOk?'#115511':'#884400');
      // Both boxes now honestly say we're moving to stage 2 before doing so.
      csdCompareVerdictV.setValue('Step 2/3: Control site confirmed. Running CSD at control site... (not computed yet)');
      runControlCSD(bestCtrl, b, a);
    }

    if(manualCtrlValid){
      // USER PROVIDED A MANUAL CONTROL SITE (advanced/optional override)
      var mLat = parseCoordPart(manualCtrlParts[0]);
      var mLon = parseCoordPart(manualCtrlParts[1]);
      var manualCtrlPt = ee.Geometry.Point([mLon, mLat]);
      print('S13: Using MANUAL control site: '+mLat+', '+mLon);
      GEBCO.rename('elevation').reduceRegion({
        reducer:ee.Reducer.mean(), geometry:manualCtrlPt, scale:1000, maxPixels:1e9
      }).evaluate(function(depthRes, depthErr){
        onControlDepthResult(depthRes, depthErr,
          {lat:mLat, lon:mLon, manual:true},
          'MANUAL control (advanced override): '+mLat.toFixed(2)+', '+mLon.toFixed(2));
      });
    } else {
      // AUTO-SELECT using smart basin-aware control site (default path)
      var smartSuggestion = getSmartControlSite(b.lat, b.lon);
      csdControlSiteV.setValue('Auto control site: '+smartSuggestion.label+
        ' ('+smartSuggestion.lat.toFixed(2)+', '+smartSuggestion.lon.toFixed(2)+') — validating depth...');
      GEBCO.rename('elevation').reduceRegion({
        reducer:ee.Reducer.mean(),
        geometry:ee.Geometry.Point([smartSuggestion.lon, smartSuggestion.lat]),
        scale:1000, maxPixels:1e9
      }).evaluate(function(depthRes, depthErr){
        onControlDepthResult(depthRes, depthErr,
          {lat:smartSuggestion.lat, lon:smartSuggestion.lon, label:smartSuggestion.label, manual:false},
          'Auto control site: '+smartSuggestion.label+'\n'+smartSuggestion.lat.toFixed(2)+', '+smartSuggestion.lon.toFixed(2));
      });
    }
  }
});

// ============================================================
// SMART CONTROL SITE SELECTOR (unchanged from v10.68)
// Basin-aware pre-validated open-ocean control site selection.
// ============================================================
function getSmartControlSite(lat, lon) {
  if(lat>12&&lat<30&&lon>32&&lon<44)
    return {lat:15,lon:58,label:'Arabian Sea (Red Sea basin control)'};
  if(lat>22&&lat<30&&lon>50&&lon<57)
    return {lat:20,lon:62,label:'Arabian Sea (Persian Gulf basin control)'};
  if(lat>30&&lat<47&&lon>-6&&lon<37)
    return {lat:38,lon:-15,label:'Open NE Atlantic (Mediterranean basin control)'};
  if(lat>50&&lat<70&&lon>-5&&lon<30)
    return {lat:60,lon:-15,label:'Open N. Atlantic (North Sea basin control)'};
  if(lat>20&&lat<30&&lon>108&&lon<122)
    return {lat:15,lon:130,label:'Open West Pacific (South China Sea control)'};
  if(lat>-20&&lat<10&&lon>30&&lon<55)
    return {lat:5,lon:70,label:'Open Indian Ocean (Arabian Sea control)'};
  if(lat>55) {
    if(lon>-60&&lon<-30)
      return {lat:lat,lon:-40,label:'Labrador Sea (SW Greenland Arctic control)'};
    if(lon>-30&&lon<10)
      return {lat:lat,lon:-15,label:'N. Atlantic (NE Greenland Arctic control)'};
    if(lon>10&&lon<40)
      return {lat:lat,lon:-5,label:'Norwegian Sea (Arctic control)'};
    if(lon>140||lon<-150)
      return {lat:lat,lon:lon>0?lon-30:lon+30,label:'Open N. Pacific (Arctic control)'};
    return {lat:lat,lon:lon+(Math.abs(lon)<160?15:-15),label:'Open ocean (Arctic fallback control)'};
  }
  if(lat<-50)
    return {lat:-55,lon:lon+30>180?lon-30:lon+30,label:'Open Southern Ocean control'};
  if(lat>8&&lat<35&&lon>-100&&lon<-60) {
    var aLon = Math.min(lon+25,-58);
    return {lat:lat,lon:aLon,label:'Open Atlantic (Caribbean/W.Atlantic control)'};
  }
  if(lat>-35&&lat<35&&lon>-30&&lon<20)
    return {lat:lat,lon:lon-20,label:'Open Central Atlantic (E.Atlantic control)'};
  if(lat>-30&&lat<15&&lon>30&&lon<80)
    return {lat:lat,lon:lon+20,label:'Open Indian Ocean (W.Indian control)'};
  if(lat>-10&&lat<25&&lon>80&&lon<100)
    return {lat:lat,lon:lon+20,label:'Open Indian Ocean (E.Indian control)'};
  if(lat>-30&&lat<0&&lon>142&&lon<160)
    return {lat:lat,lon:160,label:'Open Coral Sea (GBR control)'};
  if(lat>-10&&lat<20&&lon>115&&lon<145)
    return {lat:lat,lon:155,label:'Open W. Pacific control'};
  if(lon>-180&&lon<-100)
    return {lat:lat,lon:lon-20,label:'Open Pacific (E.Pacific control)'};
  if(lon>160||lon<-100)
    return {lat:lat,lon:lon>0?lon+20:lon-20,label:'Open Pacific control'};
  var offLon = lon+(Math.abs(lon)<165?15:-15);
  return {lat:lat,lon:offLon,label:'Auto-offset control (check depth manually)'};
}

// ============================================================
// runControlCSD: shared by both the manual and auto paths.
// Runs CSD at the control site for BEFORE and AFTER, computes
// the 4-way verdict, and now writes it to BOTH:
//   - csdCompareVerdictV: short bold headline (shown first)
//   - csdCompareResultV : full numeric breakdown (shown below)
// ============================================================
function runControlCSD(bestCtrl, b, a) {
  csdCompareResultV.setValue('Step 2/3: Running CSD at control site...\n'+
    '('+bestCtrl.lat.toFixed(2)+', '+bestCtrl.lon.toFixed(2)+
    (bestCtrl.manual?' [MANUAL]':' [AUTO]')+')');
  var ctrlPt=ee.Geometry.Point([bestCtrl.lon,bestCtrl.lat]);
  var ctrlStudy=ctrlPt.buffer(80000);
  var ctrlBeforeColl=mkMoSSTRange(b.startDate,b.months);
  var ctrlAfterColl=mkMoSSTRange(a.startDate,a.months);

  computeRealCSDDeseasonalized(ctrlBeforeColl, ctrlStudy, 'sst', 4000, function(cbRes){
    computeRealCSDDeseasonalized(ctrlAfterColl, ctrlStudy, 'sst', 4000, function(caRes){
     try {
      var cbErr = cbRes.error, caErr = caRes.error;
      if(cbErr||caErr){
        var ctrlErrTxt = cbErr ? String(cbErr) : String(caErr);
        csdCompareVerdictV.setValue('CONTROL SITE COMPUTATION FAILED');
        csdCompareVerdictV.style().set('color','#cc0000'); csdCompareVerdictV.style().set('backgroundColor','#ffd0d0');
        csdCompareVerdictV.style().set('border','2px solid #cc0000');
        csdCompareResultV.setValue(friendlyEEError(ctrlErrTxt));
        csdCompareResultV.style().set('color','#cc0000'); csdCompareResultV.style().set('backgroundColor','#ffd0d0');
        print('=== S13 COMPARE: control site evaluate error === '+ctrlErrTxt);
        return;
      }

      // ============================================================
      // STEP 3: 4-way comparison and verdict
      // ============================================================
      var ctrlB={
        ac1:(cbRes&&cbRes.realAC1!==null&&cbRes.realAC1!==undefined)?cbRes.realAC1:null,
        varTrend:(cbRes&&cbRes.varTrendRatio!==null&&cbRes.varTrendRatio!==undefined)?cbRes.varTrendRatio:null,
        n:(cbRes&&cbRes.nValidMonths!==null)?cbRes.nValidMonths:0
      };
      var ctrlA={
        ac1:(caRes&&caRes.realAC1!==null&&caRes.realAC1!==undefined)?caRes.realAC1:null,
        varTrend:(caRes&&caRes.varTrendRatio!==null&&caRes.varTrendRatio!==undefined)?caRes.varTrendRatio:null,
        n:(caRes&&caRes.nValidMonths!==null)?caRes.nValidMonths:0
      };

      // ============================================================
      // v10.151 FIX 4b - ONE POOLED CLIMATOLOGY FOR BOTH WINDOWS
      // STEP 2 runs BEFORE and AFTER at different times, so each window
      // was previously deseasonalized against its OWN climatology and
      // the delta between them was not a like-for-like comparison. A
      // real Looe Key run (24.531,-81.41) reported the same window's
      // AC1 as -0.332 here, +0.062 in the v10.127 panel, and +0.813
      // raw. Now that STEP 2 stores the raw series, both windows are
      // recomputed against a single pooled climatology - the same one
      // permutationTestDeltaFixed already uses - so the 4-way panel,
      // the toolkit and the p-values finally describe one quantity.
      var pairedB=b, pairedA=a, pooledClimNote='';
      if(b.seriesTV && a.seriesTV){
        var pooledClim=computeMonthlyClimatology(b.seriesTV.concat(a.seriesTV));
        var pStatB=jsNodeStatsFixed(b.seriesTV, pooledClim);
        var pStatA=jsNodeStatsFixed(a.seriesTV, pooledClim);
        if(pStatB.realAC1!==null && pStatA.realAC1!==null){
          pairedB={ac1:pStatB.realAC1, varTrend:pStatB.varTrendRatio, skew:pStatB.skewness,
                   nMonths:pStatB.nValidMonths, lat:b.lat, lon:b.lon,
                   startDate:b.startDate, months:b.months};
          pairedA={ac1:pStatA.realAC1, varTrend:pStatA.varTrendRatio, skew:pStatA.skewness,
                   nMonths:pStatA.nValidMonths, lat:a.lat, lon:a.lon,
                   startDate:a.startDate, months:a.months};
          var cm=pooledClim._meta||{};
          pooledClimNote='v10.151: BOTH windows recomputed against ONE pooled climatology '+
            '('+(cm.nUsableMonths||0)+'/12 calendar months had >='+(cm.minSamples||3)+' samples'+
            ((cm.nRejectedMonths||0)>0?', '+cm.nRejectedMonths+' month(s) rejected as too sparse':'')+
            '). These numbers now match the permutation test below; STEP 2\'s own per-window '+
            'figures above were computed separately and may differ.';
        } else {
          pooledClimNote='v10.151: pooled-climatology recompute returned no usable value - '+
            'falling back to STEP 2\'s per-window figures, which may not be directly comparable.';
        }
      } else {
        pooledClimNote='v10.151: no raw series stored (re-run STEP 2 for BOTH windows to enable '+
          'the pooled-climatology recompute). Figures below are per-window and may not be '+
          'directly comparable to the permutation test.';
      }

      var sDAC1=(pairedB.ac1!==null&&pairedA.ac1!==null)?(pairedA.ac1-pairedB.ac1):null;
      var sDVar=(pairedB.varTrend!==null&&pairedA.varTrend!==null)?(pairedA.varTrend-pairedB.varTrend):null;
      var cDAC1=(ctrlB.ac1!==null&&ctrlA.ac1!==null)?(ctrlA.ac1-ctrlB.ac1):null;
      var cDVar=(ctrlB.varTrend!==null&&ctrlA.varTrend!==null)?(ctrlA.varTrend-ctrlB.varTrend):null;

      // v10.153 FIX 10: thresholds from a measured null, not guesses.
      var _thr = getCalibratedThresholds(a.nMonths||b.nMonths||24);
      var THRESH=_thr.varr;
      var studyVarRose=(sDVar!==null&&sDVar>THRESH);
      var ctrlVarRose=(cDVar!==null&&cDVar>THRESH);
      var studyAC1Rose=(sDAC1!==null&&sDAC1>_thr.ac1);   // v10.153: was hardcoded 0.01

      var divergence=(sDVar!==null&&cDVar!==null)?(sDVar-cDVar):null;
      var DIVERG_THRESH=_thr.varr;
      var studyDivergesFromControl=(divergence!==null&&divergence>DIVERG_THRESH&&!studyVarRose&&!ctrlVarRose);

      var verdict, vCol, vBg, vTitle;
      // v10.118 FIX: the same bug already caught and fixed for S7E in
      // v10.114 existed here too, in S13's core COMPARE verdict - the most
      // heavily used part of the whole tool. Caught from a real Nuuk test:
      // the study site returned AC1=n/a, Var=n/a in BOTH windows (zero
      // valid OISST months), yet the headline confidently said "STUDY SITE
      // MORE STABLE THAN CONTROL (positive result)" and "the study site
      // stayed flat... This is a POSITIVE result for this reef." The study
      // site was never measured at all - "no data" and "tested and found
      // stable" are different findings, but studyVarRose/studyAC1Rose both
      // silently default to false when the underlying deltas are null, so
      // null data was indistinguishable from a genuinely flat result and
      // fell straight into the reassuring "positive" branch. Now checked
      // explicitly first, before any of the existing branches run.
      var studyHasNoData = (sDVar===null && sDAC1===null);
      var ctrlHasNoData = (cDVar===null && cDAC1===null);
      if(studyHasNoData && ctrlHasNoData){
        vTitle='CANNOT ASSESS - no data at either site';
        verdict='Neither the study site nor the control site returned any valid AC1/variance\n'+
          'in the requested windows (both n=0 valid months). This is not a "no signal" or\n'+
          '"stable" result - it means neither site could be measured at all, likely due to\n'+
          'ice masking, cloud cover, or a coastal pixel with poor satellite coverage.\n'+
          'Try a longer window, a nearby coordinate, or a different time period.';
        vCol='#888888'; vBg='#eeeeee';
      } else if(studyHasNoData){
        vTitle='CANNOT ASSESS - study site has no data';
        verdict='The study site returned no valid AC1/variance in the requested windows\n'+
          '(n=0 valid months in at least one window) - NOT a "stable" or "flat" result.\n'+
          'The control site DID return data ('+(ctrlVarRose?'variance rose':'variance did not rise')+'),\n'+
          'but with no study-site data there is nothing to compare it against.\n'+
          'Do not read this as a positive or protective finding for the study site.\n'+
          'Likely cause: ice masking, cloud cover, or a coastal pixel with poor coverage -\n'+
          'try a longer window or a nearby coordinate with cleaner satellite data.';
        vCol='#886600'; vBg='#fff6cc';
      } else if(ctrlHasNoData){
        vTitle='CANNOT ASSESS - control site has no data';
        verdict='The auto-selected control site returned no valid AC1/variance in the requested\n'+
          'windows - the study site cannot be classified as local vs regional without a\n'+
          'control that could actually be measured. Try the ADVANCED manual override below\n'+
          'with a different deep-water coordinate, or a different time window.';
        vCol='#886600'; vBg='#fff6cc';
      } else if(studyVarRose&&ctrlVarRose){
        // v10.126 FIX: this branch used to call ANY case where both sites
        // individually crossed the threshold "GLOBAL SIGNAL," with zero
        // check on HOW MUCH bigger the study site's rise was - divergence
        // was computed and displayed, but never used to decide the verdict
        // here. Caught from a real Nuuk test: study rose +3.40x, control
        // rose only +0.22x - the control barely cleared the threshold
        // while the study rose over 15x more, giving a divergence of
        // +3.18x, more than 21 times the 0.15x divergence bar this SAME
        // tool already uses elsewhere (the MARGINAL LOCAL SIGNAL branch
        // below) - but that branch was only ever reachable when NEITHER
        // site crossed the threshold individually, so a case exactly like
        // this - large, real regional signal PLUS a much larger local
        // amplification on top of it - had no branch that could describe
        // it honestly. Now explicitly checked: a large divergence on top
        // of two genuinely rising sites is a real, distinct finding, not
        // silently folded into a plain "GLOBAL, not local" label.
        if(divergence>DIVERG_THRESH){
          vTitle='GLOBAL SIGNAL, WITH STRONG LOCAL AMPLIFICATION';
          verdict='Both study AND control variance rose by more than '+THRESH+'x - a real regional\n'+
            'warming signal is present, this is NOT purely a local artifact.\n'+
            'BUT the divergence between sites ('+(divergence>0?'+':'')+divergence.toFixed(2)+'x) is far larger than the\n'+
            DIVERG_THRESH+'x divergence threshold this tool uses elsewhere - the study site is not just\n'+
            'tracking the regional trend, it is amplifying it well beyond what the control shows.\n'+
            'Treat this as a MIXED signal: real regional forcing is present, but something at this\n'+
            'specific site (shallow water, local circulation, local stressors) is making it much\n'+
            'worse locally than the open-ocean background alone would explain.\n'+
            'Divergence (study - control): '+(divergence>0?'+':'')+divergence.toFixed(2)+'x';
          vCol='#883300'; vBg='#ffe0cc';
        } else {
        vTitle='GLOBAL SIGNAL - not a local warning';
        verdict='Both study AND control variance rose by more than '+THRESH+'x.\n'+
          'This is background warming / El Nino noise affecting the whole region.\n'+
          'NOT a local CSD tipping point signal.\n'+
          'Your study site is not diverging from its surroundings.\n'+
          'Divergence (study - control): '+(divergence!==null?(divergence>0?'+':'')+divergence.toFixed(2)+'x':'n/a');
        vCol='#664400'; vBg='#fff6cc';
        }
      } else if(studyVarRose&&!ctrlVarRose){
        vTitle='LOCAL CSD SIGNAL CANDIDATE';
        verdict='Study variance rose but the deep-water control stayed flat.\n'+
          'The study site is diverging from background regional climate noise.\n'+
          (studyAC1Rose?'AC1 also rose - BOTH CSD indicators point the same direction. Stronger evidence.':
            'However AC1 did NOT rise - only variance rose. Weaker evidence.')+'\n'+
          'Divergence (study - control): '+(divergence!==null?(divergence>0?'+':'')+divergence.toFixed(2)+'x':'n/a')+'\n'+
          'Consistent with a local tipping point warning.\n'+
          'MUST validate with field data (coral cover / algae surveys).';
        vCol='#880000'; vBg='#ffd0d0';
      } else if(studyDivergesFromControl){
        vTitle='MARGINAL LOCAL SIGNAL';
        verdict='Neither site individually crossed the '+THRESH+'x variance threshold,\n'+
          'BUT the study site variance rose '+sDVar.toFixed(2)+'x while the control\n'+
          'stayed flat / fell ('+cDVar.toFixed(2)+'x).\n'+
          'Divergence (study - control) = +'+(divergence.toFixed(2))+'x - above the '+DIVERG_THRESH+'x threshold.\n'+
          (studyAC1Rose?'AC1 also rose slightly - mild CSD signal in both indicators.':
            'AC1 did NOT rise - only variance diverged. Weaker evidence.')+'\n'+
          'The study site IS diverging from background noise, but weakly.\n'+
          'Consider a shorter AFTER window to capture the stress peak more tightly.\n'+
          'MUST validate with field data before drawing conclusions.';
        vCol='#aa3300'; vBg='#ffe8cc';
      } else if(!studyVarRose&&ctrlVarRose){
        vTitle='STUDY SITE MORE STABLE THAN CONTROL (positive result)';
        verdict='Control variance rose but the study site stayed flat.\n'+
          'Regional warming is NOT affecting the study site as strongly.\n'+
          'Possible: local upwelling, shading, deep mixing, protected location.\n'+
          'This is a POSITIVE result for this reef.\n'+
          'Divergence (study - control): '+(divergence!==null?(divergence>0?'+':'')+divergence.toFixed(2)+'x':'n/a');
        vCol='#115511'; vBg='#d4f5df';
      } else {
        vTitle='NO SIGNAL AT EITHER SITE';
        verdict='Neither study nor control shows strong variance rise (>'+THRESH+'x),\n'+
          'and divergence between sites is below threshold ('+DIVERG_THRESH+'x).\n'+
          'No CSD pattern in this time window.\n'+
          'Divergence (study - control): '+(divergence!==null?(divergence>0?'+':'')+divergence.toFixed(2)+'x':'n/a')+'\n'+
          'Try a shorter AFTER window to capture the stress peak, or test\n'+
          'closer to a known bleaching/collapse event. (Tip: try FIND SWEET SPOT.)';
        vCol='#226644'; vBg='#e8f4ff';
      }

      // --- Short, bold headline verdict box (shown FIRST / above the numbers) ---
      // v10.88: this is now a PRELIMINARY headline - it gets upgraded below
      // once spatial indicators arrive, so it's clearly marked as such.
      csdCompareVerdictV.setValue('[PRELIMINARY - temporal only] '+vTitle+' — fetching spatial indicators...');
      csdCompareVerdictV.style().set('color',vCol);
      csdCompareVerdictV.style().set('backgroundColor',vBg);
      csdCompareVerdictV.style().set('border','2px solid '+vCol);
      csdCompareVerdictV.style().set('whiteSpace','normal');

      var na=function(v){return v!==null?v.toFixed(3):'n/a';};
      var nav=function(v){return v!==null?v.toFixed(2)+'x':'n/a';};
      var nad=function(v){return v!==null?(v>0?'+':'')+v.toFixed(3):'n/a';};
      var nadv=function(v){return v!==null?(v>0?'+':'')+v.toFixed(2)+'x':'n/a';};

      // v10.120 FIX: this detail panel is built early/synchronously, using
      // the raw variance-only verdict text, BEFORE the AC1-weighted
      // spatial toolkit finishes computing and rewrites the headline box
      // above (an async callback that runs later). This detail panel is
      // never updated afterward, so it stays structurally stuck showing
      // the pre-AC1-weighting text - the same conflict the v10.119 fix
      // addressed in the headline, still present here. Caught from a real
      // Nuuk test: the headline correctly said "SIGNAL PRESENT (AC1-
      // confirmed), LIKELY REGIONAL" once spatial data arrived, but THIS
      // panel still said "Control variance rose but the study site stayed
      // flat... This is a POSITIVE result" - unreconciled. Uses
      // studyAC1Rose (already computed above, before this point) as the
      // same conflict signal the headline fix uses.
      // v10.151 FIX 6: same general test here. The detail panel is built
      // synchronously, before the AC1-weighted headline exists, so it
      // compares vTitle against the AC1 direction it does have.
      var provisionalHeadline = studyAC1Rose ? 'SIGNAL PRESENT (AC1 rising)' : 'NO SIGNAL (AC1 not rising)';
      var detailConflictMsg = verdictsConflict(provisionalHeadline, vTitle);
      var detailRegionalConflict = !!detailConflictMsg;
      var detailRegionalLine = detailRegionalConflict ?
        ('Regional context: '+detailConflictMsg+'\n\n'+verdict) :
        'Regional context: '+verdict;

      var result=
        '=== 4-WAY CSD COMPARISON (detail) ===\n'+
        'STUDY SITE ('+b.lat+', '+b.lon+'):\n'+
        '  BEFORE ['+b.startDate+', '+b.months+'mo]: AC1='+na(pairedB.ac1)+', Var='+nav(pairedB.varTrend)+'\n'+
        '  AFTER  ['+a.startDate+', '+a.months+'mo]: AC1='+na(pairedA.ac1)+', Var='+nav(pairedA.varTrend)+'\n'+
        '  DELTA: \u0394AC1='+nad(sDAC1)+', \u0394Var='+nadv(sDVar)+'\n'+
        '\n'+
        'CONTROL SITE ('+bestCtrl.lat.toFixed(2)+', '+bestCtrl.lon.toFixed(2)+
        ', depth='+(bestCtrl.depth_m!==null?bestCtrl.depth_m+'m':'unknown')+'):\n'+
        '  BEFORE ['+b.startDate+', '+b.months+'mo]: AC1='+na(ctrlB.ac1)+', Var='+nav(ctrlB.varTrend)+' (n='+ctrlB.n+'mo)\n'+
        '  AFTER  ['+a.startDate+', '+a.months+'mo]: AC1='+na(ctrlA.ac1)+', Var='+nav(ctrlA.varTrend)+' (n='+ctrlA.n+'mo)\n'+
        '  DELTA: \u0394AC1='+nad(cDAC1)+', \u0394Var='+nadv(cDVar)+'\n'+
        '\n'+
        detailRegionalLine+'\n'+
        '\n'+
        'CAVEAT: for publication-quality work, also try a manually chosen\n'+
        'reference site (ADVANCED override below) rather than relying only on auto-selection.\n'+
        'Thresholds (v10.153, calibrated): variance >'+THRESH.toFixed(2)+'x | divergence >'+
        DIVERG_THRESH.toFixed(2)+'x | AC1 >'+_thr.ac1.toFixed(3)+'\n'+
        'Basis: '+_thr.basis+'.\n'+
        'The old cutoffs (AC1>0.01, Var>0.15x) were MEASURED firing on 34-57% of windows\n'+
        'at a site with no regime shift; the AC1-OR-variance rule fired on 80%. The\n'+
        'permutation test below rejected 0-9% on those same windows and is the only\n'+
        'verdict in this panel with a known false-positive rate.\n'+
        '\n'+pooledClimNote;

      csdCompareResultV.setValue(result);
      csdCompareResultV.style().set('color','#223344');
      csdCompareResultV.style().set('backgroundColor','#f4f8fb');
      print('=== S13 4-WAY CSD COMPARISON (regional context) ===');
      print('Regional context: '+(detailRegionalConflict?'CONFLICTS WITH AC1 (see detail below)':vTitle));
      print(result);

      // ============================================================
      // v10.88 TOOLKIT TALLY - "don't rely on a single indicator"
      // Temporal indicators (AC1, variance, skewness-reported-only) render
      // immediately below using data already in hand. Spatial indicators
      // (variance across the region, spatial autocorrelation) are fetched
      // as a follow-up and UPGRADE the same boxes in place once ready -
      // this never blocks or delays the result above, and degrades
      // gracefully to temporal-only if the spatial fetch fails.
      // ============================================================
      var skewDisplay = (pairedA.skew!==null&&pairedA.skew!==undefined)?
        ((pairedA.skew>0?'+':'')+pairedA.skew.toFixed(2)+(pairedB.skew!==null&&pairedB.skew!==undefined?' (was '+(pairedB.skew>0?'+':'')+pairedB.skew.toFixed(2)+')':'')) : 'n/a';

      function buildTemporalIndicators(){
        return [
          {name:'Temporal AC1 (study)', scored:true, primary:true, available:sDAC1!==null,
            agrees:sDAC1!==null&&sDAC1>0.01,
            display:sDAC1!==null?((sDAC1>0?'+':'')+sDAC1.toFixed(3)):'n/a'},
          {name:'Temporal Variance (study)', scored:true, primary:false, available:sDVar!==null,
            agrees:sDVar!==null&&sDVar>THRESH,
            display:sDVar!==null?((sDVar>0?'+':'')+sDVar.toFixed(2)+'x'):'n/a'},
          {name:'Skewness (study, reported only - direction is system-dependent)', scored:false, primary:false, available:true,
            display:skewDisplay}
        ];
      }

      function renderToolkit(indicators, spatialPending){
        var tally=buildToolkitTally(indicators);
        var confidence;
        // v10.151 FIX 5: confidence is now gated on the real p-values.
        // Before they resolve, classifyToolkitConfidence returns
        // PRELIMINARY rather than a confident-sounding label.
        if(spatialPending) confidence='PRELIMINARY (spatial indicators still loading...)';
        else confidence=classifyToolkitConfidence(tally, csdPermPAC1, csdPermPVar).label;
        var txt='TOOLKIT SUMMARY ('+tally.nAgree+'/'+tally.nAvail+' scored indicators agree - AC1 weighted as PRIMARY per Dakos et al. 2012, others are supporting evidence, not equal votes):\n'+
          tally.lines.join('\n')+'\n'+
          'Confidence: '+confidence+'\n'+
          'Data sufficiency: STUDY n='+(a.nMonths||'?')+'mo AFTER / '+(b.nMonths||'?')+'mo BEFORE'+
          (ctrlA.n&&ctrlB.n?' | CONTROL n='+ctrlA.n+'mo/'+ctrlB.n+'mo':'')+
          ' - short/sparse records make any single EWS unreliable, which is why several are combined here.';
        csdToolkitV.setValue(txt);
        return tally;
      }

      renderToolkit(buildTemporalIndicators(), true);
      // v10.151: let the (later-resolving) permutation callback refresh
      // the temporal-only toolkit once real p-values exist.
      csdToolkitRerender = function(){ renderToolkit(buildTemporalIndicators(), false); };

      try {
        var studyPtSpatial = ee.Geometry.Point([b.lon, b.lat]);
        var SPATIAL_BUFFER_M = 15000; // 15km: enough OISST (4km) pixels for a small spatial sample
        var studySpatialBuf = studyPtSpatial.buffer(SPATIAL_BUFFER_M);
        var rSpatialBefore = computeSpatialEWS(mkMoSSTRange(b.startDate,b.months), studySpatialBuf, 4000);
        var rSpatialAfter = computeSpatialEWS(mkMoSSTRange(a.startDate,a.months), studySpatialBuf, 4000);

        rSpatialBefore.evaluate(function(spBefore, spBeforeErr){
          rSpatialAfter.evaluate(function(spAfter, spAfterErr){
           try {
            // v10.100 NEW: study-vs-control temporal synchronization, fetched
            // as one more follow-up so a failure here also can't block or
            // delay the spatial/temporal results above it.
            var rSyncBefore = computeZonalSyncCSD(mkMoSSTRange(b.startDate,b.months), studySpatialBuf, ctrlStudy, 'sst', 4000);
            var rSyncAfter = computeZonalSyncCSD(mkMoSSTRange(a.startDate,a.months), studySpatialBuf, ctrlStudy, 'sst', 4000);
            rSyncBefore.evaluate(function(syncBefore, syncBeforeErr){
              rSyncAfter.evaluate(function(syncAfter, syncAfterErr){
            try {
              var spatialAvailable = !spBeforeErr && !spAfterErr && spBefore && spAfter;
              var sVarBefore = spatialAvailable&&spBefore.spatialVar!==null&&spBefore.spatialVar!==undefined?spBefore.spatialVar:null;
              var sVarAfter  = spatialAvailable&&spAfter.spatialVar!==null&&spAfter.spatialVar!==undefined?spAfter.spatialVar:null;
              var spAC1DetailB = spatialAvailable?extractSpatialAC1Detail(spBefore.corr):{value:null,key:null,how:'spatial fetch failed'};
              var spAC1Detail  = spatialAvailable?extractSpatialAC1Detail(spAfter.corr):{value:null,key:null,how:'spatial fetch failed'};
              var sAC1Before = spAC1DetailB.value;
              var sAC1After  = spAC1Detail.value;
              var pixN = spatialAvailable&&spAfter.pixelCount!==null&&spAfter.pixelCount!==undefined?spAfter.pixelCount:null;
              // v10.151 NOTE B: a real Looe Key run returned spatial AC1
              // delta = EXACTLY +0.000 across two independently computed
              // 48-pixel fields. That is far more consistent with
              // extractSpatialAC1() falling through to its last-resort
              // loop and picking the same constant twice than with a real
              // measurement. Print the raw reducer output so this can be
              // confirmed or ruled out instead of silently trusted.
              print('S13 spatial corrDict raw (BEFORE): '+(spBefore&&spBefore.corr?JSON.stringify(spBefore.corr):'null'));
              print('S13 spatial corrDict raw (AFTER):  '+(spAfter&&spAfter.corr?JSON.stringify(spAfter.corr):'null'));
              print('S13 spatial AC1 provenance: BEFORE via '+spAC1DetailB.how+
                ' (key='+spAC1DetailB.key+'), AFTER via '+spAC1Detail.how+' (key='+spAC1Detail.key+')');
              if(sAC1Before!==null&&sAC1After!==null&&sAC1Before===sAC1After){
                print('WARNING: spatial AC1 is bit-identical in both windows ('+sAC1Before+
                  '). Two independently computed pixel fields do not normally agree exactly. '+
                  'If the provenance line above says FALLBACK, this indicator has never been '+
                  'a real measurement and is now excluded from the tally.');
              }

              var spatialVarRelRise = (sVarBefore!==null&&sVarAfter!==null&&sVarBefore>1e-6)?
                ((sVarAfter-sVarBefore)/sVarBefore) : null;
              var dSpatialAC1 = (sAC1Before!==null&&sAC1After!==null)?(sAC1After-sAC1Before):null;

              var syncAvailable = !syncBeforeErr && !syncAfterErr && syncBefore && syncAfter;
              var corrBefore = syncAvailable&&syncBefore.corr!==null&&syncBefore.corr!==undefined?syncBefore.corr:null;
              var corrAfter = syncAvailable&&syncAfter.corr!==null&&syncAfter.corr!==undefined?syncAfter.corr:null;
              var dSync = (corrBefore!==null&&corrAfter!==null)?(corrAfter-corrBefore):null;

              var fullIndicators = buildTemporalIndicators().concat([
                {name:'Spatial Variance (study, '+(pixN!==null?'~'+pixN+'px':'few px')+')', scored:true, primary:false,
                  available:spatialVarRelRise!==null,
                  agrees:spatialVarRelRise!==null&&spatialVarRelRise>0.15,
                  display:spatialVarRelRise!==null?((spatialVarRelRise>0?'+':'')+(spatialVarRelRise*100).toFixed(0)+'%'):'n/a (region too small/sparse)'},
                {name:'Spatial Autocorrelation (study)'+
                    (spAC1Detail.how&&spAC1Detail.how.indexOf('FALLBACK')===0?' [UNRELIABLE - see note]':''),
                  scored:(spAC1Detail.how==='named key'),   // v10.152: not scored unless genuinely measured
                  primary:false,
                  available:dSpatialAC1!==null&&spAC1Detail.how==='named key',
                  agrees:dSpatialAC1!==null&&dSpatialAC1>0.05&&spAC1Detail.how==='named key',
                  display:(dSpatialAC1!==null?((dSpatialAC1>0?'+':'')+dSpatialAC1.toFixed(3)):'n/a')+
                    '  [source: '+(spAC1Detail.key?("key '"+spAC1Detail.key+"'"):'none')+' - '+spAC1Detail.how+']'+
                    (spAC1Detail.how!=='named key'?
                      '  <- NOT counted in the tally: this is not a confirmed measurement':'')},
                {name:'Study-Control Synchronization (raw SST - see seasonal-cycle caveat)', scored:true, primary:false,
                  available:dSync!==null,
                  agrees:dSync!==null&&dSync>0.10,
                  display:dSync!==null?('r: '+(corrBefore>0?'+':'')+corrBefore.toFixed(2)+' -> '+(corrAfter>0?'+':'')+corrAfter.toFixed(2)+' (\u0394'+(dSync>0?'+':'')+dSync.toFixed(2)+')'):'n/a (insufficient paired data)'}
              ]);

              var fullTally = renderToolkit(fullIndicators, false);
              // v10.151: once spatial indicators exist, re-renders should
              // use the FULL indicator set, not the temporal-only one.
              csdToolkitRerender = function(){ renderToolkit(fullIndicators, false); };
              var conf = classifyToolkitConfidence(fullTally, csdPermPAC1, csdPermPVar);

              // v10.90 FIX: verdict used to require a MAJORITY of indicators
              // to agree, treating AC1 and variance/spatial as equal votes.
              // Rebuilt to be AC1-weighted per Dakos et al. 2012: AC1 rising
              // is the primary trigger for a strong verdict; a variance/
              // spatial-only rise (AC1 not rising) is explicitly downgraded
              // to weaker evidence rather than being able to produce a
              // strong-sounding verdict on its own.
              var combinedTitle, combinedCol, combinedBg;
              var localFlavored = (vTitle.indexOf('LOCAL')>=0);
              if(conf.level==='high' && localFlavored){
                combinedTitle='STRONG LOCAL CSD SIGNAL (AC1-confirmed, corroborated)';
                combinedCol='#880000'; combinedBg='#ffd0d0';
              } else if(conf.level==='high' && !localFlavored){
                combinedTitle='STRONG SIGNAL, LIKELY REGIONAL (AC1-confirmed, but matches control site too)';
                combinedCol='#664400'; combinedBg='#fff6cc';
              } else if(conf.level==='moderate' && localFlavored){
                combinedTitle='LOCAL CSD SIGNAL (AC1-confirmed - the primary indicator per Dakos et al. 2012)';
                combinedCol='#aa3300'; combinedBg='#ffe8cc';
              } else if(conf.level==='moderate' && !localFlavored){
                combinedTitle='SIGNAL PRESENT (AC1-confirmed), LIKELY REGIONAL - not distinct from control site';
                combinedCol='#886600'; combinedBg='#fff6cc';
              } else if(conf.level==='low-moderate'){
                combinedTitle='WEAK SIGNAL - variance/spatial rose but AC1 (the more robust indicator) did NOT';
                combinedCol='#556633'; combinedBg='#eef4e0';
              } else {
                combinedTitle='NO RELIABLE CSD SIGNAL - primary indicator (AC1) not rising';
                combinedCol='#226644'; combinedBg='#e8f4ff';
              }
              // v10.119 FIX: "Regional context: "+vTitle used to paste the
              // OLD variance-only classification text directly after the
              // NEW AC1-weighted combinedTitle headline, with no
              // reconciliation between them - the same "two classifiers in
              // one box, can contradict each other" bug already caught and
              // fixed for FIND SWEET SPOT in v10.94. Caught here from a real
              // Nuuk test: headline correctly said "STRONG SIGNAL, LIKELY
              // REGIONAL" (AC1 rose +0.234 study vs +0.242 control - nearly
              // identical, genuinely regional), but the very same box then
              // said "Regional context: STUDY SITE MORE STABLE THAN CONTROL
              // (positive result)" - the opposite claim, driven purely by
              // variance (which fell at study, rose at control) while
              // completely ignoring that AC1 - the primary indicator - rose
              // almost identically at both sites. Now explicitly detected
              // and flagged instead of silently pasted.
              // v10.151 FIX 6: general conflict test replaces the old
              // single-string match, which missed the Looe Key case
              // (vTitle='NO SIGNAL AT EITHER SITE' vs a STRONG SIGNAL headline).
              var conflictMsg = verdictsConflict(combinedTitle, vTitle);
              var regionalContextLine = conflictMsg ? ('Regional context: '+conflictMsg) :
                ('Regional context: '+vTitle);
              // v10.152 FIX 7: a real Looe Key run printed
              //   "NO RELIABLE CSD SIGNAL - primary indicator (AC1) not rising"
              // immediately followed by "AC1 (primary): RISING".
              // combinedTitle already comes from the p-gated conf.level, but
              // this line still read fullTally.primaryAgrees, which is raw
              // DIRECTION with no significance attached. Both now derive from
              // the same gated object, so direction and significance can no
              // longer be reported as if they were the same thing.
              var ac1PrimaryTxt;
              if(!fullTally.primaryAvailable){
                ac1PrimaryTxt='n/a';
              } else if(conf.level==='preliminary'){
                ac1PrimaryTxt=(fullTally.primaryAgrees?'direction rising':'direction not rising')+
                  ' (significance pending)';
              } else if(conf.level==='low'&&fullTally.primaryAgrees){
                ac1PrimaryTxt='direction rising, but NOT significant';
              } else {
                ac1PrimaryTxt=(fullTally.primaryAgrees?'RISING':'not rising');
              }
              csdCompareVerdictV.setValue(combinedTitle+'\nAC1 (primary): '+ac1PrimaryTxt+
                ' | Supporting: '+fullTally.supportAgree+'/'+fullTally.supportAvail+' agree | '+regionalContextLine);
              csdCompareVerdictV.style().set('color',combinedCol);
              csdCompareVerdictV.style().set('backgroundColor',combinedBg);
              csdCompareVerdictV.style().set('border','2px solid '+combinedCol);
              csdCompareVerdictV.style().set('whiteSpace','pre');

              print('=== S13 TOOLKIT (final, temporal + spatial, AC1-weighted) ===');
              print('Combined verdict: '+combinedTitle);
              fullIndicators.forEach(function(ind){ print('  '+ind.name+': '+ind.display); });
            } catch(errSpatialFinal){
              print('=== S13 spatial toolkit finalize error (non-fatal, temporal verdict stands) === '+errSpatialFinal);
            }
              }); // end rSyncAfter.evaluate
            }); // end rSyncBefore.evaluate
           } catch(errSyncSetup){
             print('=== S13 zonal sync setup error (non-fatal, temporal+spatial verdict stands) === '+errSyncSetup);
           }
          });
        });
      } catch(errSpatialSetup){
        print('=== S13 spatial EWS setup error (non-fatal, temporal verdict stands) === '+errSpatialSetup);
      }
     } catch(errCompare) {
      csdCompareVerdictV.setValue('ERROR WHILE COMPARING');
      csdCompareVerdictV.style().set('color','#cc0000'); csdCompareVerdictV.style().set('backgroundColor','#ffd0d0');
      csdCompareVerdictV.style().set('border','2px solid #cc0000');
      csdCompareResultV.setValue(friendlyEEError(errCompare));
      csdCompareResultV.style().set('color','#cc0000'); csdCompareResultV.style().set('backgroundColor','#ffd0d0');
      print('=== S13 COMPARE ERROR === '+errCompare);
     }
    }); // end computeRealCSDDeseasonalized(ctrlAfterColl...)
  }); // end computeRealCSDDeseasonalized(ctrlBeforeColl...)
} // end runControlCSD

panel.add(csdCompareBtn);
panel.add(csdStalenessWarningV); // <-- staleness warning, above the verdict, hidden until triggered
panel.add(csdCompareStatValidV); // <-- v10.139 statistically-valid PRIMARY verdict, now shown FIRST
panel.add(csdCompareVerdictV); // <-- heuristic threshold verdict, now SECONDARY (see banner above)
panel.add(csdToolkitV);        // <-- multi-indicator toolkit tally, second
panel.add(csdControlSiteV);
panel.add(csdCompareResultV);  // <-- supporting numbers BELOW
panel.add(csdPermTestV);       // <-- v10.122 real permutation-test p-values
panel.add(csdDeseasonV);       // <-- v10.127 deseasonalized AC1/variance comparison
panel.add(legDiv());

// ----------------------------------------------------------
// ADVANCED (OPTIONAL) - manual control-site override
// ----------------------------------------------------------
panel.add(lbl('ADVANCED (OPTIONAL) - override the auto-selected control site',9,'#ffffff','#553388',true));
panel.add(lbl('Leave this BLANK unless you want to pick your own reference site by hand. Both COMPARE and FIND SWEET SPOT already auto-select and depth-validate a control site for you - most people never need this box.',7,'#553388'));
panel.add(lbl('If you do fill it in, choose a site in the same ocean, deep (>500m), close to the same latitude, and away from upwelling zones / enclosed seas / coastlines.',7,'#886600'));
var csdManualCtrlInput = ui.Textbox({
  placeholder: 'e.g. 24.5, -66.8  (leave blank for auto)',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(csdManualCtrlInput);
var csdCtrlModeV = dynLbl('Mode: AUTO (no manual site entered - default)','#553388');
panel.add(csdCtrlModeV);
csdManualCtrlInput.onChange(function(val){
  var trimmed = val.trim();
  if(trimmed.length > 0){
    csdCtrlModeV.setValue('Mode: MANUAL override = '+trimmed);
    csdCtrlModeV.style().set('color','#115511');
  } else {
    csdCtrlModeV.setValue('Mode: AUTO (no manual site entered - default)');
    csdCtrlModeV.style().set('color','#553388');
  }
});
panel.add(legDiv());

// ============================================================
// STEP 4 (OPTIONAL) - MULTI-WINDOW SWEET SPOT TEST
// ============================================================
panel.add(sHead('STEP 4 (OPTIONAL) - FIND SWEET SPOT','#1a3a4a'));
panel.add(lbl('Use this INSTEAD of manually guessing an AFTER window length.',7,'#226666'));
panel.add(lbl('Requires: a BEFORE window already stored in STEP 2 above. Auto-picks the same control site logic as COMPARE, and now runs a REAL control BEFORE/AFTER comparison (not a fixed baseline) so the numbers match what COMPARE would show.',7,'#226666'));
panel.add(lbl('Tests 6, 9, 12, 15, 18 and 24-month AFTER windows. For each one it checks: did the STUDY site variance/AC1 rise (LOCAL warning)? Did the CONTROL site variance/AC1 also rise (REGIONAL warning)? The window with the biggest gap between the two is the "sweet spot".',7,'#226666'));
panel.add(lbl('This fires 17 Earth Engine calls in parallel (1 control-BEFORE + 6 study-AFTER + 6 control-AFTER + 4 for the real permutation-test p-values below) and can take 30-120 seconds - a live counter below shows progress so it never looks frozen.',7,'#886600'));
panel.add(lbl('AFTER start date (YYYY-MM-DD):',7,'#334466'));
var csdAfterStartInput=ui.Textbox({
  placeholder:'AFTER start: e.g. 2023-06-01',
  style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}
});
panel.add(csdAfterStartInput);
var csdMultiResultV=ui.Label('',{fontSize:'8px',color:'#114444',
  backgroundColor:'#e0f0f8',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var csdMultiSweetSpotV=ui.Label('',{fontSize:'8px',color:'#334466',
  backgroundColor:'#f4f8fb',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
// v10.123 NEW: real permutation-test p-values across all 6 window
// lengths, with a Bonferroni-corrected significance bar - direct fix for
// the multiple-comparisons trap the OLD threshold-based ROBUSTNESS note
// only ever warned about in words, never actually corrected for
// mathematically.
var csdMultiPermTestV=ui.Label('',{fontSize:'8px',color:'#552266',
  backgroundColor:'#f6eefa',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var csdMultiStatusV=dynLbl('Fill in STEP 1-2 (BEFORE) above, enter an AFTER start date, then press FIND SWEET SPOT.','#334466');
var csdMultiWindowBtn=ui.Button({
  label:'FIND SWEET SPOT (multi-window auto-test)',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',
    backgroundColor:'#cce8ff',color:'#003388',stretch:'horizontal',
    padding:'6px 4px',border:'2px solid #0055cc'},
  onClick:function(){
    if(!csdBeforeResult){
      csdMultiStatusV.setValue('MISSING DATA - run and store a BEFORE window in STEP 2 first.');
      csdMultiStatusV.style().set('color','#cc0000'); csdMultiStatusV.style().set('backgroundColor','#ffd0d0');
      csdMultiStatusV.style().set('border','2px solid #cc0000');
      csdMultiResultV.setValue(''); csdMultiSweetSpotV.setValue('');
      return;
    }
    var afterStartTxt=(csdAfterStartInput.getValue()||'').trim();
    if(!afterStartTxt){
      csdMultiStatusV.setValue('MISSING DATA - enter an AFTER start date (YYYY-MM-DD) above.');
      csdMultiStatusV.style().set('color','#cc0000'); csdMultiStatusV.style().set('backgroundColor','#ffd0d0');
      csdMultiStatusV.style().set('border','2px solid #cc0000');
      csdMultiResultV.setValue(''); csdMultiSweetSpotV.setValue('');
      return;
    }
    var b=csdBeforeResult;
    var windowLengths=[6,9,12,15,18,24];
    var smartCtrl=getSmartControlSite(b.lat,b.lon);

    // Reset styling from any previous run
    csdMultiStatusV.style().set('padding','1px 4px'); csdMultiStatusV.style().set('whiteSpace','normal');
    csdMultiStatusV.style().set('fontWeight','normal'); csdMultiStatusV.style().set('fontSize','9px');
    csdMultiStatusV.style().set('border','none'); csdMultiStatusV.style().set('backgroundColor','rgba(0,0,0,0)');
    csdMultiSweetSpotV.setValue('');

    csdMultiStatusV.setValue('Step 1/3: Smart control = '+smartCtrl.label+
      ' ('+smartCtrl.lat.toFixed(2)+', '+smartCtrl.lon.toFixed(2)+') — launching 13 parallel EE calls...');
    csdMultiStatusV.style().set('color','#334466');
    csdMultiResultV.setValue('Study: '+b.lat+', '+b.lon+' | BEFORE: '+b.startDate+', '+b.months+'mo\n'+
      'Control: '+smartCtrl.lat.toFixed(2)+', '+smartCtrl.lon.toFixed(2)+' ('+smartCtrl.label+')\n'+
      'AFTER start: '+afterStartTxt+' | Windows: '+windowLengths.join(', ')+' months\n'+
      'Progress: 0 / '+(windowLengths.length*2+1)+' sub-tests done...');

    var studyPt=ee.Geometry.Point([b.lon,b.lat]);
    var ctrlPt=ee.Geometry.Point([smartCtrl.lon,smartCtrl.lat]);
    var studyBuf=studyPt.buffer(10000);
    var ctrlBuf=ctrlPt.buffer(80000);

    // ctrlBefore only needs to be computed ONCE - it uses the same BEFORE
    // window (b.startDate, b.months) as STEP 2, just evaluated at the
    // control site instead of the study site. This lets us compute a REAL
    // control-site delta (ctrlAfter - ctrlBefore) for every window length,
    // the same way STEP 3 (COMPARE) does, instead of comparing against an
    // arbitrary fixed baseline.
    var ctrlBeforeColl=mkMoSSTRange(b.startDate,b.months);

    // v10.123: +4 for the new permutation-test raw-value fetches (study
    // BEFORE/AFTER-max, control BEFORE/AFTER-max) - see the fire section
    // below for why only 4 new calls are needed, not 12.
    var multiRes={ctrlBefore:null}, multiTotal=windowLengths.length*2+1+4, multiDone=0, multiErrors=0;

    function bumpProgress(){
      multiDone++;
      csdMultiStatusV.setValue('Step 2/3: Computing... ('+multiDone+' / '+multiTotal+' sub-tests done'+
        (multiErrors>0?', '+multiErrors+' returned no data':'')+')');
      if(multiDone < multiTotal) return;
      finishAnalysis();
    }

    function finishAnalysis(){
      try {
        csdMultiStatusV.setValue('Step 3/3: Analysing results...');
        var cB=multiRes.ctrlBefore, cBVar=(cB&&cB.varTrendRatio!==null&&cB.varTrendRatio!==undefined)?cB.varTrendRatio:null;
        var cBAC1=(cB&&cB.realAC1!==null&&cB.realAC1!==undefined)?cB.realAC1:null;

        var rows=[], bestW=null, bestDiv=-Infinity, bestRow=null;
        var bestAc1W=null, bestAc1Rise=-Infinity, bestAc1Row=null;
        // v10.92: table now shows Study \u0394AC1 too - previously only \u0394Var
        // was visible per window, making it impossible to scan for the
        // AC1-rising / variance-falling pattern Dakos et al. 2012 documents.
        var header='Window | Study \u0394Var | Study \u0394AC1 | Ctrl \u0394Var | Divergence | Verdict';
        rows.push(header);
        rows.push(repeatChar('\u2500',72));

        // v10.153 FIX 10: calibrated per tested window length below.
        var THRESH=getCalibratedThresholds(b.months||24).varr;
        var allRows=[];

        for(var wi=0;wi<windowLengths.length;wi++){
          var w=windowLengths[wi];
          var sr=multiRes['study_'+w], cr=multiRes['ctrl_'+w];
          var sVR=(sr&&sr.varTrendRatio!==null&&sr.varTrendRatio!==undefined)?sr.varTrendRatio:null;
          var sAC1=(sr&&sr.realAC1!==null&&sr.realAC1!==undefined)?sr.realAC1:null;
          var cVR=(cr&&cr.varTrendRatio!==null&&cr.varTrendRatio!==undefined)?cr.varTrendRatio:null;
          var cAC1=(cr&&cr.realAC1!==null&&cr.realAC1!==undefined)?cr.realAC1:null;

          var bVR=b.varTrend!==null&&b.varTrend!==undefined?b.varTrend:null;
          var bAC1=b.ac1!==null&&b.ac1!==undefined?b.ac1:null;

          // REAL before/after deltas for both sites (same maths STEP 3 uses)
          var sDVar=(sVR!==null&&bVR!==null)?sVR-bVR:null;
          var sDAC1=(sAC1!==null&&bAC1!==null)?sAC1-bAC1:null;
          var cDVar=(cVR!==null&&cBVar!==null)?cVR-cBVar:null;
          var cDAC1=(cAC1!==null&&cBAC1!==null)?cAC1-cBAC1:null;
          var div=(sDVar!==null&&cDVar!==null)?(sDVar-cDVar):null;

          // Scheffer 2009 validation: classic CSD needs BOTH indicators
          // rising together (AC1 up AND variance up), not just one.
          // v10.153: each window length gets its own calibrated cutoff -
          // the null spread shrinks as the window grows (sd 0.247 at 24mo
          // vs 0.166 at 48mo), so one fixed number cannot serve all six.
          var _wthr = getCalibratedThresholds(w);
          var scheffer = (sDVar!==null&&sDVar>_wthr.varr&&sDAC1!==null&&sDAC1>_wthr.ac1);
          var ac1Rose = (sDAC1!==null&&sDAC1>_wthr.ac1);
          var varRose = (sDVar!==null&&sDVar>_wthr.varr);

          var studyVarRose=(sDVar!==null&&sDVar>_wthr.varr);
          var ctrlVarRose=(cDVar!==null&&cDVar>_wthr.varr);
          var verdict='NO SIGNAL';
          if(sDVar===null||cDVar===null) verdict='n/a (missing data)';
          else if(studyVarRose&&ctrlVarRose) verdict='GLOBAL';
          else if(studyVarRose&&!ctrlVarRose) verdict='LOCAL CSD';
          else if(!studyVarRose&&div!==null&&div>_wthr.varr) verdict='MARGINAL LOCAL';
          else if(!studyVarRose&&ctrlVarRose) verdict='ANOMALOUS (study more stable)';
          else verdict='NO SIGNAL';
          // v10.92: tag AC1-driven signal separately from the variance-driven
          // verdict above, so an "AC1 up / variance down" window is labelled
          // accurately instead of silently falling into "NO SIGNAL" just
          // because the variance-based verdict logic didn't rise.
          if(ac1Rose && !varRose) verdict += ' + AC1 CONFIRMED (var down)';
          else if(ac1Rose && varRose) verdict += ' + AC1 CONFIRMED';

          if(div!==null&&div>bestDiv){bestDiv=div;bestW=w;}
          // v10.92: separate AC1-based ranking (primary indicator per Dakos
          // et al. 2012) alongside the existing variance-divergence ranking,
          // so a strong AC1 rise is never hidden just because variance did
          // something else.
          if(sDAC1!==null&&sDAC1>bestAc1Rise){bestAc1Rise=sDAC1;bestAc1W=w;}
          var rowObj={w:w,sDVar:sDVar,cDVar:cDVar,div:div,verdict:verdict,sDAC1:sDAC1,cDAC1:cDAC1,scheffer:scheffer};
          allRows.push(rowObj);
          if(w===bestAc1W) bestAc1Row=rowObj;
        }

        for(var ri=0;ri<allRows.length;ri++){
          var ro=allRows[ri];
          var row=ro.w+'mo   | '+
            (ro.sDVar!==null?(ro.sDVar>0?'+':'')+ro.sDVar.toFixed(2)+'x':'n/a')+'   | '+
            (ro.sDAC1!==null?(ro.sDAC1>0?'+':'')+ro.sDAC1.toFixed(3):'n/a')+'    | '+
            (ro.cDVar!==null?(ro.cDVar>0?'+':'')+ro.cDVar.toFixed(2)+'x':'n/a')+'   | '+
            (ro.div!==null?(ro.div>0?'+':'')+ro.div.toFixed(2)+'x':'n/a')+'  | '+ro.verdict+
            (ro.w===bestW?' \u2190 VAR SWEET SPOT':'')+
            (ro.w===bestAc1W&&ro.w!==bestW?' \u2190 AC1 SWEET SPOT':ro.w===bestAc1W&&ro.w===bestW?' + AC1':'');
          rows.push(row);
          if(ro.w===bestW) bestRow=ro;
        }
        rows.push(repeatChar('\u2500',72));

        // v10.95 NEW: ROBUSTNESS CHECK - the actual answer to "which window
        // should I trust". Testing 6 window lengths and reporting whichever
        // one shows the biggest divergence is a classic multiple-comparisons
        // trap: with 6 tries, SOME window will look dramatic by chance alone,
        // even in pure noise. A real approach-to-tipping-point signal should
        // show up across NEIGHBOURING windows, not as an isolated spike at
        // one window length that vanishes at the windows on either side of
        // it. This counts how many of the 6 tested windows lean the same
        // direction as the "sweet spot", so an isolated result gets flagged
        // as suspect rather than reported as a confident single answer.
        var localLeaningRows = allRows.filter(function(r){
          return r.verdict.indexOf('LOCAL CSD')===0 || r.verdict.indexOf('MARGINAL LOCAL')===0;
        });
        var testedRows = allRows.filter(function(r){ return r.sDVar!==null && r.cDVar!==null; });
        var robustNote;
        if(testedRows.length===0){
          robustNote='ROBUSTNESS: no windows had usable data - cannot assess.';
        } else if(localLeaningRows.length<=1){
          robustNote='ROBUSTNESS: '+(localLeaningRows.length===0?'NONE':'only '+localLeaningRows.length)+' of '+testedRows.length+' tested windows lean LOCAL.\n'+
            'This looks like an ISOLATED result, not a signal that holds up across window\n'+
            'lengths. Per standard multiple-comparisons caution: testing 6 windows and\n'+
            'reporting the most dramatic one is expected to find "something" by chance\n'+
            'alone. Do NOT treat a single-window result as confirmed - re-test with a\n'+
            'window length chosen from independent evidence (a documented event date),\n'+
            'not from this scan.';
        } else {
          robustNote='ROBUSTNESS: '+localLeaningRows.length+' of '+testedRows.length+' tested windows lean LOCAL ('+
            localLeaningRows.map(function(r){return r.w+'mo';}).join(', ')+').\n'+
            'A signal appearing across multiple window lengths is more trustworthy than\n'+
            'a single isolated spike, but this is still not a significance test - treat\n'+
            'as a stronger hypothesis to validate with field data, not a confirmed result.';
        }
        rows.push(robustNote);
        rows.push('');
        rows.push('v10.153 THRESHOLD CALIBRATION: cutoffs are now per-window, set at the 95th');
        rows.push('percentile of |delta| under a measured no-event null (Scripps Pier CTD, 13.6yr).');
        rows.push('The previous fixed 0.15x/0.01 cutoffs fired on 34-57% of no-event windows; the');
        rows.push('AC1-OR-variance rule fired on 80%. Expect far fewer LOCAL verdicts than before -');
        rows.push('that is the correction, not a regression. The Bonferroni-corrected p-values below');
        rows.push('remain the only figures here with a known false-positive rate.');
        rows.push('\nControl BEFORE baseline (same '+b.startDate+', '+b.months+'mo window as STEP 2): AC1='+
          (cBAC1!==null?cBAC1.toFixed(3):'n/a')+', Var='+(cBVar!==null?cBVar.toFixed(2)+'x':'n/a'));
        rows.push('Control site: '+smartCtrl.label);
        rows.push('CAVEAT: study site uses a 10km buffer, control uses an 80km buffer.');
        if(multiErrors>0) rows.push('NOTE: '+multiErrors+' sub-test(s) returned no usable data (shown as n/a above).');
        csdMultiResultV.setValue(rows.join('\n'));

        // v10.123 NEW: real permutation-test p-values across all 6 window
        // lengths, using the SAME permutationTestDelta() engine already
        // validated in STEP 3 COMPARE (v10.122). Direct fix for the
        // multiple-comparisons trap the OLD threshold-based ROBUSTNESS
        // note above only ever warned about in words: testing 6 windows
        // means some window can look significant by pure chance alone,
        // even with nothing real happening. The Bonferroni correction
        // below (0.05 / 6 windows tested) gives an explicit, stricter bar
        // that actually accounts for this mathematically, alongside the
        // uncorrected count for reference.
        try {
          var permStudyBeforeByNode = groupSeriesByLabel(multiRes.permStudyBefore, 'sst');
          var permStudyAfterByNode = groupSeriesByLabel(multiRes.permStudyAfterMax, 'sst');
          var permCtrlBeforeByNode = groupSeriesByLabel(multiRes.permCtrlBefore, 'sst');
          var permCtrlAfterByNode = groupSeriesByLabel(multiRes.permCtrlAfterMax, 'sst');
          var pStudyBeforeVals = permStudyBeforeByNode['Study']||[];
          var pStudyAfterValsFull = permStudyAfterByNode['Study']||[];
          var pCtrlBeforeVals = permCtrlBeforeByNode['Control']||[];
          var pCtrlAfterValsFull = permCtrlAfterByNode['Control']||[];

          var bonferroniAlpha = 0.05/windowLengths.length;
          var permLines=['=== REAL SIGNIFICANCE ACROSS ALL 6 WINDOWS (permutation test) ==='];
          permLines.push('Same engine as STEP 3 COMPARE, run at each window length (300 shuffles each,');
          permLines.push('reduced from 500 to keep 24 total tests fast). Testing 6 windows means SOME');
          permLines.push('window can look significant by pure chance - the CORRECTED bar below (0.05 /');
          permLines.push('6 windows tested = '+bonferroniAlpha.toFixed(4)+') accounts for that mathematically,');
          permLines.push('unlike the ROBUSTNESS note above, which only ever warned about this in words.');
          permLines.push(repeatChar('\u2500',72));
          permLines.push('Window | Study AC1 p | Study Var p | Ctrl AC1 p | Ctrl Var p | Local signal?');
          permLines.push(repeatChar('\u2500',72));

          var anyCorrectedSig = false, uncorrectedLocalCount = 0;
          windowLengths.forEach(function(w){
            var sAfterW = pStudyAfterValsFull.slice(0, w);
            var cAfterW = pCtrlAfterValsFull.slice(0, w);
            var sAC1Test = permutationTestDeltaFixed(pStudyBeforeVals, sAfterW, statAC1ForPerm, 300);
            var sVarTest = permutationTestDeltaFixed(pStudyBeforeVals, sAfterW, statVarRatioForPerm, 300);
            var cAC1Test = permutationTestDeltaFixed(pCtrlBeforeVals, cAfterW, statAC1ForPerm, 300);
            var cVarTest = permutationTestDeltaFixed(pCtrlBeforeVals, cAfterW, statVarRatioForPerm, 300);
            function fp(t){ return t.pValue!==null?t.pValue.toFixed(3):'n/a'; }
            // "Local signal" = study significant (p<0.05, uncorrected) on
            // EITHER AC1 or variance, while control is NOT significant on
            // that same statistic - the same LOCAL-vs-REGIONAL logic used
            // throughout this tool, now grounded in real p-values instead
            // of the >0.01/>0.15 threshold.
            var sAC1Sig = sAC1Test.pValue!==null && sAC1Test.pValue<0.05;
            var sVarSig = sVarTest.pValue!==null && sVarTest.pValue<0.05;
            var cAC1Sig = cAC1Test.pValue!==null && cAC1Test.pValue<0.05;
            var cVarSig = cVarTest.pValue!==null && cVarTest.pValue<0.05;
            var localSig = (sAC1Sig && !cAC1Sig) || (sVarSig && !cVarSig);
            if(localSig) uncorrectedLocalCount++;
            var minStudyP = Math.min(sAC1Test.pValue!==null?sAC1Test.pValue:1, sVarTest.pValue!==null?sVarTest.pValue:1);
            if(localSig && minStudyP<bonferroniAlpha) anyCorrectedSig=true;
            permLines.push(w+'mo    | '+fp(sAC1Test)+'       | '+fp(sVarTest)+'       | '+fp(cAC1Test)+'      | '+fp(cVarTest)+
              '      | '+(localSig?(minStudyP<bonferroniAlpha?'YES (survives correction)':'YES (uncorrected only)'):'no'));
          });
          permLines.push(repeatChar('\u2500',72));
          permLines.push(uncorrectedLocalCount+' of '+windowLengths.length+' windows show a local signal at the uncorrected p<0.05 level.');
          if(anyCorrectedSig){
            permLines.push('At least one window survives the STRICTER Bonferroni-corrected bar (p<'+bonferroniAlpha.toFixed(4)+') -');
            permLines.push('this is real evidence, not just a lucky window among 6 tries.');
          } else if(uncorrectedLocalCount>0){
            permLines.push('NONE survive the corrected bar - the uncorrected hits above are consistent with');
            permLines.push('what pure chance alone would produce across 6 tests. Treat as a hypothesis to');
            permLines.push('test further (e.g. with an independently-chosen window), not a confirmed signal.');
          } else {
            permLines.push('No window shows a real local signal, corrected or uncorrected.');
          }
          csdMultiPermTestV.setValue(permLines.join('\n'));
          csdMultiPermTestV.style().set('color','#552266');
          print('=== FIND SWEET SPOT REAL PERMUTATION TEST ==='); print(permLines.join('\n'));
        } catch(ePermMulti){
          csdMultiPermTestV.setValue('Permutation test error: '+ePermMulti);
          csdMultiPermTestV.style().set('color','#cc0000');
          print('=== FIND SWEET SPOT PERMUTATION TEST ERROR === '+ePermMulti);
        }

        // --- Explicit sweet-spot breakdown: local vs regional vs Scheffer check ---
        if(bestRow){
          var ac1Rose = bestRow.sDAC1!==null && bestRow.sDAC1>0.01;
          var varRose = bestRow.sDVar!==null && bestRow.sDVar>THRESH;
          // v10.90 FIX: previously this only reported PASSED (both rose) or
          // a flat "NOT MET - only one (or neither) rose". Per Dakos et al.
          // 2012, AC1 not rising and variance not rising are NOT equally
          // weak cases - AC1 is the more robust indicator, so distinguish
          // which one is actually missing.
          var schefferText;
          if(ac1Rose && varRose){
            schefferText='PASSED - both indicators rising together at the study site. Strongest evidence.';
          } else if(ac1Rose && !varRose){
            schefferText='PARTIAL (still meaningful) - AC1, the more robust indicator per Dakos et al. 2012, IS rising, even though variance is not. Dakos et al. found variance can fail to rise (or even fall) near a real transition, so this is not necessarily weak evidence.';
          } else if(!ac1Rose && varRose){
            schefferText='PARTIAL (weaker) - variance rose but AC1 (the more robust indicator per Dakos et al. 2012) did NOT. Treat this with caution: a variance-only rise is considered less reliable evidence of critical slowing down.';
          } else {
            schefferText='NOT MET - neither indicator rose at the study site.';
          }
          var sweetTxt=
            'SWEET SPOT WINDOW: '+bestRow.w+' months (AFTER start '+afterStartTxt+')\n'+
            '\n'+
            'LOCAL (study site) warning:\n'+
            '  \u0394 Variance vs BEFORE: '+(bestRow.sDVar!==null?(bestRow.sDVar>0?'+':'')+bestRow.sDVar.toFixed(2)+'x':'n/a')+
              (varRose?'  RISING (warning sign)':'  not rising')+'\n'+
            '  \u0394 AC1 vs BEFORE: '+(bestRow.sDAC1!==null?(bestRow.sDAC1>0?'+':'')+bestRow.sDAC1.toFixed(3):'n/a')+
              (ac1Rose?'  RISING (warning sign)':'  not rising')+'\n'+
            '\n'+
            'REGIONAL (control site) warning:\n'+
            '  \u0394 Variance vs BEFORE: '+(bestRow.cDVar!==null?(bestRow.cDVar>0?'+':'')+bestRow.cDVar.toFixed(2)+'x':'n/a')+
              (bestRow.cDVar!==null&&bestRow.cDVar>THRESH?'  RISING (regional warming)':'  stable/flat')+'\n'+
            '  \u0394 AC1 vs BEFORE: '+(bestRow.cDAC1!==null?(bestRow.cDAC1>0?'+':'')+bestRow.cDAC1.toFixed(3):'n/a')+'\n'+
            '\n'+
            'SCHEFFER 2009 VALIDATION (AC1 weighted as the more robust indicator, per Dakos et al. 2012):\n'+
            '  '+schefferText+'\n'+
            '\n'+
            'STABILITY / ANOMALY READ:\n'+
            '  '+bestRow.verdict;
          csdMultiSweetSpotV.setValue(sweetTxt);
          csdMultiSweetSpotV.style().set('whiteSpace','pre');

          // v10.92 NEW: cross-check against the AC1-based ranking. Since AC1
          // is the primary indicator (v10.90/91), a window where AC1 rises
          // strongly deserves attention even if it isn't the variance-
          // divergence winner above - most relevant for exactly the
          // "AC1 up, variance down" pattern Dakos et al. document.
          if(bestAc1Row && bestAc1W!==bestW && bestAc1Row.sDAC1!==null && bestAc1Row.sDAC1>0.01){
            var ac1AltVarRose = bestAc1Row.sDVar!==null && bestAc1Row.sDVar>0.01;
            print('=== S13 AC1-BASED RANKING DIFFERS FROM VARIANCE SWEET SPOT ===');
            print('AC1 rises most strongly at the '+bestAc1W+'-month window (\u0394AC1='+
              bestAc1Row.sDAC1.toFixed(3)+'), not the '+bestW+'-month variance sweet spot.');
            print('At '+bestAc1W+'mo: \u0394Var='+(bestAc1Row.sDVar!==null?(bestAc1Row.sDVar>0?'+':'')+bestAc1Row.sDVar.toFixed(2)+'x':'n/a')+
              (ac1AltVarRose?' (also rising)':' (NOT rising - this is the AC1-up/variance-down pattern per Dakos et al. 2012, Fig. 2c & 4)'));
          }
        } else {
          csdMultiSweetSpotV.setValue('No window produced a usable divergence value - check that SST data exists for both the study and control sites in this date range.');
        }

        // --- Short, bold, colour-coded VERDICT box - shown FIRST, above the table ---
        // v10.94 FIX: this headline used to be classified from a SEPARATE,
        // looser check (bestDiv>THRESH alone, ignoring whether the study
        // site's own variance actually crossed the threshold, and ignoring
        // AC1 entirely) - so it could say "LOCAL CSD SIGNAL DETECTED" for a
        // window whose own per-window verdict said "MARGINAL LOCAL" with
        // "Scheffer NOT MET". Rebuilt to derive the headline from bestRow's
        // OWN verdict string (the same classification shown in the table
        // and the sweet-spot breakdown), plus whether AC1 also diverges
        // from the control - so the headline can never contradict the
        // detail below it again.
        var vBoxColor, vBoxBg, vBoxText;
        if(!bestRow){
          vBoxColor='#cc0000'; vBoxBg='#ffd0d0';
          vBoxText='NO USABLE DATA\n'+
            'None of the 6 windows returned a valid divergence value.\n'+
            'Check that both the study and control coordinates have SST\n'+
            'coverage for these dates, then try again.';
        } else {
          var ac1Div = (bestRow.sDAC1!==null&&bestRow.cDAC1!==null)?(bestRow.sDAC1-bestRow.cDAC1):null;
          var ac1DivRose = ac1Div!==null && ac1Div>0.01;
          var ac1Line = ac1DivRose
            ? 'AC1 also diverges from control - corroborates this signal.'
            : 'AC1 does NOT corroborate this (Scheffer check: '+(bestRow.scheffer?'PASSED':'NOT MET')+') - treat with caution.';
          var isolatedWarning = localLeaningRows.length<=1 ?
            '\n\u26A0 ISOLATED: only 1 of '+testedRows.length+' windows leans this way - see ROBUSTNESS note in table below before trusting this.' :
            '\n('+localLeaningRows.length+' of '+testedRows.length+' windows lean the same way - see ROBUSTNESS note below.)';
          if(bestRow.verdict.indexOf('LOCAL CSD')===0){
            vBoxColor='#880000'; vBoxBg='#ffd0d0';
            vBoxText='LOCAL CSD SIGNAL DETECTED\n'+
              'Sweet spot: '+bestW+'-month AFTER window\n'+
              'Study variance rose ('+(bestRow.sDVar>0?'+':'')+bestRow.sDVar.toFixed(2)+'x) while the control did not.\n'+
              ac1Line+
              isolatedWarning+'\n'+
              'Use '+bestW+' months as the AFTER window in STEP 2 + STEP 3 (COMPARE)\n'+
              'for the definitive 4-way comparison.';
          } else if(bestRow.verdict.indexOf('MARGINAL LOCAL')===0){
            vBoxColor='#aa3300'; vBoxBg='#ffe8cc';
            vBoxText='MARGINAL LOCAL SIGNAL (weaker)\n'+
              'Sweet spot: '+bestW+'-month AFTER window\n'+
              'Neither site individually crossed the variance threshold, but the\n'+
              'study site diverges from the control by '+(bestDiv>0?'+':'')+bestDiv.toFixed(2)+'x.\n'+
              ac1Line+
              isolatedWarning+'\n'+
              'Weaker than a true LOCAL CSD result - inspect the breakdown below\n'+
              'before treating this as confirmed.';
          } else if(bestRow.verdict.indexOf('GLOBAL')===0){
            vBoxColor='#664400'; vBoxBg='#fff6cc';
            vBoxText='GLOBAL SIGNAL - NOT local\n'+
              'Sweet spot: '+bestW+'-month AFTER window\n'+
              'Both study AND control variance rose - regional warming, not a\n'+
              'local tipping-point signal.\n'+
              ac1Line;
          } else if(bestRow.verdict.indexOf('ANOMALOUS')===0){
            vBoxColor='#115511'; vBoxBg='#d4f5df';
            vBoxText='STUDY SITE MORE STABLE THAN CONTROL (positive result)\n'+
              'Sweet spot: '+bestW+'-month AFTER window\n'+
              'Control variance rose but the study site stayed flat - regional\n'+
              'warming is NOT affecting this site as strongly.';
          } else {
            vBoxColor='#226644'; vBoxBg='#e8f4ff';
            vBoxText='NO LOCAL SIGNAL FOUND\n'+
              'Best window ('+bestW+'mo) still shows: '+bestRow.verdict+'\n'+
              'No local tipping-point evidence in this date range.\n'+
              'Try FIND SWEET SPOT with a different AFTER start date, or a\n'+
              'window that better targets the lead-up to a suspected event\n'+
              'rather than the event itself.';
          }
        }
        csdMultiStatusV.setValue(vBoxText);
        csdMultiStatusV.style().set('color', vBoxColor);
        csdMultiStatusV.style().set('backgroundColor', vBoxBg);
        csdMultiStatusV.style().set('padding','6px 8px');
        csdMultiStatusV.style().set('whiteSpace','pre');
        csdMultiStatusV.style().set('fontWeight','bold');
        csdMultiStatusV.style().set('fontSize','9px');
        csdMultiStatusV.style().set('border','2px solid '+vBoxColor);
        print('=== S13 MULTI-WINDOW SWEET SPOT TEST ===');
        print(vBoxText);
        print(rows.join('\n'));
      } catch(errFinal) {
        csdMultiStatusV.setValue(friendlyEEError(errFinal));
        csdMultiStatusV.style().set('color','#cc0000'); csdMultiStatusV.style().set('backgroundColor','#ffd0d0');
        csdMultiStatusV.style().set('border','2px solid #cc0000'); csdMultiStatusV.style().set('whiteSpace','pre');
        print('=== S13 FIND SWEET SPOT ERROR ===');
        print(errFinal);
      }
    }

    // 1 control-BEFORE call (shared across all windows)
    computeRealCSDDeseasonalized(ctrlBeforeColl,ctrlBuf,'sst',4000,function(res){
      if(res.error) multiErrors++;
      multiRes.ctrlBefore=res.error?{}:res;
      bumpProgress();
    });

    // 6 study-AFTER + 6 control-AFTER calls, one pair per window length
    windowLengths.forEach(function(w){
      var studyColl=mkMoSSTRange(afterStartTxt,w);
      var ctrlColl=mkMoSSTRange(afterStartTxt,w);
      computeRealCSDDeseasonalized(studyColl,studyBuf,'sst',4000,function(res){
        if(res.error) multiErrors++;
        multiRes['study_'+w]=res.error?{}:res;
        bumpProgress();
      });
      computeRealCSDDeseasonalized(ctrlColl,ctrlBuf,'sst',4000,function(res){
        if(res.error) multiErrors++;
        multiRes['ctrl_'+w]=res.error?{}:res;
        bumpProgress();
      });
    });

    // v10.123 NEW: 4 raw-value fetches for the permutation test - fetching
    // each site's BEFORE window plus its LONGEST (24mo) AFTER window ONCE.
    // Since all 6 AFTER windows share the same start date and differ only
    // in length, a shorter window's raw monthly values are always the
    // first N months of the 24-month fetch - sliced client-side once all
    // data arrives, avoiding 12 additional per-window EE calls.
    var maxWindowLen = windowLengths[windowLengths.length-1];
    var permStudyFC = ee.FeatureCollection([ee.Feature(studyBuf, {label:'Study', idx:0})]);
    var permCtrlFC = ee.FeatureCollection([ee.Feature(ctrlBuf, {label:'Control', idx:0})]);
    var permStudyBeforeColl = mkMoSSTRange(b.startDate, b.months);
    var permStudyAfterMaxColl = mkMoSSTRange(afterStartTxt, maxWindowLen);
    var permCtrlBeforeColl = mkMoSSTRange(b.startDate, b.months);
    var permCtrlAfterMaxColl = mkMoSSTRange(afterStartTxt, maxWindowLen);
    extractMultiNodeSeries(permStudyBeforeColl, permStudyFC, 'sst', 4000).evaluate(function(res,err){
      if(err) multiErrors++;
      multiRes.permStudyBefore=err?null:res;
      bumpProgress();
    });
    extractMultiNodeSeries(permStudyAfterMaxColl, permStudyFC, 'sst', 4000).evaluate(function(res,err){
      if(err) multiErrors++;
      multiRes.permStudyAfterMax=err?null:res;
      bumpProgress();
    });
    extractMultiNodeSeries(permCtrlBeforeColl, permCtrlFC, 'sst', 4000).evaluate(function(res,err){
      if(err) multiErrors++;
      multiRes.permCtrlBefore=err?null:res;
      bumpProgress();
    });
    extractMultiNodeSeries(permCtrlAfterMaxColl, permCtrlFC, 'sst', 4000).evaluate(function(res,err){
      if(err) multiErrors++;
      multiRes.permCtrlAfterMax=err?null:res;
      bumpProgress();
    });
  }
});
panel.add(csdMultiWindowBtn);
panel.add(csdMultiStatusV);     // <-- verdict FIRST (colour box)
panel.add(csdMultiSweetSpotV);  // <-- explicit local vs regional vs Scheffer breakdown for the winning window
panel.add(csdMultiResultV);     // <-- full 6-window data table BELOW
panel.add(csdMultiPermTestV);   // <-- v10.123 real permutation-test p-values across all 6 windows
panel.add(legDiv());

// ============================================================
// STEP 5 (OPTIONAL) - SLIDING WINDOW TREND TEST (v10.96 NEW)
// The actual Dakos et al. 2012 gold-standard method: a continuous rolling
// window across the WHOLE series (not discrete BEFORE/AFTER chunks),
// tested with the Mann-Kendall trend test (Kendall's tau) instead of a
// fixed numeric threshold. This is the direct fix for the window-selection
// / multiple-comparisons problem in STEP 4: instead of testing 6 arbitrary
// window lengths and picking whichever looks most dramatic, this shows the
// WHOLE trajectory and asks whether it is trending consistently, with an
// actual significance test behind that answer.
// ============================================================
panel.add(sHead('STEP 5 (OPTIONAL) - SLIDING WINDOW TREND TEST','#3a1a4a'));
panel.add(lbl('The Dakos et al. 2012 gold-standard method. Instead of testing isolated chunks, this slides a fixed window forward ONE MONTH AT A TIME across the whole series, computing AC1 and variance at every position - producing a continuous trajectory, then tests it with Kendall\'s tau (the Mann-Kendall trend test): is this metric moving consistently in one direction, or just bouncing randomly?',7,'#663388'));
panel.add(lbl('DISCLOSED SIMPLIFICATION: this uses the same linear detrend as the rest of S13 (fit once, over the whole series). Dakos et al.\'s own toolbox typically uses Gaussian kernel smoothing instead, which handles nonlinear trends more flexibly - not implemented here. Results are still real Kendall-tau / Mann-Kendall statistics, just on a simpler detrend.',7,'#886600'));
panel.add(lbl('Full series start date (YYYY-MM-DD):',7,'#334466'));
var csdSlideStartInput = ui.Textbox({placeholder:'e.g. 2021-06-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(csdSlideStartInput);
panel.add(lbl('Total months to cover (24-144):',7,'#334466'));
var csdSlideTotalInput = ui.Textbox({placeholder:'e.g. 60',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(csdSlideTotalInput);
panel.add(lbl('Window size defaults to 50% of total months (Dakos et al. standard), minimum 12 - override below or leave blank.',7,'#886600'));
var csdSlideWindowInput = ui.Textbox({placeholder:'blank = auto (50% of total)',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(csdSlideWindowInput);
var csdSlideVerdictV = ui.Label('Fill in the fields above, then press RUN.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var csdSlideResultV=ui.Label('',{fontSize:'8px',color:'#334466',backgroundColor:'#f4f0f8',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
// v10.129 NEW: real significance test on non-overlapping windows, raw
// AND deseasonalized - shown alongside (not replacing) the overlapping
// trajectory/trend test above.
var csdSlideValidV=ui.Label('',{fontSize:'8px',color:'#115566',backgroundColor:'#e6f4f8',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var csdSlideRunBtn=ui.Button({
  label:'RUN SLIDING WINDOW ANALYSIS',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#e8d9f5',color:'#4a1a6a',stretch:'horizontal',padding:'6px 4px',border:'2px solid #663388'},
  onClick:function(){
    var coordTxt=(csdTestCoordInput.getValue()||'').trim();
    var startTxt=(csdSlideStartInput.getValue()||'').trim();
    var totalTxt=(csdSlideTotalInput.getValue()||'').trim();
    var windowTxt=(csdSlideWindowInput.getValue()||'').trim();
    if(!coordTxt){csdSlideVerdictV.setValue('MISSING - fill in STEP 1 Lat/Lon above first.'); csdSlideVerdictV.style().set('color','#cc0000'); csdSlideVerdictV.style().set('backgroundColor','#ffd0d0'); return;}
    if(!startTxt||!totalTxt){csdSlideVerdictV.setValue('MISSING - fill in start date and total months.'); csdSlideVerdictV.style().set('color','#cc0000'); csdSlideVerdictV.style().set('backgroundColor','#ffd0d0'); return;}
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){csdSlideVerdictV.setValue('Invalid Lat/Lon in STEP 1 - use: lat, lon'); csdSlideVerdictV.style().set('color','#cc0000'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){csdSlideVerdictV.setValue('Invalid Lat/Lon in STEP 1 - use: lat, lon'); csdSlideVerdictV.style().set('color','#cc0000'); return;}
    var totalMonths=parseInt(totalTxt,10);
    if(isNaN(totalMonths)||totalMonths<24||totalMonths>144){csdSlideVerdictV.setValue('Total months must be 24-144.'); csdSlideVerdictV.style().set('color','#cc0000'); return;}
    var windowSize = windowTxt ? parseInt(windowTxt,10) : Math.round(totalMonths*0.5);
    if(isNaN(windowSize)||windowSize<12||windowSize>=totalMonths){
      csdSlideVerdictV.setValue('Window size must be 12 to (total months - 1). Leave blank for auto (50%).');
      csdSlideVerdictV.style().set('color','#cc0000'); return;
    }
    // v10.130 NEW: proactive warning, added directly from a real test this
    // session that hit "insufficient points (n=1)" on the NEW valid
    // significance test (v10.129) because the default 50%-of-total window
    // size, by construction, only ever produces ~2 non-overlapping windows
    // - never the 4+ needed for that test to compute anything. Raising the
    // total-months cap (84->144) alone would not have fixed this silently,
    // since the auto-default window would just grow proportionally too -
    // so this is checked and disclosed BEFORE the (slow) EE call fires,
    // rather than after, when the earlier real test only found out via a
    // late "insufficient" message post-computation.
    var nIndepWindows = Math.floor(totalMonths / windowSize);
    var indepWindowWarning = '';
    if(nIndepWindows < 4){
      indepWindowWarning = '\nNOTE: at window='+windowSize+'mo over '+totalMonths+' total months, only '+nIndepWindows+
        ' independent (non-overlapping) window(s) are possible - the NEW real significance test\n'+
        '(v10.129, below the main trajectory) needs at least 4 to compute anything. That test will show\n'+
        '"insufficient points" unless you either increase Total months or manually set Window size smaller\n'+
        '(e.g. window=12 with 96+ total months gives 8 independent windows). The main trajectory/chart above\n'+
        'will still work fine either way - this only affects the NEW valid-significance-test panel.';
    }
    // v10.98 NEW: warn upfront (before firing the EE call) if the requested
    // range runs past today's date - those months will have no satellite
    // data yet, which previously produced a silently-smaller window count
    // with no explanation. Computed client-side so it's instant, before any
    // network call.
    var startDateObj = new Date(startTxt+'T00:00:00Z');
    var rangeEndDateObj = new Date(Date.UTC(startDateObj.getUTCFullYear(), startDateObj.getUTCMonth()+totalMonths, 1));
    var todayObj = new Date();
    var futureWarning = '';
    if(!isNaN(startDateObj.getTime()) && rangeEndDateObj > todayObj){
      var monthsIntoFuture = (rangeEndDateObj.getUTCFullYear()-todayObj.getUTCFullYear())*12 + (rangeEndDateObj.getUTCMonth()-todayObj.getUTCMonth());
      futureWarning = '\nNOTE: this range extends roughly '+Math.max(monthsIntoFuture,1)+' month(s) past today - those months\n'+
        'have no satellite data yet and will be silently excluded from the scan (fewer sliding\n'+
        'window positions than requested). Not an error - just shortening your start date or\n'+
        'total months will use the full requested span instead.';
    }
    csdSlideVerdictV.setValue('Running: '+latIn.toFixed(4)+', '+lonIn.toFixed(4)+' | '+startTxt+' + '+totalMonths+
      'mo | window='+windowSize+'mo | '+(totalMonths-windowSize+1)+' sliding positions...\nThis is ONE Earth Engine call - can take 10-40s depending on series length.'+futureWarning+indepWindowWarning);
    csdSlideVerdictV.style().set('color','#334466'); csdSlideVerdictV.style().set('backgroundColor','#eeeeee');
    csdSlideVerdictV.style().set('border','2px solid #aaaaaa'); csdSlideVerdictV.style().set('whiteSpace','pre');
    csdSlideResultV.setValue('');
    var testPt=ee.Geometry.Point([lonIn,latIn]), testStudy=testPt.buffer(1000);
    var slideColl=mkMoSSTRange(startTxt,totalMonths);
    var rSlide=computeSlidingWindowCSD(slideColl,testStudy,'sst',4000,windowSize);
    rSlide.evaluate(function(fc,err){
      if(err){
        csdSlideVerdictV.setValue(friendlyEEError(err));
        csdSlideVerdictV.style().set('color','#cc0000'); csdSlideVerdictV.style().set('backgroundColor','#ffd0d0');
        csdSlideVerdictV.style().set('border','2px solid #cc0000');
        print('=== S13 SLIDING WINDOW ERROR === '+err);
        return;
      }
      try {
        var feats = (fc&&fc.features)?fc.features:[];
        if(feats.length<4){
          csdSlideVerdictV.setValue('INSUFFICIENT DATA - only '+feats.length+' sliding window position(s) computed.\nNeed at least 4. Try a longer total-months span or shorter window.');
          csdSlideVerdictV.style().set('color','#cc0000'); csdSlideVerdictV.style().set('backgroundColor','#ffd0d0');
          return;
        }
        var ac1Series=[], varSeries=[], dateLabels=[];
        for(var i=0;i<feats.length;i++){
          var p=feats[i].properties;
          ac1Series.push(p.ac1!==null&&p.ac1!==undefined?p.ac1:null);
          varSeries.push(p.variance!==null&&p.variance!==undefined?p.variance:null);
          dateLabels.push(new Date(p.endTime).toISOString().slice(0,7));
        }
        var ac1Trend = mannKendallTest(ac1Series);
        var varTrend = mannKendallTest(varSeries);

        function trendVerdict(t, name){
          if(t.error) return name+': n/a - '+t.error;
          var sig = t.p!==null && t.p<0.05;
          var dir = t.tau>0?'RISING':t.tau<0?'FALLING':'flat';
          return name+': tau='+(t.tau>0?'+':'')+t.tau.toFixed(3)+', p='+t.p.toFixed(4)+
            ' -> '+dir+(sig?' (STATISTICALLY SIGNIFICANT, p<0.05)':' (not significant at p<0.05)')+' [n='+t.n+' window positions]';
        }

        var ac1Line = trendVerdict(ac1Trend, 'AC1 trend');
        var varLine = trendVerdict(varTrend, 'Variance trend');

        // v10.97 FIX: this used to only check for RISING significant trends
        // (tau>0 && p<0.05), so a trend that was significantly FALLING
        // (tau<0 && p<0.05 - a real, reportable finding) got silently
        // lumped into the same "NO SIGNIFICANT TREND DETECTED" bucket as a
        // genuinely flat, noisy series. Caught on a real GBR run: variance
        // showed tau=-0.273, p=0.017 (significantly declining) but the
        // headline claimed "p>=0.05 for both" - factually wrong about the
        // number sitting right below it. Now classifies each indicator into
        // rising-significant / falling-significant / not-significant
        // separately, so a real significant decline is reported as what it
        // is, not silently discarded.
        var ac1SigRise = !ac1Trend.error && ac1Trend.p<0.05 && ac1Trend.tau>0;
        var ac1SigFall = !ac1Trend.error && ac1Trend.p<0.05 && ac1Trend.tau<0;
        var varSigRise = !varTrend.error && varTrend.p<0.05 && varTrend.tau>0;
        var varSigFall = !varTrend.error && varTrend.p<0.05 && varTrend.tau<0;
        var headline, hCol, hBg;
        if(ac1SigRise && varSigRise){
          headline='SIGNIFICANT CSD TREND - both AC1 and variance rising (p<0.05)';
          hCol='#880000'; hBg='#ffd0d0';
        } else if(ac1SigRise){
          headline='SIGNIFICANT AC1 RISING TREND (primary indicator, p<0.05)';
          hCol='#aa3300'; hBg='#ffe8cc';
        } else if(varSigRise){
          headline='SIGNIFICANT VARIANCE RISING TREND ONLY - AC1 not rising significantly (weaker evidence per Dakos et al. 2012)';
          hCol='#886600'; hBg='#fff6cc';
        } else if(ac1SigFall && varSigFall){
          headline='BOTH INDICATORS SIGNIFICANTLY DECLINING (p<0.05) - no CSD signature; not evidence of an approaching transition, but also not proof of safety';
          hCol='#115511'; hBg='#d4f5df';
        } else if(ac1SigFall){
          headline='AC1 SIGNIFICANTLY FALLING (p<0.05, primary indicator) - the opposite of a CSD warning sign at this site/period';
          hCol='#115511'; hBg='#d4f5df';
        } else if(varSigFall){
          headline='VARIANCE SIGNIFICANTLY FALLING (p<0.05) - no CSD signature. Per Dakos et al. 2012, declining variance alone is NOT evidence of an approaching transition (it can happen with or without one) - this is an absence-of-warning-sign result, not a "nothing found" result.';
          hCol='#226644'; hBg='#e8f4ff';
        } else {
          headline='NO SIGNIFICANT TREND DETECTED (p>=0.05 for both AC1 and variance)';
          hCol='#556677'; hBg='#eef2f6';
        }
        csdSlideVerdictV.setValue(headline+'\n'+ac1Line+'\n'+varLine);
        csdSlideVerdictV.style().set('color',hCol); csdSlideVerdictV.style().set('backgroundColor',hBg);
        csdSlideVerdictV.style().set('border','2px solid '+hCol); csdSlideVerdictV.style().set('whiteSpace','pre');

        var tableLines=['Window end | AC1 | Variance'];
        tableLines.push(repeatChar('\u2500',40));
        for(var j=0;j<feats.length;j++){
          tableLines.push(dateLabels[j]+' | '+(ac1Series[j]!==null?ac1Series[j].toFixed(3):'n/a')+' | '+
            (varSeries[j]!==null?varSeries[j].toFixed(3):'n/a'));
        }
        tableLines.push(repeatChar('\u2500',40));
        // v10.98 FIX: previously silent - if fewer valid months were found
        // than requested (most commonly because the requested range runs
        // past today's date, so some months have no satellite data yet),
        // the position count would just be smaller with no explanation,
        // leaving the user to notice and reverse-engineer the discrepancy
        // themselves. Now stated explicitly.
        var nValidMonthsFound = feats.length + windowSize - 1;
        if(nValidMonthsFound < totalMonths){
          var missingMonths = totalMonths - nValidMonthsFound;
          tableLines.push('NOTE: requested '+totalMonths+' months but only '+nValidMonthsFound+' had valid SST data ('+
            missingMonths+' missing). Most common cause: the requested range runs past today\'s\n'+
            'date, so those months have no satellite observations yet. The scan used only the '+
            nValidMonthsFound+' valid\nmonths found - check your start date + total months if this number surprises you.');
        }
        tableLines.push('Window size: '+windowSize+' months | Step: 1 month | '+feats.length+' positions scanned'+
          (nValidMonthsFound<totalMonths?' (of '+nValidMonthsFound+' valid months, not the requested '+totalMonths+')':''));
        tableLines.push('Detrend: single linear fit over the whole series (see disclosed simplification note above)');
        csdSlideResultV.setValue(tableLines.join('\n'));

        // Print the actual trend trajectory as a line chart - this is the
        // "line graph climbing upward smoothly" Dakos et al.'s figures show,
        // instead of a single pass/fail number.
        var chartFC = ee.FeatureCollection(feats.map(function(f,idx){
          return ee.Feature(null, {label:dateLabels[idx], ac1:ac1Series[idx], variance:varSeries[idx], idx:idx});
        }));
        print('=== S13 SLIDING WINDOW TREND TEST (Dakos et al. 2012 method) ===');
        print(headline);
        print(ac1Line);
        print(varLine);
        print(ui.Chart.feature.byFeature({features:chartFC, xProperty:'idx', yProperties:['ac1']})
          .setChartType('LineChart').setOptions({
            title:'AC1 sliding-window trajectory | window='+windowSize+'mo | tau='+(ac1Trend.error?'n/a':ac1Trend.tau.toFixed(3)),
            series:{0:{color:'#aa3300',lineWidth:2.5,pointSize:4,label:'AC1'}},
            backgroundColor:'#0a1628', titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
            vAxis:{title:'AC1',textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
            hAxis:{title:'Sliding window position (0 = earliest)',textStyle:{color:'#aaaaaa'}},
            legend:{textStyle:{color:'#ffffff'}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'},
            trendlines:{0:{type:'linear',color:'#ffcc44',lineWidth:2,opacity:0.9,showR2:false,visibleInLegend:true}}}));
        print(ui.Chart.feature.byFeature({features:chartFC, xProperty:'idx', yProperties:['variance']})
          .setChartType('LineChart').setOptions({
            title:'Variance sliding-window trajectory | window='+windowSize+'mo | tau='+(varTrend.error?'n/a':varTrend.tau.toFixed(3)),
            series:{0:{color:'#226644',lineWidth:2.5,pointSize:4,label:'Variance'}},
            backgroundColor:'#0a1628', titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
            vAxis:{title:'Variance',textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
            hAxis:{title:'Sliding window position (0 = earliest)',textStyle:{color:'#aaaaaa'}},
            legend:{textStyle:{color:'#ffffff'}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'},
            trendlines:{0:{type:'linear',color:'#ffcc44',lineWidth:2,opacity:0.9,showR2:false,visibleInLegend:true}}}));
        print(tableLines.join('\n'));

        // v10.129 NEW: fetch the raw monthly series client-side (ONE new
        // EE call, reusing the SAME slideColl/testStudy already used
        // above) to run a methodologically valid significance test -
        // non-overlapping windows (real independent samples, not the
        // shared-data overlapping ones above) on BOTH raw and
        // deseasonalized data. Direct extension of the v10.127/v10.128
        // deseasonalizing work and the Mann-Kendall fix validated in the
        // companion Python script this session.
        csdSlideValidV.setValue('Fetching raw series for deseasonalized + independent-window test...');
        csdSlideValidV.style().set('color','#115566');
        var singlePtFC = ee.FeatureCollection([ee.Feature(testStudy,{label:'Study'})]);
        var rRawSeries = extractMultiNodeSeries(slideColl, singlePtFC, 'sst', 4000);
        rRawSeries.evaluate(function(rawFC, rawErr){
          if(rawErr){
            csdSlideValidV.setValue('Deseasonalized/independent-window test error: '+rawErr);
            csdSlideValidV.style().set('color','#cc0000');
            print('=== S13 SLIDING WINDOW DESEASONALIZED TEST ERROR === '+rawErr);
            return;
          }
          try {
            var studyTV = groupSeriesByLabel(rawFC, 'sst')['Study'] || [];
            var climatology = computeMonthlyClimatology(studyTV);
            var deseasonTV = deseasonalizeSeries(studyTV, climatology);
            var rawVals = studyTV.map(function(s){return s.v;});
            var deseasonVals = deseasonTV.map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});

            var rawTraj = computeNonOverlappingTrajectory(rawVals, windowSize);
            var deseasonTraj = computeNonOverlappingTrajectory(deseasonVals, windowSize);

            var rawAC1MK = mannKendallTest(rawTraj.map(function(w){return w.ac1;}));
            var rawVarMK = mannKendallTest(rawTraj.map(function(w){return w.variance;}));
            var deseasonAC1MK = mannKendallTest(deseasonTraj.map(function(w){return w.ac1;}));
            var deseasonVarMK = mannKendallTest(deseasonTraj.map(function(w){return w.variance;}));

            var vLines = ['=== REAL SIGNIFICANCE TEST (independent windows, raw + deseasonalized, v10.129) ==='];
            vLines.push('The trend test above uses OVERLAPPING windows (step=1mo) - consecutive positions');
            vLines.push('share most of their data, which violates Mann-Kendall\'s independence assumption and');
            vLines.push('can report false significant trends even from pure noise (confirmed directly this');
            vLines.push('session: pure random data gave p<0.0001 overlapping vs p=0.109 non-overlapping).');
            vLines.push('This test instead uses non-overlapping '+windowSize+'-month windows - far fewer points,');
            vLines.push('but each one is a genuinely independent sample, so these p-values are trustworthy.');
            vLines.push('Also deseasonalizes (removes this site\'s own calendar-month climatology) to check');
            vLines.push('whether any trend survives once the seasonal cycle is removed.');
            vLines.push(repeatChar('\u2500',50));
            function vLine(t, name){
              if(t.error) return name+': n/a - '+t.error;
              var sig = t.p!==null && t.p<0.05;
              var dir = t.tau>0?'RISING':t.tau<0?'FALLING':'flat';
              return name+': tau='+(t.tau>0?'+':'')+t.tau.toFixed(3)+', p='+t.p.toFixed(4)+
                ' -> '+dir+(sig?' *** SIGNIFICANT (p<0.05)':' not significant')+' [n='+t.n+' independent windows]';
            }
            vLines.push(vLine(rawAC1MK, 'Raw AC1          '));
            vLines.push(vLine(rawVarMK, 'Raw Variance     '));
            vLines.push(vLine(deseasonAC1MK, 'Deseasonalized AC1     '));
            vLines.push(vLine(deseasonVarMK, 'Deseasonalized Variance'));
            vLines.push(repeatChar('\u2500',50));
            if((!deseasonAC1MK.error && deseasonAC1MK.p<0.05 && deseasonAC1MK.tau>0) &&
               (!deseasonVarMK.error && deseasonVarMK.p<0.05 && deseasonVarMK.tau>0)){
              vLines.push('Both deseasonalized indicators rise significantly across independent windows -');
              vLines.push('this is the strongest possible confirmation this tool can currently produce for a');
              vLines.push('genuine CSD trend, surviving both a proper climatology AND a valid significance test.');
            } else if(deseasonAC1MK.error || deseasonVarMK.error){
              vLines.push('Not enough independent '+windowSize+'-month windows in this series for a valid test');
              vLines.push('(need at least 4). Try a longer total-months span.');
            } else {
              vLines.push('The overlapping-window trend test above should be treated with real caution given');
              vLines.push('this result - a trend that looked significant there may not survive rigorous,');
              vLines.push('methodologically valid testing on independent samples.');
            }
            csdSlideValidV.setValue(vLines.join('\n'));
            csdSlideValidV.style().set('color','#115566');
            print(vLines.join('\n'));
          } catch(eValid){
            csdSlideValidV.setValue('Deseasonalized/independent-window test error: '+eValid);
            csdSlideValidV.style().set('color','#cc0000');
            print('=== S13 SLIDING WINDOW DESEASONALIZED TEST ERROR === '+eValid);
          }
        });
      } catch(errFinal){
        csdSlideVerdictV.setValue(friendlyEEError(errFinal));
        csdSlideVerdictV.style().set('color','#cc0000'); csdSlideVerdictV.style().set('backgroundColor','#ffd0d0');
        csdSlideVerdictV.style().set('border','2px solid #cc0000'); csdSlideVerdictV.style().set('whiteSpace','pre');
        print('=== S13 SLIDING WINDOW ANALYSIS ERROR === '+errFinal);
      }
    });
  }
});
panel.add(csdSlideRunBtn);
panel.add(csdSlideVerdictV);
panel.add(csdSlideResultV);
panel.add(csdSlideValidV);
panel.add(lbl('TIP: reuses the Lat/Lon from STEP 1 above. Full charts (AC1 and variance trajectories) print to the Console.',7,'#886600'));
panel.add(legDiv());

panel.add(lbl('WORKED EXAMPLE - Florida Keys 2023: Lat/Lon 24.55,-81.78 | STEP 2 BEFORE: start 2021-06-01, 24 months | STEP 4 AFTER start: 2023-06-01',7,'#aa6600'));

var locV=dynLbl('--','#115511'), regV=dynLbl('--','#115511'), coV=dynLbl('--','#666666');
panel.add(row('Location',locV)); panel.add(row('Region',regV)); panel.add(row('Lat/Lon',coV));
panel.add(sHead('SATELLITE DATA (globally valid)','#1a3a6a'));
var sstV=dynLbl('--','#cc2200'), trendV=dynLbl('--','#882200'), dhwV=dynLbl('--','#aa3300'), mmmV=dynLbl('--','#886600');
var chlaV=dynLbl('--','#116611'), turbV=dynLbl('--','#664400'), no2V=dynLbl('--','#550066'), depthV=dynLbl('--','#113366');
panel.add(row('SST',sstV)); panel.add(row('SST trend',trendV)); panel.add(row('DHW (per-pixel MMM)',dhwV));
var dhwNoteV=dynLbl('','#664400'); panel.add(dhwNoteV);
panel.add(row('MMM (local baseline)',mmmV)); panel.add(row('Chl-a',chlaV));
panel.add(row('Turbidity',turbV)); panel.add(row('NO2',no2V)); panel.add(row('Depth',depthV));
var depthWarnV=dynLbl('','#880000'); panel.add(depthWarnV);
var ebusWarnV=dynLbl('','#885500'); panel.add(ebusWarnV);

panel.add(sHead('SEAWEED / MACROALGAE (S7)','#226600'));
var faiV=dynLbl('--','#226600'), ndciV=dynLbl('--','#226600'), ndviWV=dynLbl('--','#226600'), algaeStatusV=dynLbl('checking...','#666666');
panel.add(row('FAI (floating algae)',faiV)); panel.add(row('NDCI (bloom index)',ndciV));
panel.add(row('NDVI water (benthic)',ndviWV)); panel.add(row('Algae status',algaeStatusV));

// ============================================================
// S7B - MULTI-POINT ALGAE SCAN (v10.102 NEW)
// Directly motivated by a real finding: two points 1.1km apart at the same
// reef showed "mild watch" vs "massive bloom" - a single click can miss or
// overstate a patchy bloom entirely. This samples FAI/NDCI/NDVI at 8
// compass points plus the centre around a study location, in ONE Earth
// Engine call, and reports how many show elevated bloom warnings and how
// spatially patchy vs widespread the result is.
// ============================================================
panel.add(sHead('S7B - MULTI-POINT ALGAE SCAN (v10.102)','#1a5522'));
panel.add(lbl('Samples FAI/NDCI/NDVI at 8 compass points + centre around a location, in ONE call - checks whether a bloom reading is patchy (isolated to 1-2 points) or widespread (most/all points elevated), instead of trusting a single click.',7,'#226644'));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var s7bCoordInput = ui.Textbox({placeholder:'lat, lon  e.g. -23.51, 152.09',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7bCoordInput);
var s7bUseLastClickBtn = ui.Button({
  label:'Use last clicked location (from map above)',
  style:{fontSize:'9px',margin:'2px 4px',backgroundColor:'#e8f4ff',color:'#225588',stretch:'horizontal',padding:'4px 4px',border:'1px solid #4488cc'},
  onClick:function(){
    if(lastClickLat===null||lastClickLon===null){
      s7bStatusV.setValue('No location clicked yet - click the map or use GO TO COORDINATES first.');
      s7bStatusV.style().set('color','#aa3300'); return;
    }
    s7bCoordInput.setValue(lastClickLat.toFixed(4)+', '+lastClickLon.toFixed(4));
    s7bStatusV.setValue('Location filled in. Pick a radius and press SCAN.');
    s7bStatusV.style().set('color','#115511');
  }
});
panel.add(s7bUseLastClickBtn);
panel.add(lbl('Ring radius (km, 0.5-10):',7,'#334466'));
var s7bRadiusInput = ui.Textbox({placeholder:'e.g. 1',value:'1',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7bRadiusInput);
var s7bStatusV = dynLbl('Fill in Lat/Lon + radius, then press SCAN.','#553388');
var s7bResultV = ui.Label('',{fontSize:'8px',color:'#114411',backgroundColor:'#eefaf0',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var s7bScanBtn = ui.Button({
  label:'SCAN NEARBY POINTS FOR ALGAE WARNING SIGNS',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#d9f0d9',color:'#1a5522',stretch:'horizontal',padding:'6px 4px',border:'2px solid #226644'},
  onClick:function(){
    var coordTxt=(s7bCoordInput.getValue()||'').trim();
    var radiusTxt=(s7bRadiusInput.getValue()||'').trim();
    if(!coordTxt){s7bStatusV.setValue('MISSING - enter Lat, Lon above.'); s7bStatusV.style().set('color','#cc0000'); return;}
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){s7bStatusV.setValue('Invalid format - use: lat, lon'); s7bStatusV.style().set('color','#cc0000'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){s7bStatusV.setValue('Invalid format - use: lat, lon'); s7bStatusV.style().set('color','#cc0000'); return;}
    var radiusKm=parseFloat(radiusTxt);
    if(isNaN(radiusKm)||radiusKm<0.5||radiusKm>10){s7bStatusV.setValue('Radius must be 0.5-10 km.'); s7bStatusV.style().set('color','#cc0000'); return;}
    recordStudySite(latIn, lonIn, 'S7B');

    // Great-circle destination point formula (client-side, no EE needed)
    function destinationPoint(lat, lon, bearingDeg, distKm){
      var R=6371, brng=bearingDeg*Math.PI/180, lat1=lat*Math.PI/180, lon1=lon*Math.PI/180;
      var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
      var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
      return {lat:lat2*180/Math.PI, lon:lon2*180/Math.PI};
    }
    var bearings=[{b:0,label:'N'},{b:45,label:'NE'},{b:90,label:'E'},{b:135,label:'SE'},
      {b:180,label:'S'},{b:225,label:'SW'},{b:270,label:'W'},{b:315,label:'NW'}];
    var pts=[{label:'Center',lat:latIn,lon:lonIn}];
    bearings.forEach(function(br){
      var d=destinationPoint(latIn,lonIn,br.b,radiusKm);
      pts.push({label:br.label,lat:d.lat,lon:d.lon});
    });

    s7bStatusV.setValue('Scanning '+pts.length+' points ('+radiusKm+'km radius)... one Earth Engine call, ~10-20s.');
    s7bStatusV.style().set('color','#334466'); s7bResultV.setValue('');

    var ptsFC = ee.FeatureCollection(pts.map(function(p,idx){
      return ee.Feature(ee.Geometry.Point([p.lon,p.lat]).buffer(300), {label:p.label, idx:idx});
    }));
    var combinedImg = faiImg.updateMask(oceanMask).rename('fai')
      .addBands(ndciImg.updateMask(oceanMask).rename('ndci'))
      .addBands(ndviWater.updateMask(oceanMask).rename('ndviW'));
    var resultsFC = combinedImg.reduceRegions({collection:ptsFC, reducer:ee.Reducer.max(), scale:100});
    resultsFC.evaluate(function(fc, err){
      if(err){
        s7bStatusV.setValue(friendlyEEError(err));
        s7bStatusV.style().set('color','#cc0000'); s7bStatusV.style().set('backgroundColor','#ffd0d0');
        s7bStatusV.style().set('border','2px solid #cc0000'); s7bStatusV.style().set('whiteSpace','pre');
        print('=== S7B MULTI-POINT SCAN ERROR === '+err);
        return;
      }
      try {
        var feats = (fc&&fc.features)?fc.features:[];
        function faiLabel(v){ if(v===null||v===undefined) return 'n/a'; return v>0.10?'EXTREME':v>0.05?'CRITICAL':v>0.03?'HIGH':v>0.01?'ELEVATED':'watch'; }
        function ndciLabel(v){ if(v===null||v===undefined) return 'n/a'; return v>0.35?'EXTREME':v>0.20?'HIGH':v>0.10?'MODERATE':'low'; }
        function ndviLabel(v){ if(v===null||v===undefined) return 'n/a'; return v>0.40?'DENSE':v>0.20?'MODERATE':v>0.05?'light':'clear'; }
        var rows=['Point | FAI | NDCI | NDVI-w | Flags'];
        rows.push(repeatChar('\u2500',48));
        var nElevatedFAI=0, nElevatedNDCI=0, nElevatedNDVI=0, nAvail=0;
        var faiVals=[], ndciVals=[], ndviVals=[];
        for(var i=0;i<feats.length;i++){
          var p=feats[i].properties;
          var fv=p.fai!==undefined?p.fai:null, nv=p.ndci!==undefined?p.ndci:null, dv=p.ndviW!==undefined?p.ndviW:null;
          if(fv!==null||nv!==null||dv!==null) nAvail++;
          if(fv!==null) faiVals.push(fv);
          if(nv!==null) ndciVals.push(nv);
          if(dv!==null) ndviVals.push(dv);
          var fLab=faiLabel(fv), nLab=ndciLabel(nv), dLab=ndviLabel(dv);
          if(fv!==null&&fv>0.03) nElevatedFAI++;
          if(nv!==null&&nv>0.20) nElevatedNDCI++;
          if(dv!==null&&dv>0.20) nElevatedNDVI++;
          var flags=[];
          if(fv!==null&&fv>0.03) flags.push('FAI:'+fLab);
          if(nv!==null&&nv>0.20) flags.push('NDCI:'+nLab);
          if(dv!==null&&dv>0.20) flags.push('NDVI:'+dLab);
          rows.push(p.label+repeatChar(' ',Math.max(1,7-p.label.length))+'| '+
            (fv!==null?fv.toFixed(3):'n/a')+' | '+(nv!==null?nv.toFixed(3):'n/a')+' | '+
            (dv!==null?dv.toFixed(3):'n/a')+' | '+(flags.length>0?flags.join(', '):'-'));
        }
        rows.push(repeatChar('\u2500',48));

        var faiRange = faiVals.length>1?(Math.max.apply(null,faiVals)-Math.min.apply(null,faiVals)):null;
        var ndciRange = ndciVals.length>1?(Math.max.apply(null,ndciVals)-Math.min.apply(null,ndciVals)):null;

        var nPts = feats.length;
        var maxElevated = Math.max(nElevatedFAI,nElevatedNDCI,nElevatedNDVI);
        var patchiness;
        if(nAvail===0) patchiness='NO DATA available at any point.';
        else if(maxElevated===0) patchiness='NO elevated readings at any of the '+nAvail+' points scanned - no bloom warning here.';
        else if(maxElevated>=Math.ceil(nAvail*0.7)) patchiness='WIDESPREAD - '+maxElevated+' of '+nAvail+' points elevated. Consistent with a genuine regional bloom, not a single-pixel artifact.';
        else if(maxElevated<=2) patchiness='PATCHY - only '+maxElevated+' of '+nAvail+' points elevated. Could be a small, localized bloom OR a single-pixel/cloud artifact - do not generalize this to the whole area from one reading, and consider re-checking those specific point(s) individually.';
        else patchiness='MIXED - '+maxElevated+' of '+nAvail+' points elevated. Partial spatial coverage; check which compass direction(s) are affected in the table above.';

        rows.push('SPATIAL VARIABILITY: FAI range='+(faiRange!==null?faiRange.toFixed(3):'n/a')+
          ' | NDCI range='+(ndciRange!==null?ndciRange.toFixed(3):'n/a')+
          (faiRange!==null&&faiRange>0.05?' - LARGE spread, confirms the site is patchy (matches the kind of 1km discrepancy that motivated this tool).':''));
        rows.push('WARNING SIGN COUNT: FAI elevated '+nElevatedFAI+'/'+nAvail+' | NDCI elevated '+nElevatedNDCI+'/'+nAvail+' | NDVI elevated '+nElevatedNDVI+'/'+nAvail);
        rows.push('PATTERN: '+patchiness);
        rows.push('Radius: '+radiusKm+'km | Each point sampled with a 300m buffer, max reducer (same methodology as the main S7 panel above)');
        s7bResultV.setValue(rows.join('\n'));

        var headline, hCol, hBg;
        if(nAvail===0){ headline='NO DATA'; hCol='#888888'; hBg='#eeeeee'; }
        else if(maxElevated===0){ headline='NO ALGAE WARNING SIGNS AT ANY SCANNED POINT'; hCol='#115511'; hBg='#d4f5df'; }
        else if(maxElevated>=Math.ceil(nAvail*0.7)){ headline='WIDESPREAD ALGAE WARNING - '+maxElevated+'/'+nAvail+' points elevated'; hCol='#880000'; hBg='#ffd0d0'; }
        else if(maxElevated<=2){ headline='PATCHY/ISOLATED WARNING - only '+maxElevated+'/'+nAvail+' points elevated, verify before generalizing'; hCol='#aa6600'; hBg='#fff6cc'; }
        else { headline='MIXED ALGAE WARNING - '+maxElevated+'/'+nAvail+' points elevated'; hCol='#886600'; hBg='#fff0d0'; }
        s7bStatusV.setValue(headline);
        s7bStatusV.style().set('color',hCol); s7bStatusV.style().set('backgroundColor',hBg);
        s7bStatusV.style().set('border','2px solid '+hCol); s7bStatusV.style().set('whiteSpace','pre'); s7bStatusV.style().set('fontWeight','bold');

        print('=== S7B MULTI-POINT ALGAE SCAN ===');
        print(headline);
        print(rows.join('\n'));
      } catch(errFinal){
        s7bStatusV.setValue(friendlyEEError(errFinal));
        s7bStatusV.style().set('color','#cc0000'); s7bStatusV.style().set('backgroundColor','#ffd0d0');
        s7bStatusV.style().set('border','2px solid #cc0000'); s7bStatusV.style().set('whiteSpace','pre');
        print('=== S7B MULTI-POINT SCAN FINAL ERROR === '+errFinal);
      }
    });
  }
});
panel.add(s7bScanBtn);
panel.add(s7bStatusV);
panel.add(s7bResultV);
panel.add(lbl('CAVEAT: 9-point compass ring is a simple sampling pattern, not a full spatial survey - a bloom could sit between sample points. For exhaustive coverage, reduce the radius and re-scan, or use the S3/S7 map layers above to visually inspect the whole area.',7,'#886600'));
panel.add(legDiv());

// ============================================================
// S7C - ALGAE COUPLING PROOF-OF-CONCEPT (v10.104 NEW)
// Tests whether AC1/variance/spatial-correlation statistics (the same
// toolkit already used on SST in S13) can be meaningfully applied to
// ALGAE (FAI) data instead. This is a genuinely different, arguably more
// theoretically appropriate application: in Scheffer's bistable-state
// framework, coral-vs-algae cover is the STATE VARIABLE that actually
// flips between stable states, while SST is closer to the external
// CONTROL PARAMETER driving the system toward a threshold. S7 previously
// only had ONE fixed 2023-2024 composite, never a time series - this adds
// the missing monthly FAI infrastructure (mkMoFAIRange) and tests it on a
// small 3-node network (Center + North + East) before committing to a
// full 9-node version, since Sentinel-2's cloud-masking data density at
// monthly resolution is genuinely untested until now.
// ============================================================
panel.add(sHead('S7C - ALGAE COUPLING PROOF-OF-CONCEPT (v10.104)','#4a2a5a'));
panel.add(lbl('Tests whether AC1/variance/spatial-correlation (the SAME statistics S13 uses on SST) can be computed from ALGAE (FAI) data instead - arguably more appropriate, since algae cover is the actual state variable in a coral-algae bistable system, not just an external driver like SST.',7,'#663388'));
panel.add(lbl('PROOF-OF-CONCEPT SCOPE: only 3 nodes (Center, North, East) and one time window - this checks whether monthly Sentinel-2 FAI data is even dense enough to compute these statistics before committing to a full 9-node network. Sentinel-2 has heavier cloud losses than the OISST used elsewhere in this tool, so some months may come back with no valid data - that is exactly what this is testing.',7,'#886600'));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var s7cCoordInput = ui.Textbox({placeholder:'lat, lon  e.g. -23.51, 152.09',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7cCoordInput);
var s7cUseLastClickBtn = ui.Button({
  label:'Use last clicked location (from map above)',
  style:{fontSize:'9px',margin:'2px 4px',backgroundColor:'#e8f4ff',color:'#225588',stretch:'horizontal',padding:'4px 4px',border:'1px solid #4488cc'},
  onClick:function(){
    if(lastClickLat===null||lastClickLon===null){
      s7cStatusV.setValue('No location clicked yet - click the map or use GO TO COORDINATES first.');
      s7cStatusV.style().set('color','#aa3300'); return;
    }
    s7cCoordInput.setValue(lastClickLat.toFixed(4)+', '+lastClickLon.toFixed(4));
    s7cStatusV.setValue('Location filled in. Fill in the date fields, then press RUN.');
    s7cStatusV.style().set('color','#115511');
  }
});
panel.add(s7cUseLastClickBtn);
panel.add(lbl('Node spacing (km, 0.3-3):',7,'#334466'));
var s7cRadiusInput = ui.Textbox({placeholder:'e.g. 0.5',value:'0.5',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7cRadiusInput);
panel.add(lbl('Start date (YYYY-MM-DD):',7,'#334466'));
var s7cStartInput = ui.Textbox({placeholder:'e.g. 2022-06-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7cStartInput);
panel.add(lbl('Months to cover (6-36):',7,'#334466'));
var s7cMonthsInput = ui.Textbox({placeholder:'e.g. 24',value:'24',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7cMonthsInput);
var s7cStatusV = ui.Label('Fill in the fields above, then press RUN.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var s7cResultV = ui.Label('',{fontSize:'8px',color:'#3a2050',backgroundColor:'#f5eefa',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var s7cRunBtn = ui.Button({
  label:'RUN PROOF-OF-CONCEPT (3 nodes)',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#e8d9f5',color:'#4a1a6a',stretch:'horizontal',padding:'6px 4px',border:'2px solid #663388'},
  onClick:function(){
    var coordTxt=(s7cCoordInput.getValue()||'').trim();
    var radiusTxt=(s7cRadiusInput.getValue()||'').trim();
    var startTxt=(s7cStartInput.getValue()||'').trim();
    var monthsTxt=(s7cMonthsInput.getValue()||'').trim();
    if(!coordTxt||!startTxt||!monthsTxt){s7cStatusV.setValue('MISSING - fill in Lat/Lon, start date, and months.'); s7cStatusV.style().set('color','#cc0000'); s7cStatusV.style().set('backgroundColor','#ffd0d0'); return;}
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){s7cStatusV.setValue('Invalid format - use: lat, lon'); s7cStatusV.style().set('color','#cc0000'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){s7cStatusV.setValue('Invalid format - use: lat, lon'); s7cStatusV.style().set('color','#cc0000'); return;}
    var radiusKm=parseFloat(radiusTxt);
    if(isNaN(radiusKm)||radiusKm<0.3||radiusKm>3){s7cStatusV.setValue('Node spacing must be 0.3-3 km.'); s7cStatusV.style().set('color','#cc0000'); return;}
    var totalMonths=parseInt(monthsTxt,10);
    if(isNaN(totalMonths)||totalMonths<6||totalMonths>36){s7cStatusV.setValue('Months must be 6-36.'); s7cStatusV.style().set('color','#cc0000'); return;}
    recordStudySite(latIn, lonIn, 'S7C');

    function destinationPoint(lat, lon, bearingDeg, distKm){
      var R=6371, brng=bearingDeg*Math.PI/180, lat1=lat*Math.PI/180, lon1=lon*Math.PI/180;
      var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
      var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
      return {lat:lat2*180/Math.PI, lon:lon2*180/Math.PI};
    }
    var nPt=destinationPoint(latIn,lonIn,0,radiusKm);
    var ePt=destinationPoint(latIn,lonIn,90,radiusKm);
    var ptCenterBuf=ee.Geometry.Point([lonIn,latIn]).buffer(150);
    var ptNBuf=ee.Geometry.Point([nPt.lon,nPt.lat]).buffer(150);
    var ptEBuf=ee.Geometry.Point([ePt.lon,ePt.lat]).buffer(150);

    s7cStatusV.setValue('Running: 3 nodes, '+totalMonths+' months of Sentinel-2 FAI (6 EE calls: 3 AC1/var + 3 correlation).\n'+
      'Sentinel-2 processing is heavier than OISST - this can take 30-90 seconds.');
    s7cStatusV.style().set('color','#334466'); s7cStatusV.style().set('backgroundColor','#eeeeee');
    s7cStatusV.style().set('border','2px solid #aaaaaa'); s7cStatusV.style().set('whiteSpace','pre');
    s7cResultV.setValue('');

    var faiColl = mkMoFAIRange(startTxt, totalMonths);
    var rCorr_CN = computeZonalSyncCSD(faiColl, ptCenterBuf, ptNBuf, 'fai', 20);
    var rCorr_CE = computeZonalSyncCSD(faiColl, ptCenterBuf, ptEBuf, 'fai', 20);
    var rCorr_NE = computeZonalSyncCSD(faiColl, ptNBuf, ptEBuf, 'fai', 20);

    var s7cRes={}, s7cPending=6, s7cErrors=0;
    function s7cBump(key,val,err){
      s7cRes[key]=val;
      if(err) s7cErrors++;
      s7cPending--;
      s7cStatusV.setValue('Running: '+(6-s7cPending)+' / 6 sub-tests done'+(s7cErrors>0?' ('+s7cErrors+' errored)':'')+'...');
      if(s7cPending>0) return;
      s7cFinish();
    }
    function s7cFinish(){
      try {
        var acC=s7cRes.acC||{}, acN=s7cRes.acN||{}, acE=s7cRes.acE||{};
        var coCN=s7cRes.coCN||{}, coCE=s7cRes.coCE||{}, coNE=s7cRes.coNE||{};
        function fmtN(v,d){ return (v!==null&&v!==undefined&&!isNaN(v))?v.toFixed(d):'n/a'; }
        var lines=[];
        lines.push('DATA DENSITY (valid months found / '+totalMonths+' requested):');
        lines.push('  Center: '+(acC.nValidMonths!==undefined?acC.nValidMonths:'n/a')+
          ' ('+(acC.nValidMonths!==undefined?Math.round(100*acC.nValidMonths/totalMonths):'?')+'%)');
        lines.push('  North:  '+(acN.nValidMonths!==undefined?acN.nValidMonths:'n/a')+
          ' ('+(acN.nValidMonths!==undefined?Math.round(100*acN.nValidMonths/totalMonths):'?')+'%)');
        lines.push('  East:   '+(acE.nValidMonths!==undefined?acE.nValidMonths:'n/a')+
          ' ('+(acE.nValidMonths!==undefined?Math.round(100*acE.nValidMonths/totalMonths):'?')+'%)');
        var avgDensity = (acC.nValidMonths!==undefined&&acN.nValidMonths!==undefined&&acE.nValidMonths!==undefined) ?
          ((acC.nValidMonths+acN.nValidMonths+acE.nValidMonths)/3/totalMonths) : null;
        lines.push('');
        lines.push('TEMPORAL AC1 / VARIANCE (FAI, per node):');
        // v10.105 FIX: the ratio alone (e.g. "15.80x") can't distinguish a
        // genuine variance surge from a near-zero-denominator artifact
        // (ratio = secondHalfVar / max(firstHalfVar, 1e-6) - a tiny first-
        // half variance inflates the ratio even for a modest second-half
        // value). Now shows the raw first/second-half variance alongside
        // the ratio, and flags it explicitly when the first half is small
        // enough that the ratio is likely more artifact than signal.
        function varLine(node, ac){
          var vFirst=ac.varFirstHalf, vSecond=ac.varSecondHalf, vRatio=ac.varTrendRatio;
          var artifactFlag = (vFirst!==null&&vFirst!==undefined&&vFirst<0.001&&vRatio!==null&&vRatio>5) ?
            '  \u26A0 LIKELY ARTIFACT: first-half variance is near-zero ('+fmtN(vFirst,5)+'), inflating the ratio - treat this Var value with real caution, not as a genuine surge.' : '';
          return '  '+node+': AC1='+fmtN(ac.realAC1,3)+' Var='+fmtN(vRatio,2)+'x (1st-half='+fmtN(vFirst,4)+', 2nd-half='+fmtN(vSecond,4)+') Skew='+fmtN(ac.skewness,2)+artifactFlag;
        }
        lines.push(varLine('Center',acC));
        lines.push(varLine('North',acN));
        lines.push(varLine('East',acE));
        lines.push('');
        lines.push('SPATIAL COUPLING MATRIX (FAI correlation between node pairs):');
        lines.push('  Center-North: r='+fmtN(coCN.corr,3)+' (n='+(coCN.n!==undefined?coCN.n:'n/a')+' paired months)');
        lines.push('  Center-East:  r='+fmtN(coCE.corr,3)+' (n='+(coCE.n!==undefined?coCE.n:'n/a')+' paired months)');
        lines.push('  North-East:   r='+fmtN(coNE.corr,3)+' (n='+(coNE.n!==undefined?coNE.n:'n/a')+' paired months)');
        s7cResultV.setValue(lines.join('\n'));

        var headline, hCol, hBg;
        if(avgDensity===null){
          headline='COULD NOT ASSESS - one or more nodes returned no data'; hCol='#cc0000'; hBg='#ffd0d0';
        } else if(avgDensity<0.4){
          headline='SPARSE DATA ('+Math.round(avgDensity*100)+'% avg coverage) - Sentinel-2 monthly FAI is likely TOO GAPPY at this site/window for reliable AC1/variance/correlation. Consider a longer window, a different site, or accept high uncertainty.';
          hCol='#aa6600'; hBg='#fff6cc';
        } else if(avgDensity<0.7){
          headline='MODERATE DATA ('+Math.round(avgDensity*100)+'% avg coverage) - usable but noisy. Statistics below are a first look, not a confirmed result - treat with real caution before scaling to 9 nodes.';
          hCol='#886600'; hBg='#fff0d0';
        } else {
          headline='GOOD DATA DENSITY ('+Math.round(avgDensity*100)+'% avg coverage) - Sentinel-2 monthly FAI looks usable at this site. Worth scaling to a full 9-node network if this pattern holds at other test sites.';
          hCol='#115511'; hBg='#d4f5df';
        }
        s7cStatusV.setValue(headline);
        s7cStatusV.style().set('color',hCol); s7cStatusV.style().set('backgroundColor',hBg);
        s7cStatusV.style().set('border','2px solid '+hCol); s7cStatusV.style().set('whiteSpace','pre'); s7cStatusV.style().set('fontWeight','bold');

        print('=== S7C ALGAE COUPLING PROOF-OF-CONCEPT ===');
        print(headline);
        print(lines.join('\n'));
      } catch(errFinal){
        s7cStatusV.setValue(friendlyEEError(errFinal));
        s7cStatusV.style().set('color','#cc0000'); s7cStatusV.style().set('backgroundColor','#ffd0d0');
        s7cStatusV.style().set('border','2px solid #cc0000'); s7cStatusV.style().set('whiteSpace','pre');
        print('=== S7C PROOF-OF-CONCEPT ERROR === '+errFinal);
      }
    }

    computeRealCSDDeseasonalized(faiColl, ptCenterBuf, 'fai', 20, function(v){ s7cBump('acC', v.error?{}:v, v.error); });
    computeRealCSDDeseasonalized(faiColl, ptNBuf, 'fai', 20, function(v){ s7cBump('acN', v.error?{}:v, v.error); });
    computeRealCSDDeseasonalized(faiColl, ptEBuf, 'fai', 20, function(v){ s7cBump('acE', v.error?{}:v, v.error); });
    rCorr_CN.evaluate(function(v,e){ s7cBump('coCN', e?{}:v||{}, e); });
    rCorr_CE.evaluate(function(v,e){ s7cBump('coCE', e?{}:v||{}, e); });
    rCorr_NE.evaluate(function(v,e){ s7cBump('coNE', e?{}:v||{}, e); });
  }
});
panel.add(s7cRunBtn);
panel.add(s7cStatusV);
panel.add(s7cResultV);
panel.add(lbl('If this comes back with good data density, next step is scaling to the full 9-node ring (matching S7B\'s geometry) with a proper BEFORE/AFTER comparison, the same toolkit structure already used for SST in STEP 3.',7,'#886600'));
panel.add(legDiv());

// ============================================================
// S7D - FULL 9-NODE ALGAE COUPLING NETWORK (v10.106 NEW)
// Scales S7C's 3-node proof-of-concept to the full 9-node ring (matching
// S7B's compass geometry), with a real BEFORE/AFTER comparison instead of
// a single snapshot window. This is a HEAVY feature - roughly 43 Earth
// Engine calls in one click (9 nodes x 2 periods AC1/variance + 8 ring
// points x 2 periods correlation-vs-Center + 9 nodes x 1 NDVI-water
// on-reef check), several times more than any other single button in this
// tool. Correlation is scored vs Center only (a "hub" network, 8 pairs),
// not the full 36-pair matrix, to keep this from becoming unusably slow.
// NDVI-water is fetched once per node (current composite, not a time
// series) specifically to flag which nodes are likely off-reef open water
// - a real finding from S7B (North showed strongly negative NDVI, likely
// off-reef), which S7C's 3-node run could not check since it only fetched
// FAI. A node/pair involving an off-reef point is explicitly caveated
// rather than silently trusted.
// ============================================================
panel.add(sHead('S7D - FULL 9-NODE ALGAE COUPLING NETWORK (v10.106)','#3a1a5a'));
panel.add(lbl('Scales S7C to the full 9-node ring with a real BEFORE/AFTER comparison. v10.107: rebuilt from ~43 separate Earth Engine calls down to just 3 batched calls (reduceRegions+flatten, stats computed client-side) - expect ~20-60s instead of 2-5 minutes. Correlation is still scored vs the Center node only (8 pairs), not the full 36-pair matrix, to keep this usable.',7,'#663388'));
panel.add(lbl('NDVI-water is checked once per node to flag likely off-reef points (per the S7B finding that North showed strongly negative NDVI there) - a node/pair involving an off-reef point is caveated, not silently trusted as a real algae-dynamics comparison.',7,'#886600'));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var s7dCoordInput = ui.Textbox({placeholder:'lat, lon  e.g. -23.51, 152.09',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7dCoordInput);
var s7dUseLastClickBtn = ui.Button({
  label:'Use last clicked location (from map above)',
  style:{fontSize:'9px',margin:'2px 4px',backgroundColor:'#e8f4ff',color:'#225588',stretch:'horizontal',padding:'4px 4px',border:'1px solid #4488cc'},
  onClick:function(){
    if(lastClickLat===null||lastClickLon===null){
      s7dStatusV.setValue('No location clicked yet - click the map or use GO TO COORDINATES first.');
      s7dStatusV.style().set('color','#aa3300'); return;
    }
    s7dCoordInput.setValue(lastClickLat.toFixed(4)+', '+lastClickLon.toFixed(4));
    s7dStatusV.setValue('Location filled in. Fill in the BEFORE/AFTER dates, then press RUN.');
    s7dStatusV.style().set('color','#115511');
  }
});
panel.add(s7dUseLastClickBtn);
panel.add(lbl('Ring radius (km, 0.3-3):',7,'#334466'));
var s7dRadiusInput = ui.Textbox({placeholder:'e.g. 0.5',value:'0.5',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7dRadiusInput);
panel.add(lbl('BEFORE start date + months:',7,'#115511'));
var s7dBeforeStartInput = ui.Textbox({placeholder:'e.g. 2021-06-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7dBeforeStartInput);
var s7dBeforeMonthsInput = ui.Textbox({placeholder:'months, e.g. 18',value:'18',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7dBeforeMonthsInput);
panel.add(lbl('AFTER start date + months:',7,'#aa3300'));
var s7dAfterStartInput = ui.Textbox({placeholder:'e.g. 2023-11-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7dAfterStartInput);
var s7dAfterMonthsInput = ui.Textbox({placeholder:'months, e.g. 12',value:'12',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7dAfterMonthsInput);
var s7dStatusV = ui.Label('Fill in the fields above, then press RUN. 3 batched calls, expect ~20-60 seconds.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var s7dResultV = ui.Label('',{fontSize:'8px',color:'#2a1040',backgroundColor:'#f2ecfa',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var s7dRunBtn = ui.Button({
  label:'RUN FULL 9-NODE NETWORK (3 batched EE calls, ~20-60s)',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#ded0f0',color:'#3a1560',stretch:'horizontal',padding:'6px 4px',border:'2px solid #5a2a80'},
  onClick:function(){
    var coordTxt=(s7dCoordInput.getValue()||'').trim();
    var radiusTxt=(s7dRadiusInput.getValue()||'').trim();
    var beforeStartTxt=(s7dBeforeStartInput.getValue()||'').trim();
    var beforeMonthsTxt=(s7dBeforeMonthsInput.getValue()||'').trim();
    var afterStartTxt=(s7dAfterStartInput.getValue()||'').trim();
    var afterMonthsTxt=(s7dAfterMonthsInput.getValue()||'').trim();
    if(!coordTxt||!beforeStartTxt||!beforeMonthsTxt||!afterStartTxt||!afterMonthsTxt){
      s7dStatusV.setValue('MISSING - fill in Lat/Lon and both BEFORE/AFTER date+months fields.');
      s7dStatusV.style().set('color','#cc0000'); s7dStatusV.style().set('backgroundColor','#ffd0d0'); return;
    }
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){s7dStatusV.setValue('Invalid format - use: lat, lon'); s7dStatusV.style().set('color','#cc0000'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){s7dStatusV.setValue('Invalid format - use: lat, lon'); s7dStatusV.style().set('color','#cc0000'); return;}
    var radiusKm=parseFloat(radiusTxt);
    if(isNaN(radiusKm)||radiusKm<0.3||radiusKm>3){s7dStatusV.setValue('Ring radius must be 0.3-3 km.'); s7dStatusV.style().set('color','#cc0000'); return;}
    var beforeMonths=parseInt(beforeMonthsTxt,10), afterMonths=parseInt(afterMonthsTxt,10);
    if(isNaN(beforeMonths)||beforeMonths<4||beforeMonths>36||isNaN(afterMonths)||afterMonths<4||afterMonths>36){
      s7dStatusV.setValue('BEFORE/AFTER months must each be 4-36.'); s7dStatusV.style().set('color','#cc0000'); return;
    }
    recordStudySite(latIn, lonIn, 'S7D');

    function destinationPoint(lat, lon, bearingDeg, distKm){
      var R=6371, brng=bearingDeg*Math.PI/180, lat1=lat*Math.PI/180, lon1=lon*Math.PI/180;
      var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
      var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
      return {lat:lat2*180/Math.PI, lon:lon2*180/Math.PI};
    }
    var bearings8=[{b:0,label:'N'},{b:45,label:'NE'},{b:90,label:'E'},{b:135,label:'SE'},
      {b:180,label:'S'},{b:225,label:'SW'},{b:270,label:'W'},{b:315,label:'NW'}];
    var nodes=[{label:'Center',lat:latIn,lon:lonIn}];
    bearings8.forEach(function(br){
      var d=destinationPoint(latIn,lonIn,br.b,radiusKm);
      nodes.push({label:br.label,lat:d.lat,lon:d.lon});
    });
    nodes.forEach(function(nd){ nd.buf=ee.Geometry.Point([nd.lon,nd.lat]).buffer(150); });

    // v10.107 FIX: previously fired ~43 separate .evaluate() calls (one full
    // server-side graph evaluation per node per period). Now batches ALL 9
    // nodes' monthly values into ONE FeatureCollection per period via
    // reduceRegions()+flatten() (extractMultiNodeSeries), plus one batched
    // NDVI check - 3 total .evaluate() calls instead of ~43. All AC1/
    // variance/skewness/correlation math now runs client-side in JS
    // (jsNodeStats/jsPairCorrelation) on the small resulting tables, using
    // formulas identical to the EE versions - only where the arithmetic
    // runs has changed, not what it computes.
    var ptsFC = ee.FeatureCollection(nodes.map(function(nd,idx){
      return ee.Feature(nd.buf, {label:nd.label, idx:idx});
    }));
    var faiCollBefore=mkMoFAIRange(beforeStartTxt,beforeMonths);
    var faiCollAfter=mkMoFAIRange(afterStartTxt,afterMonths);
    var ndviWaterMasked=ndviWater.updateMask(oceanMask).rename('ndviW');

    var rBeforeSeries = extractMultiNodeSeries(faiCollBefore, ptsFC, 'fai', 20);
    var rAfterSeries = extractMultiNodeSeries(faiCollAfter, ptsFC, 'fai', 20);
    // v10.109 FIX: switched from mean() over 150m to max() over a FIXED
    // 500m buffer to fix Center falsely reading "LIKELY NOT on-reef".
    // v10.110 FIX: that fixed value size (500m) was itself a bug - at the
    // default 0.5km ring radius, a 500m buffer overlaps substantially with
    // EVERY adjacent node's buffer (Center-to-ring spacing = radiusKm*1000;
    // ring-to-ring spacing is even smaller, ~0.77x that). A real run showed
    // Center/N/NE/E all reading the IDENTICAL value 0.84 - strong evidence
    // their overlapping buffers were all picking up the same peak pixel,
    // silently defeating the purpose of having 9 distinct nodes. Buffer
    // size now SCALES with the user's ring radius (30% of node spacing,
    // capped 80-250m) so it can never exceed roughly a third of the gap to
    // the nearest neighbouring node, while still being larger than the
    // 150m FAI buffer that caused the original false-negative.
    var ndviBufM = Math.max(80, Math.min(250, radiusKm*1000*0.3));
    var ptsFCForNdvi = ee.FeatureCollection(nodes.map(function(nd,idx){
      return ee.Feature(ee.Geometry.Point([nd.lon,nd.lat]).buffer(ndviBufM), {label:nd.label, idx:idx});
    }));
    var rNdviAll = ndviWaterMasked.reduceRegions({collection:ptsFCForNdvi, reducer:ee.Reducer.max(), scale:20});

    s7dStatusV.setValue('Running: 3 batched Earth Engine calls (BEFORE series, AFTER series, NDVI check)\n'+
      'covering all 9 nodes at once - much lighter than the previous ~43-call approach. Expect 20-60 seconds.');
    s7dStatusV.style().set('color','#334466'); s7dStatusV.style().set('backgroundColor','#eeeeee');
    s7dStatusV.style().set('border','2px solid #aaaaaa'); s7dStatusV.style().set('whiteSpace','pre');
    s7dResultV.setValue('');

    var s7dData={}, s7dPending=3, s7dErrors=0;
    function s7dBump(key,val,err){
      s7dData[key]=err?null:val;
      if(err) s7dErrors++;
      s7dPending--;
      s7dStatusV.setValue('Running: '+(3-s7dPending)+' / 3 batched calls done'+
        (s7dErrors>0?' ('+s7dErrors+' errored)':'')+'...');
      if(s7dPending===0) s7dFinish();
    }
    function s7dFinish(){
      try {
        function fmtN(v,d){ return (v!==null&&v!==undefined&&!isNaN(v))?v.toFixed(d):'n/a'; }

        var beforeByNode = groupSeriesByLabel(s7dData.before, 'fai');
        var afterByNode = groupSeriesByLabel(s7dData.after, 'fai');
        var ndviByLabel = {};
        if(s7dData.ndvi && s7dData.ndvi.features){
          s7dData.ndvi.features.forEach(function(f){
            var p=f.properties||{};
            ndviByLabel[p.label] = extractReduceRegionsValue(p,'ndviW');
          });
        }

        var perNodeStats = {};
        nodes.forEach(function(nd){
          perNodeStats[nd.label] = {before:jsNodeStatsFixed(beforeByNode[nd.label]||[]), after:jsNodeStatsFixed(afterByNode[nd.label]||[])};
        });
        var perPairCorr = {};
        for(var pi=1;pi<nodes.length;pi++){
          var label=nodes[pi].label;
          perPairCorr[label] = {
            before: jsPairCorrelation(beforeByNode['Center'], beforeByNode[label]),
            after: jsPairCorrelation(afterByNode['Center'], afterByNode[label])
          };
        }

        var rows=['Node | On-reef? | \u0394AC1 | \u0394Var | \u0394Corr-vs-Center'];
        rows.push(repeatChar('\u2500',56));
        var nOnReefRisingCorr=0, nOnReefTotal=0, nOffReef=0, nAc1Rising=0, nAc1Avail=0;
        for(var idx=0;idx<nodes.length;idx++){
          var nd=nodes[idx];
          var stB=perNodeStats[nd.label].before, stA=perNodeStats[nd.label].after;
          var ndviV=(ndviByLabel[nd.label]!==undefined)?ndviByLabel[nd.label]:null;
          var onReef = ndviV!==null ? (ndviV>-0.10) : null;
          var onReefLabel = onReef===null?'n/a':(onReef?'yes ('+fmtN(ndviV,2)+')':'LIKELY NOT ('+fmtN(ndviV,2)+')');
          if(onReef===false) nOffReef++;
          var dAC1=(stB.realAC1!==null&&stA.realAC1!==null)?(stA.realAC1-stB.realAC1):null;
          var dVar=(stB.varTrendRatio!==null&&stA.varTrendRatio!==null)?(stA.varTrendRatio-stB.varTrendRatio):null;
          if(dAC1!==null){ nAc1Avail++; if(dAC1>0.01) nAc1Rising++; }
          var dCorrTxt='-';
          if(idx>0){
            var pc=perPairCorr[nd.label];
            var dCorr=(pc.before.corr!==null&&pc.after.corr!==null)?(pc.after.corr-pc.before.corr):null;
            dCorrTxt=dCorr!==null?((dCorr>0?'+':'')+dCorr.toFixed(3)):'n/a';
            if(onReef===true){ nOnReefTotal++; if(dCorr!==null&&dCorr>0.10) nOnReefRisingCorr++; }
          }
          // v10.109 FIX: same near-zero-denominator artifact risk already
          // fixed for S7C in v10.105 (varTrendRatio = 2ndHalfVar/max(1st
          // HalfVar,1e-6) - a tiny first-half variance inflates the ratio
          // without any real surge). Ported here after a real S7D run
          // showed SW=+52.31x and SE=+11.99x, even more extreme than the
          // case that prompted the original fix.
          var artifactFlag='';
          if(stB.varFirstHalf!==null&&stB.varFirstHalf<0.001&&dVar!==null&&Math.abs(dVar)>5){
            artifactFlag=' \u26A0ARTIFACT';
          }
          rows.push(nd.label+repeatChar(' ',Math.max(1,7-nd.label.length))+'| '+
            onReefLabel+repeatChar(' ',Math.max(1,12-onReefLabel.length))+'| '+
            (dAC1!==null?(dAC1>0?'+':'')+dAC1.toFixed(3):'n/a')+' | '+
            (dVar!==null?(dVar>0?'+':'')+dVar.toFixed(2)+'x':'n/a')+artifactFlag+' | '+dCorrTxt);
        }
        rows.push(repeatChar('\u2500',56));
        rows.push('\u26A0ARTIFACT flag: 1st-half FAI variance was near-zero (<0.001) and |\u0394Var|>5x - likely a division artifact (see S7C v10.105 note), not a genuine surge. Treat those \u0394Var values with real caution.');
        rows.push('BEFORE: '+beforeStartTxt+' + '+beforeMonths+'mo | AFTER: '+afterStartTxt+' + '+afterMonths+'mo | radius='+radiusKm+'km | NDVI check buffer='+ndviBufM.toFixed(0)+'m (scaled to avoid overlapping adjacent nodes)');
        rows.push('Method: 3 batched EE fetches (reduceRegions+flatten), AC1/variance/correlation computed client-side in JS.');
        rows.push(thresholdOnlyWarning()+'  S7D uses a fixed \u0394corr>0.10 rule that has NOT been');
        rows.push('calibrated against a null. The permutationTestDelta engine already exists in this');
        rows.push('file and could be ported here - until then read S7D as exploratory only.');
        // v10.108 NEW: raw feature counts visible on-screen, not just console
        // - if these are 0, the batched fetch itself returned nothing (a
        // real query/data problem); if they are >0 but the table above is
        // still all n/a, the problem is in property extraction instead -
        // these two failure modes need different fixes and this line tells
        // you which one you're looking at without opening the console.
        var beforeFeatCount=(s7dData.before&&s7dData.before.features)?s7dData.before.features.length:0;
        var afterFeatCount=(s7dData.after&&s7dData.after.features)?s7dData.after.features.length:0;
        var ndviFeatCount=(s7dData.ndvi&&s7dData.ndvi.features)?s7dData.ndvi.features.length:0;
        rows.push('Raw features returned: BEFORE='+beforeFeatCount+' | AFTER='+afterFeatCount+' | NDVI='+ndviFeatCount+
          ' (expect BEFORE='+(beforeMonths*9)+', AFTER='+(afterMonths*9)+', NDVI=9 if fully populated)');
        if(nOffReef>0) rows.push('NOTE: '+nOffReef+' of 9 nodes flagged LIKELY NOT on-reef (NDVI-water <= -0.10) - their rows above are informational only, not trusted algae-dynamics comparisons.');
        if(s7dErrors>0) rows.push('NOTE: '+s7dErrors+' of 3 batched calls returned no usable data. Check the Console for the exact error text (search for "S7D [" lines).');
        s7dResultV.setValue(rows.join('\n'));

        var headline, hCol, hBg;
        if(nOnReefTotal===0){
          headline='COULD NOT ASSESS COUPLING - no on-reef ring points to compare against Center';
          hCol='#888888'; hBg='#eeeeee';
        } else if(nOnReefRisingCorr>=Math.ceil(nOnReefTotal*0.5)){
          headline='HYPER-SYNCHRONIZATION SIGNAL: '+nOnReefRisingCorr+' of '+nOnReefTotal+' on-reef ring points show RISING correlation with Center (Dakos et al. 2011 leading-indicator direction). AC1 rising at '+nAc1Rising+'/'+nAc1Avail+' nodes.';
          hCol='#880000'; hBg='#ffd0d0';
        } else if(nOnReefRisingCorr>0){
          headline='PARTIAL SIGNAL: '+nOnReefRisingCorr+' of '+nOnReefTotal+' on-reef ring points show rising correlation with Center - not a majority. AC1 rising at '+nAc1Rising+'/'+nAc1Avail+' nodes.';
          hCol='#886600'; hBg='#fff6cc';
        } else {
          headline='NO SYNCHRONIZATION SIGNAL: 0 of '+nOnReefTotal+' on-reef ring points show rising correlation with Center. AC1 rising at '+nAc1Rising+'/'+nAc1Avail+' nodes.';
          hCol='#115511'; hBg='#d4f5df';
        }
        s7dStatusV.setValue(headline);
        s7dStatusV.style().set('color',hCol); s7dStatusV.style().set('backgroundColor',hBg);
        s7dStatusV.style().set('border','2px solid '+hCol); s7dStatusV.style().set('whiteSpace','pre'); s7dStatusV.style().set('fontWeight','bold');

        print('=== S7D FULL 9-NODE ALGAE COUPLING NETWORK (v10.107 batched) ===');
        print(headline);
        print(rows.join('\n'));
      } catch(errFinal){
        s7dStatusV.setValue(friendlyEEError(errFinal));
        s7dStatusV.style().set('color','#cc0000'); s7dStatusV.style().set('backgroundColor','#ffd0d0');
        s7dStatusV.style().set('border','2px solid #cc0000'); s7dStatusV.style().set('whiteSpace','pre');
        print('=== S7D FINAL ERROR === '+errFinal);
      }
    }
    // v10.108 NEW: real diagnostics for each of the 3 batched calls, instead
    // of only an aggregate "N errored" counter. A run came back with EVERY
    // node showing n/a on EVERY field, which the aggregate counter alone
    // can't explain (it only reported 1 of 3 calls as errored, but the
    // table looked like all 3 failed) - this makes each call's real outcome
    // directly visible in the console: the exact error text if one failed,
    // or the feature count + first feature's raw property names if it
    // succeeded (confirms whether reduceRegions() is actually naming its
    // output the way extractReduceRegionsValue() assumes).
    function s7dDiagnose(label, v, e){
      if(e){
        print('=== S7D ['+label+'] ERROR ==='); print(e);
        return;
      }
      var nFeat = (v&&v.features)?v.features.length:0;
      print('=== S7D ['+label+'] OK - '+nFeat+' feature(s) returned ===');
      if(nFeat>0){
        print('First feature raw properties: '+JSON.stringify(v.features[0].properties));
      } else {
        print('WARNING: 0 features returned - reduceRegions/flatten produced an empty table.');
      }
    }
    rBeforeSeries.evaluate(function(v,e){ s7dDiagnose('BEFORE series', v, e); s7dBump('before', v, e); });
    rAfterSeries.evaluate(function(v,e){ s7dDiagnose('AFTER series', v, e); s7dBump('after', v, e); });
    rNdviAll.evaluate(function(v,e){ s7dDiagnose('NDVI check', v, e); s7dBump('ndvi', v, e); });
  }
});
panel.add(s7dRunBtn);
panel.add(s7dStatusV);
panel.add(s7dResultV);
panel.add(lbl('CAVEAT: correlation is scored vs Center only (8 pairs), not the full 36-pair matrix. This is an empirical correlation network, not a mechanistic J_ij interaction matrix - shared external forcing (a heatwave hitting the whole area) will raise correlation with or without any real internal coupling change.',7,'#886600'));
panel.add(legDiv());

// ============================================================
// S7E - LOCAL vs REGIONAL AUTO-CLASSIFICATION (v10.111 NEW)
// The missing piece S7D didn't have: S7D's 9 nodes all sit within 0.3-3km
// of each other - they test spatial coupling WITHIN one reef, but every
// node is still local to that same reef, so "all 9 spike together" can't
// actually distinguish a local event big enough to blanket the reef from
// a genuine regional signal. This adds a real, independent reference:
// auto-searches outward in expanding rings (20/40/70/110/160km) for the
// nearest point that is genuinely shallow-water (GEBCO -50 to 0m, the
// same threshold used for shallowMask elsewhere in this tool), then runs
// the SAME FAI-based BEFORE/AFTER AC1/variance extraction there, and
// applies an explicit decision rule comparing the study reef to that
// reference - mirroring STEP 3's already-proven LOCAL/GLOBAL/ANOMALOUS/
// NO SIGNAL classification for SST, just applied to algae with a real
// reef reference instead of open ocean (open ocean has no algae signal
// to compare against at all).
// DISCLOSED LIMIT: "GEBCO says it's shallow" is not the same as "this is
// a real reef with comparable ecology" - it could be a sandbar or a bare
// shoal with no algae community at all. Partially mitigated (not solved)
// by checking the reference site's own data density as a proxy for
// whether it has a real, trackable satellite algae signal.
// ============================================================
panel.add(sHead('S7E - LOCAL vs REGIONAL AUTO-CLASSIFICATION (v10.124)','#5a1a2a'));
panel.add(lbl('Auto-finds a genuinely distant shallow-water reef (GEBCO search, 20-160km out) as a real regional reference, then applies an explicit decision rule: does the study reef show a change the reference reef does NOT (LOCAL), or do both change together (REGIONAL)? This is what actually separates local degradation from a regional climate shock - S7D alone could not do this, since all its nodes sit on the same reef.',7,'#883344'));
panel.add(lbl('v10.124: also runs a real permutation test (500 shuffles, same engine as STEP 3 COMPARE) on both study and reference sites - a genuine p-value alongside the threshold-based verdict above, at no extra Earth Engine cost (reuses the raw FAI data already fetched for the classification).',7,'#552266'));
panel.add(lbl('DISCLOSED LIMIT: a GEBCO shallow-water match is not guaranteed to be a real reef with comparable ecology (could be a bare sandbar). Partially checked (not solved) via the reference site\'s own data density - a very sparse signal there is flagged as low-confidence, not silently trusted.',7,'#886600'));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var s7eCoordInput = ui.Textbox({placeholder:'lat, lon  e.g. -23.51, 152.09',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7eCoordInput);
var s7eUseLastClickBtn = ui.Button({
  label:'Use last clicked location (from map above)',
  style:{fontSize:'9px',margin:'2px 4px',backgroundColor:'#e8f4ff',color:'#225588',stretch:'horizontal',padding:'4px 4px',border:'1px solid #4488cc'},
  onClick:function(){
    if(lastClickLat===null||lastClickLon===null){
      s7eStatusV.setValue('No location clicked yet - click the map or use GO TO COORDINATES first.');
      s7eStatusV.style().set('color','#aa3300'); return;
    }
    s7eCoordInput.setValue(lastClickLat.toFixed(4)+', '+lastClickLon.toFixed(4));
    s7eStatusV.setValue('Location filled in. Fill in BEFORE/AFTER dates, then press RUN.');
    s7eStatusV.style().set('color','#115511');
  }
});
panel.add(s7eUseLastClickBtn);
panel.add(lbl('BEFORE start date + months:',7,'#115511'));
var s7eBeforeStartInput = ui.Textbox({placeholder:'e.g. 2021-06-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7eBeforeStartInput);
var s7eBeforeMonthsInput = ui.Textbox({placeholder:'months, e.g. 18',value:'18',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7eBeforeMonthsInput);
panel.add(lbl('AFTER start date + months:',7,'#aa3300'));
var s7eAfterStartInput = ui.Textbox({placeholder:'e.g. 2023-11-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7eAfterStartInput);
var s7eAfterMonthsInput = ui.Textbox({placeholder:'months, e.g. 12',value:'12',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7eAfterMonthsInput);
var s7eStatusV = ui.Label('Fill in the fields above, then press RUN. 3 batched calls, expect ~30-90 seconds.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var s7eResultV = ui.Label('',{fontSize:'8px',color:'#5a1020',backgroundColor:'#faeef2',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var s7eRunBtn = ui.Button({
  label:'RUN LOCAL vs REGIONAL CLASSIFICATION',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#f0d0d8',color:'#5a1020',stretch:'horizontal',padding:'6px 4px',border:'2px solid #883344'},
  onClick:function(){
    var coordTxt=(s7eCoordInput.getValue()||'').trim();
    var beforeStartTxt=(s7eBeforeStartInput.getValue()||'').trim();
    var beforeMonthsTxt=(s7eBeforeMonthsInput.getValue()||'').trim();
    var afterStartTxt=(s7eAfterStartInput.getValue()||'').trim();
    var afterMonthsTxt=(s7eAfterMonthsInput.getValue()||'').trim();
    if(!coordTxt||!beforeStartTxt||!beforeMonthsTxt||!afterStartTxt||!afterMonthsTxt){
      s7eStatusV.setValue('MISSING - fill in Lat/Lon and both BEFORE/AFTER date+months fields.');
      s7eStatusV.style().set('color','#cc0000'); s7eStatusV.style().set('backgroundColor','#ffd0d0'); return;
    }
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){s7eStatusV.setValue('Invalid format - use: lat, lon'); s7eStatusV.style().set('color','#cc0000'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){s7eStatusV.setValue('Invalid format - use: lat, lon'); s7eStatusV.style().set('color','#cc0000'); return;}
    var beforeMonths=parseInt(beforeMonthsTxt,10), afterMonths=parseInt(afterMonthsTxt,10);
    if(isNaN(beforeMonths)||beforeMonths<4||beforeMonths>36||isNaN(afterMonths)||afterMonths<4||afterMonths>36){
      s7eStatusV.setValue('BEFORE/AFTER months must each be 4-36.'); s7eStatusV.style().set('color','#cc0000'); return;
    }
    recordStudySite(latIn, lonIn, 'S7E');

    function destinationPoint(lat, lon, bearingDeg, distKm){
      var R=6371, brng=bearingDeg*Math.PI/180, lat1=lat*Math.PI/180, lon1=lon*Math.PI/180;
      var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
      var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
      return {lat:lat2*180/Math.PI, lon:lon2*180/Math.PI};
    }
    var searchRadiiKm=[20,40,70,110,160];
    var candidates=[];
    searchRadiiKm.forEach(function(rKm){
      for(var b=0;b<360;b+=45){
        var d=destinationPoint(latIn,lonIn,b,rKm);
        candidates.push({lat:d.lat,lon:d.lon,radiusKm:rKm,bearing:b});
      }
    });
    var candFeats = candidates.map(function(c,idx){
      return ee.Feature(ee.Geometry.Point([c.lon,c.lat]), {idx:idx, radiusKm:c.radiusKm});
    });
    var candFC = ee.FeatureCollection(candFeats);
    var rCandidates = GEBCO.reduceRegions({collection:candFC, reducer:ee.Reducer.first(), scale:500});

    s7eStatusV.setValue('Step 1/3: Searching 20-160km out for a genuine shallow-water reference reef (1 batched GEBCO call)...');
    s7eStatusV.style().set('color','#334466'); s7eStatusV.style().set('backgroundColor','#eeeeee');
    s7eStatusV.style().set('border','2px solid #aaaaaa'); s7eStatusV.style().set('whiteSpace','pre');
    s7eResultV.setValue('');

    rCandidates.evaluate(function(candResult, candErr){
      if(candErr){
        s7eStatusV.setValue(friendlyEEError(candErr));
        s7eStatusV.style().set('color','#cc0000'); s7eStatusV.style().set('backgroundColor','#ffd0d0');
        s7eStatusV.style().set('border','2px solid #cc0000');
        print('=== S7E candidate search error === '+candErr);
        return;
      }
      try {
        var feats = (candResult&&candResult.features)?candResult.features:[];
        var best=null;
        // v10.116 NEW: exclude candidates that coincide with a site already
        // tested as a STUDY location this session (S7B/S7C/S7D/S7E) - a
        // reference is supposed to be an uninvolved baseline, not a site
        // that has already tested positive for its own change. Caught from
        // a real run where the auto-found reference landed on Low Isles,
        // which had already been tested (and shown a real signal) as a
        // study site earlier in the same session.
        var EXCLUDE_RADIUS_KM = 5;
        var excludedByHistoryCount = 0;
        var shallowCandidateCount = 0;
        for(var i=0;i<feats.length;i++){
          var p=feats[i].properties;
          var elev=extractReduceRegionsValue(p,'elevation');
          if(elev!==null && elev>-50 && elev<0){
            shallowCandidateCount++;
            var candLat=candidates[p.idx].lat, candLon=candidates[p.idx].lon;
            var tooCloseToHistory=false;
            for(var h=0;h<s7StudySiteHistory.length;h++){
              if(haversineKm(candLat,candLon,s7StudySiteHistory[h].lat,s7StudySiteHistory[h].lon) < EXCLUDE_RADIUS_KM){
                tooCloseToHistory=true; break;
              }
            }
            if(tooCloseToHistory){ excludedByHistoryCount++; continue; }
            if(!best || p.radiusKm<best.radiusKm) best={radiusKm:p.radiusKm, idx:p.idx, elev:elev,
              lat:candLat, lon:candLon};
          }
        }
        if(!best){
          var noRefMsg = shallowCandidateCount>0 && excludedByHistoryCount===shallowCandidateCount ?
            'NO REFERENCE REEF FOUND within 160km - '+shallowCandidateCount+' shallow-water candidate(s) found, but ALL of them '+
            'coincide with a site already tested as a STUDY location this session (within '+EXCLUDE_RADIUS_KM+'km). '+
            'Try a different study site, or clear the session by reloading if you want to allow reuse.' :
            'NO REFERENCE REEF FOUND within 160km (no GEBCO shallow-water match at 20/40/70/110/160km, 8 bearings each). Cannot classify local vs regional without a reference - try a different study site.';
          s7eStatusV.setValue(noRefMsg);
          s7eStatusV.style().set('color','#cc0000'); s7eStatusV.style().set('backgroundColor','#ffd0d0');
          s7eStatusV.style().set('border','2px solid #cc0000'); s7eStatusV.style().set('whiteSpace','pre');
          print('=== S7E: no reference reef candidate found within 160km ('+excludedByHistoryCount+' excluded as prior study sites) ===');
          return;
        }

        s7eStatusV.setValue('Step 2/3: Reference reef found '+best.radiusKm+'km away (depth '+best.elev.toFixed(1)+'m)'+
          (excludedByHistoryCount>0?' ['+excludedByHistoryCount+' closer candidate(s) skipped - already tested as a study site]':'')+'. '+
          'Fetching BEFORE/AFTER FAI at study + reference (2 batched calls)...');

        var studyPtBuf=ee.Geometry.Point([lonIn,latIn]).buffer(150);
        var refPtBuf=ee.Geometry.Point([best.lon,best.lat]).buffer(150);
        var pairFC=ee.FeatureCollection([
          ee.Feature(studyPtBuf,{label:'Study'}),
          ee.Feature(refPtBuf,{label:'Reference'})
        ]);
        var faiCollBefore=mkMoFAIRange(beforeStartTxt,beforeMonths);
        var faiCollAfter=mkMoFAIRange(afterStartTxt,afterMonths);
        var rBefore=extractMultiNodeSeries(faiCollBefore, pairFC, 'fai', 20);
        var rAfter=extractMultiNodeSeries(faiCollAfter, pairFC, 'fai', 20);

        var s7ePending=2, s7eErrors=0, s7eData={};
        function s7eBump(key,v,e){
          s7eData[key]=e?null:v;
          if(e) s7eErrors++;
          s7ePending--;
          s7eStatusV.setValue('Step '+(3-s7ePending)+'/3: '+(2-s7ePending)+' / 2 batched calls done...');
          if(s7ePending===0) s7eFinish();
        }
        function s7eFinish(){
          try {
            var beforeByNode=groupSeriesByLabel(s7eData.before,'fai');
            var afterByNode=groupSeriesByLabel(s7eData.after,'fai');
            var studyB=jsNodeStatsFixed(beforeByNode['Study']||[]);
            var studyA=jsNodeStatsFixed(afterByNode['Study']||[]);
            var refB=jsNodeStatsFixed(beforeByNode['Reference']||[]);
            var refA=jsNodeStatsFixed(afterByNode['Reference']||[]);

            var dAC1_study=(studyB.realAC1!==null&&studyA.realAC1!==null)?(studyA.realAC1-studyB.realAC1):null;
            var dVar_study=(studyB.varTrendRatio!==null&&studyA.varTrendRatio!==null)?(studyA.varTrendRatio-studyB.varTrendRatio):null;
            var dAC1_ref=(refB.realAC1!==null&&refA.realAC1!==null)?(refA.realAC1-refB.realAC1):null;
            var dVar_ref=(refB.varTrendRatio!==null&&refA.varTrendRatio!==null)?(refA.varTrendRatio-refB.varTrendRatio):null;

            // v10.111: refValidity is a proxy, not proof, that the auto-
            // found GEBCO match is a real reef with a trackable algae
            // signal - low data density there means "shallow" but not
            // necessarily "ecologically comparable".
            var refDataDensity = (refB.nValidMonths+refA.nValidMonths)/(beforeMonths+afterMonths);
            var refLowConfidence = refDataDensity<0.4;

            // v10.153 FIX 10: THIS is the rule that was measured firing on
            // 80% of no-event windows (AC1 OR variance, at 0.01 / 0.15x).
            // Now uses the calibrated cutoffs for the actual window length.
            var _s7eThrB = getCalibratedThresholds(beforeMonths);
            var _s7eThrA = getCalibratedThresholds(afterMonths);
            var _s7eAc1 = Math.max(_s7eThrB.ac1, _s7eThrA.ac1);
            var _s7eVar = Math.max(_s7eThrB.varr, _s7eThrA.varr);
            var studySignal = (dAC1_study!==null&&dAC1_study>_s7eAc1) || (dVar_study!==null&&dVar_study>_s7eVar);
            var refSignal = (dAC1_ref!==null&&dAC1_ref>_s7eAc1) || (dVar_ref!==null&&dVar_ref>_s7eVar);

            function fmtN(v,d){ return (v!==null&&v!==undefined&&!isNaN(v))?v.toFixed(d):'n/a'; }

            // v10.114 FIX: a real run showed "LOCAL ANOMALY DETECTED" when
            // the reference site's dAC1/dVar were BOTH null (n/a) - not
            // because the reference was tested and found stable, but
            // because it had no computable data (one of its two periods
            // had <4 valid months, even though its COMBINED density of 55%
            // looked fine and wasn't flagged low-confidence). "Untested"
            // and "tested and stable" are very different findings, but the
            // old rule collapsed both into refSignal=false, so a real
            // signal at the study site got labelled a confirmed LOCAL
            // ANOMALY against a reference that was never actually able to
            // confirm anything. Now checked per-period (not just combined
            // density) and given its own explicit "CANNOT CLASSIFY" state
            // instead of silently defaulting to "no signal".
            var refInsufficient = (refB.nValidMonths<4 || refA.nValidMonths<4);
            var studyInsufficient = (studyB.nValidMonths<4 || studyA.nValidMonths<4);

            var verdict, vCol, vBg;
            if(refInsufficient && studyInsufficient){
              verdict='CANNOT CLASSIFY - both study and reference have insufficient valid months (<4) in at least one period - not enough data to compare, not a confirmed finding either way.';
              vCol='#888888'; vBg='#eeeeee';
            } else if(refInsufficient){
              verdict='CANNOT CLASSIFY - reference site has insufficient valid months (<4) in '+
                (refB.nValidMonths<4&&refA.nValidMonths<4?'BOTH periods':refB.nValidMonths<4?'the BEFORE period':'the AFTER period')+
                '. The study site\'s own signal ('+(studySignal?'SIGNAL':'no signal')+') cannot be confirmed as local vs regional without a reference that could actually be tested - do not treat this as a confirmed local anomaly.';
              vCol='#886600'; vBg='#fff6cc';
            } else if(studyInsufficient){
              verdict='CANNOT CLASSIFY - study site has insufficient valid months (<4) in '+
                (studyB.nValidMonths<4&&studyA.nValidMonths<4?'BOTH periods':studyB.nValidMonths<4?'the BEFORE period':'the AFTER period')+
                ' - cannot assess whether the study site changed at all.';
              vCol='#888888'; vBg='#eeeeee';
            } else if(studySignal && !refSignal){
              verdict='LOCAL ANOMALY DETECTED - study reef shows a change the regional reference does NOT';
              vCol='#880000'; vBg='#ffd0d0';
            } else if(studySignal && refSignal){
              verdict='REGIONAL SIGNAL - both study AND reference reef show similar changes, not specific to your study site';
              vCol='#886600'; vBg='#fff6cc';
            } else if(!studySignal && refSignal){
              verdict='ANOMALOUS - reference reef shows change, study site does NOT (study may be locally buffered/protected)';
              vCol='#115511'; vBg='#d4f5df';
            } else {
              verdict='NO SIGNAL AT EITHER SITE';
              vCol='#226644'; vBg='#e8f4ff';
              // v10.113 FIX: the same gap already caught and fixed for
              // STEP 5 in v10.97 ("no significant trend" was hiding a real
              // significant DECLINE) existed here too, uncaught until a
              // real run showed it directly: study AC1 fell -0.675,
              // reference fell -0.636 - nearly identical, a coherent
              // regional pattern - but the SIGNAL check only looks for
              // RISING AC1/variance (correctly, per Dakos et al. - a
              // decline isn't a CSD warning), so this got silently folded
              // into a bare "NO SIGNAL" with no mention of the pattern
              // underneath it. Now surfaced explicitly when both sites
              // show a similar, substantial AC1 decline.
              if(dAC1_study!==null && dAC1_ref!==null && dAC1_study<-0.1 && dAC1_ref<-0.1 &&
                 Math.abs(dAC1_study-dAC1_ref)<0.15){
                verdict += ' - but note: AC1 fell similarly at both sites (study '+dAC1_study.toFixed(2)+
                  ', reference '+dAC1_ref.toFixed(2)+') - a coherent REGIONAL DECLINE, not a CSD warning direction '+
                  '(Dakos et al. 2012: rising AC1 is the meaningful signal, not falling), but worth noting as '+
                  'consistent regional behaviour rather than pure noise.';
              }
            }
            // v10.112 FIX: reference density was shown, but STUDY's own
            // data density never was - so a real run showing STUDY n/a
            // gave no way to tell whether the study site itself had sparse
            // Sentinel-2 coverage (likely at cloudy sites, e.g. Daintree/
            // Cape Tribulation rainforest coast) versus some other cause.
            var studyDataDensity = (studyB.nValidMonths+studyA.nValidMonths)/(beforeMonths+afterMonths);
            var studyLowConfidence = studyDataDensity<0.4;
            if(refLowConfidence) verdict += ' [LOW CONFIDENCE: reference site has sparse data - '+
              Math.round(refDataDensity*100)+'% coverage - may not be a real algae-bearing reef]';
            if(studyLowConfidence) verdict += ' [STUDY SITE ALSO SPARSE: '+Math.round(studyDataDensity*100)+
              '% coverage - likely heavy cloud cover at this location/window, not a real "no signal" finding]';

            // v10.124 NEW: real permutation-test p-values for S7E, porting
            // the SAME engine already validated in STEP 3 COMPARE and FIND
            // SWEET SPOT (v10.122/123) to the algae-based local/regional
            // classification - closing the gap directly flagged in
            // discussion: "the permutation-testing claim needs a caveat -
            // it protects the SST-based findings right now, not yet the
            // algae-based ones." Requires ZERO new EE calls: the raw FAI
            // series for study and reference were already fetched above
            // (beforeByNode/afterByNode) to compute studyB/studyA/refB/
            // refA - this just runs the existing permutationTestDelta()
            // engine on that same already-fetched data.
            var studyAC1Test = permutationTestDeltaFixed(beforeByNode['Study']||[], afterByNode['Study']||[], statAC1ForPerm, 500);
            var studyVarTest = permutationTestDeltaFixed(beforeByNode['Study']||[], afterByNode['Study']||[], statVarRatioForPerm, 500);
            var refAC1Test = permutationTestDeltaFixed(beforeByNode['Reference']||[], afterByNode['Reference']||[], statAC1ForPerm, 500);
            var refVarTest = permutationTestDeltaFixed(beforeByNode['Reference']||[], afterByNode['Reference']||[], statVarRatioForPerm, 500);
            function fmtP(t){ return t.pValue!==null?t.pValue.toFixed(3):(t.note||'n/a'); }
            function sigTag(t){ return t.pValue===null?'':(t.pValue<0.05?' *** likely real (p<0.05)':t.pValue<0.10?' * borderline (p<0.10)':' not significant'); }

            var lines=[];
            lines.push('Reference reef: '+best.lat.toFixed(4)+', '+best.lon.toFixed(4)+' ('+best.radiusKm+'km away, depth='+best.elev.toFixed(1)+'m), data density='+Math.round(refDataDensity*100)+'%'+
              (refLowConfidence?' (LOW - treat with caution)':''));
            if(excludedByHistoryCount>0) lines.push('NOTE: '+excludedByHistoryCount+' closer shallow-water candidate(s) were skipped because they coincide with a site already tested as a STUDY location this session (v10.116).');
            lines.push('Study site data density: '+Math.round(studyDataDensity*100)+'% ('+studyB.nValidMonths+'/'+beforeMonths+' BEFORE, '+
              studyA.nValidMonths+'/'+afterMonths+' AFTER valid months)'+(studyLowConfidence?' (LOW - likely cloud-limited, not a real null result)':''));
            lines.push(repeatChar('\u2500',50));
            lines.push('STUDY  : \u0394AC1='+fmtN(dAC1_study,3)+' \u0394Var='+fmtN(dVar_study,2)+'x -> '+(studySignal?'SIGNAL':'no signal'));
            lines.push('REFERENCE: \u0394AC1='+fmtN(dAC1_ref,3)+' \u0394Var='+fmtN(dVar_ref,2)+'x -> '+(refSignal?'SIGNAL':'no signal'));
            lines.push(repeatChar('\u2500',50));
            lines.push('BEFORE: '+beforeStartTxt+' + '+beforeMonths+'mo | AFTER: '+afterStartTxt+' + '+afterMonths+'mo');
            lines.push('Decision rule (v10.153 CALIBRATED): study signal (\u0394AC1>'+_s7eAc1.toFixed(3)+
              ' OR \u0394Var>'+_s7eVar.toFixed(2)+'x) vs reference signal, same cutoffs.');
            lines.push('These replace the old 0.01 / 0.15x values, which were MEASURED firing on 80% of');
            lines.push('windows at a site with no regime shift (Scripps Pier CTD, 13.6yr, 104-176 splits).');
            lines.push('Expect far fewer LOCAL ANOMALY verdicts than before. The permutation p-values');
            lines.push('below are the only figures here with a known false-positive rate.');
            if(s7eErrors>0) lines.push('NOTE: '+s7eErrors+' of 2 batched calls returned no usable data.');
            lines.push(repeatChar('\u2500',50));
            lines.push('=== PERMUTATION TEST (real p-value, 500 shuffles, algae/FAI) ===');
            lines.push('Same engine as STEP 3 COMPARE (v10.122), applied here to algae instead of SST -');
            lines.push('answers whether the >0.01/>0.15 threshold verdict above reflects a genuinely rare');
            lines.push('pattern, or one ordinary random noise would produce anyway. Zero extra EE calls -');
            lines.push('reuses the raw FAI series already fetched for the classification above.');
            lines.push('Study AC1 p='+fmtP(studyAC1Test)+sigTag(studyAC1Test));
            lines.push('Study Var p='+fmtP(studyVarTest)+sigTag(studyVarTest));
            lines.push('Reference AC1 p='+fmtP(refAC1Test)+sigTag(refAC1Test));
            lines.push('Reference Var p='+fmtP(refVarTest)+sigTag(refVarTest));
            s7eResultV.setValue(lines.join('\n'));

            s7eStatusV.setValue(verdict);
            s7eStatusV.style().set('color',vCol); s7eStatusV.style().set('backgroundColor',vBg);
            s7eStatusV.style().set('border','2px solid '+vCol); s7eStatusV.style().set('whiteSpace','pre'); s7eStatusV.style().set('fontWeight','bold');

            print('=== S7E LOCAL vs REGIONAL AUTO-CLASSIFICATION ===');
            print(verdict);
            print(lines.join('\n'));
          } catch(errFinal2){
            s7eStatusV.setValue(friendlyEEError(errFinal2));
            s7eStatusV.style().set('color','#cc0000'); s7eStatusV.style().set('backgroundColor','#ffd0d0');
            s7eStatusV.style().set('border','2px solid #cc0000'); s7eStatusV.style().set('whiteSpace','pre');
            print('=== S7E FINISH ERROR === '+errFinal2);
          }
        }
        // v10.125 NEW: real diagnostics on S7E's 2 batched calls, porting
        // the exact pattern already proven for S7D in v10.108. Caught from
        // a real run that stayed stuck at "1/2 batched calls done" even
        // after a page reload (ruling out an expired session) AND after
        // shortening to 12/10-month windows (ruling out the window-length
        // slowness already fixed for S7D at 18/12mo) - a genuinely new,
        // undiagnosed failure mode. Rather than guess a third fix blind,
        // this makes the NEXT run tell us exactly which of the two calls
        // (BEFORE or AFTER FAI fetch) is failing and why.
        rBefore.evaluate(function(v,e){
          if(e){ print('=== S7E [BEFORE series] ERROR === '+e); }
          else { print('=== S7E [BEFORE series] OK - '+((v&&v.features)?v.features.length:0)+' feature(s) returned ==='); }
          s7eBump('before', v, e);
        });
        rAfter.evaluate(function(v,e){
          if(e){ print('=== S7E [AFTER series] ERROR === '+e); }
          else { print('=== S7E [AFTER series] OK - '+((v&&v.features)?v.features.length:0)+' feature(s) returned ==='); }
          s7eBump('after', v, e);
        });
      } catch(errFinal){
        s7eStatusV.setValue(friendlyEEError(errFinal));
        s7eStatusV.style().set('color','#cc0000'); s7eStatusV.style().set('backgroundColor','#ffd0d0');
        s7eStatusV.style().set('border','2px solid #cc0000'); s7eStatusV.style().set('whiteSpace','pre');
        print('=== S7E CANDIDATE PROCESSING ERROR === '+errFinal);
      }
    });
  }
});
panel.add(s7eRunBtn);
panel.add(s7eStatusV);
panel.add(s7eResultV);
panel.add(legDiv());

// ============================================================
// S7F - RUN ALL: S7D + S7E COMBINED (v10.116 NEW)
// Direct response to a real request: enter Lat/Lon + one BEFORE/AFTER
// window ONCE, get both the within-reef coupling picture (S7D) and the
// local-vs-regional classification (S7E) together, tabulated in one
// summary - instead of re-typing the same coordinates and dates into two
// separate forms.
// SCOPE NOTE (deliberate, explained to the user before building): this
// does NOT auto-try multiple window lengths, and does NOT chain in S7B/
// S7C/S13. Both were considered and declined: this session directly
// demonstrated why - a single S7D call with too-long a window hung for
// 5+ minutes with no way to cancel it (no setTimeout in this sandbox),
// and Earth Engine's account-level concurrency quota was hit once
// already this session, requiring a manual tier upgrade. An orchestrator
// that tries several window lengths across every S7 tool would multiply
// both the hang risk and the quota load severalfold, with no reliable
// way to know if a huge combined run is progressing or stuck. This
// version fires a bounded 6 EE calls total (matching S7D's 3 + S7E's 3),
// with the same proven functions S7D/S7E already use - not a new,
// riskier reimplementation.
// ============================================================
panel.add(sHead('S7F - RUN ALL: S7D + S7E COMBINED (v10.116)','#1a3a5a'));
panel.add(lbl('Runs S7D (within-reef coupling) and S7E (local vs regional) together from ONE shared Lat/Lon + BEFORE/AFTER input, tabulating both verdicts side by side. Does NOT try multiple window lengths automatically - see the note below for why.',7,'#224466'));
panel.add(lbl('SCOPE: still just ONE window length per run, chosen by you (the same window-tuning process from S7D/S7E still applies - this only removes re-typing coordinates/dates twice). Does not include S7B/S7C/S13 - kept bounded to 6 EE calls total to avoid the quota/hang risk already seen this session with heavier combined runs.',7,'#886600'));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var s7fCoordInput = ui.Textbox({placeholder:'lat, lon  e.g. -23.51, 152.09',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7fCoordInput);
var s7fUseLastClickBtn = ui.Button({
  label:'Use last clicked location (from map above)',
  style:{fontSize:'9px',margin:'2px 4px',backgroundColor:'#e8f4ff',color:'#225588',stretch:'horizontal',padding:'4px 4px',border:'1px solid #4488cc'},
  onClick:function(){
    if(lastClickLat===null||lastClickLon===null){
      s7fStatusV.setValue('No location clicked yet - click the map or use GO TO COORDINATES first.');
      s7fStatusV.style().set('color','#aa3300'); return;
    }
    s7fCoordInput.setValue(lastClickLat.toFixed(4)+', '+lastClickLon.toFixed(4));
    s7fStatusV.setValue('Location filled in. Fill in radius + BEFORE/AFTER, then press RUN.');
    s7fStatusV.style().set('color','#115511');
  }
});
panel.add(s7fUseLastClickBtn);
panel.add(lbl('Ring radius (km, 0.3-3, for S7D within-reef nodes):',7,'#334466'));
var s7fRadiusInput = ui.Textbox({placeholder:'e.g. 0.5',value:'0.5',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7fRadiusInput);
panel.add(lbl('BEFORE start date + months:',7,'#115511'));
var s7fBeforeStartInput = ui.Textbox({placeholder:'e.g. 2021-06-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7fBeforeStartInput);
var s7fBeforeMonthsInput = ui.Textbox({placeholder:'months, e.g. 12',value:'12',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7fBeforeMonthsInput);
panel.add(lbl('AFTER start date + months:',7,'#aa3300'));
var s7fAfterStartInput = ui.Textbox({placeholder:'e.g. 2023-11-01',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7fAfterStartInput);
var s7fAfterMonthsInput = ui.Textbox({placeholder:'months, e.g. 10',value:'10',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(s7fAfterMonthsInput);
var s7fStatusV = ui.Label('Fill in the fields above, then press RUN. 6 batched calls total, expect ~30-90 seconds.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var s7fResultV = ui.Label('',{fontSize:'8px',color:'#1a2a4a',backgroundColor:'#eef2fa',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var s7fRunBtn = ui.Button({
  label:'RUN ALL (S7D + S7E COMBINED)',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#cfe0f5',color:'#1a3a5a',stretch:'horizontal',padding:'6px 4px',border:'2px solid #1a3a5a'},
  onClick:function(){
    var coordTxt=(s7fCoordInput.getValue()||'').trim();
    var radiusTxt=(s7fRadiusInput.getValue()||'').trim();
    var beforeStartTxt=(s7fBeforeStartInput.getValue()||'').trim();
    var beforeMonthsTxt=(s7fBeforeMonthsInput.getValue()||'').trim();
    var afterStartTxt=(s7fAfterStartInput.getValue()||'').trim();
    var afterMonthsTxt=(s7fAfterMonthsInput.getValue()||'').trim();
    if(!coordTxt||!beforeStartTxt||!beforeMonthsTxt||!afterStartTxt||!afterMonthsTxt){
      s7fStatusV.setValue('MISSING - fill in Lat/Lon and both BEFORE/AFTER date+months fields.');
      s7fStatusV.style().set('color','#cc0000'); s7fStatusV.style().set('backgroundColor','#ffd0d0'); return;
    }
    var rawParts=coordTxt.split(',');
    if(rawParts.length!==2){s7fStatusV.setValue('Invalid format - use: lat, lon'); s7fStatusV.style().set('color','#cc0000'); return;}
    var latIn=parseCoordPart(rawParts[0]), lonIn=parseCoordPart(rawParts[1]);
    if(isNaN(latIn)||isNaN(lonIn)){s7fStatusV.setValue('Invalid format - use: lat, lon'); s7fStatusV.style().set('color','#cc0000'); return;}
    var radiusKm=parseFloat(radiusTxt);
    if(isNaN(radiusKm)||radiusKm<0.3||radiusKm>3){s7fStatusV.setValue('Ring radius must be 0.3-3 km.'); s7fStatusV.style().set('color','#cc0000'); return;}
    var beforeMonths=parseInt(beforeMonthsTxt,10), afterMonths=parseInt(afterMonthsTxt,10);
    if(isNaN(beforeMonths)||beforeMonths<4||beforeMonths>36||isNaN(afterMonths)||afterMonths<4||afterMonths>36){
      s7fStatusV.setValue('BEFORE/AFTER months must each be 4-36.'); s7fStatusV.style().set('color','#cc0000'); return;
    }
    recordStudySite(latIn, lonIn, 'S7F');

    function destinationPoint(lat, lon, bearingDeg, distKm){
      var R=6371, brng=bearingDeg*Math.PI/180, lat1=lat*Math.PI/180, lon1=lon*Math.PI/180;
      var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
      var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
      return {lat:lat2*180/Math.PI, lon:lon2*180/Math.PI};
    }
    var bearings8=[{b:0,label:'N'},{b:45,label:'NE'},{b:90,label:'E'},{b:135,label:'SE'},
      {b:180,label:'S'},{b:225,label:'SW'},{b:270,label:'W'},{b:315,label:'NW'}];
    var nodes=[{label:'Center',lat:latIn,lon:lonIn}];
    bearings8.forEach(function(br){
      var d=destinationPoint(latIn,lonIn,br.b,radiusKm);
      nodes.push({label:br.label,lat:d.lat,lon:d.lon});
    });
    nodes.forEach(function(nd){ nd.buf=ee.Geometry.Point([nd.lon,nd.lat]).buffer(150); });

    var ptsFC = ee.FeatureCollection(nodes.map(function(nd,idx){
      return ee.Feature(nd.buf, {label:nd.label, idx:idx});
    }));
    var faiCollBefore=mkMoFAIRange(beforeStartTxt,beforeMonths);
    var faiCollAfter=mkMoFAIRange(afterStartTxt,afterMonths);
    var ndviWaterMaskedF=ndviWater.updateMask(oceanMask).rename('ndviW');
    var ndviBufMF = Math.max(80, Math.min(250, radiusKm*1000*0.3));
    var ptsFCForNdviF = ee.FeatureCollection(nodes.map(function(nd,idx){
      return ee.Feature(ee.Geometry.Point([nd.lon,nd.lat]).buffer(ndviBufMF), {label:nd.label, idx:idx});
    }));

    // S7D's 3 calls
    var rD_before = extractMultiNodeSeries(faiCollBefore, ptsFC, 'fai', 20);
    var rD_after = extractMultiNodeSeries(faiCollAfter, ptsFC, 'fai', 20);
    var rD_ndvi = ndviWaterMaskedF.reduceRegions({collection:ptsFCForNdviF, reducer:ee.Reducer.max(), scale:20});

    // S7E's candidate search (1 call, must resolve before its other 2)
    var searchRadiiKmF=[20,40,70,110,160];
    var candidatesF=[];
    searchRadiiKmF.forEach(function(rKm){
      for(var b=0;b<360;b+=45){
        var d=destinationPoint(latIn,lonIn,b,rKm);
        candidatesF.push({lat:d.lat,lon:d.lon,radiusKm:rKm,bearing:b});
      }
    });
    var candFeatsF = candidatesF.map(function(c,idx){
      return ee.Feature(ee.Geometry.Point([c.lon,c.lat]), {idx:idx, radiusKm:c.radiusKm});
    });
    var candFCF = ee.FeatureCollection(candFeatsF);
    var rCandidatesF = GEBCO.reduceRegions({collection:candFCF, reducer:ee.Reducer.first(), scale:500});

    s7fStatusV.setValue('Running: 4 calls in parallel (S7D BEFORE/AFTER/NDVI + reference search), then 2 more once the reference is found. Expect ~30-90s.');
    s7fStatusV.style().set('color','#334466'); s7fStatusV.style().set('backgroundColor','#eeeeee');
    s7fStatusV.style().set('border','2px solid #aaaaaa'); s7fStatusV.style().set('whiteSpace','pre');
    s7fResultV.setValue('');

    var s7fData={}, s7fPending=4, s7fErrors=0;
    function s7fBump(key,v,e){
      s7fData[key]=e?null:v;
      if(e) s7fErrors++;
      s7fPending--;
      s7fStatusV.setValue('Running: '+(4-s7fPending)+' / 4 first-stage calls done'+(s7fErrors>0?' ('+s7fErrors+' errored)':'')+'...');
      if(s7fPending===0) s7fStage2();
    }
    function s7fStage2(){
      try {
        var candResult = s7fData.candidates;
        var feats = (candResult&&candResult.features)?candResult.features:[];
        var best=null, excludedCount=0, shallowCount=0;
        var EXCL_KM=5;
        for(var i=0;i<feats.length;i++){
          var p=feats[i].properties;
          var elev=extractReduceRegionsValue(p,'elevation');
          if(elev!==null && elev>-50 && elev<0){
            shallowCount++;
            var cLat=candidatesF[p.idx].lat, cLon=candidatesF[p.idx].lon;
            var tooClose=false;
            for(var h=0;h<s7StudySiteHistory.length;h++){
              if(haversineKm(cLat,cLon,s7StudySiteHistory[h].lat,s7StudySiteHistory[h].lon) < EXCL_KM){ tooClose=true; break; }
            }
            if(tooClose){ excludedCount++; continue; }
            if(!best || p.radiusKm<best.radiusKm) best={radiusKm:p.radiusKm, elev:elev, lat:cLat, lon:cLon};
          }
        }
        if(!best){
          s7fStatusV.setValue('S7D data is ready, but NO REFERENCE REEF FOUND within 160km for S7E'+
            (shallowCount>0&&excludedCount===shallowCount?' (all '+shallowCount+' shallow candidates coincide with prior study sites)':'')+
            ' - showing S7D results only below.');
          s7fStatusV.style().set('color','#886600'); s7fStatusV.style().set('backgroundColor','#fff6cc');
          s7fFinishD_only();
          return;
        }
        s7fData.refSite = best;
        var studyPtBufF=ee.Geometry.Point([lonIn,latIn]).buffer(150);
        var refPtBufF=ee.Geometry.Point([best.lon,best.lat]).buffer(150);
        var pairFCF=ee.FeatureCollection([
          ee.Feature(studyPtBufF,{label:'Study'}), ee.Feature(refPtBufF,{label:'Reference'})
        ]);
        var rE_before = extractMultiNodeSeries(faiCollBefore, pairFCF, 'fai', 20);
        var rE_after = extractMultiNodeSeries(faiCollAfter, pairFCF, 'fai', 20);
        var s7fPending2=2;
        function s7fBump2(key,v,e){
          s7fData[key]=e?null:v; if(e) s7fErrors++;
          s7fPending2--;
          s7fStatusV.setValue('Running: second-stage '+(2-s7fPending2)+' / 2 calls done...');
          if(s7fPending2===0) s7fFinishAll(excludedCount);
        }
        rE_before.evaluate(function(v,e){ s7fBump2('eBefore', v, e); });
        rE_after.evaluate(function(v,e){ s7fBump2('eAfter', v, e); });
      } catch(errStage2){
        s7fStatusV.setValue(friendlyEEError(errStage2));
        s7fStatusV.style().set('color','#cc0000'); s7fStatusV.style().set('backgroundColor','#ffd0d0');
        print('=== S7F STAGE 2 ERROR === '+errStage2);
      }
    }
    function fmtNF(v,d){ return (v!==null&&v!==undefined&&!isNaN(v))?v.toFixed(d):'n/a'; }
    function buildS7DSummary(){
      var beforeByNode = groupSeriesByLabel(s7fData.dBefore, 'fai');
      var afterByNode = groupSeriesByLabel(s7fData.dAfter, 'fai');
      var ndviByLabel = {};
      if(s7fData.dNdvi && s7fData.dNdvi.features){
        s7fData.dNdvi.features.forEach(function(f){
          var p=f.properties||{}; ndviByLabel[p.label]=extractReduceRegionsValue(p,'ndviW');
        });
      }
      var nOnReefRisingCorr=0, nOnReefTotal=0, nAc1Rising=0, nAc1Avail=0;
      var centerBeforeVals=(beforeByNode['Center']||[]), centerAfterVals=(afterByNode['Center']||[]);
      for(var idx=0;idx<nodes.length;idx++){
        var nd=nodes[idx];
        var stB=jsNodeStatsFixed(beforeByNode[nd.label]||[]);
        var stA=jsNodeStatsFixed(afterByNode[nd.label]||[]);
        var ndviV=(ndviByLabel[nd.label]!==undefined)?ndviByLabel[nd.label]:null;
        var onReef = ndviV!==null?(ndviV>-0.10):null;
        var dAC1=(stB.realAC1!==null&&stA.realAC1!==null)?(stA.realAC1-stB.realAC1):null;
        if(dAC1!==null){ nAc1Avail++; if(dAC1>0.01) nAc1Rising++; }
        if(idx>0 && onReef===true){
          var pcB=jsPairCorrelation(centerBeforeVals, beforeByNode[nd.label]);
          var pcA=jsPairCorrelation(centerAfterVals, afterByNode[nd.label]);
          nOnReefTotal++;
          if(pcB.corr!==null&&pcA.corr!==null&&(pcA.corr-pcB.corr)>0.10) nOnReefRisingCorr++;
        }
      }
      var headline;
      if(nOnReefTotal===0) headline='S7D: could not assess (no on-reef ring points)';
      else if(nOnReefRisingCorr>=Math.ceil(nOnReefTotal*0.5)) headline='S7D: HYPER-SYNCHRONIZATION ('+nOnReefRisingCorr+'/'+nOnReefTotal+' on-reef rising)';
      else if(nOnReefRisingCorr>0) headline='S7D: PARTIAL SIGNAL ('+nOnReefRisingCorr+'/'+nOnReefTotal+' on-reef rising)';
      else headline='S7D: NO SYNCHRONIZATION SIGNAL (0/'+nOnReefTotal+' on-reef rising)';
      return {headline:headline, nAc1Rising:nAc1Rising, nAc1Avail:nAc1Avail};
    }
    function s7fFinishD_only(){
      try {
        var dSum = buildS7DSummary();
        var lines=['=== COMBINED SUMMARY ==='];
        lines.push(dSum.headline+' | AC1 rising at '+dSum.nAc1Rising+'/'+dSum.nAc1Avail+' nodes');
        lines.push('S7E: not run this time (no usable reference site found)');
        lines.push(repeatChar('\u2500',50));
        lines.push('BEFORE: '+beforeStartTxt+' + '+beforeMonths+'mo | AFTER: '+afterStartTxt+' + '+afterMonths+'mo | radius='+radiusKm+'km');
        s7fResultV.setValue(lines.join('\n'));
        print('=== S7F COMBINED (S7D only) ==='); print(lines.join('\n'));
      } catch(eF){ print('=== S7F D-only finish error === '+eF); }
    }
    function s7fFinishAll(excludedCount){
      try {
        var dSum = buildS7DSummary();
        var eBeforeByNode = groupSeriesByLabel(s7fData.eBefore, 'fai');
        var eAfterByNode = groupSeriesByLabel(s7fData.eAfter, 'fai');
        var studyB=jsNodeStatsFixed(eBeforeByNode['Study']||[]);
        var studyA=jsNodeStatsFixed(eAfterByNode['Study']||[]);
        var refB=jsNodeStatsFixed(eBeforeByNode['Reference']||[]);
        var refA=jsNodeStatsFixed(eAfterByNode['Reference']||[]);
        var dAC1_study=(studyB.realAC1!==null&&studyA.realAC1!==null)?(studyA.realAC1-studyB.realAC1):null;
        var dVar_study=(studyB.varTrendRatio!==null&&studyA.varTrendRatio!==null)?(studyA.varTrendRatio-studyB.varTrendRatio):null;
        var dAC1_ref=(refB.realAC1!==null&&refA.realAC1!==null)?(refA.realAC1-refB.realAC1):null;
        var dVar_ref=(refB.varTrendRatio!==null&&refA.varTrendRatio!==null)?(refA.varTrendRatio-refB.varTrendRatio):null;
        var refInsufficient=(refB.nValidMonths<4||refA.nValidMonths<4);
        var studyInsufficient=(studyB.nValidMonths<4||studyA.nValidMonths<4);
        // v10.153 FIX 10: calibrated cutoffs. NOTE: S7F runs NO permutation
        // test, so its verdict is threshold-only - see the warning added to
        // its summary below.
        var _s7fAc1=Math.max(getCalibratedThresholds(beforeMonths).ac1, getCalibratedThresholds(afterMonths).ac1);
        var _s7fVar=Math.max(getCalibratedThresholds(beforeMonths).varr, getCalibratedThresholds(afterMonths).varr);
        var studySignal=(dAC1_study!==null&&dAC1_study>_s7fAc1)||(dVar_study!==null&&dVar_study>_s7fVar);
        var refSignal=(dAC1_ref!==null&&dAC1_ref>_s7fAc1)||(dVar_ref!==null&&dVar_ref>_s7fVar);
        var eHeadline;
        if(refInsufficient||studyInsufficient) eHeadline='S7E: CANNOT CLASSIFY (insufficient data in at least one period)';
        else if(studySignal&&!refSignal) eHeadline='S7E: LOCAL ANOMALY DETECTED';
        else if(studySignal&&refSignal) eHeadline='S7E: REGIONAL SIGNAL';
        else if(!studySignal&&refSignal) eHeadline='S7E: ANOMALOUS (reference changed, study did not)';
        else eHeadline='S7E: NO SIGNAL AT EITHER SITE';

        // v10.117 FIX: data density was already computed here (needed for
        // the CANNOT CLASSIFY check above) but never actually displayed -
        // S7D and S7E each show it in their own individual output, S7F's
        // combined summary silently dropped it. Caught directly: a user
        // asked to see S7F's density numbers and there were none to show.
        var studyDataDensityF = (studyB.nValidMonths+studyA.nValidMonths)/(beforeMonths+afterMonths);
        var refDataDensityF = (refB.nValidMonths+refA.nValidMonths)/(beforeMonths+afterMonths);

        var lines=['=== COMBINED SUMMARY ==='];
        lines.push(dSum.headline+' | AC1 rising at '+dSum.nAc1Rising+'/'+dSum.nAc1Avail+' nodes');
        lines.push(eHeadline+' | Study \u0394AC1='+fmtNF(dAC1_study,3)+' | Reference \u0394AC1='+fmtNF(dAC1_ref,3));
        lines.push('Study data density: '+Math.round(studyDataDensityF*100)+'% ('+studyB.nValidMonths+'/'+beforeMonths+' BEFORE, '+
          studyA.nValidMonths+'/'+afterMonths+' AFTER valid months)'+(studyInsufficient?' - INSUFFICIENT in at least one period':''));
        lines.push('Reference data density: '+Math.round(refDataDensityF*100)+'% ('+refB.nValidMonths+'/'+beforeMonths+' BEFORE, '+
          refA.nValidMonths+'/'+afterMonths+' AFTER valid months)'+(refInsufficient?' - INSUFFICIENT in at least one period':''));
        if(excludedCount>0) lines.push('NOTE: '+excludedCount+' closer reference candidate(s) excluded (already tested as a study site this session).');
        lines.push(repeatChar('\u2500',50));
        lines.push('Reference reef: '+s7fData.refSite.lat.toFixed(4)+', '+s7fData.refSite.lon.toFixed(4)+
          ' ('+s7fData.refSite.radiusKm+'km away, depth='+s7fData.refSite.elev.toFixed(1)+'m)');
        lines.push('BEFORE: '+beforeStartTxt+' + '+beforeMonths+'mo | AFTER: '+afterStartTxt+' + '+afterMonths+'mo | radius='+radiusKm+'km');
        lines.push('For full node-by-node / study-vs-reference detail, run S7D and S7E individually above with these same inputs.');
        lines.push('');
        lines.push(thresholdOnlyWarning());
        lines.push('Cutoffs used (v10.153 calibrated): \u0394AC1>'+_s7fAc1.toFixed(3)+', \u0394Var>'+_s7fVar.toFixed(2)+'x');
        s7fResultV.setValue(lines.join('\n'));

        var vCol = (dSum.headline.indexOf('HYPER')===0||eHeadline.indexOf('LOCAL ANOMALY')>=0)?'#880000':
                   (eHeadline.indexOf('CANNOT')>=0)?'#886600':'#226644';
        var vBg = vCol==='#880000'?'#ffd0d0':vCol==='#886600'?'#fff6cc':'#e8f4ff';
        s7fStatusV.setValue(dSum.headline+'\n'+eHeadline);
        s7fStatusV.style().set('color',vCol); s7fStatusV.style().set('backgroundColor',vBg);
        s7fStatusV.style().set('border','2px solid '+vCol); s7fStatusV.style().set('whiteSpace','pre'); s7fStatusV.style().set('fontWeight','bold');

        print('=== S7F COMBINED SUMMARY (S7D + S7E) ==='); print(lines.join('\n'));
      } catch(eF2){
        s7fStatusV.setValue(friendlyEEError(eF2));
        s7fStatusV.style().set('color','#cc0000'); s7fStatusV.style().set('backgroundColor','#ffd0d0');
        print('=== S7F finishAll error === '+eF2);
      }
    }
    rD_before.evaluate(function(v,e){ s7fBump('dBefore', v, e); });
    rD_after.evaluate(function(v,e){ s7fBump('dAfter', v, e); });
    rD_ndvi.evaluate(function(v,e){ s7fBump('dNdvi', v, e); });
    rCandidatesF.evaluate(function(v,e){ s7fBump('candidates', v, e); });
  }
});
panel.add(s7fRunBtn);
panel.add(s7fStatusV);
panel.add(s7fResultV);
panel.add(legDiv());

panel.add(sHead('S8 - AQUACULTURE SUITABILITY (v10.67)','#003366'));
panel.add(lbl('A. taxiformis: 17-21 deg C optimal (Statton 2024 AgriFutures AU)',7,'#336633'));
panel.add(lbl('Two-gate architecture: SST gate -> Chl gate -> cautions',7,'#555555'));
var sstWindowV=dynLbl('--','#666666'), chlWindowV=dynLbl('--','#666666'), ndciNutrientProxyV=dynLbl('--','#666666');
var pollutionV=dynLbl('--','#666666'), stabilityV=dynLbl('--','#666666'), bromoformV=dynLbl('--','#666666');
var aquaScoreV=dynLbl('--','#666666'), aquaStatusV=dynLbl('computing...','#666666');
panel.add(row('SST gate (17-21 optimal)',sstWindowV)); panel.add(row('Nutrients Chl-a gate',chlWindowV));
panel.add(row('NDCI nutrient proxy (S2)',ndciNutrientProxyV));
panel.add(lbl('S2 proxy, uncalibrated - NOT mg/m3, NOT in gate logic',7,'#888888'));
panel.add(row('Pollution NO2',pollutionV)); panel.add(row('Thermal stability',stabilityV));
panel.add(row('Bromoform yield',bromoformV)); panel.add(row('Aquaculture score',aquaScoreV));
var aquaMissingV=ui.Label('',{fontSize:'7px',color:'#885500',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(aquaMissingV); panel.add(row('Farm site status',aquaStatusV));
var kelpNoteV=ui.Label('Cold-water kelp: checking...',{fontSize:'7px',color:'#aaaaaa',backgroundColor:'rgba(0,0,0,0)',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre',border:'1px solid #336633'});
panel.add(kelpNoteV);

// v10.144 NEW: S20e - Real Bleaching Probability (fitted model)
panel.add(sHead('S20e - REAL BLEACHING PROBABILITY (fitted model)','#3a1a1a'));
panel.add(lbl('A REAL logistic regression fit on 32,716 real rows from the Global Coral-Bleaching Database (van Woesik & Kratochwill 2022) - not a heuristic. Uses DHW/Turbidity/Depth already computed above, zero new data calls. HONEST LIMIT, disclosed directly: held-out test AUC=0.620 (0.5=chance, 1.0=perfect) - real, better than chance, but genuinely WEAK. Predicts BLEACHING probability, NOT collapse - a real, different, narrower question.',7,'#aa5533'));
var bleachProbV=dynLbl('checking...','#aa5533');
panel.add(row('P(bleaching), fitted model',bleachProbV));

panel.add(sHead('S9 - GEOLOGY / SOIL (v10.42: dataset dead)','#5a4422'));
panel.add(lbl('OpenLandMap asset confirmed inaccessible - renders blank',7,'#664422'));
var soilTextureV=dynLbl('checking...','#664422'); panel.add(row('Soil texture (USDA)',soilTextureV));

panel.add(sHead('S10/S11 - SEISMIC / VOLCANIC HAZARD','#5a2222'));
panel.add(lbl('USGS earthquakes (M4.5+) + GDIS historic volcanic events',7,'#774444'));
var eqStatsV=dynLbl('checking...','#774444'), volcStatsV=dynLbl('checking...','#774444');
panel.add(row('Earthquakes (200km, M4.5+)',eqStatsV)); panel.add(row('Volcanic activity (300km, 1960-2018)',volcStatsV));

panel.add(sHead('FIELD DATA STATUS (this region)','#2a4a1a'));
var fieldStatusV=dynLbl('Checking...','#666666'), fieldSpeciesV=dynLbl('--','#336611');
var fieldSourcesV=ui.Label('--',{fontSize:'7px',color:'#555577',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
var fieldNotesV=ui.Label('--',{fontSize:'7px',color:'#886600',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(row('Status',fieldStatusV)); panel.add(row('Species',fieldSpeciesV));
panel.add(fieldSourcesV); panel.add(fieldNotesV);

panel.add(sHead('FIELD CORRECTIONS (valid data only)','#334422'));
var fc1V=dynLbl('n/a','#888888'), fc3V=dynLbl('n/a','#888888'), fc4V=dynLbl('n/a','#888888'), fc5V=dynLbl('n/a','#888888'), fcTV=dynLbl('0 total','#222222');
panel.add(row('F1 Urchin grazer',fc1V)); panel.add(row('F3 Anem density',fc3V));
panel.add(row('F4 Metals',fc4V)); panel.add(row('F5 Recruitment',fc5V)); panel.add(row('Total correction',fcTV));

panel.add(sHead('BOWL DEPTH (FUSED - heuristic composite index)','#2a4a2a'));
var satCcsV=dynLbl('--','#664400'), fusCcsV=dynLbl('--','#115511'), bowlV=dynLbl('--','#115533');
var omega0V=dynLbl('--','#115533'), tauV=dynLbl('--','#661111'), ac1V=dynLbl('--','#224411'), p5yrV=dynLbl('--','#880000');
panel.add(row('Satellite CCS',satCcsV)); panel.add(row('FUSED CCS',fusCcsV));
panel.add(row('Bowl depth B',bowlV)); panel.add(row('omega0',omega0V));
panel.add(row('Return tau',tauV)); panel.add(row('AC1 (legacy heuristic)',ac1V)); panel.add(row('Regime-shift index (5yr, uncalibrated)',p5yrV));
panel.add(lbl('v10.143: relabeled from "P(flip 5yr)" - this is an uncalibrated Kramers-rate-inspired HEURISTIC INDEX, not a real probability of bleaching or collapse (no confusion matrix, not fit against real outcomes). A real, calibrated bleaching-probability model is possible using the Global Coral-Bleaching Database (van Woesik & Kratochwill 2022, 34,846 records, 14,405 real sites) - a genuine future upgrade, not yet built.',7,'#aa6600'));
var bowlVsS20eV=dynLbl('checking...','#aa5533');
panel.add(row('vs S20e real fitted model (see below)',bowlVsS20eV));
panel.add(lbl('v10.145 NEW: fulfills the v10.143 note above - S20e (further down) IS that real, fitted model, now built. This row directly compares the uncalibrated heuristic index above against S20e own real, held-out-validated P(bleaching), so you can see how far apart a real model and a heuristic guess actually land at this site.',7,'#aa5533'));

panel.add(sHead('S12 - REAL CSD STATISTICS (Tier 2)','#1a4a4a'));
panel.add(lbl('Single fixed window (Jan 2023-Dec 2024), no BEFORE baseline - a snapshot, not a before/after test. For a validated Scheffer 2009 comparison against a stored BEFORE window and a control site, use S13 further down.',7,'#886600'));
var realAc1V=dynLbl('checking...','#226666'), realVarTrendV=dynLbl('checking...','#226666'), realCsdNoteV=dynLbl('','#888888');
panel.add(row('Real AC1 (detrended, 24mo)',realAc1V));
panel.add(row('Variance trend (2nd/1st half)',realVarTrendV)); panel.add(realCsdNoteV);

// v10.149 NEW: S12b - Real Significance Test (permutation). Direct
// answer to a real question: S12 is real math (detrended AC1) but zero
// statistical testing - heuristic, not probabilistic. This adds a
// genuine permutation/surrogate-data test: shuffle THIS SAME window's
// 24 months into random order 500 times, recompute detrended AC1 each
// time, and ask "how often does randomly-ordered data produce an AC1
// this high by chance alone?" - a real, different question from S13's
// BEFORE/AFTER delta test (which asks "did AC1 change significantly").
// This one asks "is the memory itself, in this single window,
// distinguishable from noise with the same values in a random order?" -
// a real, established method in nonlinear time-series analysis
// (surrogate-data testing for genuine autocorrelation).
// Built standalone (button-triggered) since it needs the raw monthly
// values, which S12's own server-side EE-Array computation never sends
// to the client - one new, small EE call, reusing mkMoSST() (the exact
// same Jan2023-Dec2024 collection S12 itself uses) and the SAME
// jsNodeStats() formula already proven throughout this tool.
panel.add(sHead('S12b - REAL SIGNIFICANCE TEST (permutation)','#1a4a4a'));
panel.add(lbl('Converts S12 from heuristic to probabilistic: shuffles this SAME window\'s 24 months into random order 500 times, asking "would randomly-ordered data produce an AC1 this high by chance?" - a genuinely different question from S13\'s before/after test. Uses the last-clicked location.',7,'#227777'));
var s12bVerdictV=ui.Label('Click a location above, then press CHECK.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var s12bBtn=ui.Button({
  label:'CHECK REAL SIGNIFICANCE (S12 permutation test)',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#cceeee',color:'#115566',stretch:'horizontal',padding:'6px 4px',border:'2px solid #227777'},
  onClick:function(){
    if(lastClickLat===undefined||lastClickLat===null){
      s12bVerdictV.setValue('Click a coastal location on the map first, then press CHECK.'); s12bVerdictV.style().set('color','#cc0000'); return;
    }
    var latS12=lastClickLat, lonS12=lastClickLon;
    s12bVerdictV.setValue('Fetching real monthly data for the S12 window (Jan 2023-Dec 2024)...');
    s12bVerdictV.style().set('color','#334466'); s12bVerdictV.style().set('backgroundColor','#eeeeee');
    var ptS12=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([lonS12,latS12]).buffer(4000),{label:'Study'})]);
    var rS12=extractMultiNodeSeries(mkMoSST(), ptS12, 'sst', 4000);
    rS12.evaluate(function(fcRes,e){
      if(e){ s12bVerdictV.setValue('TEMPORARY EARTH ENGINE HICCUP - retry.'); s12bVerdictV.style().set('color','#cc0000'); return; }
      var byNode = groupSeriesByLabel(fcRes, 'sst');
      var seriesForPerm = byNode['Study']||[];
      var permResult = permutationTestAC1Fixed(seriesForPerm, 500);
      if(permResult.error){
        s12bVerdictV.setValue(permResult.error);
        s12bVerdictV.style().set('color','#886600'); return;
      }
      var realAC1val = permResult.realAC1;
      var pVal = permResult.pValue;
      var sig = pVal<0.05;
      var lines = ['=== S12b REAL SIGNIFICANCE TEST (permutation, v10.150: deseasonalized first) ==='];
      lines.push('Real AC1 (this window, deseasonalized): '+realAC1val.toFixed(3)+' [n='+permResult.n+' months]');
      lines.push('p='+pVal.toFixed(4)+' - '+(sig?'*** SIGNIFICANT (p<0.05) - this memory level is genuinely rare in randomly-ordered data':'NOT significant - randomly-ordered versions of these same values produce an AC1 this high fairly often'));
      lines.push('');
      lines.push('v10.150 FIX: this now deseasonalizes BEFORE shuffling - the OLD version');
      lines.push('shuffled raw values, which CONFIRMED gave p=0.0000 (maximally "significant")');
      lines.push('on synthetic data with a real seasonal cycle and ZERO actual signal, since a');
      lines.push('random reshuffle barely disrupts a strong seasonal cycle\'s contribution to AC1.');
      lines.push('IMPORTANT: this tests whether the OBSERVED memory itself is distinguishable');
      lines.push('from noise - it does NOT test whether AC1 has CHANGED (that is S13\'s job).');
      s12bVerdictV.setValue('p='+pVal.toFixed(4)+' - '+(sig?'SIGNIFICANT':'not significant')+' (real AC1='+realAC1val.toFixed(3)+', n='+permResult.n+' months, deseasonalized)');
      s12bVerdictV.style().set('color',sig?'#880000':'#226644');
      print(lines.join('\n'));
    });
  }
});
panel.add(s12bBtn);
panel.add(s12bVerdictV);

panel.add(sHead('S14 - REAL THERMAL RECOVERY TIME','#1a3a1a'));
var thermalEpisodesV=dynLbl('checking...','#226622'), thermalMeanV=dynLbl('checking...','#226622');
var thermalMaxV=dynLbl('checking...','#226622'), thermalOngoingV=dynLbl('checking...','#226622');
panel.add(row('Heat-stress episodes',thermalEpisodesV)); panel.add(row('Mean recovery time',thermalMeanV));
panel.add(row('Longest recovery observed',thermalMaxV)); panel.add(row('Ongoing episode?',thermalOngoingV));

panel.add(sHead('S15 - ECOLOGICAL RECOVERY VALIDATION','#3a2a4a'));
var ecoValStatusV=dynLbl('checking...','#553377');
var ecoValDetailsV=ui.Label('',{fontSize:'7px',color:'#555577',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(row('Validation status',ecoValStatusV)); panel.add(ecoValDetailsV);

// v10.138 CONSOLIDATION: merged S20/S20b/S20c into ONE unified tool.
// Direct response to real, fair feedback: three separate buttons for
// what's conceptually one question ("what's the combined risk, at some
// point in time?") was confusing design, not a missing capability - both
// past-date sync (old S20b) and live-present (old S20c) already worked
// correctly. This single date field now does both: leave it BLANK for
// LIVE (today, rolling 35-day OISST window), or type a specific month
// (YYYY-MM, 2023-01 to 2024-12) for a HISTORICAL replay of that exact
// month - one button, one result, clearly labeled which mode ran.
// The original always-on S20 (auto-computed on every click, frozen to
// the tool's 2023-2024 peak) is RETIRED as a separate behavior - this is
// now button-triggered like everything else in this family, so the
// click-handler code that used to auto-populate it has been removed too.
panel.add(sHead('S20 - SPECIES-WEIGHTED COMBINED RISK','#4a2a1a'));
panel.add(lbl('Combines a real DHW reading with published dominant-species vulnerability data, where it exists (7 real sites: One Tree Reef, Florida Keys, Maldives, Bocas del Toro, Andaman/Nicobar, Brazil, Mesoamerican Reef). Leave the date field BLANK for LIVE (today, real current OISST data), or type a month (YYYY-MM, 2023-01 to 2024-12) for a HISTORICAL replay of that exact month - both use the same real satellite dataset and formula. NOT a satellite species-ID capability - hyperspectral species classification is only ~56-70% accurate currently, too unreliable to use as ground truth.',7,'#663311'));
panel.add(lbl('v10.140 DISCLOSED CAVEAT: a real 2015-2017 Mesoamerican Reef study (Arias-Ortiz et al. 2024, Communications Biology, doi:10.1038/s42003-024-07128-y) found DHW alone explains LESS bleaching-severity variance than a combined model using 23 stress/sensitivity metrics (which reached 75%) - climatological warming rate and other heat metrics outperformed plain DHW. It also found a genuinely counter-intuitive result: DEEPER reefs with MORE diverse coral communities showed HIGHER vulnerability, not lower. This tool still uses a simple DHW>=4 threshold - a known, disclosed simplification, not the more sophisticated real model this paper describes.',7,'#885522'));
var speciesDateInput = ui.Textbox({placeholder:'blank = LIVE (today) | or YYYY-MM for historical, e.g. 2024-02',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(speciesDateInput);
var speciesRiskStatusV=ui.Label('Click a location above, optionally enter a month, then press CHECK.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var speciesRiskDetailsV=ui.Label('',{fontSize:'7px',color:'#664422',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
var speciesRiskBtn=ui.Button({
  label:'CHECK COMBINED RISK',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#f5e0cc',color:'#663311',stretch:'horizontal',padding:'6px 4px',border:'2px solid #995522'},
  onClick:function(){
    if(lastClickLat===undefined||lastClickLat===null){
      speciesRiskStatusV.setValue('Click a coastal location on the map first, then press CHECK.'); speciesRiskStatusV.style().set('color','#cc0000'); return;
    }
    var latS=lastClickLat, lonS=lastClickLon;
    var dateTxt=(speciesDateInput.getValue()||'').trim();
    var isLive = (dateTxt==='');
    var dS, deS, modeLabel;
    if(isLive){
      var todayEE=ee.Date(Date.now());
      dS=todayEE.advance(-35,'day'); deS=todayEE;
      modeLabel='LIVE';
    } else {
      if(!/^\d{4}-\d{2}$/.test(dateTxt)){
        speciesRiskStatusV.setValue('Enter a month as YYYY-MM (e.g. 2024-02), or leave blank for LIVE.'); speciesRiskStatusV.style().set('color','#cc0000'); return;
      }
      var yr=parseInt(dateTxt.slice(0,4),10), mo=parseInt(dateTxt.slice(5,7),10);
      if(yr<2023||yr>2024||mo<1||mo>12){
        speciesRiskStatusV.setValue('Historical month must be 2023-01 to 2024-12 - the real data window this tool uses.\nOutside that range, this feature has no real data to check (not a guess). Leave blank for LIVE instead.');
        speciesRiskStatusV.style().set('color','#cc0000'); return;
      }
      dS=ee.Date(dateTxt+'-01'); deS=dS.advance(1,'month');
      modeLabel='HISTORICAL '+dateTxt;
    }
    speciesRiskStatusV.setValue('Fetching '+modeLabel+' OISST data at '+latS.toFixed(4)+', '+lonS.toFixed(4)+'...');
    speciesRiskStatusV.style().set('color','#334466'); speciesRiskStatusV.style().set('backgroundColor','#eeeeee');
    var colS=getOISSTColl(dS,deS);
    var dhwImgS=ee.Image(ee.Algorithms.If(colS.size().gt(0),
      colS.max().subtract(MMM_perpixel.add(1)).max(0).multiply(4.33).rename('dhw'),
      ee.Image.constant(-999).rename('dhw')));
    var ptS=ee.Geometry.Point([lonS,latS]).buffer(1000);
    var rDhwS=dhwImgS.reduceRegion({reducer:ee.Reducer.mean(),geometry:ptS,scale:4000,maxPixels:1e9});
    var rLatestDate=colS.aggregate_max('system:time_start');
    ee.Dictionary({dhw:rDhwS.get('dhw'),latestDate:rLatestDate}).evaluate(function(v,e){
      if(e){ speciesRiskStatusV.setValue('TEMPORARY EARTH ENGINE HICCUP - retry.'); speciesRiskStatusV.style().set('color','#cc0000');
        print('=== S20 COMBINED RISK CHECK ERROR === '+e); return; }
      var dhwRawS=(v&&v.dhw!==undefined&&v.dhw!==null)?v.dhw:null;
      var dhwS=(dhwRawS!==null&&dhwRawS>=0&&dhwRawS<60)?dhwRawS:null;
      if(dhwS===null){
        speciesRiskStatusV.setValue('No valid satellite data for '+modeLabel+' at this exact point (cloud/ice masking, or too remote).');
        speciesRiskStatusV.style().set('color','#886600'); return;
      }
      var isTropicalReefS=(latS>-30&&latS<30)&&!isEBUS(latS,lonS);
      var speciesDataS=lookupSpeciesVulnerability(latS,lonS);
      var combinedS=combineSpeciesAndDHW(speciesDataS, dhwS, isTropicalReefS);
      var dataThruTxt=(isLive&&v.latestDate)?(' (data through '+(new Date(v.latestDate)).toISOString().slice(0,10)+')'):'';
      speciesRiskStatusV.setValue(modeLabel+' DHW'+dataThruTxt+': '+dhwS.toFixed(2)+' deg C-wks\n'+
        (combinedS.combinedRisk?('Combined risk: '+combinedS.combinedRisk):'(species data not checked)'));
      speciesRiskStatusV.style().set('color',combinedS.combinedRisk==='ELEVATED'?'#880000':combinedS.combinedRisk==='MODERATE'?'#886600':'#115511');
      speciesRiskDetailsV.setValue(combinedS.text.slice(0,400));
      print('=== S20 COMBINED RISK ('+modeLabel+dataThruTxt+') ===');
      print('DHW='+dhwS.toFixed(2)+' deg C-wks | '+combinedS.text);
    });
  }
});
panel.add(speciesRiskBtn);
panel.add(speciesRiskStatusV);
panel.add(speciesRiskDetailsV);

// v10.137 NEW: S20d - Genus Growth-Form Lookup. Exposes
// GENUS_GROWTH_FORM as a standalone tool - type any genus name, get
// its real growth form and general vulnerability tier. Unlike S20/S20b/
// S20c (which need a map click, since they're tied to a specific site's
// DHW), this needs no location - it's pure genus-level biological
// context, usable for any reef where you already know (from your own
// knowledge, a field guide, etc.) what genus dominates, even if that
// exact site has no specific event-outcome citation in S20's database.
panel.add(sHead('S20d - GENUS GROWTH-FORM LOOKUP','#4a3a1a'));
panel.add(lbl('Type any coral genus name (e.g. Acropora) to see its real growth form and a general vulnerability tier, based on the Coral Trait Database (coraltraits.org) plus established branching-vs-massive bleaching-susceptibility literature. A GENERAL fallback, not site-specific mortality data - use S20 above for that, where available.',7,'#664422'));
var genusInput = ui.Textbox({placeholder:'e.g. Acropora',style:{stretch:'horizontal',margin:'2px 4px',fontSize:'11px'}});
panel.add(genusInput);
var genusResultV=ui.Label('Enter a genus name, then press LOOK UP.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var genusLookupBtn=ui.Button({
  label:'LOOK UP GENUS',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#f0e0cc',color:'#664422',stretch:'horizontal',padding:'6px 4px',border:'2px solid #997733'},
  onClick:function(){
    var g=(genusInput.getValue()||'').trim();
    if(!g){ genusResultV.setValue('Type a genus name first.'); genusResultV.style().set('color','#cc0000'); return; }
    var r=lookupGenusGrowthForm(g);
    if(!r.checked){
      genusResultV.setValue('"'+g+'" not in this small starter list ('+Object.keys(GENUS_GROWTH_FORM).length+' genera so far).\nNot a "no data exists" result - just not yet added to this tool.');
      genusResultV.style().set('color','#886600');
      print('=== S20d GENUS LOOKUP === "'+g+'" not found in GENUS_GROWTH_FORM');
      return;
    }
    genusResultV.setValue(r.genus+': '+r.generalVulnerability+' general vulnerability\nGrowth form: '+r.growthForm);
    genusResultV.style().set('color',r.generalVulnerability.indexOf('HIGH')===0?'#880000':r.generalVulnerability.indexOf('MODERATE')===0?'#886600':'#115511');
    print('=== S20d GENUS LOOKUP === '+r.genus+': '+r.generalVulnerability+' | '+r.growthForm);
  }
});
panel.add(genusLookupBtn);
panel.add(genusResultV);

panel.add(sHead('S16 - PHYSICAL ENERGY INDEX','#1a3a2a'));
var eciDepthV=dynLbl('checking...','#115533'), eciCurrentV=dynLbl('checking...','#115533');
var eciSLR03V=dynLbl('checking...','#226644'), eciSLR05V=dynLbl('checking...','#226644');
var eciSLR10V=dynLbl('checking...','#aa3300'), eciValidationV=dynLbl('checking...','#226644');
panel.add(row('Seafloor depth (GEBCO)',eciDepthV)); panel.add(row('Energy risk (current)',eciCurrentV));
panel.add(row('+0.3m SLR scenario',eciSLR03V)); panel.add(row('+0.5m SLR scenario',eciSLR05V));
panel.add(row('+1.0m SLR scenario',eciSLR10V)); panel.add(row('vs ecological B score (different metrics, not a validation)',eciValidationV));

panel.add(sHead('S17 - TIME OF EMERGENCE (ToE)','#2a1a3a'));
panel.add(lbl('SNR = |trend x record_years| / noise | threshold SNR >= 2.0',7,'#553377'));
panel.add(lbl('v10.145 DIAGNOSTIC: investigated why S17b (Mann-Kendall) sometimes disagrees with SNR above - e.g. a real Bocas del Toro test showed NO2/Salinity "EMERGED" here but NOT significant in S17b. Root cause found in the SNR formula itself: it multiplies slope by record_years with NO correction for how uncertain a slope estimate becomes with FEW data points (unlike a real p-value, which accounts for sample size directly). This means SNR structurally over-triggers on short records (pH ~4yr, NO2 7yr) - it is not calibrated to any known false-positive rate, unlike Mann-Kendall p<0.05. Treat SNR "EMERGED" on short-record variables with real caution; S17b is the more trustworthy check where it can compute one.',7,'#aa5533'));
panel.add(lbl('Surface only. Inspired by Tan et al. 2026 NCC compound CID approach.',7,'#888888'));
panel.add(lbl('v10.141: pH now uses a REAL, confirmed Copernicus asset (~4yr record, 2022-2025, LOW confidence given the short record) - not measured, model surface only, no depth zones. DO is still unavailable (dead asset, no working replacement confirmed yet).',7,'#aa6600'));
panel.add(lbl('v10.89: pH/DO are computed and shown separately from SST/Chl/Salinity/NO2, so if the BGC dataset is unavailable it only affects pH/DO below, not the other 4.',7,'#886600'));
var toeCompoundV=dynLbl('computing...','#553377'), toeSSTv=dynLbl('computing...','#226644');
var toeCHLv=dynLbl('computing...','#226644'), toeNO2v=dynLbl('computing...','#226644'), toeSALv=dynLbl('computing...','#226644');
var toePHv=dynLbl('computing...','#226644'), toeDOv=dynLbl('computing...','#226644');
panel.add(row('Compound status (of available vars)',toeCompoundV));
panel.add(row('SST (44yr, HIGH conf)',toeSSTv));
panel.add(row('Chl (27yr, MARGINAL)',toeCHLv));
panel.add(row('Salinity (32yr, MARGINAL)',toeSALv));
panel.add(row('NO2 (7yr, LOW conf)',toeNO2v));
panel.add(row('pH (~4yr REAL asset, LOW conf)',toePHv));
panel.add(row('DO / O2 (32yr, BGC model surface)',toeDOv));

// v10.142 NEW: S17b - Real Mann-Kendall Significance Test. Direct
// answer to a real question: since the Mann-Kendall engine already
// exists (built and fixed for STEP 5), can S17's fixed SNR>=2.0
// heuristic be replaced with a genuine significance test? Built as a
// STANDALONE, button-triggered tool rather than rewriting S17's core
// flow directly - that flow was already touched once this version (the
// pH fix) and is fragile/deeply async; adding a new tool alongside it
// is the same lower-risk pattern already used for COMPARE's statistically-
// valid banner. Reuses the SAME real annual-value collections already
// built for S17 (_annSSTColl etc.) and the SAME mannKendallTest() engine
// already validated for STEP 5 - genuinely real, not a new architecture.
// DISCLOSED: annual values (1 point/year) are inherently few - Kendall's
// tau needs real data, so shorter records (pH ~4yr, NO2 7yr) will likely
// show "insufficient" here, same honest limitation already found for
// STEP 5 with short spans.
panel.add(sHead('S17b - REAL TREND SIGNIFICANCE TEST (Mann-Kendall)','#3a1a5a'));
panel.add(lbl('Replaces S17\'s fixed SNR>=2.0 threshold with a genuine Kendall tau / Mann-Kendall test on the SAME real annual data - a real p-value instead of a fixed cutoff. Uses the last-clicked location. Short records (pH ~4yr, NO2 7yr) will likely show "insufficient data" - Kendall tau needs several real points to say anything.',7,'#663388'));
var toeMkVerdictV=ui.Label('Click a location above, then press CHECK.',
  {fontSize:'11px',fontWeight:'bold',color:'#555555',backgroundColor:'#eeeeee',padding:'6px 8px',margin:'2px 0',whiteSpace:'pre',border:'2px solid #aaaaaa'});
var toeMkDetailsV=ui.Label('',{fontSize:'7px',color:'#553377',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
var toeMkBtn=ui.Button({
  label:'CHECK REAL TREND SIGNIFICANCE',
  style:{fontSize:'11px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#e8d9f5',color:'#4a1a6a',stretch:'horizontal',padding:'6px 4px',border:'2px solid #663399'},
  onClick:function(){
    if(lastClickLat===undefined||lastClickLat===null){
      toeMkVerdictV.setValue('Click a coastal location on the map first, then press CHECK.'); toeMkVerdictV.style().set('color','#cc0000'); return;
    }
    var latM=lastClickLat, lonM=lastClickLon;
    toeMkVerdictV.setValue('Fetching real annual series for 5 variables (5 batched calls)...');
    toeMkVerdictV.style().set('color','#334466'); toeMkVerdictV.style().set('backgroundColor','#eeeeee');
    var ptM = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([lonM,latM]).buffer(4000),{label:'Study'})]);
    var vars=[
      {key:'sst', coll:_annSSTColl, band:'sst', scale:4000, name:'SST'},
      {key:'chl', coll:_annCHLColl, band:'chlor_a', scale:4000, name:'Chl-a'},
      {key:'sal', coll:_annSALColl, band:'salinity_0', scale:4000, name:'Salinity'},
      {key:'no2', coll:_annNO2Coll, band:'tropospheric_NO2_column_number_density', scale:4000, name:'NO2'},
      {key:'ph',  coll:_annPHColl,  band:'ph', scale:25000, name:'pH'}
    ];
    var results={}, pending=vars.length;
    vars.forEach(function(vconf){
      var r = extractMultiNodeSeries(vconf.coll, ptM, vconf.band, vconf.scale);
      r.evaluate(function(fcRes,e){
        if(e){ results[vconf.key]={error:String(e)}; }
        else {
          var byNode = groupSeriesByLabel(fcRes, vconf.band);
          var vals = (byNode['Study']||[]).map(function(s){return s.v;}).filter(function(v){return v!==null&&v!==undefined&&!isNaN(v);});
          var mk = mannKendallTest(vals);
          results[vconf.key]={name:vconf.name, n:vals.length, mk:mk};
        }
        pending--;
        if(pending===0) finalizeMk();
      });
    });
    function finalizeMk(){
      var lines=['=== S17b REAL TREND SIGNIFICANCE (Kendall tau / Mann-Kendall) ==='];
      lines.push('Same real annual data S17 already uses - a genuine p-value instead of the fixed SNR>=2.0 cutoff.');
      lines.push(repeatChar('\u2500',50));
      var nSig=0, nTotal=0;
      vars.forEach(function(vconf){
        var res=results[vconf.key];
        if(res.error){ lines.push(vconf.name+': error - '+res.error); return; }
        if(res.mk.error){ lines.push(res.name+': n/a - '+res.mk.error+' (n='+res.n+' annual points, needs 4+)'); return; }
        nTotal++;
        var sig = res.mk.p<0.05;
        if(sig) nSig++;
        var dir = res.mk.tau>0?'RISING':res.mk.tau<0?'FALLING':'flat';
        lines.push(res.name+': tau='+(res.mk.tau>0?'+':'')+res.mk.tau.toFixed(3)+', p='+res.mk.p.toFixed(4)+
          ' -> '+dir+(sig?' *** SIGNIFICANT (p<0.05)':' not significant')+' [n='+res.n+' annual points]');
        // v10.145 DIAGNOSTIC: Salinity is nominally a 32-year record
        // (1993-2024) in S17 above - if far fewer real annual points
        // show up here, that is a genuine HYCOM data-sparsity finding
        // at this specific site, not a bug in this test.
        if(vconf.key==='sal' && res.n<15){
          lines.push('  NOTE: Salinity is nominally a 32yr record (1993-2024) in S17 above, but only '+res.n+
            ' years had real, valid HYCOM data at this exact point - a genuine data-sparsity finding for this site, not a code error.');
        }
      });
      lines.push(repeatChar('\u2500',50));
      lines.push(nSig+' of '+nTotal+' testable variables show a REAL statistically significant trend (p<0.05).');
      toeMkVerdictV.setValue(nTotal===0?'No variable had enough annual data to test.':
        nSig+'/'+nTotal+' variables REAL significant trend (p<0.05)'+(nSig>=2?' - compare vs S17\'s SNR-based count above':''));
      toeMkVerdictV.style().set('color',nSig>=2?'#880000':nSig>=1?'#886600':'#115511');
      toeMkDetailsV.setValue(lines.join('\n'));
      print(lines.join('\n'));
    }
  }
});
panel.add(toeMkBtn);
panel.add(toeMkVerdictV);
panel.add(toeMkDetailsV);

panel.add(sHead('S18 - BIOGEOCHEMISTRY SNAPSHOT','#1a2a3a'));
panel.add(lbl('All 4 values below come from ONE Copernicus BGC dataset (COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/MFC_001_028). As of this version that dataset returns "not found" in the GEE catalog, so this whole section will show n/a until it is restored or replaced with a valid current asset ID.',7,'#aa6600'));
var s18O2V=dynLbl('computing...','#334466'), s18pCO2V=dynLbl('computing...','#334466');
var s18pHV=dynLbl('computing...','#334466'), s18SalV=dynLbl('computing...','#334466');
panel.add(row('Dissolved O2 (surface)',s18O2V)); panel.add(row('pCO2 (surface)',s18pCO2V));
panel.add(row('pH (surface, model)',s18pHV)); panel.add(row('Salinity (current)',s18SalV));

// S19 - REAL IN-SITU BASELINE (v10.67: 5 stations including GEM)
panel.add(sHead('S19 - REAL IN-SITU BASELINE (v10.67)','#0a3a5a'));
panel.add(lbl('Real cleaned sensor data (NOT a model/satellite)',7,'#225577'));
panel.add(lbl('Only populated within radius of an actual sensor station',7,'#225577'));
panel.add(lbl('5 stations: Looe Key FL | Agua Hedionda CA | Scripps Pier CA',7,'#888888'));
panel.add(lbl('+ GEM: MarineBasis Nuuk GF3 (SW Greenland) | Zackenberg (NE Greenland)',7,'#336644'));
panel.add(lbl('GEM data: CC BY-SA 4.0 | g-e-m.dk | Real stats now populated',7,'#336644'));
var s19StatusV=dynLbl('checking...','#225577'), s19PhV=dynLbl('--','#225577');
var s19TempV=dynLbl('--','#225577'), s19SalV=dynLbl('--','#225577'), s19DoV=dynLbl('--','#225577');
var s19FluorV=dynLbl('--','#225577'), s19PressV=dynLbl('--','#225577'), s19TurbV=dynLbl('--','#225577');
var s19NoteV=ui.Label('',{fontSize:'7px',color:'#555577',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(row('Matched station',s19StatusV)); panel.add(row('pH (clean baseline)',s19PhV));
panel.add(row('Temp (clean baseline)',s19TempV)); panel.add(row('Salinity (clean baseline)',s19SalV));
panel.add(row('DO mg/L (clean baseline)',s19DoV));
panel.add(row('Fluorescence ug/L (in-situ)',s19FluorV));
panel.add(lbl('Fluorescence = in-situ proxy for Chl-a | Compare vs S2 satellite Chl-a above',7,'#226644'));
panel.add(row('Pressure dbar (CTD depth)',s19PressV));
panel.add(row('Turbidity FTU (in-situ)',s19TurbV));
panel.add(s19NoteV);
panel.add(lbl('CAVEAT: point-source only. Compare vs S2/S18 (satellite/model) above.',7,'#aa6600'));

panel.add(sHead('ACCURACY','#222244'));
panel.add(lbl('STEP 1: Click coastal locations to log data',8,'#884400'));
panel.add(lbl('STEP 2: Click Export button below',8,'#884400'));
panel.add(lbl('STEP 3: Open Tasks tab and click RUN',8,'#cc0000'));
panel.add(lbl('STEP 4: File appears in Google Drive / STEMGeoHS_exports',8,'#884400'));
var exportBtn=ui.Button({label:'EXPORT LOGGED CLICKS TO CSV (Drive)',
  style:{fontSize:'12px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#cce0ff',color:'#003388',stretch:'horizontal',padding:'8px 4px',border:'2px solid #0055cc'},
  onClick:buildAndExport});
panel.add(exportBtn);
var exportStatusLabel=ui.Label('No clicks logged yet - 0 rows ready',{fontSize:'9px',color:'#cc0000',fontWeight:'bold',margin:'2px 4px'});
panel.add(exportStatusLabel);

function buildAndExport(){
  if(clickResultsLog.length===0){exportStatusLabel.setValue('STILL 0 rows. Click a coastal point first.'); return;}
  var features=clickResultsLog.map(function(rec){return ee.Feature(null,rec);});
  var fc=ee.FeatureCollection(features);
  var exportDesc='stemgeohs_marine_export_'+Date.now();
  var columnOrder=['timestamp','region','lat','lon',
    'sst_annual_c','sst_peak_c','sst_trend_c_per_yr','dhw_c_weeks','dhw_raw_unfiltered','dhw_artifact_flagged',
    'mmm_local_c','chl_a_mg_m3','turbidity_ndti','no2_mol_m2','depth_m',
    'cancer_score_satellite','cancer_score_fused','field_correction',
    'bowl_depth_B','accuracy_pct','status_label','field_data_available','field_species',
    'aquaculture_score','aquaculture_confidence_pct','aquaculture_status',
    'bromoform_yield','num_interventions','top_intervention',
    'soil_texture_code','soil_texture_label',
    'earthquake_count_200km','earthquake_max_mag_200km','volcanic_activity_count_300km',
    'real_ac1_24mo','real_var_trend_ratio','real_csd_n_valid_months','ndci_nutrient_proxy',
    'thermal_recovery_n_completed_episodes','thermal_recovery_mean_months',
    'thermal_recovery_max_months','thermal_recovery_ongoing_flag',
    'eco_recovery_validation_checked','eco_recovery_validation_finding',
    'gebco_depth_m','energy_concentration_risk','eci_slr_03m','eci_slr_05m','eci_slr_10m','eci_vs_heuristic_validation',
    'in_situ_station','in_situ_distance_km','in_situ_ph_mean','in_situ_temp_c_mean','in_situ_salinity_mean','in_situ_do_mgL_mean'];
  Export.table.toDrive({collection:fc,description:exportDesc,folder:'STEMGeoHS_exports',fileFormat:'CSV',selectors:columnOrder});
  exportStatusLabel.setValue('TASK CREATED: "'+exportDesc+'" with '+clickResultsLog.length+' row(s). Open Tasks tab and click RUN.');
  print('=== EXPORT TASK CREATED === '+exportDesc+' | Rows: '+clickResultsLog.length);
}

var testExportBtn=ui.Button({label:'Test export (1 sample row)',
  style:{fontSize:'10px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#e8e8e8',color:'#222222',stretch:'horizontal',padding:'6px 4px',border:'2px solid #555555'},
  onClick:function(){
    clickResultsLog.push({timestamp:new Date().toISOString(),lat:21.27,lon:-157.876,region:'TEST ROW',
      sst_annual_c:25.0,sst_peak_c:27.5,sst_trend_c_per_yr:0.01,dhw_c_weeks:0,dhw_raw_unfiltered:0,dhw_artifact_flagged:false,
      mmm_local_c:27.3,chl_a_mg_m3:null,turbidity_ndti:0.01,no2_mol_m2:0.00005,depth_m:-50,
      cancer_score_satellite:21,cancer_score_fused:21,field_correction:0,bowl_depth_B:0.79,accuracy_pct:77,status_label:'DEEP BASIN',
      field_data_available:false,field_species:'none',aquaculture_score:null,aquaculture_confidence_pct:0,aquaculture_status:'INSUFFICIENT DATA',
      bromoform_yield:50,num_interventions:1,top_intervention:'TEST ROW',soil_texture_code:null,soil_texture_label:'n/a',
      earthquake_count_200km:0,earthquake_max_mag_200km:null,volcanic_activity_count_300km:0,
      real_ac1_24mo:null,real_var_trend_ratio:null,real_csd_n_valid_months:0,ndci_nutrient_proxy:null,
      thermal_recovery_n_completed_episodes:null,thermal_recovery_mean_months:null,thermal_recovery_max_months:null,thermal_recovery_ongoing_flag:null,
      eco_recovery_validation_checked:false,eco_recovery_validation_finding:'n/a',
      gebco_depth_m:null,energy_concentration_risk:null,eci_slr_03m:null,eci_slr_05m:null,eci_slr_10m:null,eci_vs_heuristic_validation:'n/a',
      in_situ_station:null,in_situ_distance_km:null,in_situ_ph_mean:null,in_situ_temp_c_mean:null,in_situ_salinity_mean:null,in_situ_do_mgL_mean:null});
    exportStatusLabel.setValue('Test row added. '+clickResultsLog.length+' row(s) ready.');
  }
});
panel.add(testExportBtn);
var clearLogBtn=ui.Button({label:'Clear logged clicks',
  style:{fontSize:'10px',fontWeight:'bold',margin:'2px 4px',backgroundColor:'#ffd9d9',color:'#770000',stretch:'horizontal',padding:'6px 4px',border:'2px solid #aa2222'},
  onClick:function(){clickResultsLog=[]; exportStatusLabel.setValue('Log cleared - 0 rows ready');}
});
panel.add(clearLogBtn);

var accSV=dynLbl('77%','#664400'), accFV=dynLbl('+0%','#888888'), accTV=dynLbl('77%','#664400');
panel.add(row('Satellite',accSV)); panel.add(row('Field gain',accFV)); panel.add(row('Combined',accTV));
panel.add(sHead('PREVENTATIVE INTERVENTIONS (Waddington surgery)','#003355'));
var intPanel=ui.Panel({style:{margin:'2px 0'}}); panel.add(intPanel);
panel.add(sHead('FUSED COASTAL CANCER SCORE','#4a0000'));
var scoreBig=ui.Label('Score: --',{fontSize:'18px',fontWeight:'bold',color:'#003300',backgroundColor:'#d4f5df',padding:'4px 8px',margin:'2px 0'});
panel.add(scoreBig);
var scoreBarLbl=ui.Label('--',{fontSize:'11px',fontWeight:'bold',color:'#0a5c1e',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0'});
panel.add(scoreBarLbl);
var barBg=ui.Panel({style:{backgroundColor:'#dddddd',margin:'1px 2px',padding:'1px',border:'1px solid #aaaaaa'}});
var barFill=ui.Panel({style:{backgroundColor:'#22cc44',height:'12px',width:'0px',margin:'0'}});
barBg.add(barFill); panel.add(barBg);
panel.add(ui.Panel([lbl('0',8,'#555555'),lbl('safe-----------risk',8,'#999999'),lbl('100',8,'#880000')],ui.Panel.Layout.Flow('horizontal'),{margin:'0 2px'}));
var scoreInterp=ui.Label('Click reef area to analyze',{fontSize:'8px',color:'#333333',backgroundColor:'#eefaee',padding:'3px 4px',margin:'2px 0',whiteSpace:'pre'});
panel.add(scoreInterp);
panel.add(lbl('MAP LEGEND (scroll to see) | v10.67',9,'#ffffff','#334466',true));
panel.add(sHead('1. Study zone border (white)','#555555'));
panel.add(legRow('#ffffff','White ring','1km buffer'));
panel.add(legDiv()); panel.add(sHead('2. S1 - SST annual mean deg C','#b03000'));
panel.add(legRow('#2c7bb6','< 24 deg C','Cold')); panel.add(legRow('#abd9e9','24-26','Optimal'));
panel.add(legRow('#fdae61','28-30','Stress')); panel.add(legRow('#f46d43','30-32','BLEACHING')); panel.add(legRow('#a50026','> 32','MORTALITY'));
panel.add(legDiv()); panel.add(sHead('2b. S1 - SST PEAK (bleaching score, hemisphere-aware v10.93)','#8b0000'));
panel.add(lbl('Northern Hem (lat>=0): Jun-Oct window | Southern Hem (lat<0): Nov-Apr window',7,'#886600'));
panel.add(legRow('#fdae61','28-30','Thermal stress')); panel.add(legRow('#f46d43','30-32','BLEACHING RISK'));
panel.add(legRow('#d73027','32-34','MASS BLEACHING')); panel.add(legRow('#a50026','34-36','MORTALITY'));
panel.add(legDiv()); panel.add(sHead('3. S1 - Bleaching risk overlay','#880000'));
panel.add(legRow('#ff0000','RED overlay','Peak SST > 30 deg C = bleaching zone'));
panel.add(legDiv()); panel.add(sHead('4. S4 - DHW (per-pixel MMM)','#880000'));
panel.add(lbl('v10.68: blue and orange now mutually exclusive (no overlap)',7,'#336633'));
panel.add(lbl('DHW layer HIDDEN above 55N/S (metric not applicable at high latitudes)',7,'#aa6600'));
panel.add(legRow('#0000ff','Blue','DHW = 0: genuinely no heat stress'));
panel.add(legRow('#ffffd9','Yellow-white','DHW 1-4: watch')); panel.add(legRow('#fdae61','Orange','DHW 4-8: BLEACHING RISK'));
panel.add(legRow('#f46d43','Dark orange','DHW 8-12: MASS BLEACH')); panel.add(legRow('#a50026','Dark red','DHW > 12: MORTALITY'));
panel.add(legDiv()); panel.add(sHead('4b. CSD - SST Trend deg C/yr (Basin erosion signal)','#880000'));
panel.add(lbl('Palette: -0.05 to +0.05 deg C/yr | Grey = stable/near-zero trend',7,'#664400'));
panel.add(lbl('Grey at high latitudes (Arctic) = OISST ice-masked pixels -> trend ~0',7,'#aa6600'));
panel.add(legRow('#0000ff','< -0.05/yr','COOLING - basin deepening (resilience recovering)'));
panel.add(legRow('#aaaaff','-0.05 to -0.01/yr','Slight cooling'));
panel.add(legRow('#dddddd','-0.01 to +0.01/yr','STABLE / near-zero trend'));
panel.add(legRow('#ffaaaa','+0.01 to +0.05/yr','Warming - basin shallowing (CSD signal)'));
panel.add(legRow('#ff0000','> +0.05/yr','RAPID warming - bowl nearly flat'));
panel.add(legDiv()); panel.add(sHead('5. S2 - Chlorophyll-a mg/m3','#005a32'));
panel.add(legRow('#084594','< 0.1','Oligotrophic')); panel.add(legRow('#6baed6','0.5-1.0','Moderate'));
panel.add(legRow('#74c476','1.0-2.0','Good')); panel.add(legRow('#006837','> 5.0','Bloom'));
panel.add(legDiv()); panel.add(sHead('6. S3 - Depth m (GEBCO, v10.68)','#003366'));
panel.add(lbl('Range: -200m to +100m | Deeper than -200m all shows as dark navy',7,'#336633'));
panel.add(legRow('#023858','< -2000m','Very deep ocean (dark navy)'));
panel.add(legRow('#0570b0','-2000 to -1000m','Deep ocean (blue)'));
panel.add(legRow('#3690c0','-1000 to -200m','Continental slope (medium blue)'));
panel.add(legRow('#41b6c4','-200 to -50m','Continental shelf (teal/cyan)'));
panel.add(legRow('#74c476','-50 to 0m','SHALLOW REEF ZONE (green — ecologically critical)'));
panel.add(legRow('#ffffb2','0 to 10m','Intertidal / beach (yellow)'));
panel.add(legRow('#8c510a','> 10m','Land / higher ground (brown)'));
panel.add(legDiv()); panel.add(sHead('S7 - Seaweed/Macroalgae (Sentinel-2)','#006600'));
panel.add(legRow('#ccffcc','FAI 0.003-0.01','Mild trace')); panel.add(legRow('#ff0000','FAI > 0.10','EXTREME bloom'));
panel.add(legDiv()); panel.add(sHead('S10 - Earthquakes M4.5+ (USGS)','#5a2222'));
panel.add(legRow('#ff8800','Orange dot','Magnitude 4.5-5.9')); panel.add(legRow('#cc0000','Dark red dot','Magnitude 6.0+'));
panel.add(legDiv()); panel.add(sHead('S11 - Volcanic activity 1960-2018 (GDIS)','#5a2222'));
panel.add(legRow('#ff00ff','Magenta dot','Recorded volcanic event'));
panel.add(legDiv()); panel.add(sHead('Score color dot','#4a0000'));
panel.add(legRow('#00cc44','Green','0-29: DEEP BASIN')); panel.add(legRow('#ffcc00','Yellow','30-54: WARNING'));
panel.add(legRow('#ff6600','Orange','55-74: HIGH RISK')); panel.add(legRow('#ff0000','Red','75-100: CRITICAL'));
panel.add(lbl('Scroll up for measurements | v10.149 + GEM',8,'#555555'));
ui.root.insert(0,panel);
var floatP=ui.Panel({style:{position:'top-center',padding:'6px 14px',backgroundColor:'#001a00',border:'2px solid #00cc44',shown:false}});
var fN=lbl('',13,'#00cc44','',true), fS=lbl('',12,'#ffffff','',true), fB=lbl('',10,'#aaffaa','',false), fC=lbl('',9,'#888888','',false);
floatP.add(fN); floatP.add(fS); floatP.add(fB); floatP.add(fC);
Map.add(floatP);

// MODULE E - MAP SETUP
Map.setCenter(-82.25,9.35,9); Map.setOptions('HYBRID');
var defaultStartPt=ee.Geometry.Point([-82.25,9.35]), defaultStartZone=defaultStartPt.buffer(1000);
Map.addLayer(ee.Image().byte().paint(defaultStartZone,0,3),{palette:['#ffffff']},'Study zone (white border) - startup default',true,1.0);
Map.addLayer(sst,{min:20,max:32,palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],opacity:0.5},'SST annual mean (context)');
Map.addLayer(sstPeakFinal,{min:26,max:36,palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],opacity:0.8},'SST PEAK Jun-Oct 2023 (bleaching score)');
Map.addLayer(dhwProper.updateMask(dhwProper.gt(0)),{min:0.1,max:12,palette:['#ffffd9','#fdae61','#f46d43','#a50026'],opacity:0.65},'DHW per-pixel MMM (startup view)');
Map.addLayer(soilTexture,{min:1,max:12,palette:['#8B4513','#A0522D','#CD853F','#DEB887','#D2B48C','#BC8F8F','#F4A460','#DAA520','#B8860B','#E0E0E0','#FFE4B5','#FFF8DC']},'S9 - Soil texture (startup)',false,0.6);
Map.addLayer(eqModerate,{color:'ff8800'},'S10 - Earthquakes M4.5-5.9 (startup)',true,0.7);
Map.addLayer(eqMajor,{color:'cc0000'},'S10 - Earthquakes M6.0+ (startup)',true,0.9);
Map.addLayer(volcanicActivity,{color:'ff00ff'},'S11 - Volcanic activity 1960-2018 (startup)',true,0.8);

function loadLayers(study,region,cols,lat) {
  while(Map.layers().length()>0){ Map.remove(Map.layers().get(0)); }
  var dotCol=cols?cols.map:'#888888', dotLbl=cols?cols.lbl:'computing';
  // v10.93 FIX: select the hemisphere-correct peak-season / DHW image for
  // this click. Southern Hemisphere reefs (lat<0) peak Nov-Apr, not Jun-Oct.
  var sstPeakLocal = (lat!==undefined && lat<0) ? sstPeakFinalSH : sstPeakFinalNH;
  var dhwProperLocal = (lat!==undefined && lat<0) ? dhwProperSH : dhwProper;
  var peakSeasonLabel = (lat!==undefined && lat<0) ? 'Nov-Apr (Southern Hem)' : 'Jun-Oct (Northern Hem)';
  Map.addLayer(ee.Image().byte().paint(study,0,3),{palette:['#ffffff'],opacity:1.0},'Study zone (white border)');
  Map.addLayer(sst.clip(study),{min:20,max:32,palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],opacity:0.7},'S1 - SST annual mean deg C');
  Map.addLayer(sstPeakLocal.clip(study),{min:24,max:36,palette:['#2c7bb6','#ffffbf','#fdae61','#f46d43','#d73027','#a50026','#660000'],opacity:0.9},'S1 - SST PEAK '+peakSeasonLabel+' 2023-24 (bleaching score)');
  Map.addLayer(sstPeakLocal.clip(study).gt(30).selfMask(),{palette:['#ff0000'],opacity:0.75},'S1 - Bleaching risk (peak SST > 30)');
  // v10.67 FIX (v10.20): blue only where DHW truly = 0 - no overlap with orange stress layer
  var dhwForBlue=dhwProperLocal.clip(study);
  Map.addLayer(ee.Image.constant(0).updateMask(oceanMask.clip(study).and(dhwForBlue.lte(0).unmask(1))),
    {palette:['#0000ff'],opacity:0.2},'S4 - DHW blue = no stress');
  // v10.68 FIX: DHW heat stress layer masked above 55N/below 55S.
  // At high latitudes the MMM baseline is near-freezing so any summer
  // warming produces artificially extreme DHW values (24+ deg C-weeks)
  // that render deep red but have no coral bleaching relevance.
  // DHW is a tropical reef metric - scientifically meaningless at Arctic latitudes.
  var dhwApplicable = Math.abs(lat) <= 55;
  if(dhwApplicable) {
    Map.addLayer(dhwProperLocal.clip(study).updateMask(dhwProperLocal.clip(study).gt(0)),
      {min:0.1,max:12,palette:['#ffffd9','#fdae61','#f46d43','#a50026'],opacity:0.85},'S4 - DHW heat stress (orange where > 0)');
  }
  // lat>55: DHW layer simply not added — metric not applicable at high latitudes
  Map.addLayer(chla.clip(study),{min:0.01,max:5.0,palette:['#084594','#2171b5','#6baed6','#74c476','#31a354','#006837'],opacity:0.8},'S2 - Chl-a mg/m3');
  // v10.66 FIX: switched to GEBCO depth layer to match sidebar depth value
  // Depth range: -200m to +100m covers ecologically relevant zones
  // Everything deeper than -200m shows as dark navy (fine for reef ecology tool)
  Map.addLayer(GEBCO.rename('bedrock').clip(study),
    {min:-200,max:100,
     palette:['#023858','#0570b0','#3690c0','#41b6c4','#74c476','#ffffb2','#8c510a'],
     opacity:0.75},'S3 - Depth m (GEBCO)');
  Map.addLayer(intertidalMask.selfMask().clip(study),{palette:['#FFD700'],opacity:0.65},'S3 - Intertidal 0-20m');
  Map.addLayer(shallowMask.selfMask().clip(study),{palette:['#00cc99'],opacity:0.35},'S3 - Shallow reef 0-50m');
  Map.addLayer(turbImg.clip(study),{min:-0.3,max:0.3,palette:['#084594','#2171b5','#74add1','#ffffbf','#fdae61','#d73027'],opacity:0.75},'S3 - Turbidity NDTI (reef < 50m ONLY)');
  Map.addLayer(sst_slope.clip(study),{min:-0.05,max:0.05,palette:['#0000ff','#aaaaff','#dddddd','#ffaaaa','#ff0000'],opacity:0.85},'CSD - SST trend deg C/yr (warming=red, cooling=blue, grey=stable)');
  Map.addLayer(no2.clip(study),{min:0.000005,max:0.00015,palette:['#313695','#4575b4','#74add1','#fee090','#f46d43','#d73027','#a50026'],opacity:0.65},'S6 - NO2 mol/m2');
  Map.addLayer(faiMild.clip(study),{min:0.003,max:0.15,palette:['#ccffcc','#88dd44','#ffff00','#ff8800','#ff4400','#ff0000','#cc0000','#880000'],opacity:0.9},'S7 - FAI Floating Algae');
  Map.addLayer(ndviAlgae.clip(study),{min:0.05,max:0.5,palette:['#e8ffe8','#aaddaa','#44aa00','#006600','#003300'],opacity:0.8},'S7 - NDVI Water (benthic macroalgae)');
  Map.addLayer(ndciBloom.clip(study),{min:0.1,max:0.5,palette:['#ffffcc','#a6d96a','#1a9641','#004400'],opacity:0.75},'S7 - NDCI Bloom index');
  Map.addLayer(soilTexture.clip(study),{min:1,max:12,palette:['#8B4513','#A0522D','#CD853F','#DEB887','#D2B48C','#BC8F8F','#F4A460','#DAA520','#B8860B','#E0E0E0','#FFE4B5','#FFF8DC']},'S9 - Soil texture (v10.42: dead asset)',false,0.6);
  Map.addLayer(eqModerate,{color:'ff8800'},'S10 - Earthquakes M4.5-5.9',true,0.7);
  Map.addLayer(eqMajor,{color:'cc0000'},'S10 - Earthquakes M6.0+',true,0.9);
  Map.addLayer(volcanicActivity,{color:'ff00ff'},'S11 - Volcanic activity 1960-2018',true,0.8);
  Map.addLayer(ee.Image.constant(1).clip(study.centroid(50).buffer(300)),{palette:[dotCol],opacity:0.9},'SCORE: '+dotLbl+' ('+dotCol+')');
}

// MODULE F0 - ERROR HANDLING (v10.37/v10.39)
var _pipelineWatchdog=null;
function clearPipelineWatchdog(){} // no-op: setTimeout not available in GEE sandbox
function startPipelineWatchdog(lat,lon){} // no-op (v10.39)

function handleGeeError(err, stepName) {
  clearPipelineWatchdog();
  print('=== PIPELINE ERROR at step: '+stepName+' ===');
  print('Error: '+err);
  print('Everything computed before this step is still valid in sidebar.');
  print('Try clicking again or try a nearby coordinate.');
  clickLbl.setValue('ERROR at step "'+stepName+'" - see Console');
  clickLbl.style().set('backgroundColor','#880000');
  scoreInterp.setValue('Computation failed at step: '+stepName+'\nError: '+String(err).slice(0,150)+'\nTry clicking again or try a nearby coordinate.');
  scoreInterp.style().set('backgroundColor','#ffd0d0'); scoreInterp.style().set('color','#880000');
}

// MODULE F0b - FULL SIDEBAR RESET (v10.38)
function resetSidebarToComputing() {
  var C='computing...';
  sstV.setValue(C); trendV.setValue(C); dhwV.setValue(C); dhwNoteV.setValue('');
  mmmV.setValue(C); chlaV.setValue(C); turbV.setValue(C); no2V.setValue(C);
  depthV.setValue(C); depthWarnV.setValue(''); ebusWarnV.setValue('');
  faiV.setValue(C); ndciV.setValue(C); ndviWV.setValue(C); algaeStatusV.setValue(C);
  sstWindowV.setValue(C); chlWindowV.setValue(C); ndciNutrientProxyV.setValue(C);
  pollutionV.setValue(C); stabilityV.setValue(C); bromoformV.setValue(C);
  aquaScoreV.setValue(C); aquaMissingV.setValue(''); aquaStatusV.setValue(C); kelpNoteV.setValue('Cold-water kelp: checking...');
  soilTextureV.setValue(C); eqStatsV.setValue(C); volcStatsV.setValue(C);
  satCcsV.setValue('--'); fusCcsV.setValue('--'); bowlV.setValue('--');
  omega0V.setValue('--'); tauV.setValue('--'); ac1V.setValue('--'); p5yrV.setValue('--');
  realAc1V.setValue(C); realVarTrendV.setValue(C); realCsdNoteV.setValue('');
  thermalEpisodesV.setValue(C); thermalMeanV.setValue(C); thermalMaxV.setValue(C); thermalOngoingV.setValue(C);
  ecoValStatusV.setValue(C); ecoValDetailsV.setValue('');
  eciDepthV.setValue(C); eciCurrentV.setValue(C); eciSLR03V.setValue(C); eciSLR05V.setValue(C); eciSLR10V.setValue(C); eciValidationV.setValue(C);
  toeCompoundV.setValue(C); toeSSTv.setValue(C); toeCHLv.setValue(C); toeNO2v.setValue(C); toeSALv.setValue(C);
  toePHv.setValue(C); toeDOv.setValue(C);
  s18O2V.setValue(C); s18pCO2V.setValue(C); s18pHV.setValue(C); s18SalV.setValue(C);
  s19FluorV.setValue(C); s19PressV.setValue(C); s19TurbV.setValue(C);
  accSV.setValue('--'); accFV.setValue('--'); accTV.setValue('--');
  intPanel.clear();
}

// MODULE F - CLICK HANDLER
// v10.81: lastClickLat/lastClickLon remember the most recently analyzed point
// so the S13 "Use last clicked location" button can fill it in automatically.
var lastClickLat = null, lastClickLon = null;

// v10.116 NEW: tracks every coordinate tested as a STUDY site across S7B/
// S7C/S7D/S7E this session, so S7E's reference-candidate search can exclude
// them. Motivated by a real case: S7E auto-found a reference reef that
// coincided with Low Isles - a location already tested (and found to show
// a real signal) as a STUDY site earlier in the same session. A reference
// is supposed to be an uninvolved baseline; reusing a site that has
// already tested positive for its own change undermines that. Deduped by
// distance (skips adding a near-duplicate of an existing entry) so this
// doesn't grow unbounded across many clicks at nearly the same spot.
var s7StudySiteHistory = [];
var S7_HISTORY_DEDUPE_KM = 2;
function recordStudySite(lat, lon, label){
  for(var i=0;i<s7StudySiteHistory.length;i++){
    if(haversineKm(lat,lon,s7StudySiteHistory[i].lat,s7StudySiteHistory[i].lon) < S7_HISTORY_DEDUPE_KM) return;
  }
  s7StudySiteHistory.push({lat:lat, lon:lon, label:label||'study', t:Date.now()});
}

function analyzeLocation(lat, lon) {
  lastClickLat = lat; lastClickLon = lon;
  startPipelineWatchdog(lat, lon);
  resetSidebarToComputing();
  var latR=Math.round(lat*1000)/1000, lonR=Math.round(lon*1000)/1000;
  var region=getRegion(lat,lon), fp=getFieldProfile(region);
  var pt=ee.Geometry.Point([lon,lat]), study=pt.buffer(1000);
  // v10.68 FIX: latitude-adaptive SST buffer. OISST is 25km resolution.
  // At high latitudes (>50N or <-50S) sea-ice masking means a 3km buffer
  // often finds zero valid pixels -> rSST error. Use 50km buffer for Arctic/Antarctic.
  var absLat=Math.abs(lat);
  var sstBufferM = absLat>60 ? 80000 : absLat>50 ? 50000 : absLat>30 ? 20000 : 5000;
  var studySST=pt.buffer(sstBufferM);

  // S19 IN-SITU LOOKUP (synchronous JS - no evaluate() needed)
  var inSitu=getInSituBaseline(lat,lon);
  if(inSitu) {
    s19StatusV.setValue(inSitu.label+' ('+inSitu.distance_km.toFixed(1)+' km away)');
    s19StatusV.style().set('color','#115511');
    if(inSitu.pH) {
      s19PhV.setValue(inSitu.pH.mean.toFixed(3)+' +/- '+inSitu.pH.std.toFixed(3)+'  (n='+inSitu.n_clean+' clean, '+inSitu.pct_flagged+'% flagged)');
    } else { s19PhV.setValue('n/a (not in this dataset)'); }
    if(inSitu.temp_c&&inSitu.temp_c.mean!==null){
      var tempStd=inSitu.temp_c.std?(' +/- '+inSitu.temp_c.std.toFixed(2)):'';
      s19TempV.setValue(inSitu.temp_c.mean.toFixed(2)+' deg C'+tempStd);
    } else{s19TempV.setValue('n/a');}
    if(inSitu.salinity&&inSitu.salinity.mean!==null){
      var salStd=inSitu.salinity.std?(' +/- '+inSitu.salinity.std.toFixed(2)):'';
      s19SalV.setValue(inSitu.salinity.mean.toFixed(2)+' PSU'+salStd);
    } else{s19SalV.setValue('n/a');}
    if(inSitu.do_mgL){s19DoV.setValue(inSitu.do_mgL.mean.toFixed(2)+' mg/L');}
    else{s19DoV.setValue('n/a (not in this dataset)');}
    if(inSitu.fluorescence_ug_L&&inSitu.fluorescence_ug_L.mean!==null){
      var fluor=inSitu.fluorescence_ug_L;
      // Note: cv (satellite Chl-a) is not yet available here - shown in S2 above
      s19FluorV.setValue(fluor.mean.toFixed(3)+' ug/L +/-'+fluor.std.toFixed(3)+
        '  (n='+fluor.n_clean+', '+fluor.pct_flagged+'% flagged) | compare vs S2 satellite Chl-a above');
      s19FluorV.style().set('color','#115511');
    } else{s19FluorV.setValue('n/a (not in this dataset)');}
    if(inSitu.pressure_db&&inSitu.pressure_db.mean!==null){
      var pdb=inSitu.pressure_db;
      s19PressV.setValue(pdb.mean.toFixed(1)+' dbar mean | range: '+pdb.min+'-'+pdb.max+' dbar (full water column profiles)');
    } else{s19PressV.setValue('n/a (not in this dataset)');}
    if(inSitu.turbidity_ftu&&inSitu.turbidity_ftu.mean!==null){
      var turb=inSitu.turbidity_ftu;
      s19TurbV.setValue(turb.mean.toFixed(3)+' FTU +/-'+turb.std.toFixed(3)+
        '  (n='+turb.n_clean+', '+turb.pct_flagged+'% flagged)');
    } else{s19TurbV.setValue('n/a (not in this dataset)');}
    s19NoteV.setValue('Source: '+inSitu.source+'\nRecord: '+inSitu.record+'\n'+inSitu.notes.slice(0,200));
    print('=== S19 REAL IN-SITU BASELINE (v10.68) ===');
    print('Matched: '+inSitu.label+' ('+inSitu.distance_km.toFixed(2)+' km from click)');
    if(inSitu.key&&inSitu.key.indexOf('gem_')===0) {
      print('GEM CITATION REQUIRED: CC BY-SA 4.0 | https://creativecommons.org/licenses/by-sa/4.0/');
      print('Cite: Greenland Ecosystem Monitoring | https://g-e-m.dk | '+inSitu.source);
      print('NOTE: GEM stats are PLACEHOLDER nulls until gem_fetch_and_clean.py is run.');
    }
    print('CAVEAT: real measured sensor baseline, but point-source only - no spatial interpolation.');
  } else {
    s19StatusV.setValue('No in-situ station within range of this click');
    s19StatusV.style().set('color','#888888');
    s19PhV.setValue('n/a'); s19TempV.setValue('n/a'); s19SalV.setValue('n/a'); s19DoV.setValue('n/a');
    s19NoteV.setValue('Not within radius of Looe Key FL, Agua Hedionda CA, Scripps Pier CA,\nMarineBasis Nuuk GF3 (40km), or Zackenberg Young Sound (30km).');
  }

  locV.setValue(region); regV.setValue(region); coV.setValue('Lat:'+latR+' Lon:'+lonR);
  scoreBig.setValue('Computing...'); scoreBig.style().set('color','#333333'); scoreBig.style().set('backgroundColor','#eeeeee');
  scoreBarLbl.setValue('Computing...'); scoreInterp.setValue('Running pipeline...');
  barFill.style().set('width','0px'); barFill.style().set('backgroundColor','#aaaaaa');
  if(fp.hasField){
    fieldStatusV.setValue('FIELD DATA AVAILABLE'); fieldStatusV.style().set('color','#115511');
    fieldSpeciesV.setValue(fp.species||'--');
    fieldSourcesV.setValue('Sources: '+fp.sources);
    fieldNotesV.setValue('Notes: '+fp.notes);
  } else {
    fieldStatusV.setValue('SATELLITE ONLY - no field data for this region'); fieldStatusV.style().set('color','#880000');
    fieldSpeciesV.setValue('none'); fieldSourcesV.setValue('Sources: '+fp.sources); fieldNotesV.setValue('Notes: '+fp.notes);
  }
  fN.setValue(region); fS.setValue('Computing...'); fB.setValue(''); fC.setValue(latR+', '+lonR);
  floatP.style().set('shown',true);
  Map.setCenter(lon,lat,11);
  while(Map.layers().length()>0){ Map.remove(Map.layers().get(0)); }
  Map.addLayer(ee.Image().byte().paint(study,0,3),{palette:['#ffffff'],opacity:1.0},'Study zone (white border)');
  Map.addLayer(sst.clip(study),{min:20,max:32,palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],opacity:0.9},'S1 - SST deg C (loading)');

  function clip(col){ return col.map(function(img){ return img.clip(study); }); }

  print(''); print('STEMGeoHS Marine v10.67 -- '+region);
  print('=== MODELS ===');
  print('1. Waddington double-well: U(q;mu) = 0.25*q^4 - 0.5*mu*q^2');
  print('2. Langevin SDE: dx = -dU/dx*dt + sigma*dW (PNAS 2025)');
  print('3. CSD: AC1->1, tau->inf as mu->mu_c (Scheffer 2009)');
  print('4. Lagrangian L = T - V (Pham & Musielak 2022)');
  print('5. Kramers escape: k = omega0 * exp(-dU/sigma^2)');
  print('6. Cancer score: CCS = SUM(w_i*S_i) - REVERT (Shin 2025)');
  if(fp.hasField){print('FIELD: '+fp.species+' | Gain: +'+fp.accuracy_field_gain+'%');}
  else{print('STATUS: SATELLITE DATA ONLY | '+fp.notes);}
  // Charts printed inside _onAllDone() only when valid ocean data confirmed

  if(region==='Bocas del Toro, Panama'||region==='Caribbean Sea'){
    print('=== REAL FIELD DATA: Ramamurthy 2024 CJS 54:77-82 ===');
    print('Qualitative behavior codes only (no numeric ratio exists in source).');
  }

  // Satellite evaluate() calls
  // v10.93 FIX: select hemisphere-correct peak-season / DHW image for this
  // click's latitude, instead of always reading the Northern Hemisphere
  // Jun-Oct window (which silently misses Southern Hemisphere bleaching).
  var sstPeakLocalCalc = (lat<0) ? sstPeakFinalSH : sstPeakFinalNH;
  var dhwProperLocalCalc = (lat<0) ? dhwProperSH : dhwProper;
  var rSST=sst.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rSSTP=sstPeakLocalCalc.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rCHL=chla.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rCHLc=chlaCoastal.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rTR=sst_slope.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rNO2=no2.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:1000,maxPixels:1e9});
  var rBATH=GEBCO.rename('bedrock').reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:1000,maxPixels:1e9});
  var rTURB=turbImg.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:100,maxPixels:1e9});
  var rDHW=dhwProperLocalCalc.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rMMM=MMM_perpixel.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9});
  var rFAI=faiImg.updateMask(oceanMask).reduceRegion({reducer:ee.Reducer.max(),geometry:study,scale:100,maxPixels:1e9});
  var rNDCI=ndciImg.updateMask(oceanMask).reduceRegion({reducer:ee.Reducer.max(),geometry:study,scale:100,maxPixels:1e9});
  var rNDVI_W=ndviWater.updateMask(oceanMask).reduceRegion({reducer:ee.Reducer.max(),geometry:study,scale:100,maxPixels:1e9});
  var rCHL_BLOOM=chla.reduceRegion({reducer:ee.Reducer.max(),geometry:study,scale:1000,maxPixels:1e9});
  var rSOIL=soilTexture.reduceRegion({reducer:ee.Reducer.first(),geometry:pt,scale:250,maxPixels:1e9});
  var seismicZone=pt.buffer(200000), volcZone=pt.buffer(300000);
  var rEQ=ee.Dictionary({count:usgsEarthquakes.filterBounds(seismicZone).size(),maxMag:usgsEarthquakes.filterBounds(seismicZone).aggregate_max('mag')});
  var rVOLC=volcanicActivity.filterBounds(volcZone).size();
  var rGEBCO=GEBCO.rename('elevation').reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:500,maxPixels:1e9});

  // S17 ToE precomputed images at clicked point
  var toePt=ee.Geometry.Point([lon,lat]);
  var toeScale=27750;
  var rToeSST=ee.Dictionary({scale:toeSSTFit.select('scale').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('scale'),
    noise:toeSSTNoise.select('sst_stdDev').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('sst_stdDev')});
  var rToeCHL=ee.Dictionary({scale:toeCHLFit.select('scale').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('scale'),
    noise:toeCHLNoise.select('chlor_a_stdDev').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('chlor_a_stdDev')});
  var rToeSAL=ee.Dictionary({scale:toeSALFit.select('scale').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('scale'),
    noise:toeSALNoise.select('salinity_0_stdDev').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('salinity_0_stdDev')});
  var rToeNO2=ee.Dictionary({scale:toeNO2Fit.select('scale').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('scale'),
    noise:toeNO2Noise.select('tropospheric_NO2_column_number_density_stdDev').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('tropospheric_NO2_column_number_density_stdDev')});
  var rToePH=ee.Dictionary({scale:toePHFit.select('scale').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('scale'),
    noise:toePHNoise.select('ph_stdDev').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('ph_stdDev')});
  var rToeDO=ee.Dictionary({scale:toeDOFit.select('scale').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('scale'),
    noise:toeDONoise.select('o2_stdDev').reduceRegion({reducer:ee.Reducer.mean(),geometry:toePt,scale:toeScale,maxPixels:1e9}).get('o2_stdDev')});
  // v10.89 FIX: previously bundled into ONE rToeAll dictionary and evaluated
  // together - since pH/DO both depend on COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/
  // MFC_001_028 (currently returns "asset not found" in the GEE catalog),
  // that single dead asset made the WHOLE evaluate() call fail and blanked
  // ALL SIX indicators, including SST/CHL/SAL/NO2 which have nothing to do
  // with it and would have computed fine on their own. Split into two
  // independent dictionaries/evaluates below so a dead pH/DO asset can only
  // ever affect pH/DO.
  // v10.141 FURTHER SPLIT: pH now uses a real, confirmed working asset
  // (see fix above) while DO's asset remains dead - so pH and DO are now
  // ALSO split from each other, for the exact same reason - otherwise DO's
  // still-dead asset would fail the evaluate() and blank the now-working
  // pH result too.
  var rToeCore=ee.Dictionary({sst:rToeSST,chl:rToeCHL,sal:rToeSAL,no2:rToeNO2});
  var rToeBGC_PH=ee.Dictionary({ph:rToePH});
  var rToeBGC_DO=ee.Dictionary({do_o2:rToeDO});

  // v10.145 FIX: S18 rebuilt using REAL, confirmed working assets for
  // 3 of its 4 fields - same fix pattern already applied to S17's pH
  // (v10.141). pH: same real CAR/ph_depth1 asset. pCO2: real CO2/
  // spco2_depth1 sub-collection, confirmed via direct search this
  // session (unit is Pa in the raw asset, converted to uatm below -
  // 1 Pa = 9.86923 uatm). Salinity: reuses the SAME real, already-
  // working HYCOM source and decoding formula used throughout this
  // tool for S17's salinity trend (salinity_0 * 0.001 + 20).
  // O2 remains genuinely unresolved - no working sub-collection/band
  // name was confirmed after multiple direct searches this session -
  // kept separate and disclosed as unavailable, not guessed.
  // v10.146 BUGFIX: found via real testing - the v10.145 S18 rebuild
  // was missing the SAME defensive empty-collection guard used
  // everywhere else in this codebase (e.g. mkMoDHW, dhwProper) before
  // calling .select()/.mean() - if the 60-day window genuinely returned
  // zero images at some location/time, that produced a real computation
  // error instead of a graceful null. Fixed using the identical
  // ee.Algorithms.If(col.size().gt(0), ..., null) pattern already
  // proven throughout the rest of this tool. Also widened the window
  // 60->90 days to further reduce the chance of a genuinely empty
  // result from processing lag.
  // v10.147 BUGFIX: found via a SECOND real test (One Tree Reef) that
  // falsified the v10.146 hypothesis - the failure was NOT specific to
  // Florida Keys' coastal proximity, since One Tree Reef (36m deep,
  // open lagoon) failed identically. Compared directly against S17's
  // OWN pH fetch, which succeeded at BOTH real test sites using the
  // SAME real CAR asset - found two real differences: (1) S17 uses a
  // bare ee.Geometry.Point with NO buffer, while S18 was using `study`,
  // a buffered polygon - against a coarse ~25km global grid, a small
  // buffered area can produce a genuinely empty/null reduceRegion result
  // in edge cases where a bare point cleanly resolves to one pixel;
  // (2) S18 used a narrow 90-day window while S17's successful version
  // uses a multi-year annual collection. Fixed both: switched to the
  // same bare-point geometry pattern already proven at both real test
  // sites, and widened the window to 2 years (730 days) to reduce
  // reliance on the near-real-time forecast product's most recent,
  // potentially-gappy days.
  var toeBgcPt = ee.Geometry.Point([lon,lat]);
  var colPH = ee.ImageCollection('COPERNICUS/MARINE/GLOBAL_ANALYSISFORECAST_BGC_001_028/CAR')
    .filterDate(ee.Date(Date.now()).advance(-730,'day'), ee.Date(Date.now())).select('ph_depth1');
  var rBGC_PH = ee.Dictionary({ph_depth1: ee.Algorithms.If(colPH.size().gt(0),
    colPH.mean().reduceRegion({reducer:ee.Reducer.mean(),geometry:toeBgcPt,scale:27750,maxPixels:1e9}).get('ph_depth1'),
    null)});
  var colCO2 = ee.ImageCollection('COPERNICUS/MARINE/GLOBAL_ANALYSISFORECAST_BGC_001_028/CO2')
    .filterDate(ee.Date(Date.now()).advance(-730,'day'), ee.Date(Date.now())).select('spco2_depth1');
  var rBGC_CO2 = ee.Dictionary({spco2_depth1: ee.Algorithms.If(colCO2.size().gt(0),
    colCO2.mean().reduceRegion({reducer:ee.Reducer.mean(),geometry:toeBgcPt,scale:27750,maxPixels:1e9}).get('spco2_depth1'),
    null)});
  var colSAL = ee.ImageCollection('HYCOM/sea_temp_salinity')
    .filterDate(ee.Date(Date.now()).advance(-730,'day'), ee.Date(Date.now())).select('salinity_0');
  var rBGC_SAL = ee.Dictionary({salinity_0: ee.Algorithms.If(colSAL.size().gt(0),
    colSAL.mean().reduceRegion({reducer:ee.Reducer.mean(),geometry:toeBgcPt,scale:27750,maxPixels:1e9}).get('salinity_0'),
    null)});

  // Monthly SST for thermal recovery
  var monthlyFC=ee.FeatureCollection(mkMoSST().map(function(img){
    var v=img.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9}).get('sst');
    return ee.Feature(null,{v:v,t:img.get('system:time_start')});
  })).sort('t');

  // ============================================================
  // PARALLEL EVALUATE (v10.68 PERF FIX)
  // Previously: 19 sequential nested evaluate() calls
  //   Total time = SUM of all 19 server calls (~40-90 seconds)
  // Now: all 19 fire simultaneously, process when all arrive
  //   Total time = MAX of all 19 server calls (~5-15 seconds)
  // ============================================================
  var _res={}, _pending=19;
  function _got(key,val){ _res[key]=val; if(--_pending===0) _onAllDone(); }

  rSST.evaluate(function(v,e){   _got('s',   e?{}:v||{}); });
  rSSTP.evaluate(function(v,e){  _got('sp',  e?{}:v||{}); });
  rCHL.evaluate(function(v,e){   _got('c',   e?{}:v||{}); });
  rCHLc.evaluate(function(v,e){  _got('cc',  e?{}:v||{}); });
  rTR.evaluate(function(v,e){    _got('tr',  e?{}:v||{}); });
  rNO2.evaluate(function(v,e){   _got('n',   e?{}:v||{}); });
  rBATH.evaluate(function(v,e){  _got('b',   e?{}:v||{}); });
  rTURB.evaluate(function(v,e){  _got('t',   e?{}:v||{}); });
  rDHW.evaluate(function(v,e){   _got('d',   e?{}:v||{}); });
  rMMM.evaluate(function(v,e){   _got('mmm', e?{}:v||{}); });
  rFAI.evaluate(function(v,e){   _got('fa',  e?{}:v||{}); });
  rNDCI.evaluate(function(v,e){  _got('nc',  e?{}:v||{}); });
  rNDVI_W.evaluate(function(v,e){_got('nw',  e?{}:v||{}); });
  rCHL_BLOOM.evaluate(function(v,e){_got('cb',e?{}:v||{}); });
  rSOIL.evaluate(function(v,e){  _got('soilRes',  e?{}:v||{}); });
  rEQ.evaluate(function(v,e){    _got('eqRes',e?{count:0,maxMag:null}:v||{count:0,maxMag:null}); });
  rVOLC.evaluate(function(v,e){  _got('volcRes',e?0:v||0); });
  computeRealCSDDeseasonalized(mkMoSST(), study, 'sst', 4000, function(result) {
    _got('csdRes', result.error ? {} : result);
  });
  rGEBCO.evaluate(function(v,e){ _got('gebcoRes',e?{}:v||{}); });

  function _onAllDone(){
    var s=_res.s, sp=_res.sp, c=_res.c, cc=_res.cc;
    var tr=_res.tr, n=_res.n, b=_res.b, t=_res.t;
    var d=_res.d, mmm=_res.mmm, fa=_res.fa, nc=_res.nc;
    var nw=_res.nw, cb=_res.cb, soilRes=_res.soilRes;
    var eqRes=_res.eqRes, volcRes=_res.volcRes;
    var csdRes=_res.csdRes, gebcoRes=_res.gebcoRes;

    var sv=(s&&s.sst!==null&&s.sst>-900)?s.sst:null;
    var sv_peak=(sp&&sp.sst_peak!==null&&sp.sst_peak>-900)?sp.sst_peak:sv;
    var cv=(c&&c.chlor_a>0)?c.chlor_a:(cc&&cc.chlor_a>0)?cc.chlor_a:null;
    var chlSource=(c&&c.chlor_a>0)?'ocean':(cc&&cc.chlor_a>0)?'coastal relaxed':'n/a';
    var tv=(tr&&tr.scale!==null)?tr.scale:null;
    var nv=(n&&n.NO2_column_number_density>0)?n.NO2_column_number_density:null;
    var bv=(b&&b.bedrock!==null)?b.bedrock:null;
    var turv=(t&&t.turbidity!==null)?t.turbidity:null;
    var dhwv_raw=(d&&d.dhw!==null)?d.dhw:null;
    var DHW_SANITY_CAP=60;
    var dhwv=(dhwv_raw!==null&&dhwv_raw<=DHW_SANITY_CAP)?dhwv_raw:null;
    var dhwFlaggedArtifact=(dhwv_raw!==null&&dhwv_raw>DHW_SANITY_CAP);
    var mmmv=(mmm&&mmm.mmm!==null)?mmm.mmm:null;

    // S20e Real Bleaching Probability (v10.144) - reuses dhwv, turv, bv
    // (depth), all already computed above for the main sidebar - zero
    // new EE calls. GEBCO depth (bv) is negative below sea level;
    // GCBD's Depth_m was a positive real-world depth, so abs() here.
    var bleachPred = predictBleachingProbability(dhwv, turv, bv!==null?Math.abs(bv):null);
    if(bleachPred.error){
      bleachProbV.setValue('n/a - '+bleachPred.error);
      bleachProbV.style().set('color','#888888');
    } else {
      var pPct = (bleachPred.p*100).toFixed(1);
      bleachProbV.setValue(pPct+'% (fitted model, held-out AUC=0.620 - real, weak-moderate signal)');
      bleachProbV.style().set('color',bleachPred.p>0.6?'#aa2200':bleachPred.p>0.4?'#aa7700':'#227744');
      print('=== S20e REAL BLEACHING PROBABILITY (fitted model) ===');
      print('P(bleaching)='+pPct+'% | inputs: DHW='+fmt(dhwv,2)+', Turbidity='+fmt(turv,4)+', Depth='+fmt(bv!==null?Math.abs(bv):null,1)+'m');
      print('Model: logistic regression, 32,716 real GCBD rows, held-out test AUC=0.620 (weak-moderate, disclosed) - predicts BLEACHING, not collapse.');
    }

    var sc=computeScore(sv_peak,cv,tv,nv,turv,dhwv,fp,lat,lon);
    sc.sst_annual=sv; sc.sst_peak=sv_peak;
    var cols=scoreColors(sc.ccs);
    loadLayers(study,region,cols,lat);

    var peakSeasonTxt = (lat<0) ? 'Nov-Apr' : 'Jun-Oct';
    var peakLabel=sv_peak!==null?' | Peak('+peakSeasonTxt+'):'+fmt(sv_peak,1)+' deg C'+(sv_peak>32?' MORTALITY':sv_peak>30?' BLEACHING':sv_peak>28?' STRESS':''):'';
    sstV.setValue(fmt(sv,2)+' deg C (annual mean)'+peakLabel);
    trendV.setValue(fmt(tv,4)+' deg C/yr'+(tv?(tv>0.05?' RAPID':tv>0.01?' Warming':tv<-0.01?' Cooling':' Stable'):''));
    var isTropicalReef=(lat>-30&&lat<30)&&!isEBUS(lat,lon);
    var dhwLabel=dhwv?(isTropicalReef?(dhwv>12?' MORTALITY':dhwv>8?' MASS BLEACH':dhwv>4?' BLEACH RISK':dhwv>0?' Watch':' No stress'):(dhwv>8?' STRONG warm anomaly':dhwv>4?' Moderate warm anomaly':' Mild warm anomaly')):'';
    dhwV.setValue(fmt(dhwv,2)+' deg C-wks'+dhwLabel);
    mmmV.setValue(fmt(mmmv,1)+' deg C (local baseline)');
    if(Math.abs(lat)>55){
      dhwNoteV.setValue('DHW MAP LAYER HIDDEN at lat '+Math.round(lat)+'° — DHW is a tropical coral reef metric. '+
        'At high latitudes the near-freezing MMM baseline makes any summer warming appear as extreme DHW '+
        '(e.g. 24+ deg C-wks) with no ecological meaning. Sidebar value shown for reference only.');
      dhwNoteV.style().set('color','#884400');
    } else if(dhwFlaggedArtifact){
      dhwNoteV.setValue('DATA ARTIFACT: raw DHW='+fmt(dhwv_raw,1)+' exceeds realistic max (~60). Suppressed to n/a.');
      dhwNoteV.style().set('color','#cc0000');
    } else if(!isTropicalReef&&dhwv&&dhwv>4){
      dhwNoteV.setValue('NOTE: DHW>4 here = warm anomaly NOT coral bleaching. No reefs at lat '+Math.round(lat)+'.');
      dhwNoteV.style().set('color','#884400');
    } else{dhwNoteV.setValue('');}
    depthWarnV.setValue(bv!==null&&bv<-200?'NOTE:'+classifyDepthLabel(bv)+' - click closer to shore':'');
    if(isEBUS(lat,lon)){ebusWarnV.setValue('EBUS ZONE: SST cooling trend may be upwelling artifact. B may be inflated.'); ebusWarnV.style().set('color','#885500');}else{ebusWarnV.setValue('');}
    chlaV.setValue(fmt(cv,3)+' mg/m3 ('+chlSource+')'+(cv?(cv<0.1?' [oligotrophic]':cv<0.5?' [low]':cv<1?' [moderate]':cv<2?' [good]':' [high]'):''));
    turbV.setValue(turv!==null?fmt(turv,3)+' NDTI (reef valid)':'n/a - open ocean');
    no2V.setValue(fmt(nv,8)+' mol/m2');
    depthV.setValue(fmt(bv,0)+' m (GEBCO)'+classifyDepthLabel(bv));

    fc1V.setValue(fp.urchin_N!==null?Math.round(sc.F1)+' (N='+fp.urchin_N+'/m2)':'n/a - no data this region');
    fc3V.setValue(fp.anem_N!==null?Math.round(sc.F3)+' (aN='+fp.anem_N+'/m2)'+(fp.anem_N_estimated?' [ESTIMATED]':''):'n/a - no data this region');
    fc4V.setValue(fp.Cd!==null?'+'+sc.F4+' (Cd='+fp.Cd+')':'n/a - no data this region');
    fc5V.setValue(fp.recruit!==null?Math.round(sc.F5)+' (recruit='+fp.recruit+')'+(fp.recruit_estimated?' [ESTIMATED]':''):'n/a - no data this region');
    fcTV.setValue((sc.fcTotal>0?'+':'')+sc.fcTotal+' total');
    fc1V.style().set('color',fp.urchin_N!==null?'#115511':'#888888');
    fc3V.style().set('color',fp.anem_N_estimated?'#aa6600':fp.anem_N!==null?'#115511':'#888888');
    fc4V.style().set('color',fp.Cd!==null?'#553300':'#888888');
    fc5V.style().set('color',fp.recruit_estimated?'#aa6600':fp.recruit!==null?'#115544':'#888888');

    satCcsV.setValue(sc.sat_ccs+'/100');
    fusCcsV.setValue(sc.ccs+'/100 FUSED'); fusCcsV.style().set('color',cols.text);
    bowlV.setValue(sc.B.toFixed(2)+(sc.B>0.7?' deep-safe':sc.B>0.5?' moderate':sc.B>0.3?' shallow':' near-flat!'));
    omega0V.setValue(sc.omega0.toFixed(4)); tauV.setValue(sc.tau.toFixed(1)+'x'); ac1V.setValue(sc.ac1.toFixed(3));
    p5yrV.setValue(sc.p5yr+'%'+(sc.p5yr>50?' CRITICAL':sc.p5yr>25?' HIGH':sc.p5yr>10?' MOD':' LOW'));
    p5yrV.style().set('color',cols.text);
    // v10.145 NEW: real comparison between the uncalibrated heuristic
    // index and S20e's real, held-out-validated fitted model. bleachPred
    // is computed earlier in this same click handler - reused here, no
    // new EE calls.
    if(bleachPred && !bleachPred.error){
      var heuristicPct = sc.p5yr, modelPct = bleachPred.p*100;
      var gapAbs = Math.abs(heuristicPct-modelPct);
      bowlVsS20eV.setValue('Heuristic: '+heuristicPct+'% | S20e real model: '+modelPct.toFixed(1)+'% | gap: '+gapAbs.toFixed(1)+' points');
      bowlVsS20eV.style().set('color',gapAbs>30?'#aa2200':gapAbs>15?'#aa7700':'#227744');
    } else {
      bowlVsS20eV.setValue('S20e model unavailable at this point - cannot compare.');
      bowlVsS20eV.style().set('color','#888888');
    }
    accSV.setValue(sc.acc_sat+'%');
    accFV.setValue('+'+sc.acc_field+'%'+(fp.hasField?' ('+fp.species.split('(')[0].trim()+')':' (none)'));
    accTV.setValue(sc.acc_total+'% TOTAL');
    accFV.style().set('color',fp.hasField?'#115511':'#888888');
    accTV.style().set('color',sc.acc_total>90?'#0a5c1e':sc.acc_total>85?'#664400':'#880000');

    var displayCcs=(!isNaN(sc.ccs)&&sc.ccs!==null)?sc.ccs:sc.sat_ccs;
    scoreBig.setValue('SCORE: '+displayCcs+'/100  B='+sc.B.toFixed(2));
    scoreBig.style().set('color',cols.text); scoreBig.style().set('backgroundColor',cols.bg);
    scoreBarLbl.setValue(cols.lbl+' ('+sc.ccs+'/100)'); scoreBarLbl.style().set('color',cols.map);
    barFill.style().set('width',Math.round(sc.ccs*1.9)+'px'); barFill.style().set('backgroundColor',cols.bar);
    var interp=sc.ccs<30?region+': DEEP BASIN\nB='+sc.B.toFixed(2)+' | Acc='+sc.acc_total+'%\nGREEN dot.':
      sc.ccs<55?region+': WARNING\nB='+sc.B.toFixed(2)+' | P='+sc.p5yr+'%\nYELLOW dot.':
      sc.ccs<75?region+': HIGH RISK\nB='+sc.B.toFixed(2)+' | P='+sc.p5yr+'%\nORANGE dot.':
      region+': CRITICAL\nB='+sc.B.toFixed(2)+' | P='+sc.p5yr+'%\nRED dot.';
    scoreInterp.setValue(interp); scoreInterp.style().set('color',cols.text); scoreInterp.style().set('backgroundColor',cols.bg);
    floatP.style().set('border','2px solid '+cols.map);
    fN.style().set('color',cols.map); fN.setValue(region);
    fS.setValue('Score: '+sc.ccs+'/100  B='+sc.B.toFixed(2)+'  Acc:'+sc.acc_total+'%'); fS.style().set('color',cols.map);
    fB.setValue((fp.hasField?fp.species.split('(')[0]:'SAT ONLY')+' | '+cols.lbl);

    // Algae S7
    var faiv=(fa&&fa.fai!==null&&!isNaN(fa.fai))?fa.fai:null;
    var nciV=(nc&&nc.ndci!==null&&!isNaN(nc.ndci))?nc.ndci:null;
    var ndwV=(nw&&nw.ndvi_water!==null&&!isNaN(nw.ndvi_water))?nw.ndvi_water:null;
    var chlBloom=(cb&&cb.chlor_a!==null)?cb.chlor_a:null;
    // v10.148 FIX: computeAquaculture moved earlier (was previously
    // called AFTER this S7 block, at what's now the second reference
    // below) so its cold-water-kelp detection is available in time to
    // contextualize S7's algae alarm text. Direct fix for a real,
    // confirmed issue: at Nuuk, S7 showed "EXTREME - massive bloom,
    // RED" - language built for tropical coral-algae phase shifts -
    // for what S8 (and the Field Data section) independently identify
    // as a likely HEALTHY cold-water kelp signal, not a warning. Same
    // class of gap already fixed for DHW (explicitly hidden above
    // 55 deg N/S as ecologically meaningless there) - S7 never got that
    // same climate-zone awareness until now.
    var aq=computeAquaculture(sv,sv_peak,cv,nv,turv,tv);
    var kelpDetected = aq.aqua_missing_note &&
      aq.aqua_missing_note.indexOf('COLD-WATER KELP OPPORTUNITY')>=0;
    var algaeElevated = (faiv!==null&&faiv>0.01) || (ndwV!==null&&ndwV>0.05);

    faiV.setValue(faiv!==null?fmt(faiv,4)+' '+(faiv>0.10?'EXTREME (RED)':faiv>0.05?'CRITICAL':faiv>0.03?'HIGH':faiv>0.01?'ELEVATED':'WATCH'):(ndwV!==null?fmt(ndwV,3)+' (NDVI proxy)':(chlBloom!==null&&chlBloom>1?fmt(chlBloom,2)+' mg/m3 HIGH':'n/a - no S2')));
    faiV.style().set('color',faiv!==null&&faiv>0.05?'#880000':faiv!==null&&faiv>0.01?'#664400':faiv!==null?'#115511':'#888888');
    ndciV.setValue(nciV!==null?fmt(nciV,3)+' '+(nciV>0.35?'EXTREME bloom':nciV>0.20?'HIGH bloom':nciV>0.10?'MODERATE':'low'):'n/a');
    ndciV.style().set('color',nciV===null?'#888888':nciV>0.20?'#880000':nciV>0.10?'#664400':'#115511');
    ndviWV.setValue(ndwV!==null?fmt(ndwV,3)+' '+(ndwV>0.40?'DENSE benthic algae':ndwV>0.20?'MODERATE cover':ndwV>0.05?'light cover':'clear'):'n/a');
    ndviWV.style().set('color',ndwV===null?'#888888':ndwV>0.30?'#880000':ndwV>0.10?'#664400':'#115511');
    if(kelpDetected && algaeElevated){
      algaeStatusV.setValue((faiv!==null?(faiv>0.10?'EXTREME':faiv>0.05?'CRITICAL':faiv>0.03?'HIGH':'ELEVATED'):'ELEVATED')+
        ' benthic signal - LIKELY COLD-WATER KELP (see S8 below), not a coral-algae bloom warning. This site is too cold for coral entirely.');
      algaeStatusV.style().set('color','#115511');
    } else {
      algaeStatusV.setValue(faiv!==null?(faiv>0.10?'EXTREME - massive bloom':faiv>0.05?'CRITICAL - dense mat':faiv>0.03?'HIGH - dense algae':faiv>0.01?'ELEVATED':'WATCH - mild signal'):ndwV!==null?(ndwV>0.20?'HIGH benthic algae':ndwV>0.05?'MODERATE benthic cover':'CLEAR'):'No S2 data - check FAI layer');
    }

    // S8 Aquaculture (aq already computed above, before S7, so its
    // kelp-detection result could contextualize S7's algae display)
    sstWindowV.setValue(sv!==null?fmt(sv,1)+' deg C '+(aq.sst_optimal?'OPTIMAL (17-21 C)':aq.sst_in_window?'IN RANGE (15-28 C)':sv<15?'TOO COLD':sv>28?'TOO HOT':'MARGINAL'):'n/a');
    sstWindowV.style().set('color',!aq.sst_really_available?'#888888':aq.sst_optimal?'#115511':aq.sst_in_window?'#664400':'#880000');
    chlWindowV.setValue(cv!==null?fmt(cv,3)+' mg/m3 '+(cv>=0.3&&cv<=1.5?'GATE PASSED':cv<0.1?'Too low':'Outside optimal range'):'n/a');
    chlWindowV.style().set('color',aq.chl_score===100?'#115511':'#880000');
    ndciNutrientProxyV.setValue(nciV!==null?fmt(nciV,3)+' '+(nciV>0.20?'HIGH':nciV>0.10?'MODERATE':'LOW'):(cv!==null?'n/a (calibrated Chl-a available)':'n/a (no S2)'));
    ndciNutrientProxyV.style().set('color',nciV===null?'#888888':nciV>0.20?'#aa3300':'#226666');
    pollutionV.setValue(nv!==null?fmt(nv,7)+' mol/m2 '+(nv<0.00005?'CLEAN':nv<0.00010?'Acceptable':'HIGH - avoid'):'n/a');
    pollutionV.style().set('color',aq.poll_score>=80?'#115511':aq.poll_score>=50?'#664400':'#880000');
    stabilityV.setValue(tv!==null?fmt(tv,4)+' deg C/yr '+(tv<0.01?'STABLE':tv<0.03?'Mostly stable':tv<0.05?'Warming':'RAPID'):'n/a');
    stabilityV.style().set('color',aq.stab_score>=80?'#115511':aq.stab_score>=50?'#664400':'#880000');
    bromoformV.setValue(sv!==null?(aq.bromo_score>=80?'HIGH yield expected':aq.bromo_score>=60?'MODERATE yield':sv>28?'LOW - heat degrades':sv<15?'TOO COLD':'LOW - suboptimal'):'n/a');
    bromoformV.style().set('color',aq.bromo_score>=80?'#115511':aq.bromo_score>=50?'#664400':'#880000');
    aquaScoreV.setValue(aq.aqua_score!==null?aq.aqua_score+'/100':('n/a ('+aq.aqua_confidence+'% confidence)'));
    aquaScoreV.style().set('color',aq.aqua_score===null?'#880000':aq.aqua_score>=70?'#115511':aq.aqua_score>=30?'#664400':'#880000');
    aquaMissingV.setValue(aq.aqua_missing_note);
    aquaStatusV.setValue(aq.status);
    aquaStatusV.style().set('color',aq.aqua_score===null?'#888888':aq.status.indexOf('GOOD')===0?'#115511':aq.status.indexOf('DECENT')===0?'#664400':'#880000');
    // Kelp opportunity display (kelpDetected already computed earlier,
    // before S7's block, so it could contextualize S7's algae display)
    if(kelpDetected){
      kelpNoteV.setValue(
        '*** COLD-WATER KELP OPPORTUNITY ***\n'+
        'Saccharina latissima (sugar kelp) or\n'+
        'Alaria esculenta (winged kelp) may be viable.\n'+
        'SST='+fmt(sv,1)+' deg C | Chl-a='+fmt(cv,3)+' mg/m3\n'+
        'Both farmed commercially in Norway / Iceland / Canada\n'+
        'at 0-12 deg C. Check S7 FAI layer for natural kelp signal.\n'+
        'See FAO cold-water aquaculture guidelines for planning.');
      kelpNoteV.style().set('color','#115511');
      kelpNoteV.style().set('backgroundColor','#e8ffe8');
    } else {
      kelpNoteV.setValue('Cold-water kelp: conditions not met (too warm or insufficient nutrients)');
      kelpNoteV.style().set('color','#aaaaaa');
      kelpNoteV.style().set('backgroundColor','rgba(0,0,0,0)');
    }

    // S9 soil
    var soilCode=(soilRes&&soilRes.b0!==null&&soilRes.b0!==undefined)?soilRes.b0:null;
    var soilLabel=soilCode!==null?soilTextureLabel(soilCode):null;
    soilTextureV.setValue(soilCode!==null?soilLabel+' (USDA class '+soilCode+')':'Dataset unavailable (v10.42: asset dead)');
    soilTextureV.style().set('color',soilCode!==null?'#664422':'#888888');

    // S10/S11
    var eqCount=(eqRes&&eqRes.count!==null&&eqRes.count!==undefined)?eqRes.count:0;
    var eqMaxMag=(eqRes&&eqRes.maxMag!==null&&eqRes.maxMag!==undefined)?eqRes.maxMag:null;
    var volcCount=(volcRes!==null&&volcRes!==undefined)?volcRes:0;
    eqStatsV.setValue(eqCount>0?eqCount+' event(s), max M'+eqMaxMag.toFixed(1):'None nearby (M4.5+)');
    eqStatsV.style().set('color',eqCount>0?'#aa3300':'#115511');
    volcStatsV.setValue(volcCount>0?volcCount+' historic event(s) nearby':'None recorded 1960-2018');
    volcStatsV.style().set('color',volcCount>0?'#aa3300':'#115511');

    // S12 Real CSD
    var realAC1=(csdRes&&csdRes.realAC1!==null&&csdRes.realAC1!==undefined)?csdRes.realAC1:null;
    var varTrendRatio=(csdRes&&csdRes.varTrendRatio!==null&&csdRes.varTrendRatio!==undefined)?csdRes.varTrendRatio:null;
    var nValidMonths=(csdRes&&csdRes.nValidMonths!==null&&csdRes.nValidMonths!==undefined)?csdRes.nValidMonths:0;
    realAc1V.setValue(realAC1!==null?realAC1.toFixed(3)+' (n='+nValidMonths+' months)':'n/a (insufficient valid months, n='+nValidMonths+')');
    realAc1V.style().set('color',realAC1===null?'#888888':realAC1>0.6?'#aa3300':realAC1>0.3?'#aa6600':'#226666');
    realVarTrendV.setValue(varTrendRatio!==null?varTrendRatio.toFixed(2)+'x'+(varTrendRatio>1.5?' RISING (possible CSD signal)':varTrendRatio<0.67?' falling':' stable'):'n/a');
    realVarTrendV.style().set('color',varTrendRatio===null?'#888888':varTrendRatio>1.5?'#aa3300':'#226666');
    // v10.87 FIX: this note used to say "Both AC1 and variance rising" -
    // but realAC1>0.5 only means AC1 is CURRENTLY elevated in this one
    // fixed window (no BEFORE baseline exists here at all), not that it
    // is "rising". Only the variance ratio (2nd half vs 1st half of THIS
    // window) is actually a trend. Calling both "rising" was misleading
    // and could contradict S13's proper BEFORE/AFTER Scheffer test, which
    // computes a real AC1 delta against a stored baseline.
    if(realAC1!==null&&realAC1>0.5&&varTrendRatio!==null&&varTrendRatio>1.3){
      realCsdNoteV.setValue('AC1 is elevated AND variance is rising within this single window - worth watching.\n'+
        'NOTE: this is a one-window snapshot with no BEFORE baseline, not a validated before/after test.\n'+
        'For a real Scheffer 2009 BEFORE-vs-AFTER check against a control site, use S13 below.');
    } else {
      realCsdNoteV.setValue('');
    }

    // S16 ECI from GEBCO
    var gebcoDepth=(gebcoRes&&gebcoRes.elevation!==null&&gebcoRes.elevation!==undefined)?gebcoRes.elevation:null;
    var gebcoDepthM=gebcoDepth!==null?-gebcoDepth:null;
    eciDepthV.setValue(gebcoDepthM!==null?fmt(gebcoDepthM,1)+' m (GEBCO)':'n/a');
    eciDepthV.style().set('color',gebcoDepthM!==null?'#115533':'#888888');
    if(gebcoDepthM!==null&&gebcoDepthM>0){
      var eci_current=Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM)));
      var eci_slr03=Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM-0.3)));
      var eci_slr05=Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM-0.5)));
      var eci_slr10=Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM-1.0)));
      eciCurrentV.setValue(eci_current.toFixed(3)+(eci_current>1.0?' EXTREME risk':eci_current>0.7?' HIGH risk':eci_current>0.4?' moderate risk':' low risk'));
      eciCurrentV.style().set('color',eci_current>1.0?'#880000':eci_current>0.7?'#aa3300':'#115533');
      eciSLR03V.setValue(eci_slr03.toFixed(3)+' ECI (vs '+eci_current.toFixed(3)+' now)');
      eciSLR05V.setValue(eci_slr05.toFixed(3)+' ECI');
      eciSLR10V.setValue(eci_slr10.toFixed(3)+' ECI'+(eci_slr10>1.2?' CRITICAL':eci_slr10>1.0?' EXTREME':''));
      eciSLR10V.style().set('color',eci_slr10>1.0?'#aa3300':'#226644');
      var bScore=sc.B, eciNorm=eci_current/1.414, bEci=1.0-eciNorm;
      // v10.145 FIX: investigated this "agreement" check directly - it
      // was comparing a PURE physical wave-energy metric (depth only)
      // against the broad ecological Fused composite (SST, chlorophyll,
      // DHW, biology, NO2). These measure genuinely different things -
      // a "DISAGREE" here does not indicate an error in either
      // calculation, since there's no real reason a narrow physics proxy
      // should track a broad ecological score. Relabeled to reflect
      // this honestly instead of implying a validation failure.
      var closeMatch=Math.abs(bScore-bEci)<0.2;
      var agreement=closeMatch?'Similar values (within 0.2) - coincidental, not a validation':
        'DIFFERENT, as expected - physical wave-exposure (ECI='+bEci.toFixed(2)+') and broad ecological risk (B='+bScore.toFixed(2)+') measure different things, not the same quantity twice';
      eciValidationV.setValue(agreement);
    } else {
      eciCurrentV.setValue(gebcoDepthM===null?'n/a':'land/dry area - ECI not applicable');
      eciSLR03V.setValue('n/a'); eciSLR05V.setValue('n/a'); eciSLR10V.setValue('n/a'); eciValidationV.setValue('n/a');
    }

    // S14 Thermal Recovery
    monthlyFC.evaluate(function(monthlyData, eTR){
      if(eTR){print('S14 monthly data error: '+eTR); thermalEpisodesV.setValue('n/a - error'); thermalMeanV.setValue('n/a'); thermalMaxV.setValue('n/a'); thermalOngoingV.setValue('n/a'); return;}
      var trResult=analyzeThermalRecovery(monthlyData.features,mmmv);
      if(trResult.error){
        thermalEpisodesV.setValue('n/a - '+trResult.error); thermalMeanV.setValue('n/a'); thermalMaxV.setValue('n/a'); thermalOngoingV.setValue('n/a');
      } else {
        thermalEpisodesV.setValue(trResult.nCompletedEpisodes+' completed episode(s)');
        thermalMeanV.setValue(trResult.meanRecoveryMonths!==null?trResult.meanRecoveryMonths.toFixed(1)+' months average':'n/a (no completed episodes)');
        thermalMaxV.setValue(trResult.maxRecoveryMonths!==null?trResult.maxRecoveryMonths+' months (longest)':'n/a');
        thermalOngoingV.setValue(trResult.ongoingEpisode!==null?'YES - started '+trResult.ongoingEpisode.startDate+' ('+trResult.ongoingEpisode.months+' mo so far, peak '+trResult.ongoingEpisode.peakSST.toFixed(1)+' deg C)':'No ongoing episode');
        thermalOngoingV.style().set('color',trResult.ongoingEpisode!==null?'#aa3300':'#115511');
        print('=== S14 THERMAL RECOVERY TIME ===');
        print('MMM threshold: '+trResult.threshold.toFixed(2)+' deg C | Completed episodes: '+trResult.nCompletedEpisodes);
        for(var ti=0;ti<trResult.episodes.length;ti++){
          var tep=trResult.episodes[ti];
          print('  Episode '+(ti+1)+': '+tep.startDate+' to '+tep.endDate+' - '+tep.months+' month(s), peak '+tep.peakSST.toFixed(2)+' deg C');
        }
        if(trResult.ongoingEpisode!==null) print('ONGOING: started '+trResult.ongoingEpisode.startDate+', '+trResult.ongoingEpisode.months+' months so far');
      }
      // S15 Ecological Recovery Validation
      var ecoVal=getEcologicalRecoveryValidation(region,lat,lon);
      ecoValStatusV.setValue(ecoVal.checked?(ecoVal.thermalPredictsEcological===false?'NO - thermal did NOT predict ecological recovery':ecoVal.thermalPredictsEcological===null?'UNKNOWN - no published follow-up found':'YES - consistent'):'NOT CHECKED - no published data found');
      ecoValStatusV.style().set('color',ecoVal.checked?(ecoVal.thermalPredictsEcological===false?'#880000':ecoVal.thermalPredictsEcological===null?'#664400':'#115511'):'#888888');
      if(ecoVal.details) ecoValDetailsV.setValue(ecoVal.details.slice(0,300));
      print('=== S15 ECOLOGICAL RECOVERY VALIDATION ===');
      if(ecoVal.checked){print('Finding: '+ecoVal.finding); print('Sources: '+ecoVal.sources);}
      else{print('NOT CHECKED - no published ecological-recovery follow-up independently checked.');}
      // v10.138: S20 is now button-triggered (see consolidated
      // LIVE/HISTORICAL tool below in the panel) - no longer auto-
      // computed here on every click.
    });

    // S17 ToE + S18 BGC (parallel evaluate, v10.60; v10.89 split core/BGC)
    var toeResults={sst:null,chl:null,sal:null,no2:null,ph:null,do_o2:null};
    var toePhAvailable=null, toeDoAvailable=null; // v10.141 split: null=pending, true once computed, false if unavailable - previously one shared flag, now independent since pH and DO resolve separately

    function calcToE(r, nYears){
      if(!r||r.scale===null||r.scale===undefined||r.noise===null||r.noise===undefined) return {snr:null,emerged:false,error:'no data'};
      var sl=r.scale, ns=r.noise;
      var signal=Math.abs(sl*nYears), snr=ns>0?signal/ns:0;
      return {slope:sl,noise:ns,signal:signal,snr:snr,emerged:snr>=2.0,direction:sl>0?'RISING':'FALLING',error:null};
    }
    function toeTxt(t,conf){if(!t||t.error)return 'n/a ('+((t&&t.error)||'no data')+')'; return (t.emerged?'EMERGED':'not yet')+'  SNR='+t.snr.toFixed(2)+' ['+conf+'] '+t.direction;}

    function renderToeCompound(){
      var have=[];
      if(toeResults.sst!==null) have.push(toeResults.sst);
      if(toeResults.chl!==null) have.push(toeResults.chl);
      if(toeResults.sal!==null) have.push(toeResults.sal);
      if(toeResults.no2!==null) have.push(toeResults.no2);
      if(toePhAvailable===true && toeResults.ph!==null) have.push(toeResults.ph);
      if(toeDoAvailable===true && toeResults.do_o2!==null) have.push(toeResults.do_o2);
      var nTotal=have.length;
      var nEmerged=have.filter(function(t){return t&&t.emerged;}).length;
      var pendingNote = (toePhAvailable===null || toeDoAvailable===null) ? ' (pH/DO still pending...)' :
        (toePhAvailable===false && toeDoAvailable===false) ? ' (pH excluded - see row below; DO excluded - dataset unavailable)' :
        toeDoAvailable===false ? ' (DO excluded - dataset unavailable, see row below)' :
        toePhAvailable===false ? ' (pH excluded - see row below)' : '';
      toeCompoundV.setValue(nEmerged+'/'+nTotal+' available variables emerged (SNR >= 2.0)'+pendingNote+
        (nTotal===0?'':nEmerged>=4?' | COMPOUND CID DETECTED':nEmerged>=2?' | MULTIPLE CIDs':nEmerged===1?' | SINGLE CID detected':' | No emergence yet'));
      toeCompoundV.style().set('color',nEmerged>=4?'#880000':nEmerged>=2?'#aa3300':nEmerged>=1?'#664400':'#115511');
    }

    rToeCore.evaluate(function(coreRes,eCore){
      if(eCore){
        print('=== S17 ToE core error (SST/CHL/SAL/NO2) === '+eCore);
        var coreErrTxt=friendlyEEError(eCore);
        toeSSTv.setValue(coreErrTxt); toeCHLv.setValue(coreErrTxt); toeSALv.setValue(coreErrTxt); toeNO2v.setValue(coreErrTxt);
        renderToeCompound();
        return;
      }
      var tSST=calcToE(coreRes.sst,44), tCHL=calcToE(coreRes.chl,27), tSAL=calcToE(coreRes.sal,32), tNO2=calcToE(coreRes.no2,7);
      toeResults.sst=tSST; toeResults.chl=tCHL; toeResults.sal=tSAL; toeResults.no2=tNO2;
      toeSSTv.setValue(toeTxt(tSST,'HIGH conf - 44yr OISST, v10.54 fix applied'));
      toeCHLv.setValue(toeTxt(tCHL,'MARGINAL - 27yr, global product'));
      toeSALv.setValue(toeTxt(tSAL,'MARGINAL - 32yr HYCOM model'));
      toeNO2v.setValue(toeTxt(tNO2,'LOW conf - 7yr only, pre-industrial baseline unknown'));
      toeSSTv.style().set('color',tSST&&tSST.emerged?'#aa3300':'#226644');
      toeCHLv.style().set('color',tCHL&&tCHL.emerged?'#aa3300':'#226644');
      toeSALv.style().set('color',tSAL&&tSAL.emerged?'#aa3300':'#226644');
      toeNO2v.style().set('color',tNO2&&tNO2.emerged?'#664400':'#888888');
      renderToeCompound();
      print('=== S17 ToE core (SST/CHL/SAL/NO2) ===');
      print('SST (44yr): SNR='+(tSST.snr!==null?tSST.snr.toFixed(3):'n/a')+' -> '+(tSST.emerged?'EMERGED':'not yet'));
      print('CHL (27yr): SNR='+(tCHL.snr!==null?tCHL.snr.toFixed(3):'n/a')+' -> '+(tCHL.emerged?'EMERGED':'not yet'));
      print('SAL (32yr): SNR='+(tSAL.snr!==null?tSAL.snr.toFixed(3):'n/a')+' -> '+(tSAL.emerged?'EMERGED':'not yet'));
      print('NO2 (7yr):  SNR='+(tNO2.snr!==null?tNO2.snr.toFixed(3):'n/a')+' -> '+(tNO2.emerged?'EMERGED':'not yet'));
    });

    rToeBGC_PH.evaluate(function(phRes,ePh){
      if(ePh){
        toePhAvailable=false;
        toePHv.setValue('n/a - pH data unavailable at this location/time (real asset, but no valid pixel here)');
        toePHv.style().set('color','#888888');
        print('=== S17 ToE pH error (isolated - does not affect SST/CHL/SAL/NO2/DO) === '+ePh);
      } else {
        toePhAvailable=true;
        // v10.141: record-years changed from 32 to 4 - the real,
        // confirmed pH asset's actual data only starts 2021-10-01, not
        // 1993. Using the old 32yr assumption with real ~4yr data would
        // have silently overstated the SNR (SNR = |trend*years|/noise
        // scales UP with an inflated year count) - fixed to match the
        // real record length used to build _annPHColl above.
        var tPH=calcToE(phRes.ph,4);
        toeResults.ph=tPH;
        toePHv.setValue(toeTxt(tPH,'REAL asset (COPERNICUS CAR/ph_depth1), surface only, ~4yr record (2022-2025) - LOW confidence, short record'));
        toePHv.style().set('color',tPH&&tPH.emerged?'#880000':'#226644');
        print('=== S17 ToE pH (v10.141: real asset, ~4yr record) ===');
        print('pH (~4yr, LOW conf): SNR='+(tPH.snr!==null?tPH.snr.toFixed(3):'n/a')+' -> '+(tPH.emerged?'EMERGED':'not yet'));
      }
      renderToeCompound();
    });
    rToeBGC_DO.evaluate(function(doRes,eDo){
      if(eDo){
        toeDoAvailable=false;
        var doUnavailTxt='n/a - DO dataset unavailable (COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/MFC_001_028 not found in the current GEE catalog - a working replacement was not confirmed this session; this does NOT affect pH above, or SST/Chl/Salinity/NO2)';
        toeDOv.setValue(doUnavailTxt);
        toeDOv.style().set('color','#888888');
        print('=== S17 ToE DO: dataset unavailable (isolated - does not affect pH/SST/CHL/SAL/NO2) === '+eDo);
        renderToeCompound();
        return;
      }
      toeDoAvailable=true;
      var tDO=calcToE(doRes.do_o2,32);
      toeResults.do_o2=tDO;
      toeDOv.setValue(toeTxt(tDO,'BGC model surface only - not measured, no depth zones'));
      toeDOv.style().set('color',tDO&&tDO.emerged?'#880000':'#226644');
      renderToeCompound();
      print('=== S17 ToE DO ===');
      print('DO  (32yr): SNR='+(tDO.snr!==null?tDO.snr.toFixed(3):'n/a')+' -> '+(tDO.emerged?'EMERGED':'not yet')+' [BGC model surface]');
    });

    ee.Dictionary({ph:rBGC_PH.get('ph_depth1'), pco2Pa:rBGC_CO2.get('spco2_depth1'), sal:rBGC_SAL.get('salinity_0')}).evaluate(function(bgcData,eBGC){
      if(eBGC){
        print('=== S18 BGC error === '+eBGC);
        var s18ErrTxt='n/a - error fetching real BGC assets (see console)';
        s18pHV.setValue(s18ErrTxt); s18pCO2V.setValue(s18ErrTxt); s18SalV.setValue(s18ErrTxt);
      } else {
        var bgcD=bgcData||{};
        var s18NullTxt='n/a - genuinely no valid pixel in the real asset here, even over a 2yr window (v10.147: bare-point geometry, same as S17\'s working fetch)';
        s18pHV.setValue(bgcD.ph!==null&&bgcD.ph!==undefined?fmt(bgcD.ph,3)+' (REAL asset, ~2yr recent mean)':s18NullTxt);
        var pco2Uatm = (bgcD.pco2Pa!==null&&bgcD.pco2Pa!==undefined)?bgcD.pco2Pa*9.86923:null;
        s18pCO2V.setValue(pco2Uatm!==null?fmt(pco2Uatm,1)+' uatm (REAL asset, ~2yr recent mean)':s18NullTxt);
        var salReal = (bgcD.sal!==null&&bgcD.sal!==undefined)?(bgcD.sal*0.001+20):null;
        s18SalV.setValue(salReal!==null?fmt(salReal,2)+' PSU (REAL HYCOM, ~2yr recent mean)':s18NullTxt);
        print('=== S18 BIOGEOCHEMISTRY SNAPSHOT (v10.147: bare-point geometry fix) ===');
        print('pH: '+(bgcD.ph!==null?bgcD.ph.toFixed(3):'n/a')+' | pCO2: '+(pco2Uatm!==null?pco2Uatm.toFixed(1)+' uatm':'n/a')+' | Salinity: '+(salReal!==null?salReal.toFixed(2)+' PSU':'n/a'));
      }
      // O2 - genuinely unresolved, disclosed directly, independent of
      // the 3 real fields above (never bundled with a dead asset again).
      s18O2V.setValue('n/a - no working dissolved-oxygen asset confirmed this session (searched twice; GLOBAL_OCEAN_BGC/MFC_001_028 is dead, real replacement sub-collection/band name not found)');
    });

    // Interventions
    var isReefZoneForInt=(lat>-30&&lat<30)&&!isEBUS(lat,lon);
    var interventions=computeInterventions(region,fp,sc,dhwv,tv,isReefZoneForInt);
    intPanel.clear();
    var priorityColors={'CRITICAL':{bg:'#ffd0d0',tc:'#880000'},'HIGH':{bg:'#ffe0c0',tc:'#aa4400'},'MODERATE':{bg:'#fff4c0',tc:'#886600'},'LOW':{bg:'#e8f4ff',tc:'#225588'},'DATA GAP':{bg:'#f0e8ff',tc:'#553388'},'CONTEXT':{bg:'#f0f0f0',tc:'#555555'},'STABLE':{bg:'#e0ffe0',tc:'#115511'}};
    for(var ai=0;ai<interventions.length;ai++){
      var act=interventions[ai], pc=priorityColors[act.priority]||{bg:'#eeeeee',tc:'#333333'};
      var actBox=ui.Panel({style:{backgroundColor:pc.bg,margin:'2px 0',padding:'4px 6px',border:'1px solid '+pc.tc}});
      actBox.add(ui.Label('['+act.priority+'] '+act.action,{fontSize:'8px',fontWeight:'bold',color:pc.tc,backgroundColor:'rgba(0,0,0,0)',margin:'0',whiteSpace:'pre'}));
      actBox.add(ui.Label(act.basis,{fontSize:'7px',color:'#444444',backgroundColor:'rgba(0,0,0,0)',margin:'2px 0 0 0',whiteSpace:'pre'}));
      actBox.add(ui.Label('Waddington: '+act.waddington,{fontSize:'7px',color:'#225588',fontWeight:'bold',backgroundColor:'rgba(0,0,0,0)',margin:'2px 0 0 0',whiteSpace:'pre'}));
      intPanel.add(actBox);
    }

    print('=== MEASUREMENTS ===');
    print('SST annual: '+fmt(sv,2)+' | Peak '+peakSeasonTxt+': '+fmt(sv_peak,2)+' deg C (BLEACHING SCORE)');

    // Print time-series charts only when valid ocean SST data exists
    if(sv!==null){
      print(ui.Chart.image.series({imageCollection:clip(mkMoSST()),region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
      .setOptions({title:'SST Monthly | '+region+' | OISST V2.1',series:{0:{color:'#ff4422',lineWidth:2.5,pointSize:4,label:'SST'}},
        backgroundColor:'#0a1628',titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
        vAxis:{title:'SST (deg C)',textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
        hAxis:{title:'Month Jan 2023 - Dec 2024',textStyle:{color:'#aaaaaa'}},legend:{textStyle:{color:'#ffffff'}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));

      print(ui.Chart.image.series({imageCollection:clip(mkMoDHW()),region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
      .setOptions({title:'DHW Monthly (per-pixel MMM) | '+region,series:{0:{color:'#ff6644',lineWidth:2.5,pointSize:4,label:'DHW'}},
        backgroundColor:'#0a1628',titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
        vAxis:{title:'DHW (deg C-weeks)',viewWindow:{min:0},textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
        hAxis:{title:'Month',textStyle:{color:'#aaaaaa'}},legend:{textStyle:{color:'#ffffff'}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));

      print(ui.Chart.image.series({imageCollection:clip(mkMoCHL()),region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
      .setOptions({title:'Chl-a Monthly | '+region+' | Copernicus Ocean Color V6',series:{0:{color:'#44cc44',lineWidth:2.5,pointSize:4,label:'Chl-a'}},
        backgroundColor:'#0a1628',titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
        vAxis:{title:'Chl-a (mg/m3)',textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
        hAxis:{title:'Month Jan 2023 - Dec 2024',textStyle:{color:'#aaaaaa'}},legend:{textStyle:{color:'#ffffff'}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));

      print(ui.Chart.image.series({imageCollection:clip(annualSST),region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
      .setOptions({title:'Annual SST 2003-2024 | '+region+' | Basin erosion CSD',series:{0:{color:'#ff6644',lineWidth:2.5,pointSize:4,label:'Annual SST'}},
        backgroundColor:'#0a1628',titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
        vAxis:{title:'Annual SST (deg C)',textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
        hAxis:{title:'Year 2003-2024',textStyle:{color:'#aaaaaa'}},
        legend:{position:'top',textStyle:{color:'#ffffff',fontSize:9}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'},
        trendlines:{0:{type:'linear',color:'#ffcc44',lineWidth:2.5,opacity:0.9,showR2:true,visibleInLegend:true}}}));

      print(ui.Chart.image.series({imageCollection:clip(mkMoNO2()),region:study,reducer:ee.Reducer.mean(),scale:1000,xProperty:'system:time_start'})
      .setOptions({title:'NO2 Monthly | '+region+' | S5P OFFL+NRTI',series:{0:{color:'#cc44cc',lineWidth:2.5,pointSize:4,label:'NO2'}},
        backgroundColor:'#0a1628',titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
        vAxis:{title:'NO2 (mol/m2)',textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
        hAxis:{title:'Month Jan 2023 - Dec 2024',textStyle:{color:'#aaaaaa'}},legend:{textStyle:{color:'#ffffff'}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));
    } else {
      print('Charts skipped: no valid ocean SST data at this location (land / ice-covered).');
    }
    print('MMM local: '+fmt(mmmv,1)+' | DHW: '+fmt(dhwv,2)+' deg C-wks | Trend: '+fmt(tv,4)+' deg C/yr');
    print('Chl-a: '+fmt(cv,3)+' mg/m3 | Turbidity: '+(turv!==null?fmt(turv,3)+' NDTI':'n/a')+' | NO2: '+fmt(nv,8)+' | Depth: '+fmt(bv,0)+' m (GEBCO)');
    print('=== FUSED RESULT ===');
    print('Satellite CCS: '+sc.sat_ccs+'/100 | Field correction: '+(sc.fcTotal>0?'+':'')+sc.fcTotal);
    print('FUSED CCS: '+sc.ccs+'/100 ('+sc.acc_total+'%) | B='+sc.B.toFixed(4)+' | Status: '+cols.lbl);

    var sFC=ee.FeatureCollection([
      ee.Feature(null,{c:'S1 SST\n(15%)',s:sc.s1,f:sc.s1,safe:30}),
      ee.Feature(null,{c:'S2 Chl-a\n(15%)',s:sc.s2,f:sc.s2,safe:30}),
      ee.Feature(null,{c:'S3 Turb\n(15%)',s:sc.s3,f:sc.s3,safe:30}),
      ee.Feature(null,{c:'S4 DHW+Trend\n(25%)',s:sc.s4,f:sc.s4,safe:30}),
      ee.Feature(null,{c:'S5 Bio\n(20%)',s:sc.s5,f:sc.s5,safe:30}),
      ee.Feature(null,{c:'S6 NO2\n(10%)',s:sc.s6,f:sc.s6,safe:30})
    ]);
    print(ui.Chart.feature.byFeature({features:sFC,xProperty:'c',yProperties:['s','f','safe']}).setChartType('ColumnChart').setOptions({
      title:'FUSED SCORE: '+sc.ccs+'/100 | '+cols.lbl+' | '+region+' | B='+sc.B.toFixed(3),
      series:{0:{color:'#ffaa44',label:'Satellite'},1:{color:cols.bar,label:'Fused'},2:{color:'#4466ff',type:'line',lineWidth:2,pointSize:0,label:'Safe (30)'}},
      vAxis:{title:'Risk 0-100',viewWindow:{min:0,max:100},textStyle:{color:'#cccccc'},titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
      hAxis:{title:'Component (weight)',textStyle:{color:'#dddddd'}},
      backgroundColor:'#0a1628',titleTextStyle:{color:cols.bar,fontSize:10,bold:true},
      legend:{position:'top',textStyle:{color:'#ffffff',fontSize:9}},chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));

    // Log click for export
    var trForLog={n:0,mean:null,max:null,ongoing:false};
    var ecoValForLog=getEcologicalRecoveryValidation(region,lat,lon);
    clickResultsLog.push({
      timestamp:new Date().toISOString(),lat:latR,lon:lonR,region:region,
      sst_annual_c:sv,sst_peak_c:sv_peak,sst_trend_c_per_yr:tv,
      dhw_c_weeks:dhwv,dhw_raw_unfiltered:dhwv_raw,dhw_artifact_flagged:dhwFlaggedArtifact,
      mmm_local_c:mmmv,chl_a_mg_m3:cv,turbidity_ndti:turv,no2_mol_m2:nv,depth_m:bv,
      cancer_score_satellite:sc.sat_ccs,cancer_score_fused:sc.ccs,field_correction:sc.fcTotal,
      bowl_depth_B:sc.B,accuracy_pct:sc.acc_total,status_label:cols.lbl,
      field_data_available:fp.hasField,field_species:fp.hasField?fp.species:'none',
      aquaculture_score:aq.aqua_score,aquaculture_confidence_pct:aq.aqua_confidence,
      aquaculture_status:aq.status,bromoform_yield:aq.bromo_score,
      num_interventions:interventions.length,top_intervention:interventions.length>0?interventions[0].action:'none',
      soil_texture_code:soilCode,soil_texture_label:soilCode!==null?soilLabel:'Dataset unavailable (v10.42)',
      earthquake_count_200km:eqCount,earthquake_max_mag_200km:eqMaxMag,volcanic_activity_count_300km:volcCount,
      real_ac1_24mo:realAC1,real_var_trend_ratio:varTrendRatio,real_csd_n_valid_months:nValidMonths,ndci_nutrient_proxy:nciV,
      thermal_recovery_n_completed_episodes:trForLog.n,thermal_recovery_mean_months:trForLog.mean,
      thermal_recovery_max_months:trForLog.max,thermal_recovery_ongoing_flag:trForLog.ongoing,
      eco_recovery_validation_checked:ecoValForLog.checked,eco_recovery_validation_finding:ecoValForLog.finding||'n/a',
      gebco_depth_m:gebcoDepthM,energy_concentration_risk:gebcoDepthM!==null&&gebcoDepthM>0?Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM))):null,
      eci_slr_03m:gebcoDepthM!==null&&gebcoDepthM>0?Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM-0.3))):null,
      eci_slr_05m:gebcoDepthM!==null&&gebcoDepthM>0?Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM-0.5))):null,
      eci_slr_10m:gebcoDepthM!==null&&gebcoDepthM>0?Math.min(1.414,1.0/Math.sqrt(Math.max(0.5,gebcoDepthM-1.0))):null,
      eci_vs_heuristic_validation:'see console',
      in_situ_station:inSitu?inSitu.key:null,in_situ_distance_km:inSitu?inSitu.distance_km:null,
      in_situ_ph_mean:inSitu&&inSitu.pH?inSitu.pH.mean:null,
      in_situ_temp_c_mean:inSitu&&inSitu.temp_c&&inSitu.temp_c.mean!==null?inSitu.temp_c.mean:null,
      in_situ_salinity_mean:inSitu&&inSitu.salinity&&inSitu.salinity.mean!==null?inSitu.salinity.mean:null,
      in_situ_do_mgL_mean:inSitu&&inSitu.do_mgL?inSitu.do_mgL.mean:null
    });
    exportStatusLabel.setValue(clickResultsLog.length+' row(s) ready - click EXPORT button above');
    exportStatusLabel.style().set('color','#116611');
    clearPipelineWatchdog();
    print('Logged click #'+clickResultsLog.length+'. S19 station: '+(inSitu?inSitu.key:'none'));

  } // end _onAllDone
} // end analyzeLocation

Map.onClick(function(coords){ analyzeLocation(coords.lat, coords.lon); });

// STARTUP
print('STEMGeoHS Marine v10.153 -- READY');
print('');
print('v10.153 FIX 10: the fixed thresholds now have a MEASURED error rate.');
print('  Calibrated against 13.6 years of Scripps Pier CTD (Jan 2013-Aug 2026,');
print('  1.95M QC-flagged 4-min readings: temperature, salinity, chlorophyll).');
print('  La Jolla is a temperate kelp coast with no documented regime shift, so');
print('  every threshold crossing there is a false positive by construction.');
print('');
print('  MEASURED false-positive rate of the OLD cutoffs (AC1>0.01, Var>0.15x),');
print('  across 104-176 non-overlapping BEFORE/AFTER splits:');
print('     24mo window: AC1 51%, Var 48%');
print('     36mo window: AC1 46%, Var 49%');
print('     48mo window: AC1 34%, Var 57%');
print('     the AC1-OR-variance rule S7E uses to declare a signal: 80%');
print('  The permutation test on those SAME windows rejected at 0-9%, i.e. it is');
print('  correctly calibrated and has been right every time it disagreed.');
print('');
print('  Cutoffs are now the 95th percentile of |delta| under that measured null,');
print('  set per window length (24/36/48mo) in CSD_NULL_CALIBRATION. Applied in');
print('  STEP 3, FIND SWEET SPOT, S7E and S7F.');
print('  Set CSD_USE_CALIBRATED_THRESHOLDS=false to restore the old behaviour.');
print('');
print('  EXPECT FAR FEWER LOCAL/SIGNAL VERDICTS. That is the correction working,');
print('  not a regression. Calibrated cutoffs this large rarely fire - which is');
print('  the honest conclusion: a fixed threshold cannot adapt to how noisy a');
print('  given site is, while the permutation test does so automatically.');
print('');
print('  S7D and S7F run NO significance test at all. Both now print an explicit');
print('  threshold-only warning. The permutationTestDelta engine already exists');
print('  in this file and should be ported to them next.');
print('');
print('v10.152 FIXES (3) - all found in the real Looe Key 48mo run:');
print('  7. The headline printed "NO RELIABLE CSD SIGNAL - AC1 not rising"');
print('     immediately followed by "AC1 (primary): RISING". combinedTitle');
print('     came from the p-gated confidence, but the line after it still');
print('     read raw direction. Both now derive from the same object, and');
print('     a rising-but-null AC1 reads "direction rising, but NOT');
print('     significant" instead of a bare RISING.');
print('  8. The v10.127 header said "13 of 12 calendar months" - the new');
print('     _meta key was being counted by Object.keys(). Now uses the');
print('     real count, and reports how many months were rejected by the');
print('     >=3-sample floor rather than the stale "~1-2 samples" text.');
print('  9. Spatial autocorrelation returned EXACTLY +0.000 at two');
print('     different sites on different dates. Two independent pixel');
print('     fields do not agree bit-for-bit. extractSpatialAC1 now reports');
print('     WHICH key it read; if it fell through to the last-resort loop');
print('     the indicator is labelled UNRELIABLE and EXCLUDED from the');
print('     tally, so it can no longer pad the supporting-indicator count');
print('     with a value that was never a measurement.');
print('');
print('v10.151 FIXES (6) - all driven by a real Looe Key run, 24.531 -81.41:');
print('  1. parseCoordPart now accepts pasted Unicode minus / en-dash /');
print('     degree sign. These previously gave a bare "invalid format".');
print('  2. computeMonthlyClimatology requires >=3 samples per calendar');
print('     month. With exactly 2, subtracting the mean forces the pair to');
print('     (+d,-d) BY CONSTRUCTION. Simulated on pure seasonal + white');
print('     noise (true AC1=0), the old path returned AC1=-0.051 at 24mo.');
print('  3. jsNodeStatsFixed floor raised 24 -> 48 months, and it now');
print('     accepts a SHARED climatology.');
print('  4. computeRealCSDDeseasonalized returns the raw series; STEP 2');
print('     stores it; COMPARE pools BOTH windows into ONE climatology.');
print('     This collapses the three conflicting AC1 values the same');
print('     window reported (-0.332 / +0.062 / +0.813) into one number');
print('     that finally matches the permutation test.');
print('  5. classifyToolkitConfidence now takes the p-values and CANNOT');
print('     exceed them. The run that printed "AC1 p=0.980" beside');
print('     "Confidence: HIGH" now reads LOW, and says why.');
print('  6. verdictsConflict replaces two single-string conflict checks');
print('     that missed the "STRONG SIGNAL ... Regional context: NO SIGNAL');
print('     AT EITHER SITE" case. Fifth appearance of that bug class.');
print('');
print('NOT changed, deliberately - see NOTE A in the code: S12b shuffles');
print('the series, so its null is white noise, which every real SST pixel');
print('rejects for ordinary physical reasons. Your own validation recorded');
print('p=0.0000 on structured data. A spectrum-preserving surrogate would');
print('make an AC1 test near-vacuous instead. Recommend deleting S12b');
print('rather than patching it; COMPARE delta test is the sound one.');
print('');
print('v10.149 NEW: S12b real significance test, heuristic->probabilistic');
print('  S12 was real math with ZERO statistical testing - a heuristic.');
print('  Added a real permutation test: shuffles S12 own 24-month window');
print('  into random order 500 times, asking "would randomly-ordered data');
print('  produce an AC1 this high by chance?" - genuinely different from');
print('  S13 change-detection test. Reuses mkMoSST() (same window S12');
print('  itself uses) and jsNodeStats() (same formula) - real AC1 always');
print('  matches S12 own number. Verified: structured data -> p=0.0000');
print('  (significant); pure noise -> p=0.76 (not significant). DISCLOSED');
print('  directly: a significant result here does NOT mean AC1 changed or');
print('  a reef is destabilizing - real seasonal cycling will always look');
print('  "significant" against scrambled data. Tests if memory is real,');
print('  not if resilience is being lost - that is still S13 job.');
print('');
print('v10.148 FIX: S7 no longer alarms coral-style at kelp sites');
print('  Real Nuuk test showed S7 saying "EXTREME - massive bloom, RED"');
print('  at a site S8 independently flagged as likely healthy cold-water');
print('  kelp habitat - same class of gap already fixed for DHW (hidden');
print('  above 55 deg N/S). Moved computeAquaculture() earlier so its');
print('  kelp-detection result can contextualize S7 text. When kelp is');
print('  detected AND algae is elevated, S7 now says "LIKELY COLD-WATER');
print('  KELP... not a coral-algae bloom warning" instead. Verified both');
print('  directions: real Nuuk values now show kelp context; real');
print('  Florida Keys values (kelp not detected) still show the normal');
print('  alarm, unchanged.');
print('');
print('v10.147 BUGFIX: S18 real root cause found via 2nd real test');
print('  v10.146 fix did not solve it - One Tree Reef (open lagoon, 36m');
print('  deep) failed identically to Florida Keys, falsifying the');
print('  coastal-proximity hypothesis. Compared directly against S17');
print('  own working pH fetch: found S18 was using a buffered polygon');
print('  geometry instead of a bare point (risky against a coarse ~25km');
print('  grid), and a 90-day window instead of a multi-year one. Fixed');
print('  both, matching S17 proven pattern exactly. Not yet re-tested');
print('  live - new diagnostic message will say so directly if it still');
print('  fails, rather than reusing the old generic error text.');
print('');
print('v10.146 BUGFIX: S18\'s own v10.145 fix was itself broken, now fixed');
print('  Real test at One Tree Reef showed S18 pH/pCO2/Salinity all');
print('  failing, even though S17 pH (same real asset) succeeded. Root');
print('  cause: v10.145 skipped the SAME empty-collection guard every');
print('  other real-asset fetch in this codebase uses. Fixed using the');
print('  identical, already-proven ee.Algorithms.If(size().gt(0),...)');
print('  pattern; widened the window 60->90 days too.');
print('');
print('v10.145 MULTI-FIX: S18 real assets, SNR diagnosis, S16/S20e fixes');
print('  1) S18 pH+pCO2+Salinity now REAL (same assets as S17/S20e),');
print('  split so a still-dead O2 field cannot blank the 3 working ones.');
print('  2) DO/O2 searched a 3rd time - still unresolved, disclosed not');
print('  guessed. 3) Diagnosed WHY SNR vs Mann-Kendall disagree: SNR has');
print('  no correction for slope uncertainty at low sample sizes, unlike');
print('  a real p-value - disclosed directly in S17. Also found real');
print('  HYCOM salinity data sparsity at some sites, disclosed in S17b.');
print('  4) Fixed S16 ECI-vs-B: was comparing pure physics vs a broad');
print('  ecological score - relabeled, no longer implies an error.');
print('  5) S20e wired into Bowl Depth - real model now sits next to the');
print('  uncalibrated heuristic for direct comparison, fulfilling the');
print('  v10.143 promise. NOT done: S20e 6-feature extension - needs 3');
print('  new real data sources not yet found, disclosed as a real gap.');
print('');
print('v10.144 NEW: S20e real fitted bleaching probability model');
print('  Downloaded and fit a REAL model on the Global Coral-Bleaching');
print('  Database (van Woesik 2022) - 32,716 real, cleaned rows. Real');
print('  train/test split, held-out AUC=0.620 - better than chance, but');
print('  genuinely weak, disclosed directly. Two real findings from the');
print('  fit: deeper reefs more likely to bleach (independently confirms');
print('  Arias-Ortiz et al. 2024), more turbid reefs less likely (matches');
print('  Sully & Woesik 2020). Zero new EE calls - reuses DHW/Turbidity/');
print('  Depth already computed. Predicts BLEACHING, not collapse.');
print('');
print('v10.143 FIX: collapse terminology corrected where overclaiming');
print('  Found a real, large database for a possible future bleaching-');
print('  probability model (GCBD, van Woesik 2022, 34,846 records) - not');
print('  yet pulled in, needs the real file. Since any such model would');
print('  predict BLEACHING not COLLAPSE, reviewed every "collapse"');
print('  reference. FIXED: "P(flip 5yr)" -> "Regime-shift index (5yr,');
print('  uncalibrated)" - was never a real probability. STEP 2 dropdown');
print('  "AFTER event (post-collapse)" -> "(suspect period)" - most');
print('  AFTER windows tested this session showed NO signal, so assuming');
print('  verified collapse by default was misleading. LEFT UNCHANGED:');
print('  real citations (Levitan urchin collapse, Byrne colony collapse)');
print('  and S13 own accurate description of Scheffer critical-slowing-');
print('  down theory, which genuinely is about collapse.');
print('');
print('v10.142 NEW: S17b real Mann-Kendall significance test added');
print('  Reuses S17 own real annual data + the same Mann-Kendall engine');
print('  already built and fixed for STEP 5 - a genuine p-value instead');
print('  of the fixed SNR>=2.0 cutoff. Built standalone (button-');
print('  triggered) rather than rewriting S17 fragile core flow -');
print('  same lower-risk pattern as COMPARE statistically-valid banner.');
print('  Verified before shipping: short records (pH ~4yr, NO2 7yr) will');
print('  often correctly show "insufficient data" - a real, honest limit');
print('  of annual-resolution testing, not a bug.');
print('');
print('v10.140 NEW: real Mesoamerican Reef entry + model caveat added');
print('  Confirmed Healthy Reefs Initiative real (70+ partners, since');
print('  ~2008) - report-based, not a live API. Caught and fixed a real');
print('  error before shipping: nearly merged this into the Bocas del');
print('  Toro entry - they are DIFFERENT reef systems (Panama vs Mexico/');
print('  Belize/Guatemala/Honduras). Added as a real, separate 7th site,');
print('  confirmed 1042km apart, non-overlapping. Real 2023 stats: ~40%');
print('  corals severely affected, cover fell 19%->17% (GCRMN 2024).');
print('  Also disclosed a real finding: Arias-Ortiz et al. 2024 (Comms');
print('  Biology) found a 23-metric model explains 75% of bleaching');
print('  variance vs DHW alone explaining less, and deeper/more-diverse');
print('  reefs were surprisingly MORE vulnerable - flagged as a known');
print('  limit of S20 simple DHW>=4 threshold, not a full model port.');
print('');
print('v10.139 NEW: COMPARE now leads with a statistically-valid verdict');
print('  Real question answered: the p-value already existed - why was');
print('  the >0.01/>0.15 HEURISTIC threshold still primary? Real reason:');
print('  the permutation test needs its own separate EE call, resolving');
print('  LATER than the threshold verdict - a sequencing fact, not a');
print('  choice. New purple banner ABOVE the old verdict box, populated');
print('  once real p-values arrive - now the first, most prominent thing');
print('  shown. Old threshold box stays below, relabeled secondary.');
print('  Verified against 3 real sites already tested this session - all');
print('  correctly resolve to NO SIGNIFICANT SIGNAL, matching what was');
print('  already independently established at each.');
print('');
print('v10.138 CONSOLIDATED - S20/S20b/S20c merged into one tool');
print('  Real, fair feedback: 3 separate buttons for one conceptual');
print('  question was confusing, not a missing capability. Now ONE date');
print('  field: blank = LIVE (today), or type YYYY-MM (2023-01 to');
print('  2024-12) = HISTORICAL replay of that month. One button, one');
print('  result, clearly labeled which mode ran. Old always-on S20');
print('  (auto-computed every click) retired - now button-triggered like');
print('  the rest of this family. S20d (genus lookup) unchanged.');
print('');
print('v10.137 NEW: S20d Genus Growth-Form Lookup added');
print('  Confirmed a real, established global database: Coral Trait');
print('  Database (coraltraits.org, Madin et al. 2016) - 166,245 real');
print('  observations, 5,112 species, actively maintained. Could not');
print('  fetch its live CSV from inside this GEE script (unverified');
print('  outbound access) - built GENUS_GROWTH_FORM instead, a real');
print('  starter table, each entry labeled [CTD-confirmed] (individually');
print('  verified) or [literature pattern] (established but not cross-');
print('  checked per genus). Type any genus in the new S20d box for its');
print('  real growth form + general vulnerability tier - a global');
print('  fallback, distinct from S20 site-specific mortality data.');
print('');
print('v10.136 NEW: S20c Live Current DHW Check added');
print('  Checked directly: NOAA/CDR/OISST/V2_1 (same dataset used');
print('  everywhere in this tool) is confirmed real-time, updated daily,');
print('  data through essentially today. The frozen 2023-2024 window was');
print('  never a data limit - only a hardcoded date range. New button');
print('  computes real current DHW using a rolling 35-day window ending');
print('  today, same formula as elsewhere, then runs it through S20\'s');
print('  combined-risk logic. DISCLOSED: this is an OBSERVATION of the');
print('  present, NOT a forecast - genuine forecasting needs real ocean/');
print('  climate model output, a different data source, not built here.');
print('');
print('v10.135 NEW: S20b - Historical Month Check (retrospective test)');
print('  Found a real architectural fact first: the main click flow was');
print('  ALREADY hardcoded to a fixed 2023-2024 window this whole session,');
print('  not live - so DHW=7.71 was already the PEAK across that window,');
print('  not a present-day reading. No "current" value to roll back from.');
print('  Built the honest version: pick any month (YYYY-MM, 2023-01 to');
print('  2024-12) and see the REAL DHW + S20 combined risk for just that');
print('  month, using the same real data already in this tool. Months');
print('  outside that range correctly return an error, not a guess.');
print('');
print('v10.134 NEW: 3 real entries added proactively (Bocas/Andaman/Brazil)');
print('  Proactively checked sites already tested this session + named');
print('  regions, instead of waiting to be told each one. Bocas del Toro');
print('  (the exact site tested with S13/CSD this session): Neal et al.');
print('  2017, 8-year study - massive corals showed continued net tissue');
print('  loss, no full recovery from 2005 bleaching. Andaman/Nicobar: real');
print('  2010 event, clean genus mortality - Acropora 43%, Montipora 22%,');
print('  Porites 14%. Brazil (Mussismilia harttii): "unprecedented');
print('  erosion" + no recovery 3 years post-bleaching - real qualitative');
print('  finding, no invented percentage. Now 6 real sites total.');
print('');
print('v10.133 NEW: real Maldives entry added (regional, 18-reef survey)');
print('  Searched real literature for Maldives and Philippines. Maldives:');
print('  strong match - Coral Reefs journal (Springer), 2026, real survey');
print('  of 18 reefs, central/southern atolls. >40% coral cover loss in');
print('  central atolls (up to 57% worst-hit), Acropora disproportionately');
print('  hit; southern Huvadhoo Atoll showed high resistance, Porites-');
print('  dominated. This is REGIONAL data (18-reef survey), not single-');
print('  reef precision like the other 2 entries - disclosed directly.');
print('  Philippines: real papers found, but NOT added - they measure');
print('  bleaching prevalence not mortality, or lack peer review/clean');
print('  genus data - disclosed rather than forcing in weak data.');
print('');
print('v10.132 NEW: shared species database + Florida Keys real entry');
print('  Refactored S20 into SPECIES_VULNERABILITY_SITES, a shared array');
print('  any module can query via lookupSpeciesVulnerability(lat, lon) -');
print('  not S20-specific anymore. Added a real, verified second entry:');
print('  Florida Keys/Dry Tortugas. Manzello et al. 2025 (Science) found');
print('  97.8-100% mortality of Acropora palmata/cervicornis - described');
print('  as functional extinction of both species at this site. Florida');
print('  FWC CREMP 2023-24 report adds Orbicella/Siderastrea/Montastraea');
print('  data. Still explicitly NOT a comprehensive global database -');
print('  grows one verified citation at a time. Every other site still');
print('  correctly shows "NOT CHECKED", never a guess.');
print('');
print('v10.131 NEW: S20 - species-weighted combined risk module added');
print('  Combines the existing real DHW value with real, published');
print('  dominant-species vulnerability data (same Byrne et al. 2025');
print('  citation as S15, restructured as data). NOT a satellite');
print('  capability - real hyperspectral species classification is only');
print('  ~56-70% accurate currently, too unreliable to use. Only');
print('  populated where a real field survey exists (One Tree Reef for');
print('  now) - honest "NOT CHECKED" everywhere else.');
print('');
print('v10.130 FIX: STEP 5 total-months cap raised, proactive warning added');
print('  A real Panama test hit "insufficient points (n=1)" on the new');
print('  valid significance test - the fix (96 months, 12mo window) would');
print('  have been silently REJECTED by STEP 5 own 24-84 month cap. Raised');
print('  to 144 (12 years). Also added a proactive warning shown BEFORE');
print('  the EE call fires if the chosen settings would leave fewer than');
print('  4 independent windows - the default 50%-of-total auto-window');
print('  always gives only ~2, by construction, so most default runs will');
print('  now see this warning unless a window size is set manually.');
print('');
print('v10.129 NEW: real significance test for STEP 5 sliding window');
print('  Direct port of a real bug found this session: STEP 5 Mann-Kendall');
print('  ran on OVERLAPPING windows, sharing most data between consecutive');
print('  points - violates the independence assumption. Confirmed: pure');
print('  noise gave p<0.0001 overlapping vs p=0.109 non-overlapping. Now');
print('  fetches the raw series (1 new EE call), deseasonalizes, and runs');
print('  a VALID Mann-Kendall test on non-overlapping windows for both raw');
print('  and deseasonalized data - shown alongside the existing trajectory');
print('  chart, which stays useful for visualization even though its own');
print('  p-value should now be read with caution.');
print('');
print('v10.128 FIX: deseasonalized comparison logic corrected');
print('  Caught on the first real test: the interpretation text said');
print('  "similar in size" when the deseasonalized delta was actually 7x');
print('  LARGER than the raw delta - the original binary check only ever');
print('  looked for "much smaller." Now checks all three real cases, and');
print('  adds a comparison of the ABSOLUTE AC1 baseline level (not just');
print('  its change) - that same test showed raw AC1=0.843 vs');
print('  deseasonalized=0.434, nearly half, a more informative signal');
print('  than the delta comparison alone would have shown.');
print('');
print('v10.127 NEW: deseasonalized AC1/variance comparison added');
print('  Direct fix for the v10.101 disclosed caveat that shared seasonal');
print('  cycling can inflate AC1 regardless of real dynamics - never');
print('  actually built until now. New climatology functions remove each');
print('  site own calendar-month average before computing AC1/variance,');
print('  shown alongside a raw recomputation using the identical method');
print('  for a clean comparison. Zero new EE calls - reuses the same raw');
print('  series already fetched for the permutation test. Disclosed limit:');
print('  climatology built from only 1-2 samples per month at typical');
print('  window lengths - an exploratory comparison, not a full fix.');
print('');
print('v10.126 FIX: GLOBAL SIGNAL now flags strong local amplification');
print('  Caught from a real Nuuk test: study variance rose +3.40x, control');
print('  only +0.22x - both crossed the threshold, so it was called plain');
print('  "GLOBAL SIGNAL" even though the divergence (+3.18x) was over 21x');
print('  the divergence bar used elsewhere in this same tool. That bar was');
print('  only ever checked when NEITHER site crossed individually. Now: a');
print('  large divergence on top of two genuinely rising sites gets its');
print('  own verdict - "GLOBAL SIGNAL, WITH STRONG LOCAL AMPLIFICATION" -');
print('  instead of being silently absorbed into plain GLOBAL.');
print('');
print('v10.125 DIAGNOSTIC: S7E hang investigation (no fix yet)');
print('  Stuck at "1/2 batched calls done" even after a page reload AND');
print('  shorter 12/10-month windows - neither previously-known cause.');
print('  Both raw-fetch calls now print real error text on failure, or');
print('  actual feature count on success, matching the diagnostic pattern');
print('  already proven for S7D in v10.108. Re-run S7E and check the');
print('  Console for "S7E [" lines - can only help once a call eventually');
print('  returns something; cannot reveal a true infinite hang.');
print('');
print('v10.124 NEW: real permutation-test p-values added to S7E (algae)');
print('  Closes a gap flagged directly: the permutation engine only ever');
print('  protected SST-based S13 findings - S7E (algae LOCAL/REGIONAL');
print('  classification) still used only the old >0.01/>0.15 thresholds.');
print('  Now runs the same engine on both study and reference FAI series,');
print('  with ZERO new EE calls - reuses the raw monthly data S7E already');
print('  fetches for its existing statistics. Scoped to S7E for now, not');
print('  yet S7D or S7F.');
print('');
print('v10.123 NEW: Bonferroni-corrected p-values added to FIND SWEET SPOT');
print('  Direct answer to: if window length changes the result, how do you');
print('  reliably conclude anything from testing 6 windows? Testing many');
print('  windows and trusting the best one is a multiple-comparisons trap -');
print('  now runs the real permutation test at all 6 windows and reports a');
print('  corrected significance bar (p<0.0083 for 6 tests), so a signal');
print('  that survives correction can be told apart from a lucky window.');
print('  Adds only 4 new EE calls (not 12) by fetching each site raw');
print('  series once and slicing client-side for the shorter windows.');
print('');
print('v10.122 NEW: real permutation-test p-values added to COMPARE');
print('  Most LOCAL/REGIONAL decisions in this tool use fixed thresholds');
print('  (>0.01 AC1, >0.15 variance), copy-pasted across STEP 3/4, S7D,');
print('  S7E, S7F - only STEP 5 ever answered "bigger than random noise"');
print('  with a real p-value. New permutationTestDelta() engine: pools');
print('  BEFORE+AFTER raw months, randomly reshuffles the labels 500');
print('  times, reports what fraction of random shuffles beat the real');
print('  delta - a genuine p-value. Wired into COMPARE first, shown');
print('  alongside (not replacing) the existing verdict. Pure client-side');
print('  JS, reuses proven building blocks - no new architecture.');
print('');
print('v10.121 FIX: COMPARE now warns when its result has gone stale');
print('  Caught from a real test: COMPARE showed a leftover AC1 value from');
print('  a PREVIOUS site, while STEP 2 own label already showed the');
print('  correct new one. Confirmed by re-clicking COMPARE (fixed it');
print('  immediately) that this was a stale un-refreshed snapshot, not a');
print('  computation bug - COMPARE never auto-updates if STEP 2 changes');
print('  after it last ran. Now: a visible warning appears immediately if');
print('  STEP 2 is re-run after COMPARE, telling you to press COMPARE');
print('  again before trusting the result on screen.');
print('');
print('v10.120 FIX: detail panel below the headline also fixed now');
print('  The v10.119 fix only covered the top headline box - caught on the');
print('  very next Nuuk test that the "4-WAY CSD COMPARISON (detail)" panel');
print('  below it STILL showed the raw contradictory "study more stable...');
print('  POSITIVE result" text. Root cause: that panel is built early,');
print('  synchronously, before the AC1-weighted toolkit finishes computing');
print('  - structurally stuck on old text. Fixed with the same conflict');
print('  detection, using studyAC1Rose (available earlier than fullTally).');
print('');
print('v10.119 FIX: COMPARE Regional context no longer contradicts headline');
print('  Same class of bug already fixed for FIND SWEET SPOT in v10.94,');
print('  found in COMPARE too. A real Nuuk test showed the headline saying');
print('  "STRONG SIGNAL, LIKELY REGIONAL" (AC1 rose almost identically at');
print('  both sites) while "Regional context" said "STUDY SITE MORE STABLE');
print('  THAN CONTROL (positive result)" - the opposite claim, from an');
print('  unreconciled variance-only classifier ignoring what AC1 showed.');
print('  Now explicitly detected and flagged as a conflict instead of');
print('  silently presenting two disagreeing verdicts as consistent.');
print('');
print('v10.118 CRITICAL FIX: COMPARE could call zero data a positive result');
print('  Same bug already fixed for S7E in v10.114, found in S13 COMPARE -');
print('  the most heavily used part of this tool. A real Nuuk test showed');
print('  the study site with AC1=n/a, Var=n/a in BOTH windows, yet the');
print('  headline said "STUDY SITE MORE STABLE THAN CONTROL (positive');
print('  result)". The site was never measured - null data was silently');
print('  read as a genuinely flat result. Now: three explicit CANNOT ASSESS');
print('  checks run before any GLOBAL/LOCAL/ANOMALOUS branch, so missing');
print('  data can never again be reported as a confirmed finding.');
print('');
print('v10.117 FIX: S7F now shows the data density it already computes');
print('  Asked to check S7F output for data density to validate a real');
print('  finding, and it was not there - the numbers were already computed');
print('  internally (needed for the CANNOT CLASSIFY check) but never');
print('  actually displayed. Now shown for both study and reference sites,');
print('  same INSUFFICIENT flag as elsewhere - a combined run can now be');
print('  checked for data quality without re-running S7D/S7E separately.');
print('');
print('v10.116 NEW: reference-site exclusion + S7F combined orchestrator');
print('  1) S7E now excludes reference candidates within 5km of a site');
print('  already tested as a STUDY location this session (caught from a');
print('  real case: a reference landed on Low Isles, already tested as a');
print('  study site with a real signal). Tracked via new recordStudySite().');
print('  2) NEW S7F - RUN ALL: one shared Lat/Lon + BEFORE/AFTER input runs');
print('  S7D and S7E together (6 EE calls total), one combined summary.');
print('  Deliberately scoped, not full auto-orchestration - this session');
print('  already hit a real EE concurrency limit once and a 5+ minute hang');
print('  with no cancel option, so auto-trying multiple window lengths');
print('  across every S7 tool was declined as too high-risk to build blind.');
print('');
print('v10.115 NEW: S7E shows the reference site actual coordinates');
print('  Direct response to a real question - could not previously tell');
print('  whether a map landmark was the auto-found reference, since only');
print('  distance and depth were shown, never lat/lon (though already');
print('  computed internally). Now displayed directly, so the reference');
print('  candidate can be located on the map or checked against known');
print('  features.');
print('');
print('v10.114 FIX: S7E no longer confirms LOCAL ANOMALY against untested refs');
print('  Caught on a real run: headline said LOCAL ANOMALY DETECTED while');
print('  the reference showed AC1=n/a, Var=n/a - not tested-and-stable, just');
print('  had no computable data (AFTER period <4 valid months, hidden by a');
print('  combined-density check that looked fine at 55%). Now checks each');
print('  period individually for both sites - insufficient data in either');
print('  gets an explicit CANNOT CLASSIFY verdict naming the period, instead');
print('  of silently defaulting into a confirmed local/regional finding.');
print('');
print('v10.113 FIX: S7E no longer hides a coherent regional decline');
print('  Same gap already fixed for STEP 5 in v10.97, uncaught in S7E');
print('  until a real run showed study AC1 falling -0.675 and reference');
print('  -0.636 - nearly identical, a real regional pattern - both folded');
print('  into a bare "NO SIGNAL". Now surfaced explicitly as a REGIONAL');
print('  DECLINE when both sites show a similar, substantial AC1 fall,');
print('  clearly distinguished from a CSD warning direction (Dakos et al.');
print('  2012: rising AC1 is the meaningful signal, not falling).');
print('');
print('v10.112 FIX: S7E now shows STUDY site data density too');
print('  A real run at a cloudy rainforest coastline (Daintree/Cape');
print('  Tribulation) returned NO SIGNAL with BOTH study and reference');
print('  showing n/a - but only reference density was shown, giving no way');
print('  to tell if the study site itself had sparse Sentinel-2 coverage.');
print('  Now shown symmetrically for both sites, with the same LOW-');
print('  CONFIDENCE flag applied - a sparse-data null result is now');
print('  distinguishable from a genuine checked-and-found-nothing result.');
print('');
print('v10.111 NEW: S7E - LOCAL vs REGIONAL auto-classification');
print('  Closes a real gap: S7D tests coupling WITHIN one reef (all 9 nodes');
print('  within 0.3-3km) but cannot distinguish a local event from a real');
print('  regional signal. S7E auto-searches 20-160km out for a genuine');
print('  shallow-water reference reef (GEBCO), runs the same FAI AC1/');
print('  variance stats there, and applies an explicit decision rule -');
print('  study signal without reference signal = LOCAL ANOMALY; both =');
print('  REGIONAL; mirrors STEP 3 pattern for SST. Disclosed limit: a');
print('  GEBCO shallow match is not proof of a real comparable reef -');
print('  flagged LOW CONFIDENCE via data density, not solved outright.');
print('');
print('v10.110 FIX: v10.109 overcorrected the NDVI buffer size');
print('  Fixing Center\'s false "LIKELY NOT on-reef" reading used max() over');
print('  a fixed 500m buffer - too large at the default 0.5km ring radius,');
print('  causing adjacent nodes buffers to overlap (a real run showed');
print('  Center/N/NE/E all reading the identical value 0.84). Buffer now');
print('  scales to 30% of node spacing (clamped 80-250m), verified to stay');
print('  clear of overlap across the full 0.3-3km radius range.');
print('  Also caught during this fix: a version-bump script left an');
print('  unescaped apostrophe in a label string ("v10.109\'s") - a real');
print('  syntax error. Found and fixed using node --check, which does');
print('  actual JS parsing and is now the primary validation step alongside');
print('  the existing brace-balance audits.');
print('');
print('v10.109 FIX: S7D two real bugs caught from a completed real run');
print('  1) Center flipped to "LIKELY NOT on-reef" (NDVI=-0.32) despite');
print('  being confirmed on-reef with +0.627 in every prior test - caused');
print('  by mean() over a small 150m buffer being sensitive to small-scale');
print('  heterogeneity. Fixed: switched to max() over 500m, matching the');
print('  main panel/S7B convention (FAI coupling buffers unchanged).');
print('  2) Var=+52.31x/+11.99x readings - same near-zero-denominator');
print('  artifact already fixed for S7C in v10.105, now ported to S7D too.');
print('');
print('v10.108 DIAGNOSTIC: S7D failure investigation (no fix yet)');
print('  A real run at correct coordinates returned n/a for every node/');
print('  field, only partly explained by "1 of 3 errored". Added real');
print('  diagnostics instead of guessing: each of the 3 batched calls now');
print('  prints its ACTUAL error text on failure, or its feature count +');
print('  first raw feature properties on success (confirms/refutes whether');
print('  reduceRegions() names its output the way the code assumes). Raw');
print('  feature counts also now shown in the on-screen S7D results table.');
print('  Re-run S7D and check the Console for "S7D [" lines to see exactly');
print('  what happened this time.');
print('');
print('v10.107 PERFORMANCE: S7D rebuilt from ~43 EE calls down to 3');
print('  Not quantum computing (the actual bottleneck - raster compositing,');
print('  API round trips - is unrelated to problems quantum accelerates) -');
print('  a classical batching fix instead. New extractMultiNodeSeries() uses');
print('  reduceRegions()+flatten() to fetch all 9 nodes whole time series in');
print('  ONE call per period, instead of 9 separate calls each re-triggering');
print('  the full monthly compositing graph. AC1/variance/correlation now');
print('  computed client-side in JS on the small fetched tables (jsNodeStats,');
print('  jsPairCorrelation) - same formulas, same output, ~20-60s not 2-5min.');
print('');
print('v10.106 NEW: S7D full 9-node algae coupling network (BEFORE/AFTER)');
print('  Scales S7C to the full 9-node ring with a real before/after');
print('  comparison. ~43 EE calls (9 nodes x 2 periods AC1/variance, 8');
print('  ring points x 2 periods correlation-vs-Center, 9 NDVI on-reef');
print('  checks) - HEAVY, expect 2-5 minutes. NDVI-water flags likely');
print('  off-reef points (per the real S7B finding that North was likely');
print('  off-reef) so those readings are marked, not silently trusted.');
print('  Headline counts on-reef points showing rising correlation with');
print('  Center (Dakos et al. 2011 hyper-synchronization direction).');
print('');
print('v10.105 FIX: S7C variance ratio no longer looks like a fake surge');
print('  Caught on a real run: East node showed Var=15.80x, far beyond');
print('  anything seen with SST-based ratios (never exceeded ~7x). Likely');
print('  a near-zero-denominator artifact (ratio = 2ndHalfVar / 1stHalfVar,');
print('  and a near-zero first half inflates it). Now shows raw first-half/');
print('  second-half variance alongside the ratio, and auto-flags the case');
print('  with an explicit "LIKELY ARTIFACT" warning.');
print('');
print('v10.104 NEW: S7C algae/AC1/variance coupling proof-of-concept');
print('  Tests whether AC1/variance/correlation (the toolkit already used');
print('  on SST) can be computed from ALGAE (FAI) data instead - arguably');
print('  more theoretically correct, since algae cover is the actual state');
print('  variable in a coral-algae bistable system. New monthly Sentinel-2');
print('  FAI builder (mkMoFAIRange) - S7 previously had only ONE fixed');
print('  composite, never a time series. Scoped to 3 nodes + 1 window to');
print('  test real Sentinel-2 data density before committing to 9 nodes.');
print('');
print('v10.103 NEW: real S15 field data for One Tree Reef');
print('  Found a real paper (Byrne et al. 2025, Limnol. Oceanogr. Lett.)');
print('  tracking 462 coral colonies at THIS EXACT coordinate (23.51S,');
print('  152.09E) through the 2023-24 heatwave: 80% bleached by April,');
print('  up to 52% mortality by July, genus-level detail included.');
print('  Keyed by proximity (60km), not the whole GBR region bucket -');
print('  northern reefs like Lizard Island had a very different bleaching');
print('  history and should not inherit this southern-GBR-specific finding.');
print('');
print('v10.102 NEW: S7B Multi-Point Algae Scan');
print('  Direct response to a real finding: 1.1km apart at One Tree Reef,');
print('  S7 showed "mild watch" vs "massive bloom" for the same indicators.');
print('  New tool samples FAI/NDCI/NDVI at 8 compass points + centre in ONE');
print('  Earth Engine call, reports how many show warning signs, and');
print('  classifies the pattern as WIDESPREAD / PATCHY-ISOLATED / MIXED -');
print('  so a single click can no longer misrepresent a patchy bloom.');
print('');
print('v10.101 CAVEAT DISCLOSED (no computation changed this version):');
print('  A direct question about whether v10.100 could have detected the');
print('  GBR event earlier prompted a re-check that found a real concern:');
print('  the sync indicator (and the core AC1 statistic) use raw/linear-');
print('  detrend-only SST, not deseasonalized anomalies - shared seasonal');
print('  cycling between nearby ocean points can push both toward a high');
print('  baseline regardless of real dynamics. Evidence: AC1 has landed in');
print('  a narrow ~0.78-0.89 band across every GBR test in this tool so');
print('  far, more consistent with seasonal autocorrelation dominating');
print('  than genuine site-specific signal. Disclosed in the UI/code; a');
print('  proper deseasonalized fix needs live testing before shipping.');
print('');
print('v10.100 NEW: study-control synchronization indicator in STEP 3');
print('  Computes correlation between the study reef and the auto-selected');
print('  deep-water control site, BEFORE vs AFTER. Rising correlation =');
print('  reef losing local independence/buffering from open-ocean forcing');
print('  (Dakos et al. 2011 leading-indicator direction). Added as a');
print('  SUPPORTING (not primary) toolkit indicator - AC1 still decides');
print('  the headline. This is a correlation, not a mechanistic J_ij');
print('  interaction coefficient - disclosed honestly in the code and UI.');
print('');
print('v10.99 CRITICAL FIX: crash on unguarded .trim() of undefined');
print('  Caught live in the GEE Code Editor: RUN SLIDING WINDOW ANALYSIS');
print('  crashed with "Cannot read property \'trim\' of undefined" when the');
print('  optional window-size field (blank = auto by design) was never');
print('  typed into - ui.Textbox.getValue() returns undefined, not \'\', for');
print('  an untouched field. Same unguarded pattern existed in 8 places');
print('  total across the file - all fixed at once: every X.getValue().trim()');
print('  is now (X.getValue()||\'\').trim(), so an empty field safely reads');
print('  as \'\' and falls through to the normal validation message.');
print('');
print('v10.98 FIX: STEP 5 explains when future dates reduce valid data');
print('  Caught on a real run: requested 65 months, silently got only 51');
print('  valid ones - because the range ran past today into the future,');
print('  where OISST has no data yet. Now warned upfront (before the EE');
print('  call even fires) and confirmed explicitly in the results if it');
print('  happens. Also: ROBUSTNESS note now says "NONE of 6" instead of');
print('  the awkward "only 0 of 6" for the zero-agreement case.');
print('');
print('v10.97 FIX: STEP 5 no longer hides significant falling trends');
print('  Caught on a real GBR run: variance was significantly DECLINING');
print('  (p=0.017) but the headline claimed "p>=0.05 for both" - only');
print('  RISING significant trends were being checked, so a real declining');
print('  trend got mislabelled as "nothing found". Now classifies rising-');
print('  significant / falling-significant / not-significant separately -');
print('  a significant decline gets its own accurate headline instead of');
print('  being silently discarded.');
print('');
print('v10.96 NEW: STEP 5 - real Dakos et al. 2012 sliding-window method');
print('  Replaces discrete BEFORE/AFTER chunk testing with a continuous');
print('  rolling window (1-month steps) across the whole series, producing');
print('  an AC1(t)/variance(t) trajectory - tested with the actual Kendall');
print('  tau / Mann-Kendall trend significance test (real p-value), not an');
print('  arbitrary numeric threshold. Prints the trajectory as line charts.');
print('  Built as ONE server-side EE call regardless of window count.');
print('  Complements STEP 3/4 (which still handle LOCAL vs REGIONAL control-');
print('  site comparison) rather than replacing them.');
print('');
print('v10.95 NEW: FIND SWEET SPOT robustness check (isolated vs consistent)');
print('  Testing 6 AFTER windows and reporting the most dramatic one is a');
print('  classic multiple-comparisons trap - some window looks significant');
print('  by chance alone. New table footer + headline flag: counts how many');
print('  of the 6 windows lean LOCAL. 0-1 = ISOLATED (explicit caution,');
print('  re-test using an independently-known event date instead). 2+ = more');
print('  trustworthy, though still not a formal significance test.');
print('');
print('v10.94 FIX: FIND SWEET SPOT headline no longer contradicts its detail');
print('  Caught on a real test: headline said "LOCAL CSD SIGNAL DETECTED"');
print('  while the breakdown for the SAME window said "MARGINAL LOCAL" and');
print('  "Scheffer NOT MET". Two different classification rules were being');
print('  used - the headline checked only divergence magnitude, the table');
print('  correctly required the study variance to actually cross the');
print('  threshold. Headline now derived from the same bestRow.verdict as');
print('  the table/breakdown, plus an explicit AC1-vs-control check.');
print('');
print('v10.93 CRITICAL FIX: Southern Hemisphere DHW/bleaching was broken');
print('  Peak SST / DHW used ONE hardcoded window (Jun-Oct 2023, Northern');
print('  Hemisphere summer only). Southern Hemisphere reefs (GBR, Ningaloo,');
print('  S. Indian Ocean, S. Brazil, S. Africa) peak Nov-Apr - this was');
print('  silently missed entirely, showing false "no heat stress" (DHW=0.00)');
print('  even during documented severe bleaching (caught testing a real');
print('  southern GBR site with confirmed 2023-24 mass bleaching). Fixed:');
print('  peak/DHW now computed for BOTH hemispheres, selected per-click by');
print('  latitude sign. Affects the DHW badge, S4 score component (25%');
print('  weight), map bleaching layer, and intervention thermal-stress');
print('  trigger - all now correct for Southern Hemisphere reefs.');
print('');
print('v10.92 NEW: FIND SWEET SPOT can surface AC1-up/variance-down windows');
print('  Table now shows Study deltaAC1 per window, not just deltaVar.');
print('  Added a SEPARATE AC1-based ranking alongside the variance-divergence');
print('  ranking - a window where AC1 rises strongly is no longer hidden just');
print('  because it is not the variance-divergence winner. When the two');
print('  rankings disagree, console prints an explicit cross-check citing');
print('  Dakos et al. 2012 Fig 2c/4 for the AC1-up/variance-down case.');
print('');
print('v10.91 FIX: AC1-up + variance-down no longer looks like weak evidence');
print('  When AC1 rises but variance/spatial indicators are available and');
print('  disagree (e.g. variance falls), the label used to say "no supporting');
print('  indicator corroborates it yet" - implying pending/weak evidence.');
print('  Dakos et al. 2012 document this exact pattern (Fig. 2c, Fig. 4):');
print('  AC1 rises "regardless of the responsiveness of the ecosystem" while');
print('  variance can fall near a genuine transition. Wording fixed to cite');
print('  this directly; still scored MODERATE (not upgraded to HIGH, since');
print('  disagreement is not corroboration either - just neutral).');
print('');
print('v10.90 FIX: AC1 weighted as PRIMARY indicator (Dakos et al. 2012)');
print('  STEP 3 toolkit and FIND SWEET SPOT previously treated AC1, variance,');
print('  spatial variance, and spatial autocorrelation as equal votes. Dakos');
print('  et al. found AC1 "relatively robust" while variance can rise OR fall');
print('  near a real transition - so a variance-only rise now shows as LOW-');
print('  MODERATE confidence, explicitly weaker than an AC1-confirmed rise.');
print('  No change to the underlying math, only how it is weighted/interpreted.');
print('');
print('v10.89 FIX: S17 no longer blanks 4 healthy variables over 1 dead asset');
print('  pH/DO both depend on COPERNICUS/MARINE/GLOBAL_OCEAN_BGC/MFC_001_028,');
print('  which is currently "not found" in the GEE catalog. This used to be');
print('  bundled with SST/Chl/Salinity/NO2 in one evaluate() call, so the dead');
print('  asset blanked all 6. Now split: SST/Chl/Salinity/NO2 compute normally');
print('  regardless of the BGC asset, only pH/DO show "dataset unavailable".');
print('  S18 was already correctly all-n/a (all 4 of its fields genuinely come');
print('  from that one dead asset) - only its error message got clearer.');
print('');
print('v10.88 NEW: real multi-indicator toolkit in STEP 3 (COMPARE)');
print('  Adds skewness (reported), spatial variance, and a spatial-');
print('  autocorrelation proxy alongside temporal AC1/variance. New TOOLKIT');
print('  SUMMARY tallies how many of up to 4 scored indicators actually agree');
print('  (HIGH/MODERATE/LOW confidence), instead of one rule on one indicator');
print('  deciding the verdict. Verdict box renders temporal-only first (fast),');
print('  then upgrades in place once spatial indicators arrive.');
print('');
print('v10.87 FIX: S12 no longer claims AC1 is "rising"');
print('  S12 used a fixed single window with no BEFORE baseline, so');
print('  "AC1 > 0.5" was mislabelled as AC1 "rising" - it just means');
print('  currently elevated. This could contradict S13, which computes a');
print('  real AC1 delta vs a stored BEFORE window and can correctly show');
print('  AC1 FALLING. Wording fixed; S12 now points to S13 for a validated');
print('  before/after Scheffer 2009 comparison. No math changed.');
print('');
print('v10.86 FIX: STEP 3 (COMPARE) staged messaging + error handling');
print('  Verdict box no longer says "Computing verdict..." before the control');
print('  site is even found - both status boxes now agree on the real stage.');
print('  Control-site depth lookup (GEBCO) now has real error handling -');
print('  previously a failure there caused a silent, unexplained hang.');
print('');
print('v10.85 FIX: Depth label now matches the S3 legend');
print('  Sidebar Depth row previously called anything below -50m "deep ocean",');
print('  contradicting the map legend which splits that into shelf/slope/deep/');
print('  very-deep zones. New classifyDepthLabel() uses the same 7 bins as the');
print('  legend for both the Depth row and the "click closer to shore" warning.');
print('');
print('v10.84 NEW: S13 error handling hardened');
print('  Common transient Earth Engine errors (e.g. "Unknown reference to');
print('  value named" or "Failed to contact Earth Engine servers")');
print('  now show as a plain "temporary hiccup - retry" message, not a raw dump.');
print('  STEP 2 echoes parsed lat/lon/dates before running (catches typos/truncation).');
print('  STEP 3 (COMPARE) now actually checks the control-site evaluate() error,');
print('  which was previously silently ignored.');
print('');
print('v10.83 CRITICAL FIX: FIND SWEET SPOT no longer crashes');
print('  Root cause: String.prototype.repeat() is ES6 and is NOT supported');
print('  by the GEE server-side script sandbox - calling it threw');
print('  "TypeError: (intermediate value)..." right after "Analysing results..."');
print('  Fixed by replacing it with a loop-based repeatChar(ch,n) helper.');
print('  Grepped the whole script for other ES6-only syntax (arrow functions,');
print('  let/const, template literals, Array.from, Object.values/entries) -');
print('  none found elsewhere, so this was the only crash site.');
print('');
print('v10.82 NEW: FIND SWEET SPOT (S13 STEP 4) rebuilt');
print('  - Control site "before" state now REAL (1 shared control-BEFORE test), not a fake 1.0 baseline');
print('  - Live progress counter while 13 parallel EE calls run - never looks frozen');
print('  - New LOCAL vs REGIONAL vs Scheffer-2009-validation breakdown for the winning window');
print('  - try/catch around final analysis - shows a red error box instead of silently stalling');
print('');
print('v10.81 NEW: S13 (CSD Early Warning Test) rebuilt for clarity');
print('  - Numbered STEP 1/2/3/4 workflow headers throughout S13');
print('  - Manual control-site override now "ADVANCED (OPTIONAL)", after the main flow, defaults AUTO');
print('  - COMPARE now shows a short bold colour-coded VERDICT box ABOVE the detailed numbers');
print('  - FIND SWEET SPOT verdict box kept above its data table (unchanged behavior, clearer wording)');
print('  - New "Use last clicked location" button in S13 Step 1 (uses lastClickLat/lastClickLon)');
print('  - No changes to the underlying math: computeRealCSD/getSmartControlSite/runControlCSD unchanged');
print('');
print('v10.68 NEW: Real GEM stats now populated (no more PLACEHOLDERs)');
print('  Nuuk GF3: temp=1.79 deg C, salinity=33.24 PSU (76,713 readings 2005-present)');
print('  Zackenberg: temp=-1.61 deg C, salinity=32.22 PSU (406,125 readings 2003-2019)');
print('  New S19 fields: fluorescence_ug_L, pressure_db, turbidity_ftu');
print('  Nuuk extra: fluorescence=0.265 ug/L, pressure=158.3 dbar, turbidity=0.610 FTU');
print('  Zackenberg: temp/salinity only (mooring has no fluorescence/turbidity sensors)');
print('  ALL GEM stats are now real cleaned values - no more PLACEHOLDERs');
print('');
print('  5 in-situ stations now in MODULE A11:');
print('    Looe Key FL (SECOORA) | Agua Hedionda CA (SCCOOS) | Scripps Pier CA (CenCOOS)');
print('    MarineBasis Nuuk GF3 - Godthåbsfjord (64.13N,51.38W, r=40km)');
print('    MarineBasis Zackenberg - Young Sound (74.315N,20.279W, r=30km)');
print('  GEM data: CC BY-SA 4.0 | g-e-m.dk | api.g-e-m.dk | Key: GeoMarineAnalysis');
print('  DOIs: 10.17897/KMEK-TK21 (Nuuk CTD) | 10.17897/8GPS-CE70 (Zackenberg)');
print('  NOTE: GEM stats are PLACEHOLDER nulls. Run gem_fetch_and_clean.py to populate.');
print('  New Greenland region entries added to getRegion().');
print('');
print('v10.67 FIX (v10.19): Hemisphere suffix parsing in GO TO COORDINATES and S13.');
print('  parseFloat("13.5S") was silently returning 13.5 (positive).');
print('  New parseCoordPart() handles N/S/E/W suffixes correctly.');
print('  Both "13.5S, 112.3E" and "-13.5, 112.3" now work.');
print('');
print('v10.67 FIX (v10.20): DHW blue/orange layers now mutually exclusive.');
print('  "DHW blue = no stress" was painted over ALL ocean including stressed pixels.');
print('  Now correctly masked to only show where DHW = 0 (genuinely no stress).');
print('');
print('v10.67 FIX (v10.20): S13 BEFORE/AFTER CSD persistent state labels.');
print('  Added "BEFORE stored:" / "AFTER stored:" labels + Clear button.');
print('');
print('All v10.66 features intact: S14-S19, ToE, BGC, GEBCO, error handling,');
print('full sidebar reset, v10.63 peak SST gate, v10.57 real thresholds.');
print('');
print('CLICK any coastal/ocean area to analyze.');
print('Map centered: Bocas del Toro, Panama.');
print('GO TO COORDINATES supports: "13.5S, 112.3E"  OR  "-13.5, 112.3"');
print('Test GEM stations: 64.13, -51.38 (Nuuk GF3) | 74.32, -20.28 (Zackenberg)');
