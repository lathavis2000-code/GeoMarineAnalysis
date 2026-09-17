"""proxy_skill.py - PROTOTYPE skill evaluation for the SB02 vessel proxy.

Joins the per-scene SAR table to the acoustic ship-detection record and asks the
study's actual question: does satellite-detected vessel presence discriminate
between times when a ship WAS acoustically detected and times when one was not?

-------------------------------------------------------------------------------
THIS IS THE FIRST SCRIPT IN THE STUDY THAT COMPUTES A STATISTIC
-------------------------------------------------------------------------------
Everything upstream extracts data and refuses to interpret it. This one
interprets, so it carries the obligations that come with that:

  - It reports n alongside every number. An AUC on 40 scenes is not a finding.
  - It reports the raw estimate AND the wind-stratified one, always, together.
  - It refuses to label a scene from a gap in the acoustic record.
  - It prints what it excluded and why, not just what survived.

It is a PROTOTYPE. It runs on the p002 table: persistence masking OFF and radii
to 10 km only. Both limits are stated in the header of every run so a number
lifted out of this output carries them.

-------------------------------------------------------------------------------
THE CONFOUND THAT WOULD MAKE THIS LIE, AND WHY THE STRATIFIED AUC IS NOT OPTIONAL
-------------------------------------------------------------------------------
Wind drives BOTH sides of this comparison:

  wind -> sea clutter -> fewer/more SAR detections at a fixed threshold
  wind -> ambient noise -> ship detection harder acoustically

So a healthy-looking raw correlation between SAR counts and acoustic presence
can be wind moving both, with no causal path between them at all. That is not a
hypothetical - it is the most likely way this study produces a confident wrong
answer, and it is why ERA5 wind is in the export.

The mitigation here is stratification: split the scenes into wind terciles and
recompute the AUC WITHIN each. If the discrimination is real it survives; if it
was wind, the within-stratum AUCs collapse toward 0.5 while the pooled one stays
high. Both are printed side by side, and a self-test below constructs exactly
that failure to prove the check can detect it.

Stratification is a weaker control than a regression with wind as a covariate.
It is chosen because it needs no fitting, no link function and no library, so
the prototype has nothing to get subtly wrong. A logistic model is the right
tool for the full study, not for the first look.

-------------------------------------------------------------------------------
NO THIRD-PARTY DEPENDENCIES
-------------------------------------------------------------------------------
Pure stdlib, like the other scripts here, so CI stays a zero-install job and the
statistics are readable rather than delegated. Run: python3 proxy_skill.py
"""
import bisect
import csv
import datetime
import math
import os
import sys


# =============================================================================
# 1. ACOUSTIC LABELS, WITH AN EXPLICIT "UNKNOWN"
# =============================================================================
#
# acoustic_label.ship_present() answers True/False. That is the right contract
# for it and the wrong one here, because of a distinction it does not make: the
# SanctSound ships product is a series of DEPLOYMENTS, and between them the
# recorder is not listening. A SAR scene that falls in such a gap gets False -
# "no ship" - when the truth is "nobody was listening".
#
# Labelling unknowns as negatives is not a small bias. Those scenes keep their
# SAR counts, so they land in the negative class carrying whatever vessels the
# satellite saw, which drags the negative-class score distribution UP and pushes
# the AUC toward 0.5. It costs sensitivity rather than inventing a result, but a
# study that concludes "no skill" because it scored its own blind spots is no
# better than one that concludes the opposite.
#
# Recorder-on periods cannot be recovered exactly from detection intervals
# alone - a long quiet stretch with the recorder ON looks the same as one with
# it OFF. The heuristic is a gap threshold, and it is a heuristic: a gap longer
# than GAP_HOURS between consecutive detections is treated as no-coverage.
# GAP_HOURS is swept in the output for exactly that reason.

def coverage_gaps(merged, gap_hours):
    """Gaps longer than gap_hours between consecutive merged intervals.

    Returns [(gap_start, gap_end)], each half-open. merged must be sorted and
    non-overlapping - load_intervals() guarantees both.
    """
    gap = datetime.timedelta(hours=gap_hours)
    out = []
    for i in range(len(merged) - 1):
        end = merged[i][1]
        nxt = merged[i + 1][0]
        if nxt - end > gap:
            out.append((end, nxt))
    return out


def label_scene(t, merged, gaps, tol_min):
    """True (ship), False (no ship), or None (not covered).

    None is returned for a time outside the record entirely or inside a
    no-coverage gap. Callers must DROP those rows, not coerce them.
    """
    if not merged:
        return None
    if t < merged[0][0] or t >= merged[-1][1]:
        return None
    starts = [g[0] for g in gaps]
    i = bisect.bisect_right(starts, t) - 1
    if 0 <= i < len(gaps) and gaps[i][0] <= t < gaps[i][1]:
        return None
    tol = datetime.timedelta(minutes=tol_min)
    istart = [a for a, _ in merged]
    j = bisect.bisect_right(istart, t + tol) - 1
    for k in (j, j + 1):
        if 0 <= k < len(merged):
            a, b = merged[k]
            if (a - tol) <= t < (b + tol):
                return True
    return False


# =============================================================================
# 2. STATISTICS
# =============================================================================

def _ranks(xs):
    """Ranks with ties averaged. 1-based, as the Spearman formula expects."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    out = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def auc(pos, neg):
    """Area under the ROC curve, via the Mann-Whitney identity.

    Ties contribute 0.5, which is what makes this correct for a COUNT score:
    'zero vessels' is by far the most common value on both sides, so a tie rule
    that rounded either way would bias the result by a lot. Returns None when
    either class is empty - an AUC needs both.
    """
    if not pos or not neg:
        return None
    allv = list(pos) + list(neg)
    r = _ranks(allv)
    rsum = sum(r[:len(pos)])
    n1, n2 = len(pos), len(neg)
    u = rsum - n1 * (n1 + 1) / 2.0
    return u / (n1 * n2)


def auc_se(a, n1, n2):
    """Hanley-McNeil standard error. Approximate, and labelled as such wherever
    it is printed: it assumes an exponential score distribution, which counts
    with a large spike at zero are not. It is here to stop a reader treating
    AUC = 0.61 on n = 30 as the same object as AUC = 0.61 on n = 300."""
    if a is None or n1 == 0 or n2 == 0:
        return None
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    v = (a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2)
    return math.sqrt(v) if v > 0 else 0.0


def spearman(xs, ys):
    """Rank correlation, ties averaged. None if either side is constant."""
    if len(xs) < 3:
        return None
    rx, ry = _ranks(xs), _ranks(ys)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    sxx = sum((v - mx) ** 2 for v in rx)
    syy = sum((v - my) ** 2 for v in ry)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def terciles(values):
    """Cut points at roughly the 1/3 and 2/3 quantiles, as (lo, hi).

    Strata are then w <= lo, lo < w <= hi, w > hi.

    THE INDEX IS (n-1)//3, NOT n//3, AND THE DIFFERENCE IS NOT COSMETIC. With
    n = 120 values in three tied groups of 40, v[n//3] = v[40] is the FIRST
    element of the second group, so "w <= lo" swallows groups one and two, the
    third stratum gets everything else and the top stratum is empty. The
    stratified check then silently degrades to a two-way split - which is how
    the confound self-test below first reported n=80 / n=40 / n=0.
    v[(n-1)//3] lands on the LAST element of the first group instead, and the
    three strata come out 40 / 40 / 40.

    Heavy ties can still collapse the strata - three distinct wind speeds
    cannot be cut into four. That is reported through the per-stratum n rather
    than hidden by forcing equal bins, because forcing equal bins would split
    identical wind speeds across strata, which is a worse lie than an uneven
    split. Returns (None, None) below 3 values.
    """
    v = sorted(values)
    n = len(v)
    if n < 3:
        return (None, None)
    return (v[(n - 1) // 3], v[2 * (n - 1) // 3])


# =============================================================================
# 3. THE JOIN
# =============================================================================

def read_scene_table(path):
    """Per-scene CSV -> list of dicts, numeric fields coerced, bad rows dropped.

    Unparseable numerics are dropped rather than defaulted: a row whose vessel
    count will not parse is a row whose count is unknown, and a silent 0 there
    is a fabricated negative.
    """
    rows = []
    with open(path, newline='') as fh:
        for raw in csv.DictReader(fh):
            try:
                raw['_t'] = datetime.datetime.strptime(
                    raw['timestamp_iso'], '%Y-%m-%dT%H:%M:%SZ')
            except (ValueError, KeyError, TypeError):
                continue
            rows.append(raw)
    return rows


def _num(row, key):
    try:
        return float(row[key])
    except (ValueError, KeyError, TypeError):
        return None


def evaluate(rows, tag, merged, gaps, tol_min, require_full_cov=True):
    """Skill of vessels_<tag> against the acoustic label. Returns a dict.

    Exclusions are COUNTED, not just applied, and every count is returned. A
    filter that quietly removes most of the data is the difference between a
    result and an artefact, and the caller cannot see that from the AUC alone.
    """
    kept, drop_cov, drop_label, drop_val, drop_wind = [], 0, 0, 0, 0
    for r in rows:
        count = _num(r, 'vessels_' + tag)
        if count is None:
            drop_val += 1
            continue
        if require_full_cov and _num(r, 'full_cov_' + tag) != 1:
            # A scene clipping the disc sees less water and UNDERCOUNTS. Keeping
            # it adds a low count that has nothing to do with vessel presence.
            drop_cov += 1
            continue
        lab = label_scene(r['_t'], merged, gaps, tol_min)
        if lab is None:
            drop_label += 1
            continue
        wind = _num(r, 'wind_speed_ms')
        if _num(r, 'era5_matched') != 1:
            # wind_speed_ms is 0 on an unmatched row, and 0 does not mean calm.
            wind = None
            drop_wind += 1
        kept.append((count, lab, wind))

    pos = [c for c, l, _ in kept if l]
    neg = [c for c, l, _ in kept if not l]
    a = auc(pos, neg)
    res = {
        'tag': tag, 'n': len(kept), 'n_pos': len(pos), 'n_neg': len(neg),
        'base_rate': (len(pos) / len(kept)) if kept else None,
        'auc': a, 'auc_se': auc_se(a, len(pos), len(neg)),
        'spearman': spearman([c for c, _, _ in kept],
                             [1.0 if l else 0.0 for _, l, _ in kept]),
        'dropped_coverage': drop_cov, 'dropped_unlabelled': drop_label,
        'dropped_unparseable': drop_val, 'no_wind': drop_wind,
        'strata': [],
    }

    withwind = [(c, l, w) for c, l, w in kept if w is not None]
    lo, hi = terciles([w for _, _, w in withwind])
    if lo is not None:
        for name, sel in (
            ('wind low ', lambda w: w <= lo),
            ('wind mid ', lambda w: lo < w <= hi),
            ('wind high', lambda w: w > hi),
        ):
            sp = [c for c, l, w in withwind if sel(w) and l]
            sn = [c for c, l, w in withwind if sel(w) and not l]
            sa = auc(sp, sn)
            res['strata'].append({
                'name': name, 'n': len(sp) + len(sn), 'n_pos': len(sp),
                'auc': sa, 'auc_se': auc_se(sa, len(sp), len(sn)),
            })
    return res


def format_result(res):
    def f(x, nd=3):
        return 'n/a' if x is None else ('%.*f' % (nd, x))
    out = []
    out.append('  radius %-5s n=%-4d (%d ship / %d no-ship, base rate %s)'
               % (res['tag'], res['n'], res['n_pos'], res['n_neg'],
                  f(res['base_rate'], 3)))
    se = res['auc_se']
    ci = ('' if (res['auc'] is None or se is None)
          else '  [%s, %s]' % (f(res['auc'] - 1.96 * se), f(res['auc'] + 1.96 * se)))
    out.append('    pooled    AUC %s%s   Spearman rho %s'
               % (f(res['auc']), ci, f(res['spearman'])))
    for s in res['strata']:
        out.append('    %s AUC %s   (n=%d, %d ship)'
                   % (s['name'], f(s['auc']), s['n'], s['n_pos']))
    out.append('    excluded: %d partial-coverage, %d unlabelled (record gap), '
               '%d unparseable; %d rows lack wind'
               % (res['dropped_coverage'], res['dropped_unlabelled'],
                  res['dropped_unparseable'], res['no_wind']))
    return '\n'.join(out)


# =============================================================================
# 4. SELF-TESTS
# =============================================================================
# The important one is test_wind_confound. Every other test here checks that a
# statistic is computed correctly; that one checks that the DESIGN can catch the
# specific way this study would fool itself. A confound control nobody has seen
# fire is not a control.

def _selftests():
    fails = []

    def check(name, got, want, eps=1e-9):
        ok = (got is None and want is None) or (
            got is not None and want is not None and abs(got - want) < eps)
        print('  %-46s %-8s expected %-8s %s'
              % (name, 'None' if got is None else '%.4f' % got,
                 'None' if want is None else '%.4f' % want,
                 'PASS' if ok else 'FAIL'))
        if not ok:
            fails.append(name)

    check('AUC, perfectly separated', auc([3, 4, 5], [0, 1, 2]), 1.0)
    check('AUC, perfectly inverted', auc([0, 1, 2], [3, 4, 5]), 0.0)
    check('AUC, identical distributions', auc([1, 2, 3], [1, 2, 3]), 0.5)
    check('AUC, all values tied', auc([0, 0, 0], [0, 0, 0]), 0.5)
    # The tie rule matters most where real count data lives: a big spike at 0.
    check('AUC, one class all zero, other mixed',
          auc([0, 0, 1, 1], [0, 0, 0, 0]), 0.75)
    check('AUC, empty positive class', auc([], [1, 2]), None)
    check('Spearman, monotone increasing', spearman([1, 2, 3, 4], [1, 2, 3, 4]), 1.0)
    check('Spearman, monotone decreasing', spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
    # Non-linear but monotone: Pearson would not return 1.0 here, Spearman must.
    check('Spearman, monotone but not linear',
          spearman([1, 2, 3, 4], [1, 4, 9, 16]), 1.0)
    check('Spearman, constant input', spearman([1, 1, 1, 1], [1, 2, 3, 4]), None)

    # --- unknowns are not negatives -----------------------------------------
    base = datetime.datetime(2019, 1, 1)
    merged = [(base, base + datetime.timedelta(hours=1)),
              (base + datetime.timedelta(days=10),
               base + datetime.timedelta(days=10, hours=1))]
    gaps = coverage_gaps(merged, gap_hours=24)
    print('\n  coverage gaps found: %d (expect 1)' % len(gaps))
    if len(gaps) != 1:
        fails.append('coverage_gaps count')
    cases = [
        ('inside first interval', base + datetime.timedelta(minutes=30), True),
        ('inside the 10-day gap', base + datetime.timedelta(days=5), None),
        ('before the record', base - datetime.timedelta(days=1), None),
        ('after the record', base + datetime.timedelta(days=30), None),
    ]
    for name, t, want in cases:
        got = label_scene(t, merged, gaps, 0)
        ok = got is want
        print('  %-46s %-8s expected %-8s %s'
              % ('label: ' + name, str(got), str(want), 'PASS' if ok else 'FAIL'))
        if not ok:
            fails.append('label ' + name)

    # --- the one that matters ------------------------------------------------
    # A pure wind confound: SAR count and acoustic label are both driven by
    # wind and are conditionally INDEPENDENT given it. Pooled AUC must look
    # strong; every within-stratum AUC must sit at 0.5. If the stratified check
    # cannot separate these two, it is decoration.
    print('\n  wind confound, constructed so SAR and acoustics share only wind:')
    rows, merged2 = _confounded_rows()
    gaps2 = coverage_gaps(merged2, gap_hours=24)
    res = evaluate(rows, '10km', merged2, gaps2, tol_min=0)
    print(format_result(res))
    pooled_strong = res['auc'] is not None and res['auc'] > 0.75
    strat = [s['auc'] for s in res['strata'] if s['auc'] is not None]
    strat_flat = len(strat) == 3 and all(abs(a - 0.5) < 0.12 for a in strat)
    print('  pooled AUC > 0.75 ....................... %s'
          % ('PASS' if pooled_strong else 'FAIL (%s)' % res['auc']))
    print('  all three stratified AUCs within 0.12 of 0.5 ... %s'
          % ('PASS' if strat_flat else 'FAIL (%s)' % strat))
    if not pooled_strong:
        fails.append('confound: pooled AUC not strong')
    if not strat_flat:
        fails.append('confound: stratification did not flatten it')

    # --- and the mirror: real signal must SURVIVE stratification -------------
    # Without this, a control that always returned 0.5 would pass the test above
    # and destroy every true result in the study.
    print('\n  real signal, independent of wind:')
    rows3, merged3 = _signal_rows()
    gaps3 = coverage_gaps(merged3, gap_hours=24)
    res3 = evaluate(rows3, '10km', merged3, gaps3, tol_min=0)
    print(format_result(res3))
    strat3 = [s['auc'] for s in res3['strata'] if s['auc'] is not None]
    survives = len(strat3) == 3 and all(a > 0.75 for a in strat3)
    print('  all three stratified AUCs > 0.75 ........ %s'
          % ('PASS' if survives else 'FAIL (%s)' % strat3))
    if not survives:
        fails.append('signal: stratification destroyed a real effect')

    return fails


def _row(t, count, wind, cov=1):
    return {'_t': t, 'vessels_10km': str(count), 'full_cov_10km': str(cov),
            'wind_speed_ms': str(wind), 'era5_matched': '1',
            'timestamp_iso': t.strftime('%Y-%m-%dT%H:%M:%SZ')}


def _confounded_rows():
    """SAR count and acoustic presence both determined by wind, nothing else.

    Within a wind stratum the count carries no information about the label, so
    the within-stratum AUC is 0.5 by construction and the pooled one is not.

    TWO THINGS THIS FIXTURE HAS TO GET RIGHT, both of which it got wrong first:

    1. Positives must be INTERLEAVED in time, not blocked. A contiguous run of
       negatives longer than the gap threshold becomes a no-coverage gap, and
       coverage_gaps() then - correctly - unlabels every negative in it. The
       first version blocked them and dropped all 60 negatives, which made the
       fixture untestable rather than the code wrong.
    2. The count jitter must be LABEL-BLIND. Jitter of i % 2 against a label of
       i % 2 makes the count predict the label perfectly inside a stratum, which
       is the opposite of what this fixture is for. (i // 2) % 2 is independent
       of i % 2 and of i % 5, the two label rules used below.
    """
    base = datetime.datetime(2019, 1, 1)
    rows, positives = [], []
    t = base
    #                wind   label rule -> positive rate
    for level, (wind, rule) in enumerate((
            (2.0, lambda i: i % 5 == 0),   # 20% present
            (6.0, lambda i: i % 2 == 0),   # 50% present
            (12.0, lambda i: i % 5 != 0),  # 80% present
    )):
        for i in range(40):
            t = t + datetime.timedelta(hours=3)
            count = 2 + 6 * level + ((i // 2) % 2)
            rows.append(_row(t, count, wind))
            if rule(i):
                positives.append(t)
    return rows, _intervals_around(positives)


def _signal_rows():
    """Count predicts the label within every wind stratum."""
    base = datetime.datetime(2020, 1, 1)
    rows, positives = [], []
    t = base
    for wind in (2.0, 6.0, 12.0):
        for i in range(40):
            t = t + datetime.timedelta(hours=3)
            present = (i % 2 == 0)
            count = (5 if present else 0) + (i % 3)
            rows.append(_row(t, count, wind))
            if present:
                positives.append(t)
    return rows, _intervals_around(positives)


def _intervals_around(times):
    """One 20-minute interval per positive time, merged the way the real
    loader merges - so the synthetic record goes through the same shape of
    data structure the live one does."""
    iv = sorted((t - datetime.timedelta(minutes=10),
                 t + datetime.timedelta(minutes=10)) for t in times)
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


# =============================================================================
# 5. ENTRY POINT
# =============================================================================

def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    csv_path = argv[1] if len(argv) > 1 else None
    ships_path = os.path.join(here, 'sb02_ships.json')

    print('=== SB02 vessel-proxy PROTOTYPE skill evaluation ===')
    print('PROTOTYPE. The p002 table is built with persistence masking OFF and')
    print('radii to 10 km only. Fixed objects (buoys, wrecks, the mooring) are')
    print('therefore still counted. That adds a near-constant offset per orbit,')
    print('which does not move an AUC or a rank correlation - but it does mean')
    print('absolute counts and densities here are NOT vessel counts.')
    print('No conclusion about proxy skill should be drawn from a single radius,')
    print('a single tolerance, or a pooled AUC read without its stratified row.')
    print('')

    if not csv_path or not os.path.exists(csv_path):
        print('No per-scene CSV given (or not found), so running SELF-TESTS only.')
        print('Usage: python3 proxy_skill.py <sb02_s1_vessel_detections_*.csv>')
        print('')
        fails = _selftests()
        print('\n%s' % ('all self-tests passed' if not fails
                        else '%d SELF-TEST(S) FAILED: %s' % (len(fails), fails)))
        return 0 if not fails else 1

    if not os.path.exists(ships_path):
        print('sb02_ships.json not found - run fetch_ships.py first.')
        return 1

    sys.path.insert(0, here)
    import acoustic_label
    merged = acoustic_label.load_intervals(ships_path)
    rows = read_scene_table(csv_path)
    print('scenes in table: %d' % len(rows))
    print('merged acoustic intervals: %d' % len(merged))

    modes = set(r.get('p_persistence_mode', '?') for r in rows)
    print('persistence mode recorded in the table: %s' % ', '.join(sorted(modes)))
    if modes - {'off'}:
        print('  NOTE: this table was not built with persistence off; the')
        print('  prototype caveat above may not apply as written.')

    tags = sorted(set(k[len('vessels_'):] for r in rows for k in r
                      if k.startswith('vessels_')),
                  key=lambda t: int(t.replace('km', '')))
    for gap_hours in (6, 24):
        gaps = coverage_gaps(merged, gap_hours)
        print('\n--- no-coverage gap threshold: %d h (%d gaps found) ---'
              % (gap_hours, len(gaps)))
        for tol in (0, 30, 60):
            print('\n  acoustic tolerance +/-%d min' % tol)
            for tag in tags:
                print(format_result(evaluate(rows, tag, merged, gaps, tol)))

    print('\nREAD THE STRATIFIED ROWS. A pooled AUC well above the three')
    print('wind-stratified ones is the signature of wind driving both sides,')
    print('not of a satellite proxy for acoustic ship presence.')
    print('')
    fails = _selftests()
    print('\n%s' % ('all self-tests passed' if not fails
                    else '%d SELF-TEST(S) FAILED: %s' % (len(fails), fails)))
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
