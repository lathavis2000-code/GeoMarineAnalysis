"""Acoustic ground-truth labelling for the SB02 vessel-proxy test.

Given a Sentinel-1 acquisition time, answer: was a ship acoustically detected
at that instant (within +/- tolerance)? Independent of the SAR side entirely -
it reads only the SanctSound ships detection intervals.

The intervals are half-open [start, end). Tolerance widens each interval
symmetrically, which is the right shape here: a vessel audible 10 minutes after
a SAR pass was very likely within the scene at pass time, and vice versa.
"""
import bisect, datetime, json


def load_intervals(path='sb02_ships.json'):
    raw = json.load(open(path))
    iv = sorted((datetime.datetime.fromisoformat(a),
                 datetime.datetime.fromisoformat(b)) for a, b in raw)
    # Merge overlapping/adjacent intervals so a timestamp cannot match twice and
    # so coverage arithmetic stays honest.
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def ship_present(t, merged, tol_min=0):
    """True if t falls inside any interval widened by tol_min on both sides."""
    tol = datetime.timedelta(minutes=tol_min)
    starts = [a for a, _ in merged]
    i = bisect.bisect_right(starts, t + tol) - 1
    # widening means an interval starting slightly after t can still match,
    # so check the candidate and its neighbour
    for j in (i, i + 1):
        if 0 <= j < len(merged):
            a, b = merged[j]
            if (a - tol) <= t < (b + tol):
                return True
    return False


def label_all(times, merged, tol_min=0):
    return [ship_present(t, merged, tol_min) for t in times]


def _synthetic_intervals():
    """Deterministic stand-in so the property tests run without the 715 KB data
    file, which is fetched (see fetch_ships.py) rather than committed. Shape is
    chosen to mimic the real record: frequent short events at roughly the
    observed duty cycle, so the tolerance behaviour is exercised realistically."""
    base = datetime.datetime(2019, 1, 1)
    out = []
    t = base
    for i in range(400):
        a = t + datetime.timedelta(minutes=(i * 7) % 23)
        b = a + datetime.timedelta(minutes=31)
        out.append((a, b))
        t = t + datetime.timedelta(minutes=68)
    return out


if __name__ == '__main__':
    import os
    _path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sb02_ships.json')
    if os.path.exists(_path):
        m = load_intervals(_path)
        print('using real SanctSound intervals')
    else:
        # Merge the synthetic set through the same code path the real data takes.
        _raw = sorted(_synthetic_intervals())
        m = []
        for a, b in _raw:
            if m and a <= m[-1][1]:
                if b > m[-1][1]:
                    m[-1] = (m[-1][0], b)
            else:
                m.append((a, b))
        print('sb02_ships.json not present - running property tests on synthetic '
              'intervals (run fetch_ships.py for the real record)')
    print('merged intervals: %d' % len(m))
    span = (m[-1][1] - m[0][0]).total_seconds() / 60
    cov = sum((b - a).total_seconds() / 60 for a, b in m)
    print('merged coverage: %.1f%% of span' % (100 * cov / span))

    # Property tests
    a0, b0 = m[0]
    checks = [
        ('inside an interval', ship_present(a0 + (b0 - a0) / 2, m, 0), True),
        ('exactly at start', ship_present(a0, m, 0), True),
        ('just before start, tol=0', ship_present(a0 - datetime.timedelta(minutes=5), m, 0), False),
        ('just before start, tol=15', ship_present(a0 - datetime.timedelta(minutes=5), m, 15), True),
        ('far before record', ship_present(m[0][0] - datetime.timedelta(days=30), m, 30), False),
        ('far after record', ship_present(m[-1][1] + datetime.timedelta(days=30), m, 30), False),
    ]
    bad = 0
    for name, got, want in checks:
        ok = got == want
        bad += (not ok)
        print('  %-28s %-5s expected %-5s %s' % (name, got, want, 'PASS' if ok else 'FAIL'))

    # Base rate under each tolerance, sampled on a regular grid across the record
    print('\nbase rate vs tolerance (uniform 1-h grid over the record):')
    t = m[0][0]
    grid = []
    while t < m[-1][1]:
        grid.append(t); t += datetime.timedelta(hours=1)
    for tol in (0, 15, 30, 60):
        lab = label_all(grid, m, tol)
        print('  tol +/-%2d min -> %.1f%% present  (n=%d)' % (tol, 100 * sum(lab) / len(lab), len(lab)))
    print('\n%s' % ('all property tests passed' if bad == 0 else '%d PROPERTY TEST(S) FAILED' % bad))
