"""Azimuth-ambiguity flagging for the SB02 SAR detection table.

WHY THIS IS OFFLINE. Ambiguity detection is inherently pairwise: for each
detection, is there a much brighter detection displaced along the azimuth axis
by the ambiguity distance? That is natural on a table and painful in Earth
Engine image space, which is why the GEE script emits per-detection rows and
stops there.

WHAT AN AZIMUTH AMBIGUITY IS. SAR PRF sampling aliases energy from targets
outside the unambiguous azimuth swath, so a very bright scatterer produces
fainter replicas displaced ALONG-TRACK. They look like real vessels and sit at
a characteristic distance from their parent.

FLAG, DO NOT DELETE. The spec is right about this: the displacement band is
approximate and a real vessel can legitimately sit 5 km from a bright ship in a
shipping lane. Flagged rows stay in the table with a reason, so the analysis can
be run with and without them and the difference reported.
"""
import csv, math, sys
from collections import defaultdict

# S1 orbit geometry. The azimuth axis is the ground-track heading; ambiguities
# are displaced along it in EITHER direction, so the test is on a line (mod 180).
S1_INCLINATION_DEG = 98.18


def azimuth_axis_deg(lat_deg, orbit_pass):
    """Ground-track heading (mod 180) at a given latitude, for one orbit branch.

    orbit_pass MUST be given. asin() returns only the ASCENDING branch, and for
    a long time this function returned it unconditionally: 168.88 deg at SB02.
    The descending branch is 180 - that = 11.12 deg, which is 22.24 deg away
    (mod 180) - outside the 20 deg bearing tolerance. 185 of the 278 scenes in
    this study are DESCENDING, so real ghosts on two thirds of the data could
    never be flagged, while pairs happening to lie along the ascending line were
    flagged instead.

    That is worse than missing flags. separation_histogram() below is the
    empirical check the spec makes mandatory before the flag is trusted at all,
    and an absent 4.5-6.5 km spike is documented as meaning "ambiguities are not
    a material contaminant". With the wrong axis on two thirds of the scenes the
    spike is suppressed, and the wrong conclusion is reached confidently.

    There is deliberately NO default. A row without orbit_pass raises rather than
    silently assuming ascending, which is the exact failure being fixed.
    """
    s = math.cos(math.radians(S1_INCLINATION_DEG)) / math.cos(math.radians(lat_deg))
    s = max(-1.0, min(1.0, s))
    ascending = math.degrees(math.asin(s)) % 180.0
    p = str(orbit_pass or '').strip().upper()
    if p.startswith('ASC'):
        return ascending
    if p.startswith('DESC'):
        return (180.0 - ascending) % 180.0
    raise ValueError(
        'orbit_pass must be ASCENDING or DESCENDING, got %r. The detection table '
        'exports it as orbit_pass; without it the axis cannot be chosen and a '
        'guess would silently reproduce the bug this argument exists to fix.'
        % (orbit_pass,))


def _row_axis(row):
    """Axis for a detection row, with a readable error if the column is absent."""
    if 'orbit_pass' not in row:
        raise KeyError(
            'detection row has no orbit_pass column; re-export with '
            'sb02_sar_export.js, which emits it for every detection')
    return azimuth_axis_deg(float(row['lat']), row['orbit_pass'])


def haversine_m(a, b):
    R = 6371008.8
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dp = p2 - p1
    dl = math.radians(b[1] - a[1])
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(math.sqrt(h))


def bearing_deg(a, b):
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dl = math.radians(b[1] - a[1])
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x)) % 360.0


def angular_sep_mod180(a, b):
    d = abs((a - b) % 180.0)
    return min(d, 180.0 - d)


def flag(rows, min_km=3.5, max_km=7.5, parent_db=15.0, bearing_tol_deg=20.0):
    """Adds ambig_flag / ambig_parent / ambig_range_m / ambig_bearing_off_deg."""
    by_scene = defaultdict(list)
    for r in rows:
        by_scene[r['scene_id']].append(r)

    for r in rows:
        r['ambig_flag'] = 0
        r['ambig_parent'] = ''
        r['ambig_range_m'] = ''
        r['ambig_bearing_off_deg'] = ''

    for scene, dets in by_scene.items():
        for i, d in enumerate(dets):
            dpos = (float(d['lat']), float(d['lon']))
            axis = _row_axis(d)
            best = None
            for j, pnt in enumerate(dets):
                if i == j:
                    continue
                # parent must be markedly brighter
                if float(pnt['peak_db']) - float(d['peak_db']) < parent_db:
                    continue
                ppos = (float(pnt['lat']), float(pnt['lon']))
                dist = haversine_m(dpos, ppos)
                if not (min_km * 1000 <= dist <= max_km * 1000):
                    continue
                off = angular_sep_mod180(bearing_deg(ppos, dpos), axis)
                if off > bearing_tol_deg:
                    continue
                if best is None or off < best[2]:
                    best = (pnt, dist, off)
            if best is not None:
                d['ambig_flag'] = 1
                d['ambig_parent'] = best[0].get('scene_id', '') + ':' + \
                                    ('%.5f,%.5f' % (float(best[0]['lat']), float(best[0]['lon'])))
                d['ambig_range_m'] = '%.1f' % best[1]
                d['ambig_bearing_off_deg'] = '%.1f' % best[2]
    return rows


def separation_histogram(rows, parent_db=15.0, bins_km=1.0, max_km=15.0):
    """The EMPIRICAL check the spec asks for: confirm a bump in the 3.5-7.5 km
    band before believing the flag. If the histogram is featureless there, the
    ambiguity band is not present in this data and flagging is doing nothing."""
    by_scene = defaultdict(list)
    for r in rows:
        by_scene[r['scene_id']].append(r)
    hist = defaultdict(int)
    for dets in by_scene.values():
        for i, d in enumerate(dets):
            dpos = (float(d['lat']), float(d['lon']))
            axis = _row_axis(d)
            for j, pnt in enumerate(dets):
                if i == j or float(pnt['peak_db']) - float(d['peak_db']) < parent_db:
                    continue
                ppos = (float(pnt['lat']), float(pnt['lon']))
                dist = haversine_m(dpos, ppos) / 1000.0
                if dist > max_km:
                    continue
                if angular_sep_mod180(bearing_deg(ppos, dpos), axis) > 20.0:
                    continue
                hist[int(dist / bins_km)] += 1
    return hist


if __name__ == '__main__':
    _asc = azimuth_axis_deg(42.470793, 'ASCENDING')
    _desc = azimuth_axis_deg(42.470793, 'DESCENDING')
    _sep = abs(_asc - _desc)
    _sep = min(_sep, 180.0 - _sep)
    print('azimuth axis at SB02 (42.47 N), mod 180:')
    print('  ASCENDING  %.2f deg   DESCENDING %.2f deg   separation %.2f deg'
          % (_asc, _desc, _sep))
    print('  -> the two branches are %.2f deg apart, wider than the 20 deg bearing'
          % _sep)
    print('     tolerance, so ONE axis cannot serve both. 185 of 278 scenes here')
    print('     are DESCENDING.\n')

    if len(sys.argv) > 1:
        rows = list(csv.DictReader(open(sys.argv[1])))
        print('loaded %d detections' % len(rows))
        h = separation_histogram(rows)
        print('\nalong-azimuth separation histogram (brighter-parent pairs):')
        for k in sorted(h):
            print('  %4.1f-%4.1f km: %s (%d)' % (k, k + 1, '#' * min(60, h[k]), h[k]))
        rows = flag(rows)
        n = sum(r['ambig_flag'] for r in rows)
        print('\nflagged %d of %d detections (%.1f%%) as likely ambiguities'
              % (n, len(rows), 100.0 * n / max(1, len(rows))))
        out = sys.argv[1].replace('.csv', '_flagged.csv')
        w = csv.DictWriter(open(out, 'w', newline=''), fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
        print('wrote %s' % out)
    else:
        print('no CSV given - running self-test on synthetic geometry\n')
        asc_axis = azimuth_axis_deg(42.470793, 'ASCENDING')
        desc_axis = azimuth_axis_deg(42.470793, 'DESCENDING')
        lat0, lon0 = 42.470793, -70.24294

        def offset(lat, lon, bearing, dist_m):
            R = 6371008.8
            b = math.radians(bearing); d = dist_m / R
            p1 = math.radians(lat); l1 = math.radians(lon)
            p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(b))
            l2 = l1 + math.atan2(math.sin(b) * math.sin(d) * math.cos(p1),
                                 math.cos(d) - math.sin(p1) * math.sin(p2))
            return math.degrees(p2), math.degrees(l2)

        ghost_lat, ghost_lon = offset(lat0, lon0, asc_axis, 5000)      # 5 km along ASC axis
        cross_lat, cross_lon = offset(lat0, lon0, asc_axis + 90, 5000)  # 5 km across it
        far_lat, far_lon = offset(lat0, lon0, asc_axis, 12000)          # 12 km along it
        # Placed along the DESCENDING track at a HARD-CODED bearing, deliberately
        # not via azimuth_axis_deg(). This fixture's weakness was that it sited
        # every ghost with the same function it then tested, so a wrong axis moved
        # the ghost too and the test stayed self-consistent: with the ascending-only
        # bug restored, case E still passed. A literal breaks that circle.
        #   derivation: asin(cos(98.18) / cos(42.470793)) = 11.122 deg from north,
        #   the descending ground track; the ascending branch is 180 - 11.122.
        DESC_AXIS_LITERAL = 11.122
        assert abs(desc_axis - DESC_AXIS_LITERAL) < 0.01, (
            'azimuth_axis_deg disagrees with the independently derived descending '
            'axis: %.3f vs %.3f' % (desc_axis, DESC_AXIS_LITERAL))
        dghost_lat, dghost_lon = offset(lat0, lon0, DESC_AXIS_LITERAL, 5000)

        # EACH CASE GETS ITS OWN scene_id. The first draft of this fixture put
        # them all in one scene and the "12 km, out of band" case came back
        # flagged - correctly, because another fixture row happened to sit 7 km
        # away and 15 dB brighter, which is a textbook ambiguity pair. The code
        # was right and the test was wrong. Isolating scenes keeps each case
        # testing only the condition it names.
        A = {'orbit_pass': 'ASCENDING'}
        D = {'orbit_pass': 'DESCENDING'}

        def row(scene, lat, lon, db, pas):
            r = {'scene_id': scene, 'lat': lat, 'lon': lon, 'peak_db': db}
            r.update(pas)
            return r

        rows = [
            row('A', lat0, lon0, 5.0, A),
            row('A', ghost_lat, ghost_lon, -12.0, A),

            row('B', lat0, lon0, 5.0, A),
            row('B', cross_lat, cross_lon, -12.0, A),

            row('C', lat0, lon0, 5.0, A),
            row('C', far_lat, far_lon, -12.0, A),

            row('D', lat0, lon0, 5.0, A),
            row('D', ghost_lat, ghost_lon, 3.0, A),

            # E - a DESCENDING ghost on the DESCENDING track. Must be flagged.
            #     Under the old single-axis code this returned 0: the descending
            #     track is 22.24 deg from the ascending one, outside the 20 deg
            #     tolerance. Two thirds of this study's scenes are descending, so
            #     this single row is the bug, reproduced.
            row('E', lat0, lon0, 5.0, D),
            row('E', dghost_lat, dghost_lon, -12.0, D),

            # F - a pair lying along the ASCENDING track but on a DESCENDING
            #     scene. Must NOT be flagged. This is the case the fixture could
            #     not state before, and the one that does not depend on which
            #     function placed the point: the geometry is fixed and only the
            #     label changes, so a function ignoring orbit_pass fails it.
            row('F', lat0, lon0, 5.0, D),
            row('F', ghost_lat, ghost_lon, -12.0, D),
        ]
        flag(rows)
        names = ['A parent (bright, never a ghost)',
                 'A ghost 5km ALONG asc axis, -17 dB   [ASC]',
                 'B parent',
                 'B target 5km ACROSS asc axis, -17 dB [ASC]',
                 'C parent',
                 'C target 12km along axis (out of band)',
                 'D parent',
                 'D target 5km along axis, only -2 dB fainter',
                 'E parent                              [DESC]',
                 'E ghost 5km along DESC axis, -17 dB   [DESC]',
                 'F parent                              [DESC]',
                 'F target on ASC axis but DESC scene   [DESC]']
        want = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
        bad = 0
        for nme, r, w in zip(names, rows, want):
            ok = r['ambig_flag'] == w
            bad += (not ok)
            print('  %-46s flag=%d expected=%d  %s'
                  % (nme, r['ambig_flag'], w, 'PASS' if ok else 'FAIL'))
        print('\n%s' % ('all self-tests passed' if bad == 0 else '%d FAILED' % bad))
