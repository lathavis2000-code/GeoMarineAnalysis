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


def azimuth_axis_deg(lat_deg):
    """Ground-track heading (mod 180) at a given latitude."""
    s = math.cos(math.radians(S1_INCLINATION_DEG)) / math.cos(math.radians(lat_deg))
    s = max(-1.0, min(1.0, s))
    return math.degrees(math.asin(s)) % 180.0


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
            axis = azimuth_axis_deg(dpos[0])
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
            axis = azimuth_axis_deg(dpos[0])
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
    print('azimuth axis at SB02 (42.47 N): %.2f deg (mod 180)' % azimuth_axis_deg(42.470793))
    _ax = azimuth_axis_deg(42.470793)
    print('  -> the axis is a LINE, so %.2f deg is %.2f deg off the meridian;'
          % (_ax, min(_ax, 180.0 - _ax)))
    print('     ambiguities are displaced roughly N-S, tilted ~%.1f deg.\n'
          % min(_ax, 180.0 - _ax))

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
        axis = azimuth_axis_deg(42.470793)
        lat0, lon0 = 42.470793, -70.24294

        def offset(lat, lon, bearing, dist_m):
            R = 6371008.8
            b = math.radians(bearing); d = dist_m / R
            p1 = math.radians(lat); l1 = math.radians(lon)
            p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(b))
            l2 = l1 + math.atan2(math.sin(b) * math.sin(d) * math.cos(p1),
                                 math.cos(d) - math.sin(p1) * math.sin(p2))
            return math.degrees(p2), math.degrees(l2)

        ghost_lat, ghost_lon = offset(lat0, lon0, axis, 5000)      # 5 km along axis
        cross_lat, cross_lon = offset(lat0, lon0, axis + 90, 5000)  # 5 km across axis
        far_lat, far_lon = offset(lat0, lon0, axis, 12000)          # 12 km along axis

        # EACH CASE GETS ITS OWN scene_id. The first draft of this fixture put
        # them all in one scene and the "12 km, out of band" case came back
        # flagged - correctly, because another fixture row happened to sit 7 km
        # away and 15 dB brighter, which is a textbook ambiguity pair. The code
        # was right and the test was wrong. Isolating scenes keeps each case
        # testing only the condition it names.
        rows = [
            {'scene_id': 'A', 'lat': lat0, 'lon': lon0, 'peak_db': 5.0},
            {'scene_id': 'A', 'lat': ghost_lat, 'lon': ghost_lon, 'peak_db': -12.0},

            {'scene_id': 'B', 'lat': lat0, 'lon': lon0, 'peak_db': 5.0},
            {'scene_id': 'B', 'lat': cross_lat, 'lon': cross_lon, 'peak_db': -12.0},

            {'scene_id': 'C', 'lat': lat0, 'lon': lon0, 'peak_db': 5.0},
            {'scene_id': 'C', 'lat': far_lat, 'lon': far_lon, 'peak_db': -12.0},

            {'scene_id': 'D', 'lat': lat0, 'lon': lon0, 'peak_db': 5.0},
            {'scene_id': 'D', 'lat': ghost_lat, 'lon': ghost_lon, 'peak_db': 3.0},
        ]
        flag(rows)
        names = ['A parent (bright, never a ghost)',
                 'A ghost 5km ALONG axis, -17 dB',
                 'B parent',
                 'B target 5km ACROSS axis, -17 dB',
                 'C parent',
                 'C target 12km along axis (out of band)',
                 'D parent',
                 'D target 5km along axis, only -2 dB fainter']
        want = [0, 1, 0, 0, 0, 0, 0, 0]
        bad = 0
        for nme, r, w in zip(names, rows, want):
            ok = r['ambig_flag'] == w
            bad += (not ok)
            print('  %-46s flag=%d expected=%d  %s'
                  % (nme, r['ambig_flag'], w, 'PASS' if ok else 'FAIL'))
        print('\n%s' % ('all self-tests passed' if bad == 0 else '%d FAILED' % bad))
