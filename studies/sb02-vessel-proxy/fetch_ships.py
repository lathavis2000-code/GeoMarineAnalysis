"""Fetch the SanctSound ship-detection ground truth for SB02.

The detection intervals are NOT committed to this repository. They are ~715 KB
of derived data whose source of truth is NOAA's public archive, and a fetch
script keeps the pipeline reproducible without carrying a blob that can drift
silently out of step with the archive.

Source: gs://noaa-passive-bioacoustic/sanctsound/products/detections/sb02/
Public, unauthenticated. Writes sb02_ships.json next to this script.

VERIFIED on 2026-09-16: 18 ship-detection CSVs, 15,255 events,
2018-11-12 .. 2021-12-15, 45.7% ship-present duty cycle. If the counts below
differ when you run it, the archive has changed - investigate before using the
result, do not just accept the new number.
"""
import csv, datetime, io, json, os, sys, urllib.request

API = ('https://storage.googleapis.com/storage/v1/b/noaa-passive-bioacoustic/o'
       '?prefix=sanctsound/products/detections/sb02/&maxResults=2000')
BASE = 'https://storage.googleapis.com/noaa-passive-bioacoustic/'
EXPECT_FILES, EXPECT_EVENTS, EXPECT_DUTY = 18, 15255, 45.7


def fetch(timeout=90):
    listing = json.loads(urllib.request.urlopen(API, timeout=timeout).read())
    ships = sorted(
        (i for i in listing.get('items', [])
         if '_ships' in i['name'] and i['name'].endswith('.csv')),
        key=lambda x: x['name'])
    events = []
    for obj in ships:
        txt = urllib.request.urlopen(BASE + obj['name'], timeout=timeout)
        txt = txt.read().decode('utf-8', 'replace')
        rdr = csv.reader(io.StringIO(txt))
        next(rdr, None)                      # ISOStartTime,ISOEndTime,Labels
        for row in rdr:
            if len(row) < 3 or not row[0].strip():
                continue
            try:
                a = datetime.datetime.strptime(row[0].strip()[:19], '%Y-%m-%dT%H:%M:%S')
                b = datetime.datetime.strptime(row[1].strip()[:19], '%Y-%m-%dT%H:%M:%S')
            except ValueError:
                continue                     # malformed row: skip, do not guess
            events.append((a, b))
    events.sort()
    return len(ships), events


def main():
    nfiles, ev = fetch()
    span = (ev[-1][1] - ev[0][0]).total_seconds() / 60
    dur = sum((b - a).total_seconds() / 60 for a, b in ev)
    duty = 100 * dur / span
    print('files : %d   (expected %d)' % (nfiles, EXPECT_FILES))
    print('events: %d   (expected %d)' % (len(ev), EXPECT_EVENTS))
    print('span  : %s .. %s' % (ev[0][0].date(), ev[-1][1].date()))
    print('duty  : %.1f%%  (expected %.1f%%)' % (duty, EXPECT_DUTY))

    drift = (nfiles != EXPECT_FILES or len(ev) != EXPECT_EVENTS
             or abs(duty - EXPECT_DUTY) > 0.2)
    if drift:
        print('\nWARNING: the archive no longer matches the verified snapshot.')
        print('Every number derived downstream was computed against the old one.')

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sb02_ships.json')
    json.dump([[a.isoformat(), b.isoformat()] for a, b in ev], open(out, 'w'))
    print('\nwrote %s' % out)
    return 1 if drift else 0


if __name__ == '__main__':
    sys.exit(main())
