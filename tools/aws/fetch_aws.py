"""기상청 API허브 AWS 매분자료(nph-aws2_min) → 10분 간격 자료로 줄여 저장.
1개 지점은 호출당 하루까지만 조회되므로 지점·일 단위로 반복 호출한다.
키는 환경변수 KMA_APIHUB_KEY에서만 읽고, URL·로그에 출력하지 않는다."""
import os, sys, time, gzip, csv, io, datetime as dt, urllib.request, urllib.parse
COLS = ['TM','STN','WD1','WS1','WDS','WSS','WD10','WS10','TA','RE','RN15','RN60','RN12H','RNDAY','HM','PA','PS','TD']
OUT = ['TM','STN','WD10','WS10','WS1_MEAN','WSS_MAX','WDS_AT_MAX','TA','RN15','HM','N1MIN']
URL = 'https://apihub.kma.go.kr/api/typ01/cgi-bin/url/nph-aws2_min'

def parse(text):
    rows = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        parts = [p.strip() for p in s.replace(',', ' ').split()]
        if len(parts) < len(COLS) or not parts[0].isdigit():
            continue
        r = dict(zip(COLS, parts[:len(COLS)]))
        for k in COLS[2:]:
            try:
                v = float(r[k]); r[k] = None if v <= -50 else v
            except ValueError:
                r[k] = None
        rows.append(r)
    return rows

def to_10min(rows):
    """10분 구간(끝 시각 기준, 예: 09:01~09:10 → 09:10)으로 묶음."""
    bins = {}
    for r in rows:
        t = dt.datetime.strptime(r['TM'], '%Y%m%d%H%M')
        end = t + dt.timedelta(minutes=(-t.minute) % 10)
        bins.setdefault((r['STN'], end), []).append(r)
    out = []
    for (stn, end), g in sorted(bins.items()):
        last = [r for r in g if dt.datetime.strptime(r['TM'], '%Y%m%d%H%M') == end]
        ref = last[0] if last else None
        ws1 = [r['WS1'] for r in g if r['WS1'] is not None]
        gust = [(r['WSS'], r['WDS']) for r in g if r['WSS'] is not None]
        mx = max(gust) if gust else (None, None)
        out.append({'TM': end.strftime('%Y%m%d%H%M'), 'STN': stn,
            'WD10': ref and ref['WD10'], 'WS10': ref and ref['WS10'],
            'WS1_MEAN': round(sum(ws1)/len(ws1), 2) if ws1 else None,
            'WSS_MAX': mx[0], 'WDS_AT_MAX': mx[1], 'TA': ref and ref['TA'],
            'RN15': ref and ref['RN15'], 'HM': ref and ref['HM'], 'N1MIN': len(g)})
    return out

def fetch_day(stn, day, key, tries=3):
    q = urllib.parse.urlencode({'tm1': day.strftime('%Y%m%d') + '0001',
        'tm2': (day + dt.timedelta(days=1)).strftime('%Y%m%d') + '0000',
        'stn': stn, 'disp': 1, 'help': 0, 'authKey': key})
    for i in range(tries):
        try:
            with urllib.request.urlopen(URL + '?' + q, timeout=60) as f:
                b = f.read()
            for enc in ('utf-8', 'euc-kr'):
                try:
                    return b.decode(enc)
                except UnicodeDecodeError:
                    pass
            return b.decode('utf-8', 'replace')
        except Exception as e:
            print(f'  retry {stn} {day:%Y-%m-%d}: {type(e).__name__}', flush=True)
            time.sleep(5 * (i + 1))
    return ''

def main(stations, year, start, end, outdir):
    key = os.environ['KMA_APIHUB_KEY']
    s = max(dt.date(year, 1, 1), start); e = min(dt.date(year, 12, 31), end)
    os.makedirs(outdir, exist_ok=True)
    for stn in stations:
        path = os.path.join(outdir, f'aws10_{stn}_{year}.csv.gz')
        recs, day, nodata = [], s, 0
        while day <= e:
            rows = parse(fetch_day(stn, day, key))
            if not rows: nodata += 1
            recs += to_10min(rows); day += dt.timedelta(days=1)
            time.sleep(0.2)
        with gzip.open(path, 'wt', newline='') as f:
            w = csv.DictWriter(f, fieldnames=OUT); w.writeheader(); w.writerows(recs)
        print(f'{stn} {year}: {len(recs)} rows, days without data {nodata}', flush=True)

if __name__ == '__main__':
    st = [x for x in sys.argv[1].split(',') if x]
    main(st, int(sys.argv[2]), dt.date.fromisoformat(sys.argv[3]), dt.date.fromisoformat(sys.argv[4]), sys.argv[5])
