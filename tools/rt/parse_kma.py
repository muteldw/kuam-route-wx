"""기상청 API허브 텍스트 응답 파서 (stn_inf, nph-aws2_min) — 영역 내 AWS 목록·최신 관측을 JSON으로."""
import json, re, sys

BBOX = dict(lon0=126.66, lon1=127.01, lat0=37.47, lat1=37.73)   # 콘솔 영역(EXT) + 약 2 km 여유


def header_cols(lines):
    cols = None
    for ln in lines:
        if ln.startswith('#'):
            toks = re.findall(r'[A-Za-z][A-Za-z0-9_\-]*', ln)
            if any(t.upper() == 'STN' for t in toks) and any(t.upper().startswith('LON') for t in toks):
                cols = toks
    return cols


def stations(path):
    lines = open(path, encoding='euc-kr', errors='replace').read().splitlines()
    cols = header_cols(lines) or []
    up = [c.upper() for c in cols]
    out = []
    for ln in lines:
        if not ln.strip() or ln.startswith('#'):
            continue
        t = ln.split()
        try:
            rec = {'stn': int(t[0]), 'lon': float(t[1]), 'lat': float(t[2])}
        except (ValueError, IndexError):
            continue
        for k, v in zip(up, t):
            if k.startswith('HT'):
                try: rec[k.lower()] = float(v)
                except ValueError: pass
        ko = [x for x in t if re.search('[가-힣]', x)]
        rec['name'] = ko[0] if ko else str(rec['stn'])
        if BBOX['lon0'] <= rec['lon'] <= BBOX['lon1'] and BBOX['lat0'] <= rec['lat'] <= BBOX['lat1']:
            out.append(rec)
    return {'columns': cols, 'stations': out}


def aws_min(path, keep):
    lines = open(path, encoding='euc-kr', errors='replace').read().splitlines()
    cols = None
    for ln in lines:
        if ln.startswith('#') and 'WS' in ln and 'STN' in ln:
            cols = [c.strip() for c in ln.lstrip('#').replace(',', ' ').split()]
    out = []
    for ln in lines:
        if not ln.strip() or ln.startswith('#'):
            continue
        t = [x.strip() for x in ln.replace(',', ' ').split()]
        if len(t) < 5:
            continue
        try:
            stn = int(t[1])
        except ValueError:
            continue
        if keep and stn not in keep:
            continue
        rec = {'tm': t[0], 'stn': stn}
        if cols:
            for k, v in zip(cols, t):
                if k in ('YYMMDDHHMI', 'TM', 'STN'):
                    continue
                try:
                    f = float(v); rec[k] = None if f <= -50 else f
                except ValueError:
                    pass
        else:
            rec['raw'] = t
        out.append(rec)
    return {'columns': cols, 'obs': out}


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'stn':
        print(json.dumps(stations(sys.argv[2]), ensure_ascii=False, indent=1))
    else:
        keep = set(s['stn'] for s in json.load(open(sys.argv[3]))['stations']) if len(sys.argv) > 3 else set()
        print(json.dumps(aws_min(sys.argv[2], keep), ensure_ascii=False, indent=1))
