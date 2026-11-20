"""K-UAM 회랑 기상 콘솔 준실시간 수집기 (GitHub Actions에서 10분마다 실행)

- 기상청 API허브: 콘솔 영역 안 AWS 전 지점 1분 자료, KIM 국지(KIML 1.5 km) U10·V10·T2 0~48시간
- Open-Meteo: ECMWF IFS 0.25° 지상·등압면(925·850 hPa) 시계열 48시간
- 결과: out/latest.json (rt-data 브랜치에 올라가고, 웹 콘솔이 읽음)

환경변수 HUB에 API허브 authKey(또는 authKey가 들어간 주소 전체), PREV에 직전 latest.json 경로(선택).
키는 어떤 출력에도 쓰지 않는다.
"""
import datetime as D, json, os, re, sys, time, urllib.request, urllib.parse

KST = D.timezone(D.timedelta(hours=9))
BBOX = dict(lon0=126.66, lon1=127.01, lat0=37.47, lat1=37.73)
V5 = (126.7456, 37.6688)
HF = list(range(0, 13)) + list(range(15, 49, 3))
STRIDE = 2   # KIML 격자 솎기 (1.5 km → 약 3 km)
LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a); LOG.append(s); print(s, flush=True)


def key():
    k = os.environ.get('HUB', '')
    m = re.search(r'authKey=([A-Za-z0-9_\-]+)', k)
    return m.group(1) if m else re.sub(r'\s', '', k)


def get(url, timeout=90, tries=2):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa
            err = type(e).__name__
            time.sleep(2)
    raise RuntimeError(err)


def hub(path, **q):
    q['authKey'] = KEY
    return get('https://apihub.kma.go.kr/api/' + path + '?' + urllib.parse.urlencode(q, safe=','))


def num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f <= -50 else f


# ---------- AWS ----------
def stations():
    txt = hub('typ01/url/stn_inf.php', inf='AWS').decode('euc-kr', 'replace')
    out = {}
    for ln in txt.splitlines():
        if not ln.strip() or ln.startswith('#'):
            continue
        t = ln.split()
        try:
            stn, lon, lat = int(t[0]), float(t[1]), float(t[2])
        except (ValueError, IndexError):
            continue
        if not (BBOX['lon0'] <= lon <= BBOX['lon1'] and BBOX['lat0'] <= lat <= BBOX['lat1']):
            continue
        ko = [x for x in t if re.search('[가-힣]', x)]
        rec = {'stn': stn, 'name': ko[0] if ko else str(stn), 'lon': lon, 'lat': lat}
        try:
            rec['ht'] = float(t[4]); rec['hwd'] = float(t[5])
        except (ValueError, IndexError):
            pass
        out[stn] = rec   # 같은 번호가 여러 줄이면 마지막(최신) 정보
    return out


def aws_now(stn):
    now = D.datetime.now(KST) - D.timedelta(minutes=2)
    for back in (0, 1, 2, 5, 10):
        tm = (now - D.timedelta(minutes=back)).strftime('%Y%m%d%H%M')
        txt = hub('typ01/cgi-bin/url/nph-aws2_min', tm2=tm, stn=0, disp=1, help=0).decode('euc-kr', 'replace')
        rows = []
        for ln in txt.splitlines():
            if not ln.strip() or ln.startswith('#'):
                continue
            t = [x.strip() for x in ln.split(',')]
            if len(t) < 18:
                continue
            try:
                s = int(t[1])
            except ValueError:
                continue
            if s not in stn:
                continue
            r = dict(stn[s])
            r.update(tm=t[0], wd1=num(t[2]), ws1=num(t[3]), wdg=num(t[4]), gust=num(t[5]), wd=num(t[6]), ws=num(t[7]),
                     ta=num(t[8]), rn15=num(t[10]), rn60=num(t[11]), rnday=num(t[13]), hm=num(t[14]), pa=num(t[15]), td=num(t[17]))
            if r['wd'] is not None and r['ws'] is not None and (r['ws'] < 0 or r['wd'] > 360):
                r['wd'] = r['ws'] = None
            rows.append(r)
        if len([r for r in rows if r['ta'] is not None]) >= 5:
            return tm, rows
    return None, []


# ---------- KIM 국지 ----------
def kiml_csv(var, tmfc, hf):
    b = hub('typ06/url/nwp_nc_vars_down.php', group='KIML', data='U', vars=var, tmfc=tmfc, hf=hf,
            lon1=BBOX['lon0'], lat1=BBOX['lat0'], lon2=BBOX['lon1'], lat2=BBOX['lat1'], dataType='CSV').decode('euc-kr', 'replace')
    if not b.startswith('ix,iy'):
        return None
    pts = []
    for ln in b.splitlines()[1:]:
        t = ln.split(',')
        if len(t) == 5:
            pts.append((int(t[0]), int(t[1]), float(t[2]), float(t[3]), float(t[4])))
    return pts


def kiml_latest(prev):
    now = D.datetime.now(D.timezone.utc)
    for back in range(3, 30, 6):
        t = now - D.timedelta(hours=back)
        tmfc = t.strftime('%Y%m%d') + '%02d' % (t.hour // 6 * 6)
        if prev and prev.get('tmfc') == tmfc and len(prev.get('hf', [])) == len(HF):
            log('KIML', tmfc, '직전 자료 재사용'); return prev
        p0 = kiml_csv('T2', tmfc, 0)
        if p0:
            break
    else:
        log('KIML 자료 없음'); return prev
    ix0, iy0 = min(p[0] for p in p0), min(p[1] for p in p0)
    keep = [i for i, p in enumerate(p0) if (p[0] - ix0) % STRIDE == 0 and (p[1] - iy0) % STRIDE == 0]
    out = {'tmfc': tmfc, 'hf': [], 'lon': [round(p0[i][2], 4) for i in keep], 'lat': [round(p0[i][3], 4) for i in keep],
           'u10': [], 'v10': [], 't2': []}
    for hf in HF:
        rows = {}
        for var in ('U10', 'V10', 'T2'):
            p = p0 if (var == 'T2' and hf == 0) else kiml_csv(var, tmfc, hf)
            if not p or len(p) != len(p0):
                rows = None; break
            rows[var] = [p[i][4] for i in keep]
        if rows is None:
            log('KIML', tmfc, 'hf', hf, '없음 — 여기까지'); break
        out['hf'].append(hf)
        out['u10'].append([round(v, 2) for v in rows['U10']])
        out['v10'].append([round(v, 2) for v in rows['V10']])
        out['t2'].append([round(v - 273.15, 2) for v in rows['T2']])
    log('KIML', tmfc, '시간', len(out['hf']), '격자점', len(keep))
    return out if out['hf'] else prev


# ---------- ECMWF (Open-Meteo) ----------
def ecmwf():
    hv = ['temperature_2m', 'dew_point_2m', 'precipitation', 'snowfall', 'cloud_cover', 'cloud_cover_low',
          'wind_speed_10m', 'wind_direction_10m', 'wind_gusts_10m', 'cape',
          'wind_speed_925hPa', 'wind_direction_925hPa', 'temperature_925hPa', 'relative_humidity_925hPa',
          'wind_speed_850hPa', 'wind_direction_850hPa', 'temperature_850hPa']
    q = dict(latitude=37.60, longitude=126.84, hourly=','.join(hv), wind_speed_unit='ms', models='ecmwf_ifs025',
             timezone='Asia/Seoul', past_hours=3, forecast_hours=50)
    d = json.loads(get('https://api.open-meteo.com/v1/forecast?' + urllib.parse.urlencode(q)))
    h = d['hourly']
    out = {'model': 'ECMWF IFS 0.25° (Open-Meteo)', 'time': h['time']}
    short = {'temperature_2m': 't2', 'dew_point_2m': 'td2', 'precipitation': 'pr', 'snowfall': 'sn', 'cloud_cover': 'cc',
             'cloud_cover_low': 'ccl', 'wind_speed_10m': 'ws10', 'wind_direction_10m': 'wd10', 'wind_gusts_10m': 'g10', 'cape': 'cape',
             'wind_speed_925hPa': 'ws925', 'wind_direction_925hPa': 'wd925', 'temperature_925hPa': 't925', 'relative_humidity_925hPa': 'rh925',
             'wind_speed_850hPa': 'ws850', 'wind_direction_850hPa': 'wd850', 'temperature_850hPa': 't850'}
    for k, s in short.items():
        out[s] = h.get(k)
    log('ECMWF', len(out['time']), '시간', 'ws925 첫값', (out['ws925'] or [None])[0])
    return out


def main():
    prev = {}
    if os.environ.get('PREV') and os.path.exists(os.environ['PREV']):
        try:
            prev = json.load(open(os.environ['PREV']))
        except Exception:
            prev = {}
    res = {'v': 1, 'made': D.datetime.now(KST).strftime('%Y-%m-%dT%H:%M+09:00'), 'bbox': BBOX}
    try:
        stn = stations(); log('영역 AWS 지점정보', len(stn))
        tm, rows = aws_now(stn)
        res['obs'] = {'tm': tm, 'src': '기상청 API허브 AWS 1분 자료', 'stations': rows}
        log('AWS', tm, '관측', len(rows), '바람', sum(1 for r in rows if r['ws'] is not None))
    except Exception as e:
        log('AWS 실패', type(e).__name__); res['obs'] = prev.get('obs')
    try:
        res['kiml'] = kiml_latest(prev.get('kiml'))
    except Exception as e:
        log('KIML 실패', type(e).__name__); res['kiml'] = prev.get('kiml')
    try:
        res['ecmwf'] = ecmwf()
    except Exception as e:
        log('ECMWF 실패', type(e).__name__); res['ecmwf'] = prev.get('ecmwf')
    res['log'] = LOG
    os.makedirs('out', exist_ok=True)
    json.dump(res, open('out/latest.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    log('저장', os.path.getsize('out/latest.json'), 'bytes')


KEY = key()
if __name__ == '__main__':
    if not KEY:
        sys.exit('HUB 키 없음')
    main()
