"""웹 지도용 벡터 기본도(basemap.json) 생성.

자료: Natural Earth 10m (public domain) — 한강 중심선·주요도로·도시지역
      southkorea/southkorea-maps kostat 2018 — 시군구·읍면동 경계
사용: python tools/build_basemap.py <natural_earth_dir> <southkorea-maps_dir> web/basemap.json
"""
import json, math, sys
from shapely.geometry import shape, box, mapping, LineString, Polygon, MultiPolygon
from shapely.ops import unary_union

NE, SK, OUT = sys.argv[1:4]
LON0, LAT0 = 126.7456, 37.6688          # 원점: 킨텍스 V5 (근사)
KX, KY = 111.32 * math.cos(math.radians(LAT0)), 110.95
BB = box(126.60, 37.47, 127.00, 37.76)

def to_km(lon, lat): return ((lon - LON0) * KX, (lat - LAT0) * KY)
def to_ll(x, y): return (LON0 + x / KX, LAT0 + y / KY)
def rnd(geom, nd=5):
    g = mapping(geom)
    def r(c):
        return [round(c[0], nd), round(c[1], nd)] if isinstance(c[0], (int, float)) else [r(i) for i in c]
    g['coordinates'] = r(g['coordinates']); return g
def chaikin(pts, n=3):
    for _ in range(n):
        out = [pts[0]]
        for a, b in zip(pts[:-1], pts[1:]):
            out += [(0.75*a[0]+0.25*b[0], 0.75*a[1]+0.25*b[1]), (0.25*a[0]+0.75*b[0], 0.25*a[1]+0.75*b[1])]
        out.append(pts[-1]); pts = out
    return pts

# 한강: NE 중심선 + 하류(일산~오두산) 근사 연장
han = None
for ft in json.load(open(f'{NE}/ne_10m_rivers_lake_centerlines.geojson'))['features']:
    if ft['properties'].get('name') == 'Han':
        g = ft['geometry']; ls = g['coordinates'] if g['type'] == 'MultiLineString' else [g['coordinates']]
        han = [p for l in ls for p in l if 126.4 < p[0] < 127.2]
han = [tuple(p) for p in han] + [(126.735, 37.650), (126.715, 37.660), (126.695, 37.672), (126.681, 37.700),
                                  (126.676, 37.740), (126.676, 37.770)]   # 일산대교 하류~오두산: 근사
han.reverse()                                        # 하류(서) → 상류(동)
han = [p for p in han if p[0] < 127.03]
han_km = chaikin([to_km(*p) for p in han])
def width(x):                                       # 하류일수록 넓게 (km, 근사)
    return 1.15 if x < 4 else 1.15 - 0.35 * min(1, (x - 4) / 12)
parts = [LineString(han_km[i:i+2]).buffer(width(han_km[i][0]) / 2, cap_style=1) for i in range(len(han_km) - 1)]
river_km = unary_union(parts).simplify(0.03)
river = Polygon([to_ll(*c) for c in river_km.exterior.coords]).intersection(BB)

layers = {'river': rnd(river), 'river_km': [[round(x, 3), round(y, 3)] for x, y in han_km[::3]]}

def clip(path, keep, tol):
    out = []
    for ft in json.load(open(path))['features']:
        g = shape(ft['geometry'])
        if not g.intersects(BB): continue
        g = g.intersection(BB).simplify(tol, preserve_topology=True)
        if g.is_empty: continue
        out.append({'type': 'Feature', 'properties': {k: ft['properties'].get(k) for k in keep}, 'geometry': rnd(g)})
    return {'type': 'FeatureCollection', 'features': out}

K = f'{SK}/kostat/2018/json'
layers['sigungu'] = clip(f'{K}/skorea-municipalities-2018-geo.json', ['name'], 0.0006)
layers['dong'] = clip(f'{K}/skorea-submunicipalities-2018-geo.json', ['name'], 0.0006)
layers['roads'] = clip(f'{NE}/ne_10m_roads.geojson', ['type'], 0.0005)
layers['urban'] = clip(f'{NE}/ne_10m_urban_areas.geojson', [], 0.001)
for fc in ('sigungu',):
    for ft in layers[fc]['features']:
        c = shape(ft['geometry']).representative_point(); ft['properties']['c'] = [round(c.x, 4), round(c.y, 4)]
json.dump(layers, open(OUT, 'w'), ensure_ascii=False, separators=(',', ':'))
print({k: (len(v['features']) if isinstance(v, dict) and 'features' in v else '') for k, v in layers.items()})
