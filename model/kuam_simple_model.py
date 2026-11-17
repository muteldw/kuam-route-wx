"""K-UAM 경로 기상 위험 · 배터리 · 항로 추천 단순 진단모델 (참조 구현).

작성자: 이동원 (Dong-Won Yi) <yidongwonyi@gmail.com> · 버전: 1.2.0

웹 화면(web/index.html)의 JavaScript 모델과 같은 식을 쓴다. 기체 제원은 data/fleet.json.
실측 연결 시 met() 안의 각 항을 LDAPS/ECMWF/관측 자료로 교체한다.

사용 예:
    python kuam_simple_model.py --from "킨텍스 V5" --to 여의도 --aircraft joby --adult 2 --child 1 --scenario conv
    python kuam_simple_model.py --list-aircraft
    python kuam_simple_model.py --from 마곡 --to "상암 DMC" --aircraft ehang --adult 1 --json > plan.json
"""
from __future__ import annotations

import argparse
import heapq
import json
import math
from dataclasses import dataclass, replace
from pathlib import Path

G, KAPPA = 9.81, 0.4
ROOT = Path(__file__).resolve().parent.parent

# 좌표: km, 원점 킨텍스 V5 (126.7456E, 37.6688N), x 동쪽 / y 북쪽 (근사)
RIVER = [(-6.133, 11.228), (-6.133, 10.916), (-6.133, 10.136), (-6.126, 8.922), (-6.064, 7.518), (-5.94, 5.958), (-5.738, 4.356), (-5.424, 2.942), (-4.998, 1.714), (-4.486, 0.716), (-3.9, -0.033), (-3.247, -0.557), (-2.586, -1.025), (-1.925, -1.462), (-1.263, -1.889), (-0.594, -2.349), (0.08, -2.842), (0.731, -3.341), (1.344, -3.831), (2.03, -4.411), (3.558, -5.762), (6.04, -7.982), (8.824, -10.467), (10.613, -12.014), (11.404, -12.621), (11.983, -12.906), (12.741, -13.178), (13.673, -13.457), (14.754, -13.894), (15.98, -14.512), (17.156, -15.175), (17.89, -15.612), (18.182, -15.823), (18.33, -15.998), (18.48, -16.233), (18.635, -16.52), (18.807, -16.812), (18.999, -17.1), (19.214, -17.364), (19.46, -17.563), (19.737, -17.697), (20.023, -17.752), (20.305, -17.725), (20.582, -17.617), (20.83, -17.464), (21.048, -17.268), (21.242, -17.04), (21.432, -16.796), (21.616, -16.537), (21.891, -16.204), (22.302, -15.769), (22.822, -15.257), (23.254, -14.854), (23.569, -14.586), (23.808, -14.421), (24.054, -14.293), (24.307, -14.202), (24.498, -14.144), (24.594, -14.115)]
AIRPORT = (3.97, -12.26)
YEO = (15.78, -15.9)
EXT = dict(x0=-5, x1=21, y0=-19.5, y1=4)
LANDMARKS = {"킨텍스 V5": (0, 0), "일산대교": (1.45, -3.94), "행주대교": (5.67, -7.68), "방화대교": (7.79, -9.57),
             "가양대교": (10.26, -11.78), "성산대교": (12.64, -13.15), "여의도": (15.76, -15.89), "김포공항": (3.97, -12.26),
             "상암 DMC": (12.63, -10.07), "마곡": (7.87, -12.06), "용산": (19.3, -15.4), "고양 창릉": (10.5, -3.2),
             "일산 호수공원": (1.71, -1.2), "화정": (7.4, -3.7), "계양": (-1.4, -14.3), "부천": (1.8, -18.3),
             "은평": (16.2, -7.6), "행주산성": (7.43, -7.86)}
VERTIPORTS = ["킨텍스 V5", "여의도", "김포공항", "상암 DMC", "마곡", "용산", "고양 창릉", "일산 호수공원"]
MASS = dict(pilot=85, adult=80, child=35)


def clamp(v, a, b):
    return max(a, min(b, v))


def nearest_name(x, y, max_km=2.5):
    n, d = min(((k, math.hypot(v[0] - x, v[1] - y)) for k, v in LANDMARKS.items()), key=lambda kv: kv[1])
    return n if d <= max_km else None


@dataclass
class Scenario:
    label: str
    issueH: int
    U: float
    dir: float
    T: float
    dep: float
    pr: float = 0.0
    llj: float = 0.0
    conv: bool = False
    front: dict | None = None


SCENARIOS = {
    "calm": Scenario("가을 맑음 · 약풍", 14, 6, 270, 17, 7),
    "nw": Scenario("겨울 강한 북서풍", 9, 17, 315, 1, 9, front={"h": 10, "dU": 5, "dDir": 15}),
    "llj": Scenario("야간 저층제트", 21, 9, 235, 13, 4, llj=9),
    "conv": Scenario("여름 대류성 뇌우·돌풍전선", 15, 7, 250, 28, 4, conv=True),
    "rain": Scenario("장마 비 · 남서풍", 11, 11, 225, 23, 2.2, pr=9),
    "snow": Scenario("겨울 눈 · 약한 북서풍", 7, 8, 320, -2, 3.2, pr=1.0),
    "fog": Scenario("새벽 한강 복사안개", 4, 2, 80, 9, 0.6),
}
CONV = dict(t0=-0.8, t1=3.5, x0=-1, y0=-14, vx=5.5, vy=3, R=3, W=7, O=13, P=45)


def _poly_dist(x, y, P):
    best, side = 1e9, 1
    for (ax, ay), (bx, by) in zip(P[:-1], P[1:]):
        dx, dy = bx - ax, by - ay
        t = clamp(((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy), 0.0, 1.0)
        d = math.hypot(x - (ax + t * dx), y - (ay + t * dy))
        if d < best:
            best, side = d, (-1 if dx * (y - ay) - dy * (x - ax) < 0 else 1)
    return best, side


def _smooth(e0, e1, v):
    t = clamp((v - e0) / (e1 - e0), 0.0, 1.0)
    return t * t * (3 - 2 * t)


_LZ = {}


def logz0(x, y):
    key = (round(x * 20), round(y * 20))
    if key in _LZ:
        return _LZ[key]
    d, side = _poly_dist(x, y, RIVER)
    land = -0.3
    if side < 0:
        land = -0.3 + (-0.7) * (1 - _smooth(7, 10, x))
        if x > 9:
            land = -0.2
    elif x > 9:
        land = -0.1
    ra = math.hypot(x - AIRPORT[0], y - AIRPORT[1])
    land += (-1.55 - land) * math.exp(-(ra / 1.7) ** 4)
    land += 0.45 * math.exp(-(math.hypot(x - YEO[0], y - YEO[1]) / 1.1) ** 2)
    wr = 1 / (1 + math.exp((d - 0.55) / 0.12))
    _LZ[key] = (land * (1 - wr) - 3.7 * wr, d)
    return _LZ[key]


def met(x, y, z, t, sc: Scenario):
    """지점 (x, y km), 고도 z (m), 발표 후 t 시간의 진단 기상값."""
    h = (sc.issueH + t) % 24
    S = math.sin(math.pi * (h - 6) / 12) if 6 < h < 18 else -0.8
    Us = sc.U * (1 + 0.12 * math.sin(2 * math.pi * t / 30 + 0.7))
    wdir = sc.dir + 8 * math.sin(2 * math.pi * t / 20)
    if sc.front and t > sc.front["h"]:
        k = 1 - math.exp(-(t - sc.front["h"]) / 2)
        Us += sc.front["dU"] * k
        wdir += sc.front["dDir"] * k
    Seff = S * math.exp(-Us / 12)
    zz = max(z, 10.0)
    lz, dRiver = logz0(x, y)
    lzE = lz + (-0.6 - lz) * min(1, zz / 400)
    z0e = 10 ** lzE
    a = clamp((0.096 * lzE + 0.016 * lzE * lzE + 0.24) * (1 - 0.45 * Seff), 0.06, 0.5)
    spd = Us * (zz / 600) ** a if zz <= 600 else Us * (1 + 0.0003 * (zz - 600))
    jet = dJ = 0.0
    if S < 0 and sc.llj > 0:
        jet = sc.llj * min(1, -S / 0.8) * math.exp(-((zz - 320) / 110) ** 2)
        dJ = jet * (-2 * (zz - 320) / 110 ** 2)
    r5 = math.hypot(x, y)
    A = 1.58 + 0.45 * math.cos(math.radians(wdir))
    amp = 1 + (A - 1) * math.exp(-(r5 / 1.2) ** 2) * math.exp(-zz / 80)
    bw = math.exp(-(math.hypot(x - YEO[0], y - YEO[1]) / 0.9) ** 2) * math.exp(-zz / 150)
    base = spd * (1 - 0.25 * bw) + jet
    spd = spd * amp * (1 - 0.25 * bw) + jet
    veer = (25 if Seff < 0 else 10) * (1 - min(zz, 600) / 600)
    dr = math.radians(wdir - veer)
    u, v = -spd * math.sin(dr), -spd * math.cos(dr)
    w = Pc = epsC = tiC = ts = 0.0
    sd = math.sin(2 * math.pi * (h - 9) / 24)
    T0 = sc.T + 5 * sd
    if sc.conv and CONV["t0"] <= t <= CONV["t1"]:
        tau = t - CONV["t0"]
        life = math.sin(math.pi * tau / (CONV["t1"] - CONV["t0"]))
        cx, cy = CONV["x0"] + CONV["vx"] * tau, CONV["y0"] + CONV["vy"] * tau
        r, R = max(0.01, math.hypot(x - cx, y - cy)), CONV["R"]
        O = CONV["O"] * life * (r / R) * math.exp(1 - r / R) * math.exp(-zz / 450)
        u += O * (x - cx) / r
        v += O * (y - cy) / r
        w = -CONV["W"] * life * math.exp(-(r / (0.6 * R)) ** 2) * min(1, zz / 250)
        Pc = CONV["P"] * life * math.exp(-(r / (0.7 * R)) ** 2)
        epsC = 0.04 * life * math.exp(-(r / R) ** 2)
        tiC = 0.15 * life * math.exp(-(r / R) ** 2)
        if T0 > 15:
            ts = life * math.exp(-(r / (1.3 * R)) ** 2)
    sH = math.hypot(u, v)
    lnz = math.log(zz / z0e)
    TI = (1 / lnz) * (1 + 0.5 * max(Seff, 0)) + 0.10 * bw \
        + 0.08 * math.exp(-(r5 / 0.8) ** 2) * math.exp(-zz / 100) + tiC
    us = KAPPA * (base + (max(0, sH - spd) if sc.conv else 0)) / lnz
    eps = us ** 3 / (KAPPA * max(zz, 50)) * (1 + bw) + 0.0015 * max(Seff, 0) * max(0, 1 - zz / 1500) \
        + 1.5 * dJ * dJ + epsC
    dep = max(0.15, sc.dep * (1 + 0.75 * sd))
    es = lambda T: math.exp(17.625 * T / (243.04 + T))
    RH0 = 100 * es(T0 - dep) / es(T0)
    inv = 3 * min(1, -S / 0.8) * min(zz, 150) / 150 if S < 0 else 0
    Tz = T0 - 0.0065 * zz + inv
    Tdz = (T0 - dep) - 0.0018 * zz
    RHz = min(100, 100 * es(min(Tdz, Tz)) / es(Tz))
    ceil = max(30, 125 * dep) if (dep < 2.5 or sc.pr > 0) else 3000
    in_cloud = zz > ceil
    cc = 1 if in_cloud else (min(0.9, 1 - math.sqrt((100 - RHz) / 20)) if RHz >= 80 else 0)
    f = clamp((RH0 - 95) / 4, 0, 1) * clamp((5 - 0.35 * Us) / 3, 0, 1) * (1 if S < 0.25 else max(0, 1 - (S - 0.25) * 3))
    f = min(1, f * (1 + 0.4 * math.exp(-dRiver ** 2)))
    fog_top = 60 + 200 * f
    prS = sc.pr * max(0.15, 0.6 + 0.4 * math.sin(0.28 * x + 0.17 * y - 0.9 * t)
                      + 0.25 * math.sin(0.11 * x - 0.23 * y + 0.5 * t)) if sc.pr > 0 else 0
    P = prS + Pc
    snowy = Tz <= 0.5
    rain, snow = (0, P) if snowy else (P, 0)
    vis = 10000 * math.exp(-4.6 * f) if zz < fog_top else 10000
    if in_cloud:
        vis = min(vis, 300)
    vis = min(vis, 12000 / (1 + 0.3 * rain))
    if snow > 0.05:
        vis = min(vis, 2500 / (0.3 + snow))
    vis = max(80, vis)
    icing = -12 <= Tz <= 0 and (in_cloud or (f > 0.3 and zz < fog_top) or rain > 0.2 or (snow > 0.3 and Tz > -3 and RHz > 95))
    return dict(u=u, v=v, w=w, spd=sH, gust=sH * (1 + 2.2 * TI), edr=eps ** (1 / 3), vis=vis, ceil=ceil, T=Tz, rh=RHz,
                cc=cc, in_cloud=in_cloud, icing=icing, P=P, rain=rain, snow=snow, ts=ts, z0=z0e)


# ---------- risk ----------
LV = ["정상", "주의", "경고", "불가"]
FTH, HSH, VSH = (0.07, 0.105), (5, 7.5), (1.5, 2.0)


def ri_of(m):
    if m >= 0.2:
        return 40 * clamp((1 - m) / 0.8, 0, 1)
    if m >= 0:
        return 40 + 20 * (0.2 - m) / 0.2
    if m >= -0.25:
        return 60 + 20 * (-m / 0.25)
    return min(100, 80 + 20 * min(1, (-0.25 - m) / 0.75))


def cat_ri(ri):
    return 0 if ri < 40 else 1 if ri < 60 else 2 if ri < 80 else 3


def margin_items(m, z, lim, low):
    it = [("이착륙 풍속", 1 - m["spd"] / lim["xw"]) if low else ("풍속", 1 - m["spd"] / lim["wmax"]),
          ("돌풍", 1 - m["gust"] / lim["gust"]), ("난류", 1 - m["edr"] / lim["edr"]), ("시정", 1 - lim["vis"] / m["vis"]),
          ("운고", 1 - (z + 60) / m["ceil"]), ("강우", 1 - m["rain"] / lim["rain"]), ("강설", 1 - m["snow"] / lim["snow"]),
          ("뇌우", 1 - m["ts"] / 0.25)]
    if m["icing"] and not lim["icing"]:
        it.append(("착빙", -0.5))
    return [(k, clamp(v, -1, 1)) for k, v in it]


# ---------- fleet & energy ----------
def load_fleet():
    fl = json.loads((ROOT / "data" / "fleet.json").read_text(encoding="utf-8"))["aircraft"]
    return {a["id"]: a for a in fl}


def rho_at(T):
    return 1.225 * 288.15 / (273.15 + T)


def hover_power(a, mass, T):
    W, DL = mass * G, a["DL"] * mass / a["mtow"]
    return W * math.sqrt(DL / (2 * rho_at(T))) / a["FM"] / 0.9


def level_power(a, mass, V, T):
    Ph, W, vc = hover_power(a, mass, T), mass * G, a["vc"] / 3.6
    if a["cfg"] == "multi":
        return Ph * (0.75 + 0.35 * (V / vc) ** 2)
    Vt, wing = 0.45 * vc, W * V / (a["LD"] * 0.8 * 0.9)
    if V < Vt:
        return Ph * (1 - 0.55 * V / Vt) + wing * (V / Vt)
    return wing * (0.75 + 0.25 * (V / vc) ** 2) + 0.04 * Ph


def aux_power(T):
    return 3000 + (4000 if T < 5 else 0) + (3000 if T > 28 else 0)


def batt_usable(T):
    return clamp(1 - 0.01 * max(0, 10 - T), 0.75, 1)


@dataclass
class Flight:
    air: dict
    lim: dict
    vto: float
    vc: float
    vld: float
    alt: float
    climb: float
    desc: float
    mass: float
    corr: float = 1.0


def make_flight(air, adult=0, child=0, bag=10, **over):
    d = air["def"]
    people = air["pilot"] * MASS["pilot"] + adult * MASS["adult"] + child * MASS["child"]
    f = Flight(air=air, lim=dict(air["lim"]), vto=d["vto"], vc=d["vc"], vld=d["vld"], alt=d["alt"], climb=d["climb"],
               desc=d["desc"], mass=air["empty"] + people + (adult + child) * bag)
    return replace(f, **{k: v for k, v in over.items() if v is not None})


# ---------- route ----------
DS = 0.1


def build_route(poly):
    xs, ys, tx, ty, s, acc = [], [], [], [], [], 0.0
    for (ax, ay), (bx, by) in zip(poly[:-1], poly[1:]):
        L = math.hypot(bx - ax, by - ay)
        if L < 1e-6:
            continue
        n = max(1, round(L / DS))
        for k in range(n):
            f = k / n
            xs.append(ax + (bx - ax) * f); ys.append(ay + (by - ay) * f)
            tx.append((bx - ax) / L); ty.append((by - ay) / L); s.append(acc + L * f)
        acc += L
    xs.append(poly[-1][0]); ys.append(poly[-1][1]); tx.append(tx[-1] if tx else 1); ty.append(ty[-1] if ty else 0); s.append(acc)
    return dict(xs=xs, ys=ys, tx=tx, ty=ty, s=s, n=len(xs), L=acc, poly=poly)


def evaluate(t0, fl: Flight, sc, R, fast=False):
    n, L, H0, air = R["n"], R["L"], 30.0, fl.air
    top = max(fl.alt, H0 + 20)
    dc, dd = (top - H0) / fl.climb * fl.vto / 3.6 / 1000, (top - H0) / fl.desc * fl.vld / 3.6 / 1000
    if dc + dd > 0.9 * L:
        k = 0.9 * L / (dc + dd); dc *= k; dd *= k; top = H0 + (top - H0) * k
    alt, va, ph = [], [], []
    for s in R["s"]:
        if s < dc:
            alt.append(H0 + (top - H0) * s / dc); va.append(fl.vto / 3.6); ph.append("climb")
        elif s > L - dd:
            alt.append(H0 + (top - H0) * (L - s) / dd); va.append(fl.vld / 3.6); ph.append("desc")
        else:
            alt.append(top); va.append(fl.vc / 3.6); ph.append("cruise")
    tUp, tDn = H0 / fl.climb, H0 / fl.desc

    def timing(Wx):
        gs = [max(4, va[i] + (Wx[i] if Wx else 0)) for i in range(n)]
        cum = [tUp]
        for i in range(1, n):
            cum.append(cum[-1] + (R["s"][i] - R["s"][i - 1]) * 1000 / ((gs[i] + gs[i - 1]) / 2))
        return gs, cum

    gs, cum = timing(None)
    pts, Wx, W = [], [], []
    for i in range(n):
        x, y, z, eta = R["xs"][i], R["ys"][i], alt[i], t0 + cum[i] / 3600
        m = met(x, y, z, eta, sc)
        vsh = 0.0
        if not fast:
            m1, m0 = met(x, y, z + 15, eta, sc), met(x, y, max(2, z - 15), eta, sc)
            vsh = math.hypot(m1["u"] - m0["u"], m1["v"] - m0["v"]) * 30 / (z + 15 - max(2, z - 15))
        Wx.append(m["u"] * R["tx"][i] + m["v"] * R["ty"][i]); W.append(m["w"])
        pts.append(dict(m=m, z=z, vsh=vsh))
    gs, cum = timing(Wx)
    Fr = []
    for i in range(n):
        a, b = max(0, i - 3), min(n - 1, i + 3)
        ds = (R["s"][b] - R["s"][a]) * 1000 or 100
        Fr.append((Wx[b] - Wx[a]) / ds * gs[i] / G - W[i] / max(15, va[i]))
        pts[i]["hsh"] = abs(Wx[b] - Wx[a]) * 0.5 / ((R["s"][b] - R["s"][a]) or 0.5)
    m0, mN = pts[0]["m"], pts[-1]["m"]
    eUp = ((hover_power(air, fl.mass, m0["T"]) + fl.mass * G * fl.climb / 0.72) * (1 + 0.8 * m0["edr"]) + aux_power(m0["T"])) * tUp / 3.6e6
    eDn = (0.75 * hover_power(air, fl.mass, mN["T"]) * (1 + 0.8 * mN["edr"]) + aux_power(mN["T"])) * tDn / 3.6e6
    E, phase_e, ri_sum, ri_max, worst, kw_prev = eUp + eDn, dict(vert=eUp + eDn, climb=0, cruise=0, desc=0), 0, 0, (9, 0, ""), None
    for i, p in enumerate(pts):
        win = Fr[max(0, i - 5): min(n, i + 6)]
        F = sum(win) / len(win)
        m = p["m"]
        items = margin_items(m, p["z"], fl.lim, p["z"] < 150 and ph[i] != "cruise")
        items += [("급변풍", clamp(1 - F / FTH[1], -1, 1) if F > 0 else 1), ("수평시어", clamp(1 - p["hsh"] / HSH[1], -1, 1))]
        if not fast:
            items.append(("연직시어", clamp(1 - p["vsh"] / VSH[1], -1, 1)))
        f, mg = min(items, key=lambda kv: kv[1])
        ri = ri_of(mg)
        P = level_power(air, fl.mass, va[i], m["T"])
        if ph[i] == "climb":
            P += fl.mass * G * fl.climb / 0.72
        elif ph[i] == "desc":
            P = max(0.3 * hover_power(air, fl.mass, m["T"]), P - 0.5 * fl.mass * G * fl.desc * 0.9)
        kw = (P * (1 + 0.8 * m["edr"]) + aux_power(m["T"])) / 1000
        if i:
            e = (kw + kw_prev) / 2 * (cum[i] - cum[i - 1]) / 3600
            E += e; phase_e[ph[i]] += e
            ri_sum += (ri + pts[i - 1]["ri"]) / 2 * (R["s"][i] - R["s"][i - 1])
        kw_prev = kw
        p.update(F=F, Wx=Wx[i], ri=ri, c=cat_ri(ri), limiting=f, kw=kw, eta=t0 + cum[i] / 3600)
        ri_max = max(ri_max, ri)
        if mg < worst[0]:
            worst = (mg, i, f)
    T_mean = sum(p["m"]["T"] for p in pts) / n
    usable = air["battKWh"] * batt_usable(min(m0["T"], T_mean))
    used = E * fl.corr
    pct = used / usable * 100
    c, limit = cat_ri(ri_max), worst[2]
    bc = 3 if pct > 100 else 2 if pct > 80 else 1 if pct > 65 else 0
    if bc > c:
        c, limit = bc, "배터리 여유"
    if fl.mass > air["mtow"]:
        c, limit = 3, "중량 초과"
    return dict(pts=pts, phase=ph, top=top, total_s=cum[-1] + tDn, ri_max=ri_max, ri_mean=ri_sum / (L or 1), c=c,
                limit=limit, worst_i=worst[1], energy=dict(kwh=used, raw=E, pct=pct, usable=usable, phase=phase_e,
                                                           kwh_per_km=used / (L or 1)))


# ---------- route planner ----------
PROFILES = [dict(key="fast", name="최단 시간", wR=1.2, river=1.0, lam=0.4),
            dict(key="safe", name="최소 위험", wR=10, river=1.0, lam=3.0),
            dict(key="river", name="한강 회랑", wR=4, river=0.55, lam=1.5)]
STEP = 0.4


def cost_grid(z, t, sc, lim):
    NX, NY = int((EXT["x1"] - EXT["x0"]) / STEP) + 1, int((EXT["y1"] - EXT["y0"]) / STEP) + 1
    ri, river, apt = [], [], []
    for j in range(NY):
        for i in range(NX):
            x, y = EXT["x0"] + i * STEP, EXT["y0"] + j * STEP
            m = met(x, y, z, t, sc)
            ri.append(ri_of(min(v for _, v in margin_items(m, z, lim, False))))
            river.append(logz0(x, y)[1] < 0.6)
            apt.append(math.hypot(x - AIRPORT[0], y - AIRPORT[1]))
    return dict(NX=NX, NY=NY, ri=ri, river=river, apt=apt)


def cell(G_, x, y):
    return clamp(round((x - EXT["x0"]) / STEP), 0, G_["NX"] - 1), clamp(round((y - EXT["y0"]) / STEP), 0, G_["NY"] - 1)


def mult_grid(G_, prof, A, B, avoid):
    apt_ok = math.hypot(A[0] - AIRPORT[0], A[1] - AIRPORT[1]) < 3.2 or math.hypot(B[0] - AIRPORT[0], B[1] - AIRPORT[1]) < 3.2
    pts = [(R["xs"][i], R["ys"][i]) for R in avoid for i in range(0, R["n"], 4)]
    out = []
    for k, r in enumerate(G_["ri"]):
        p = 0 if r < 30 else (r - 30) / 70
        mu = (1 + prof["wR"] * p + (30 if r >= 80 else 0)) * (prof["river"] if G_["river"][k] else 1)
        if not apt_ok and G_["apt"][k] < 2.6:
            mu = math.inf
        if pts and math.isfinite(mu):
            x, y = EXT["x0"] + (k % G_["NX"]) * STEP, EXT["y0"] + (k // G_["NX"]) * STEP
            d = min(abs(q[0] - x) + abs(q[1] - y) for q in pts)
            if d < 1.0:
                mu += 3 * (1 - d)
        out.append(mu)
    return out


def dijkstra(G_, mult, s, e):
    NX, NY = G_["NX"], G_["NY"]
    sk, ek = s[1] * NX + s[0], e[1] * NX + e[0]
    mult = list(mult); mult[sk] = mult[sk] if math.isfinite(mult[sk]) else 1; mult[ek] = mult[ek] if math.isfinite(mult[ek]) else 1
    dist, prev, done = {sk: 0.0}, {}, set()
    hq = [(0.0, sk)]
    while hq:
        d, k = heapq.heappop(hq)
        if k in done:
            continue
        done.add(k)
        if k == ek:
            break
        i, j = k % NX, k // NX
        for di, dj, l in ((1, 0, 1), (-1, 0, 1), (0, 1, 1), (0, -1, 1), (1, 1, 1.4142), (1, -1, 1.4142), (-1, 1, 1.4142), (-1, -1, 1.4142)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < NX and 0 <= jj < NY):
                continue
            kk = jj * NX + ii
            if kk in done or not math.isfinite(mult[kk]):
                continue
            nd = d + l * STEP * (mult[k] + mult[kk]) / 2
            if nd < dist.get(kk, math.inf):
                dist[kk], prev[kk] = nd, k
                heapq.heappush(hq, (nd, kk))
    if ek not in dist:
        return None, mult
    path, k = [], ek
    while True:
        path.append((EXT["x0"] + (k % NX) * STEP, EXT["y0"] + (k // NX) * STEP))
        if k == sk:
            break
        k = prev[k]
    return path[::-1], mult


def seg_cost(G_, mult, a, b):
    L = math.hypot(b[0] - a[0], b[1] - a[1]); n = max(1, math.ceil(L / 0.2)); tot = 0
    for t in range(n + 1):
        i, j = cell(G_, a[0] + (b[0] - a[0]) * t / n, a[1] + (b[1] - a[1]) * t / n)
        mu = mult[j * G_["NX"] + i]
        if not math.isfinite(mu):
            return math.inf
        tot += mu
    return tot / (n + 1) * L


def solve(G_, mult, A, B):
    path, mult = dijkstra(G_, mult, cell(G_, *A), cell(G_, *B))
    if not path:
        return None
    path[0], path[-1] = tuple(A), tuple(B)
    cum = [0.0]
    for p, q in zip(path[:-1], path[1:]):
        cum.append(cum[-1] + seg_cost(G_, mult, p, q))
    out, i = [path[0]], 0
    while i < len(path) - 1:
        j = len(path) - 1
        while j > i + 1 and seg_cost(G_, mult, path[i], path[j]) > (cum[j] - cum[i]) * 1.02:
            j -= 1
        out.append(path[j]); i = j
    return out


def similar(poly, R):
    Q = build_route(poly)
    ds = [min(math.hypot(R["xs"][k] - Q["xs"][i], R["ys"][k] - Q["ys"][i]) for k in range(0, R["n"], 3)) for i in range(0, Q["n"], 5)]
    return sum(ds) / len(ds) < 0.45


def plan_routes(A, B, t0, sc, fl: Flight, auto_alt=True):
    top = fl.alt
    alts = sorted({a for a in (top, top - 150, top - 300) if a >= 150}, reverse=True) if auto_alt else [top]
    grids = [(z, cost_grid(z, t0 + 0.1, sc, fl.lim)) for z in alts]
    picks = []
    for prof in PROFILES:
        best = None
        for z, G_ in grids:
            poly = solve(G_, mult_grid(G_, prof, A, B, []), A, B)
            if poly and any(similar(poly, p["R"]) for p in picks):
                poly = solve(G_, mult_grid(G_, prof, A, B, [p["R"] for p in picks]), A, B) or poly
            if not poly:
                continue
            R = build_route(poly)
            r = evaluate(t0, replace(fl, alt=z), sc, R, fast=True)
            J = r["total_s"] / 60 + prof["lam"] * r["ri_mean"] / 10 + (60 if r["ri_max"] >= 80 else 0) + r["energy"]["pct"] / 20 + (4 if z < 300 else 0)
            if not best or J < best["J"]:
                best = dict(prof=prof, R=R, alt=z, J=J, r=r)
        if best:
            picks.append(best)
    for p in picks:
        r = p["r"]
        p["S"] = r["total_s"] / 60 + 2.5 * r["ri_mean"] / 10 + (50 if r["ri_max"] >= 80 else 10 if r["ri_max"] >= 60 else 0) + r["energy"]["pct"] / 15 + (40 if r["c"] == 3 else 0)
    return sorted(picks, key=lambda p: p["S"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="src", default="킨텍스 V5", choices=VERTIPORTS)
    ap.add_argument("--to", dest="dst", default="여의도", choices=VERTIPORTS)
    ap.add_argument("--scenario", default="conv", choices=SCENARIOS)
    ap.add_argument("--aircraft", default="joby")
    ap.add_argument("--adult", type=int, default=2)
    ap.add_argument("--child", type=int, default=0)
    ap.add_argument("--bag", type=float, default=10, help="1인당 수하물 kg")
    ap.add_argument("--lead", type=float, default=0.0, help="출발 선행시간 (시간)")
    ap.add_argument("--vto", type=float); ap.add_argument("--vc", type=float); ap.add_argument("--vld", type=float)
    ap.add_argument("--alt", type=float, help="순항고도 상한 m")
    ap.add_argument("--fixed-alt", action="store_true", help="순항고도 자동 추천 끄기")
    ap.add_argument("--list-aircraft", action="store_true")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    fleet = load_fleet()
    if a.list_aircraft:
        for f in fleet.values():
            print(f"{f['id']:<11} {f['country']:<4} {f['name']:<22} 승객 {f['seats']} · 순항 {f['vc']} km/h · 항속 {f['range']} km · MTOW {f['mtow']} kg · 배터리 {f['battKWh']} kWh(추정)")
        return
    air = fleet[a.aircraft]
    if a.adult + a.child > air["seats"]:
        raise SystemExit(f"{air['name']} 승객 좌석은 {air['seats']}석입니다.")
    fl = make_flight(air, a.adult, a.child, a.bag, vto=a.vto, vc=a.vc, vld=a.vld, alt=a.alt)
    sc = SCENARIOS[a.scenario]
    picks = plan_routes(LANDMARKS[a.src], LANDMARKS[a.dst], a.lead, sc, fl, not a.fixed_alt)
    out = []
    for k, p in enumerate(picks, 1):
        r = evaluate(a.lead, replace(fl, alt=p["alt"]), sc, p["R"])
        w = r["pts"][r["worst_i"]]
        out.append(dict(rank=k, name=p["prof"]["name"], verdict=LV[r["c"]], limiting=r["limit"], dist_km=round(p["R"]["L"], 2),
                        time_min=round(r["total_s"] / 60, 1), cruise_alt=round(r["top"]), ri_max=round(r["ri_max"]),
                        ri_mean=round(r["ri_mean"]), kwh=round(r["energy"]["kwh"], 1), batt_pct=round(r["energy"]["pct"], 1),
                        worst_at=nearest_name(p["R"]["xs"][r["worst_i"]], p["R"]["ys"][r["worst_i"]]), worst_alt=round(w["z"]),
                        waypoints=[[round(x, 2), round(y, 2)] for x, y in p["R"]["poly"]]))
    if a.json:
        print(json.dumps(dict(aircraft=air["name"], mass_kg=round(fl.mass), scenario=sc.label, routes=out), ensure_ascii=False, indent=1))
    else:
        print(f"{air['name']} · 비행중량 {fl.mass:.0f} kg · {sc.label} · {a.src} → {a.dst} (+{a.lead} h)")
        for o in out:
            print(f" 추천 {o['rank']} {o['name']:<5} {o['verdict']} ({o['limiting']}) · {o['dist_km']} km · {o['time_min']} 분 · "
                  f"순항 {o['cruise_alt']} m · 위험도 최대 {o['ri_max']}/평균 {o['ri_mean']} · {o['kwh']} kWh ({o['batt_pct']}%)")


if __name__ == "__main__":
    main()
