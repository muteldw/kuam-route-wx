"""K-UAM 경로 기상 위험 단순 진단모델 (참조 구현).

웹 초안(web/index.html)의 JavaScript 모델과 같은 식을 쓴다.
실측 연결 시 met() 안의 각 항을 LDAPS/ECMWF/관측 자료로 교체한다.

사용 예:
    python kuam_simple_model.py --scenario conv --aircraft lift --lead 0.5 > route.json
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field, asdict

G, KAPPA = 9.81, 0.4

# 좌표: km, 원점 킨텍스 V5 (126.7456E, 37.6688N), x 동쪽 / y 북쪽. 교량 위치는 근사
WP = [("킨텍스 V5", 0, 0), ("일산대교", 1.45, -3.94), ("행주대교", 5.67, -7.68),
      ("방화대교", 7.79, -9.57), ("가양대교", 10.26, -11.78), ("성산대교", 12.64, -13.15),
      ("여의도", 15.76, -15.89)]
# 한강 중심선 (Natural Earth 10m + 하류 근사 연장, tools/build_basemap.py 와 같은 값)
RIVER = [(-6.133, 11.228), (-6.133, 10.916), (-6.133, 10.136), (-6.126, 8.922), (-6.064, 7.518), (-5.94, 5.958), (-5.738, 4.356), (-5.424, 2.942), (-4.998, 1.714), (-4.486, 0.716), (-3.9, -0.033), (-3.247, -0.557), (-2.586, -1.025), (-1.925, -1.462), (-1.263, -1.889), (-0.594, -2.349), (0.08, -2.842), (0.731, -3.341), (1.344, -3.831), (2.03, -4.411), (3.558, -5.762), (6.04, -7.982), (8.824, -10.467), (10.613, -12.014), (11.404, -12.621), (11.983, -12.906), (12.741, -13.178), (13.673, -13.457), (14.754, -13.894), (15.98, -14.512), (17.156, -15.175), (17.89, -15.612), (18.182, -15.823), (18.33, -15.998), (18.48, -16.233), (18.635, -16.52), (18.807, -16.812), (18.999, -17.1), (19.214, -17.364), (19.46, -17.563), (19.737, -17.697), (20.023, -17.752), (20.305, -17.725), (20.582, -17.617), (20.83, -17.464), (21.048, -17.268), (21.242, -17.04), (21.432, -16.796), (21.616, -16.537), (21.891, -16.204), (22.302, -15.769), (22.822, -15.257), (23.254, -14.854), (23.569, -14.586), (23.808, -14.421), (24.054, -14.293), (24.307, -14.202), (24.498, -14.144), (24.594, -14.115)]
AIRPORT = (3.97, -12.26)
YEO = (15.78, -15.9)


@dataclass
class Scenario:
    issueH: int = 15
    U: float = 7.0       # 600 m 일반풍 (m/s)
    dir: float = 250.0   # 풍향 (deg, 불어오는 방향)
    T: float = 28.0      # 일평균 기온 (degC)
    dep: float = 4.0     # 평균 이슬점차 (degC)
    llj: float = 0.0     # 야간 저층제트 세기 (m/s)
    conv: bool = False
    front: dict | None = None


SCENARIOS = {
    "calm": Scenario(14, 6, 270, 17, 7),
    "nw": Scenario(9, 17, 315, 1, 9, front={"h": 10, "dU": 5, "dDir": 15}),
    "llj": Scenario(21, 9, 235, 13, 4, llj=9),
    "conv": Scenario(15, 7, 250, 28, 4, conv=True),
    "fog": Scenario(4, 2, 80, 9, 0.6),
}
CONV = dict(t0=-0.8, t1=3.5, x0=-1, y0=-14, vx=5.5, vy=3, R=3, W=7, O=13, P=45)


@dataclass
class Aircraft:
    vc: float; alt: float; climb: float; desc: float; xw: float; tail: float
    wmax: float; gust: float; edr: float; vis: float; rain: float; icing: bool = False


AIRCRAFT = {  # 가상 기체 — 실제 제원으로 교체할 것
    "multi": Aircraft(100, 300, 3, 2.5, 8, 4, 12, 13, 0.20, 1600, 5),
    "lift": Aircraft(180, 450, 4, 3, 10, 5, 18, 16, 0.25, 1600, 10),
    "tilt": Aircraft(240, 600, 5, 4, 12, 6, 22, 19, 0.30, 1600, 15),
}


def _poly_dist(x, y, P):
    best, side = 1e9, 1
    for (ax, ay), (bx, by) in zip(P[:-1], P[1:]):
        dx, dy = bx - ax, by - ay
        t = max(0.0, min(1.0, ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy)))
        d = math.hypot(x - (ax + t * dx), y - (ay + t * dy))
        if d < best:
            best, side = d, (-1 if dx * (y - ay) - dy * (x - ax) < 0 else 1)
    return best, side


def _smooth(e0, e1, v):
    t = max(0.0, min(1.0, (v - e0) / (e1 - e0)))
    return t * t * (3 - 2 * t)


def logz0(x, y):
    """log10 z0 (m): 도시 0.5, 김포 농경지 0.1, 공항 0.03, 한강 0.0002, 여의도 고층 상향."""
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
    return land * (1 - wr) - 3.7 * wr, d


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
    lzE = lz + (-0.6 - lz) * min(1, zz / 400)          # 높이에 따른 내부경계층 혼합
    z0e = 10 ** lzE
    a = (0.096 * lzE + 0.016 * lzE * lzE + 0.24) * (1 - 0.45 * Seff)
    a = max(0.06, min(0.5, a))
    spd = Us * (zz / 600) ** a if zz <= 600 else Us * (1 + 0.0003 * (zz - 600))
    jet = dJ = 0.0
    if S < 0 and sc.llj > 0:
        jet = sc.llj * min(1, -S / 0.8) * math.exp(-((zz - 320) / 110) ** 2)
        dJ = jet * (-2 * (zz - 320) / 110 ** 2)
    r5 = math.hypot(x, y)
    A = 1.58 + 0.45 * math.cos(math.radians(wdir))       # MUKLIMO_3 V5 민감도 근사
    amp = 1 + (A - 1) * math.exp(-(r5 / 1.2) ** 2) * math.exp(-zz / 80)
    bw = math.exp(-(math.hypot(x - YEO[0], y - YEO[1]) / 0.9) ** 2) * math.exp(-zz / 150)
    base = spd * (1 - 0.25 * bw) + jet          # 난류용 기본풍 (국지 증폭 제외)
    spd = spd * amp * (1 - 0.25 * bw) + jet
    veer = (25 if Seff < 0 else 10) * (1 - min(zz, 600) / 600)
    dr = math.radians(wdir - veer)
    u, v = -spd * math.sin(dr), -spd * math.cos(dr)
    w = P = epsC = tiC = 0.0
    if sc.conv and CONV["t0"] <= t <= CONV["t1"]:
        tau = t - CONV["t0"]
        life = math.sin(math.pi * tau / (CONV["t1"] - CONV["t0"]))
        cx, cy = CONV["x0"] + CONV["vx"] * tau, CONV["y0"] + CONV["vy"] * tau
        r, R = max(0.01, math.hypot(x - cx, y - cy)), CONV["R"]
        O = CONV["O"] * life * (r / R) * math.exp(1 - r / R) * math.exp(-zz / 450)
        u += O * (x - cx) / r
        v += O * (y - cy) / r
        w = -CONV["W"] * life * math.exp(-(r / (0.6 * R)) ** 2) * min(1, zz / 250)
        P = CONV["P"] * life * math.exp(-(r / (0.7 * R)) ** 2)
        epsC = 0.04 * life * math.exp(-(r / R) ** 2)
        tiC = 0.15 * life * math.exp(-(r / R) ** 2)
    sH = math.hypot(u, v)
    lnz = math.log(zz / z0e)
    TI = (1 / lnz) * (1 + 0.5 * max(Seff, 0)) + 0.10 * bw \
        + 0.08 * math.exp(-(r5 / 0.8) ** 2) * math.exp(-zz / 100) + tiC
    us = KAPPA * (base + (sH - spd if sc.conv else 0)) / lnz
    eps = us ** 3 / (KAPPA * max(zz, 50)) * (1 + bw) + 0.0015 * max(Seff, 0) * max(0, 1 - zz / 1500) \
        + 1.5 * dJ * dJ + epsC
    sd = math.sin(2 * math.pi * (h - 9) / 24)
    T0 = sc.T + 5 * sd
    dep = max(0.15, sc.dep * (1 + 0.75 * sd))
    es = lambda T: math.exp(17.625 * T / (243.04 + T))
    RH = 100 * es(T0 - dep) / es(T0)
    ceil = max(30, 125 * dep) if dep < 2.5 else 3000  # 습윤할 때만 낮은 구름
    f = max(0, min(1, (RH - 95) / 4)) * max(0, min(1, (5 - 0.35 * Us) / 3)) \
        * (1 if S < 0.25 else max(0, 1 - (S - 0.25) * 3))
    f = min(1, f * (1 + 0.4 * math.exp(-dRiver ** 2)))
    fog_top = 60 + 200 * f
    vis = 10000 * math.exp(-4.6 * f) if zz < fog_top else 10000
    in_cloud = zz > ceil
    if in_cloud:
        vis = min(vis, 300)
    vis = max(80, min(vis, 12000 / (1 + 0.3 * P)))
    Tz = T0 - 0.0065 * zz
    icing = -12 <= Tz <= 0 and (in_cloud or (f > 0.3 and zz < fog_top) or P > 0.5)
    return dict(u=u, v=v, w=w, spd=sH, gust=sH * (1 + 2.2 * TI), edr=eps ** (1 / 3), vis=vis,
                ceil=ceil, T=Tz, rh=RH, icing=icing, P=P, z0=z0e)


def build_route(ds=0.1):
    pts, acc = [], 0.0
    for (n0, ax, ay), (_, bx, by) in zip(WP[:-1], WP[1:]):
        L = math.hypot(bx - ax, by - ay)
        n = max(1, round(L / ds))
        for k in range(n):
            f = k / n
            pts.append(dict(s=acc + L * f, x=ax + (bx - ax) * f, y=ay + (by - ay) * f,
                            tx=(bx - ax) / L, ty=(by - ay) / L, seg=n0))
        acc += L
    pts.append(dict(s=acc, x=WP[-1][1], y=WP[-1][2], tx=pts[-1]["tx"], ty=pts[-1]["ty"], seg=WP[-2][0]))
    return pts, acc


LV = ["정상", "주의", "경고", "불가"]


def cat(m):
    return 0 if m >= 0.2 else 1 if m >= 0 else 2 if m >= -0.25 else 3


def evaluate(lead, ac: Aircraft, sc: Scenario):
    pts, L = build_route()
    Vc = ac.vc / 3.6
    Vcl = 0.55 * Vc
    dc, dd, top = ac.alt / ac.climb * Vcl / 1000, ac.alt / ac.desc * Vcl / 1000, ac.alt
    if dc + dd > 0.9 * L:
        k = 0.9 * L / (dc + dd)
        dc, dd, top = dc * k, dd * k, top * k
    cum = 0.0
    for i, p in enumerate(pts):
        s = p["s"]
        if s < dc:
            p.update(z=top * s / dc, vg=Vcl, phase="climb")
        elif s > L - dd:
            p.update(z=top * (L - s) / dd, vg=Vcl, phase="desc")
        else:
            p.update(z=top, vg=Vc, phase="cruise")
        p["z"] = max(p["z"], 5)
        if i:
            cum += (s - pts[i - 1]["s"]) * 1000 / p["vg"]
        p["eta"] = lead + cum / 3600
        m = met(p["x"], p["y"], p["z"], p["eta"], sc)
        zt, zb = p["z"] + 15, max(2, p["z"] - 15)
        m1, m0 = met(p["x"], p["y"], zt, p["eta"], sc), met(p["x"], p["y"], zb, p["eta"], sc)
        p.update(m=m, Wx=m["u"] * p["tx"] + m["v"] * p["ty"],
                 cross=abs(-m["u"] * p["ty"] + m["v"] * p["tx"]),
                 vsh=math.hypot(m1["u"] - m0["u"], m1["v"] - m0["v"]) * 30 / (zt - zb))
    n = len(pts)
    Fr = []
    for i, p in enumerate(pts):
        a, b = max(0, i - 3), min(n - 1, i + 3)
        ds = (pts[b]["s"] - pts[a]["s"]) * 1000
        Va = max(15, p["vg"] - p["Wx"])
        Fr.append((pts[b]["Wx"] - pts[a]["Wx"]) / ds * p["vg"] / G - p["m"]["w"] / Va)
        p["hsh"] = abs(pts[b]["Wx"] - pts[a]["Wx"]) * 0.5 / (ds / 1000)
    out = []
    for i, p in enumerate(pts):
        win = Fr[max(0, i - 5): min(n, i + 6)]
        F = sum(win) / len(win)
        m = p["m"]
        low = p["z"] < 150 and p["phase"] != "cruise"
        # eVTOL은 바람을 향해 이착륙할 수 있으므로 저고도 구간은 총풍속으로 판정
        items = [("이착륙 풍속", 1 - m["spd"] / ac.xw)] if low else [("풍속", 1 - m["spd"] / ac.wmax)]
        items += [("돌풍", 1 - m["gust"] / ac.gust), ("난류", 1 - m["edr"] / ac.edr),
                  ("시정", 1 - ac.vis / m["vis"]), ("운고", 1 - (p["z"] + 100) / m["ceil"]),
                  ("급변풍", 1 - F / 0.105 if F > 0 else 1), ("수평시어", 1 - p["hsh"] / 7.5),
                  ("연직시어", 1 - p["vsh"] / 2.0), ("강수", 1 - m["P"] / ac.rain)]
        if m["icing"] and not ac.icing:
            items.append(("착빙", -0.5))
        items = [(k, max(-1.0, min(1.0, v))) for k, v in items]
        f, mg = min(items, key=lambda kv: kv[1])
        out.append(dict(s=round(p["s"], 2), x=round(p["x"], 3), y=round(p["y"], 3), alt=round(p["z"]),
                        eta_h=round(p["eta"], 4), spd=round(m["spd"], 2), gust=round(m["gust"], 2),
                        edr=round(m["edr"], 3), vis=round(m["vis"]), F=round(F, 3),
                        hshear=round(p["hsh"], 2), vshear30=round(p["vsh"], 2),
                        margin=round(mg, 3), limiting=f, level=LV[cat(mg)]))
    worst = min(out, key=lambda r: r["margin"])
    return dict(lead_h=lead, aircraft=asdict(ac), route_km=round(L, 2),
                verdict=LV[cat(worst["margin"])], worst=worst, points=out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scenario", default="conv", choices=SCENARIOS)
    ap.add_argument("--aircraft", default="lift", choices=AIRCRAFT)
    ap.add_argument("--lead", type=float, default=0.0, help="출발 선행시간 (시간)")
    ap.add_argument("--summary", action="store_true", help="요약만 출력")
    a = ap.parse_args()
    res = evaluate(a.lead, AIRCRAFT[a.aircraft], SCENARIOS[a.scenario])
    if a.summary:
        w = res["worst"]
        print(f"{a.scenario}/{a.aircraft} +{a.lead}h: {res['verdict']} "
              f"(제한 {w['limiting']}, 여유율 {w['margin']:.0%}, {w['s']} km, {w['alt']} m)")
    else:
        print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
