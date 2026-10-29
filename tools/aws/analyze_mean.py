"""5개 AWS 10분 자료: 품질검사 → 유사성 검사 → 지점 평균 풍속 계열 생성.
풍속 중심. 풍향은 벡터평균(참고용)으로만 계산하고, 사용 여부는 풍향 일치도 지표를 보고 판단."""
import sys, glob, json, numpy as np, pandas as pd
def load(d):
    fs = sorted(glob.glob(f'{d}/aws10_*.csv.gz'))
    df = pd.concat([pd.read_csv(f, dtype={'TM': str, 'STN': str}) for f in fs])
    df['time'] = pd.to_datetime(df.TM, format='%Y%m%d%H%M')
    df = df[(df.N1MIN >= 8)]                       # 10분 중 8분 이상 관측
    return df
def hourly(df, col):
    return df.pivot_table(index='time', columns='STN', values=col).resample('1h').mean()
def screen(df, r_min=0.8, ratio_lo=0.8, ratio_hi=1.25, avail_min=0.9):
    H = hourly(df, 'WS10'); n_exp = len(pd.date_range(H.index.min(), H.index.max(), freq='1h'))
    rep = pd.DataFrame({'가용률': H.notna().sum() / n_exp, '평균풍속': H.mean()})
    C = H.corr(); rep['상관(타지점 중앙값)'] = [(C[s].drop(s)).median() for s in C]
    med = H.median(axis=1); rep['중앙값 대비 풍속비'] = (H.div(med, axis=0)).median()
    # 풍향 일치도: 3 m/s 이상일 때 지점 풍향과 다지점 벡터평균 풍향의 차이 중앙값
    D = hourly(df, 'WD10').reindex_like(H); u = -H * np.sin(np.radians(D)); v = -H * np.cos(np.radians(D))
    mdir = (np.degrees(np.arctan2(-u.mean(axis=1), -v.mean(axis=1))) % 360)
    dd = ((D.sub(mdir, axis=0) + 180) % 360 - 180).abs().where(H >= 3)
    rep['풍향차 중앙값(°)'] = dd.median()
    rep['사용'] = (rep['가용률'] >= avail_min) & (rep['상관(타지점 중앙값)'] >= r_min) & rep['중앙값 대비 풍속비'].between(ratio_lo, ratio_hi)
    return rep.round(3), C.round(3)
def mean_series(df, use):
    d = df[df.STN.isin(use)]
    W = d.pivot_table(index='time', columns='STN', values='WS10')
    G = d.pivot_table(index='time', columns='STN', values='WSS_MAX')
    D = d.pivot_table(index='time', columns='STN', values='WD10').reindex_like(W)
    u = (-W * np.sin(np.radians(D))).mean(axis=1); v = (-W * np.cos(np.radians(D))).mean(axis=1)
    out = pd.DataFrame({'ws_mean': W.mean(axis=1), 'ws_n': W.notna().sum(axis=1), 'ws_std': W.std(axis=1),
        'gust_max': G.max(axis=1), 'wd_vec': (np.degrees(np.arctan2(-u, -v)) % 360).round(1),
        'wd_steadiness': (np.hypot(u, v) / W.mean(axis=1)).round(3)})
    return out[out.ws_n >= max(3, len(use) - 1)]
if __name__ == '__main__':
    df = load(sys.argv[1]); rep, C = screen(df)
    print(rep.to_string()); print(C.to_string())
    use = rep.index[rep['사용']].tolist(); ms = mean_series(df, use)
    ms.to_csv(sys.argv[2]); print('사용 지점', use, '평균 계열', len(ms), '건')
