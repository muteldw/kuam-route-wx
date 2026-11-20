"""rt-data/latest.json을 web/index.html 안에 실황 스냅샷(RT_SNAP)으로 넣는다.
claude.ai 아티팩트처럼 외부 자료를 못 읽는 화면에서도 실황 모드로 계산하게 하기 위함.
사용: python tools/rt/embed_snapshot.py latest.json [web/index.html]
"""
import json, re, sys

src = sys.argv[1]; html = sys.argv[2] if len(sys.argv) > 2 else 'web/index.html'
d = json.load(open(src, encoding='utf-8'))
d.pop('log', None)
blob = json.dumps(d, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
block = f'<!--RT_SNAP--><script>const RT_SNAP = {blob};</script><!--/RT_SNAP-->'
s = open(html, encoding='utf-8').read()
if '<!--RT_SNAP-->' in s:
    s = re.sub(r'<!--RT_SNAP-->.*?<!--/RT_SNAP-->', lambda m: block, s, flags=re.S)
else:
    i = s.index('<script>\n(() => {')
    s = s[:i] + block + '\n' + s[i:]
open(html, 'w', encoding='utf-8').write(s)
print('embedded', d['obs']['tm'], len(blob), 'bytes')
