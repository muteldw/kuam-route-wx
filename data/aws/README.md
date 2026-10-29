# 킨텍스 버티포트 주변 AWS 최근 5년 자료 수집

이 작업 환경에서는 기상청·공공데이터포털이 막혀 있어서, GitHub Actions(접속 확인됨)가 자료를 받아 이 브랜치에 커밋합니다.

## 준비 (1회)
1. 기상청 API허브(apihub.kma.go.kr) 가입 → 지상관측 'AWS 매분자료' API 활용신청 → 인증키
2. (선택) 공공데이터포털 '기상청_지상(종관, ASOS) 시간자료 조회서비스' 활용신청 → 일반 인증키(Decoding)
3. GitHub 저장소 Settings → Secrets and variables → Actions → New repository secret
   - `KMA_APIHUB_KEY` = API허브 인증키
   - `DATA_GO_KR_KEY` = 공공데이터포털 인증키(Decoding)

## 단계
1. discovery: 지점 목록, 자료 형식 확인 → `data/aws/discovery/`
2. fetch: 킨텍스(37.668N, 126.745E) 반경 약 15km AWS·ASOS 지점, 2021.11~2026.10 수집(연·월 단위 압축)
3. analyze: 풍향·풍속 빈도, 강풍·돌풍 사례 추출 → 민감도 실험 조견표로 버티포트 풍속 환산
