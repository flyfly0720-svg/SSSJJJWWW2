# 비만 → YAP/TAZ 과활성화 → 암 발병 가설 검증 대시보드

에너지 과잉 상태 → 지방세포 크기 증가 → 세포골격(F-actin) 긴장 상승 → LATS 억제 →
YAP/TAZ 핵 내 과활성화 → 비만 연관 암종의 발암 촉진이라는 가설을, **두 축**으로 나누어
검증하는 Streamlit 웹 앱입니다.

- **축 1 (비만 → YAP/TAZ)**: Emont Lab 인체 지방세포 단일세포 데이터(`human_adipocytes_lite.rds`)에서 비만군의 YAP/TAZ 표적 유전자 점수가 정상군보다 높은지 확인
- **축 2 (YAP/TAZ → 암)**: GEPIA2 TCGA 범암 데이터에서 CCN2(CTGF)·CCN1(CYR61) 발현이 높은 암종일수록 역학적으로 비만과 강하게 연관된 암종인지 상관관계 분석

## 실행 방법

```bash
pip install -r requirements.txt
streamlit run app.py
```

## 폴더 구성

| 파일 | 역할 |
|---|---|
| `app.py` | Streamlit 메인 앱 (4개 탭: 가설 개요 / 지방세포 / TCGA 범암 / 상관관계) |
| `data_reference.py` | 실제 문헌·공개 데이터 기반 참고 테이블 (아래 "데이터 출처" 참고) |
| `export_from_seurat.R` | `.rds`(Seurat) 파일에서 앱이 읽을 수 있는 CSV를 만드는 R 스크립트 |
| `requirements.txt` | Python 의존 패키지 |

## 왜 PDF/rds 파일을 앱이 직접 읽지 않는가

1. **GEPIA2 boxplot PDF (CTGF, CYR61)**: 두 PDF는 수천 개의 점을 그린 산점도 이미지입니다.
   PDF에서 텍스트를 추출하면 점의 좌표(실제 발현값)는 사라지고 글자로 남아 있는 정보,
   즉 각 암종의 표본 수(`num(T)`, `num(N)`)만 정확히 읽을 수 있습니다. 이 값은
   `data_reference.py`의 `GEPIA2_COHORT`에 그대로 반영해 두었습니다.
   **실제 발현값·log2FC 수치는 임의로 만들어 넣지 않았으며**, GEPIA2 웹사이트
   (http://gepia2.cancer-pku.cn → Expression DIY → Box Plots → Download)에서 직접
   내려받아 ③ 탭의 템플릿에 채워 업로드해야 합니다.
2. **`human_adipocytes_lite.rds`**: Seurat 객체는 R의 S4 클래스 구조라서 Python
   (pandas/pyreadr 등)으로 안전하게 파싱할 수 없습니다. `export_from_seurat.R`을
   R/RStudio에서 먼저 실행해 CSV로 변환한 뒤 ② 탭에 업로드하세요.

## 앱 사용 순서

1. **① 가설·데이터 개요**: 기전 설명, 사용 유전자 목록, 참고문헌 확인
2. **② 지방세포 탭**: `export_from_seurat.R`로 만든 CSV 업로드 → 비만군/정상군 간
   YAP/TAZ 점수를 Plotly box plot으로 비교, Mann-Whitney U 검정으로 유의성 확인
3. **③ TCGA 범암 탭**: 템플릿 CSV를 내려받아 GEPIA2에서 얻은 CCN2/CCN1 log2FC 값을
   채운 뒤 업로드 → 암종별 발현 변화를 막대그래프로 확인 (IARC 비만 연관 암종은 빨간색)
4. **④ 상관관계 탭**: ③에서 올린 데이터가 자동으로 연결되어, 암종별 비만-암 상대위험도
   (RR)와 YAP/TAZ 표적 유전자 log2FC 사이의 Pearson/Spearman 상관계수·산점도·회귀선을
   보여주고, IARC 비만 연관군 vs 미확립군 간 log2FC 차이를 Mann-Whitney U 검정으로 비교

## 데이터 출처 (`data_reference.py`)

- **GEPIA2 코호트 표본 수**: 업로드한 CTGF/CYR61 boxplot PDF에서 실제로 읽은 값
- **비만-암 연관성 분류·상대위험도(RR)**:
  - Lauby-Secretan B, et al. *Body Fatness and Cancer — Viewpoint of the IARC
    Working Group.* N Engl J Med. 2016;375:794-798. (비만 연관 확정 13개 암종)
  - Renehan AG, et al. *Body-mass index and incidence of cancer: a systematic
    review and meta-analysis of prospective observational studies.* Lancet.
    2008;371:569-578. (BMI +5 kg/m² 당 상대위험도, 성별 구분)
  - Bhaskaran K, et al. *Body-mass index and risk of 22 specific cancers: a
    population-based cohort study of 5.24 million UK adults.* Lancet.
    2014;384:755-765.
- **GEPIA2**: Tang Z, et al. *GEPIA2: a web server for large-scale expression
  profiling and interactive analysis.* Nucleic Acids Res. 2019;47:W556-W560.
- **지방세포 단일세포 아틀라스**: Emont MP, et al. *A single-cell atlas of human
  and mouse white adipose tissue.* Nature. 2022;603:926-933.

TCGA 암종 코드가 위 역학 연구의 암종 분류와 조직학적으로 정확히 일치하지 않는 경우
(예: `CHOL`=담관암 vs 문헌상 담낭암, `ESCA`/`STAD`는 선암·편평상피암 혼합)는
`data_reference.py`의 `note` 필드에 한계를 명시했습니다.

## 해석상 주의 (제한점)

- ④ 탭의 상관관계는 **암종 단위의 집단 통계(ecological correlation)**이며, 개별 환자의
  BMI와 유전자 발현을 직접 짝지은 값이 아닙니다.
- TCGA의 Tumor vs Normal 비교는 '비만 여부'가 아니라 '종양 조직 여부'의 비교이므로,
  이 결과만으로 "비만이 YAP/TAZ를 통해 암을 유발한다"는 **인과관계를 증명할 수 없습니다.**
  가설을 뒷받침하는 정황 증거(association) 수준으로 해석해야 합니다.
- ② 탭의 "원자료 업로드" 모드에서 계산하는 YAP/TAZ 점수는 유전자별 z-score의 단순 평균으로,
  Seurat `AddModuleScore()`의 대조유전자(control gene) 보정 방식과는 다른 **단순화된
  근사치**입니다. 정확한 값이 필요하면 `export_from_seurat.R`로 Seurat에서 직접 계산한
  점수를 사용하세요.
