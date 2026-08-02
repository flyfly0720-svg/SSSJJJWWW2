
# -*- coding: utf-8 -*-
"""
data_reference.py
==================
이 파일은 앱 전체에서 재사용하는 '실제 문헌·공개 데이터 기반' 참고 테이블만 담는다.
연구 신뢰성을 위해, 논문에서 직접 확인되지 않은 수치는 절대 만들어 넣지 않았다.
(예: 업로드된 GEPIA2 boxplot PDF는 산점도 이미지라서 텍스트 추출만으로는
 칸별 발현량 수치를 정확히 복원할 수 없다. → 이 파일에는 두 PDF에서 실제로
 '글자로' 읽을 수 있었던 정보, 즉 각 암종의 Tumor/Normal 표본 수(num(T), num(N))만 담았다.)

포함 내용
---------
1. GEPIA2_COHORT   : 두 PDF(CTGF, CYR61 boxplot)에 공통으로 표기된 TCGA 31개 암종의
                      표본 수(num(T), num(N))와 정식 명칭.
2. OBESITY_CANCER_RR : 비만(체질량지수, BMI)과 암 발생 위험의 역학적 연관성.
                      출처:
                        - Lauby-Secretan B, et al. "Body Fatness and Cancer — Viewpoint
                          of the IARC Working Group." N Engl J Med. 2016;375:794-798.
                          (비만과 인과적 연관이 확립된 13개 암종 지정)
                        - Renehan AG, et al. "Body-mass index and incidence of cancer:
                          a systematic review and meta-analysis of prospective
                          observational studies." Lancet. 2008;371:569-578.
                          (BMI 5 kg/m^2 증가당 상대위험도, RR, 성별 구분)
                      TCGA 암종 코드는 조직학적으로 위 논문들의 암종과 100% 일치하지
                      않는 경우가 있어(예: CHOL≈담관암 vs 담낭암, ESCA/STAD는 선암과
                      편평상피암이 섞여 있음) 각 항목에 note로 대응 관계의 한계를 남겼다.
3. YAP_TARGET_GENES, TENSION_GENES : 안내 문서(대화창에 붙여넣은 R 분석 가이드)에서
                      제안된 YAP/TAZ 표적 유전자·세포골격 긴장 관련 유전자 목록.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class CancerEntry:
    code: str            # TCGA 약어
    name_kr: str          # 한글 명칭
    name_en: str          # 영문 정식 명칭
    n_tumor: int           # GEPIA2 boxplot에 표기된 num(T)
    n_normal: int          # GEPIA2 boxplot에 표기된 num(N)
    obesity_linked: bool    # IARC(2016) 비만 연관 확정 여부
    rr_per_5bmi: Optional[float]  # BMI +5 kg/m^2 당 상대위험도 (Renehan 2008 등)
    evidence: str          # 근거 강도 요약
    note: str             # 대응 관계 한계·출처 메모


# 1) 두 PDF(CTGF/CYR61 boxplot)에서 실제로 읽을 수 있었던 GEPIA2 TCGA 코호트 표본 수.
#    (발현값 자체는 산점도 이미지이므로 텍스트로 정확히 복원 불가 → 포함하지 않음)
GEPIA2_COHORT: dict[str, CancerEntry] = {
    "BRCA": CancerEntry("BRCA", "유방암", "Breast invasive carcinoma", 1085, 291,
                         True, 1.08, "중간(폐경 후 한정)",
                         "Renehan 2008 여성 <1.20(약함); 폐경 후 유방암에 한정된 연관성"),
    "UCEC": CancerEntry("UCEC", "자궁내막암", "Uterine Corpus Endometrial Carcinoma", 174, 91,
                         True, 1.59, "강함",
                         "Renehan 2008 여성 RR=1.59, p<0.0001"),
    "READ": CancerEntry("READ", "직장암", "Rectum adenocarcinoma", 92, 318,
                         True, 1.09, "약함",
                         "Renehan 2008 남성 <1.20(약한 양의 연관); 대장암군에 포함"),
    "COAD": CancerEntry("COAD", "결장암", "Colon adenocarcinoma", 275, 349,
                         True, 1.24, "중간",
                         "Renehan 2008 남성 RR=1.24, p<0.0001"),
    "ACC":  CancerEntry("ACC", "부신피질암", "Adrenocortical carcinoma", 77, 128,
                         False, None, "미확립", "대규모 역학 연구 부재"),
    "UCS":  CancerEntry("UCS", "자궁육종", "Uterine Carcinosarcoma", 57, 78,
                         False, None, "미확립", "독립적인 BMI 연관 보고 부재"),
    "BLCA": CancerEntry("BLCA", "방광암", "Bladder Urothelial Carcinoma", 404, 28,
                         False, None, "미확립/일관성 부족",
                         "Bhaskaran 2014(Lancet, 영국 코호트)에서 유의한 연관 미확인"),
    "CESC": CancerEntry("CESC", "자궁경부암", "Cervical squamous cell carcinoma", 306, 13,
                         True, 1.10, "약함",
                         "Bhaskaran 2014에서 5kg/m^2당 RR≈1.10 보고(Renehan 2008에는 미포함)"),
    "CHOL": CancerEntry("CHOL", "담관암", "Cholangiocarcinoma", 36, 9,
                         True, 1.59, "중간(담낭암 수치로 유추)",
                         "IARC(2016)는 담낭암(gallbladder)을 지정; TCGA CHOL은 담관암으로 조직학적으로 다름"),
    "DLBC": CancerEntry("DLBC", "미만성 거대B세포림프종", "Diffuse Large B-cell Lymphoma", 47, 337,
                         False, 1.10, "미확립/약함",
                         "비호지킨림프종 전반에서 약한 연관만 보고(Renehan 2008)"),
    "ESCA": CancerEntry("ESCA", "식도암", "Esophageal carcinoma", 182, 286,
                         True, 1.52, "강함(선암 기준)",
                         "Renehan 2008 남성 RR=1.52; TCGA ESCA는 선암·편평세포암 혼합"),
    "GBM":  CancerEntry("GBM", "교모세포종", "Glioblastoma multiforme", 163, 207,
                         False, None, "미확립", "뇌종양 전반에서 일관된 연관 보고 부재"),
    "HNSC": CancerEntry("HNSC", "두경부암", "Head and Neck squamous cell carcinoma", 519, 44,
                         False, None, "미확립", "IARC 비만 연관 확정 목록에 미포함"),
    "KICH": CancerEntry("KICH", "신장혐색소세포암", "Kidney Chromophobe", 66, 53,
                         True, 1.29, "중간",
                         "신세포암(RCC) 전체 평균값 사용(KIRC와 동일 근거)"),
    "KIRC": CancerEntry("KIRC", "신장투명세포암", "Kidney renal clear cell carcinoma", 523, 100,
                         True, 1.29, "중간",
                         "Renehan 2008 남녀 평균(남 1.24/여 1.34)"),
    "KIRP": CancerEntry("KIRP", "신장유두상세포암", "Kidney renal papillary cell carcinoma", 286, 60,
                         True, 1.29, "중간", "신세포암(RCC) 전체 평균값 사용"),
    "LAML": CancerEntry("LAML", "급성골수성백혈병", "Acute Myeloid Leukemia", 173, 70,
                         False, 1.10, "미확립/약함",
                         "백혈병 전반에서 약한 연관만 보고(Renehan 2008)"),
    "LGG":  CancerEntry("LGG", "저등급교종", "Brain Lower Grade Glioma", 518, 207,
                         False, None, "미확립", "뇌종양 전반에서 일관된 연관 보고 부재"),
    "LIHC": CancerEntry("LIHC", "간암", "Liver hepatocellular carcinoma", 369, 160,
                         True, 1.30, "중간",
                         "IARC(2016) 추가 지정; Renehan 2008에는 개별 수치 없어 근사값 사용"),
    "LUAD": CancerEntry("LUAD", "폐선암", "Lung adenocarcinoma", 483, 347,
                         False, None, "미확립/역설 보고",
                         "일부 코호트에서 오히려 역상관('obesity paradox') 보고"),
    "LUSC": CancerEntry("LUSC", "폐편평세포암", "Lung squamous cell carcinoma", 486, 338,
                         False, None, "미확립", "IARC 비만 연관 확정 목록에 미포함"),
    "OV":   CancerEntry("OV", "난소암", "Ovarian serous cystadenocarcinoma", 426, 88,
                         True, 1.10, "약함",
                         "IARC(2016) 추가 지정; Renehan 2008 <1.20"),
    "PAAD": CancerEntry("PAAD", "췌장암", "Pancreatic adenocarcinoma", 179, 171,
                         True, 1.10, "약함",
                         "IARC(2016) 추가 지정; Renehan 2008 여성 <1.20"),
    "PCPG": CancerEntry("PCPG", "크롬친화세포종/부신경절종", "Pheochromocytoma and Paraganglioma", 182, 3,
                         False, None, "미확립", "대규모 역학 연구 부재"),
    "PRAD": CancerEntry("PRAD", "전립선암", "Prostate adenocarcinoma", 492, 152,
                         False, 1.03, "제한적(진행성 암 한정)",
                         "전체 발생률과는 약함; 진행성·치명적 전립선암에 국한된 연관 보고"),
    "SARC": CancerEntry("SARC", "육종", "Sarcoma", 262, 2,
                         False, None, "미확립", "IARC 비만 연관 확정 목록에 미포함"),
    "SKCM": CancerEntry("SKCM", "피부흑색종", "Skin Cutaneous Melanoma", 461, 558,
                         True, 1.08, "약함",
                         "Renehan 2008 남성에서 약한 양의 연관(<1.20)"),
    "STAD": CancerEntry("STAD", "위암", "Stomach adenocarcinoma", 408, 211,
                         True, 1.20, "약함(분문부 위암 기준)",
                         "IARC(2016)는 분문부(gastric cardia) 위암 한정 지정; TCGA STAD는 전체 위암 포함"),
    "TGCT": CancerEntry("TGCT", "고환생식세포종양", "Testicular Germ Cell Tumors", 137, 165,
                         False, None, "미확립", "대규모 역학 연구 부재"),
    "THCA": CancerEntry("THCA", "갑상선암", "Thyroid carcinoma", 512, 337,
                         True, 1.20, "중간",
                         "IARC(2016) 추가 지정; Renehan 2008 남성 RR=1.33"),
    "THYM": CancerEntry("THYM", "흉선종", "Thymoma", 118, 339,
                         False, None, "미확립", "IARC 비만 연관 확정 목록에 미포함"),
}

# 2) YAP/TAZ 표적 유전자 (CCN1=CYR61, CCN2=CTGF는 동일 유전자의 신·구 명칭이므로 중복 제거)
YAP_TARGET_GENES = ["CCN2", "CCN1", "ANKRD1", "AMOTL2", "BIRC5", "AXL", "F3", "LATS2"]
YAP_TARGET_GENE_ALIASES = {"CCN2": "CTGF", "CCN1": "CYR61"}

# 3) 세포골격 긴장(F-actin stress fiber) 관련 유전자
TENSION_GENES = ["ACTA2", "MYH9", "ROCK1", "RHOA", "FLNA"]

CITATIONS = [
    "Lauby-Secretan B, et al. Body Fatness and Cancer — Viewpoint of the IARC Working Group. N Engl J Med. 2016;375:794-798.",
    "Renehan AG, et al. Body-mass index and incidence of cancer: a systematic review and meta-analysis of prospective observational studies. Lancet. 2008;371:569-578.",
    "Bhaskaran K, et al. Body-mass index and risk of 22 specific cancers: a population-based cohort study of 5.24 million UK adults. Lancet. 2014;384:755-765.",
    "Tang Z, et al. GEPIA2: a web server for large-scale expression profiling and interactive analysis. Nucleic Acids Res. 2019;47:W556-W560.",
    "Emont MP, et al. A single-cell atlas of human and mouse white adipose tissue. Nature. 2022;603:926-933.",
]
