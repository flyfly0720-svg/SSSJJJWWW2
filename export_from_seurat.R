# export_from_seurat.R
# ---------------------------------------------------------------------------
# human_adipocytes_lite.rds (Emont Lab, Seurat 객체)에서
# Streamlit 앱(app.py, ② 탭)이 바로 읽을 수 있는 CSV 두 세트를 만든다.
# Python은 R의 S4 객체(Seurat)를 직접 읽지 못하므로, 이 단계는 R(RStudio 등)에서
# 먼저 실행해야 한다.
#
# 실행 전 준비: install.packages(c("Seurat", "dplyr"))
# ---------------------------------------------------------------------------

library(Seurat)
library(dplyr)

adipo <- readRDS("human_adipocytes_lite.rds")

# 1) 메타데이터 컬럼을 먼저 확인해서, 비만 여부를 나타내는 실제 컬럼명을 찾는다.
#    (예: BMI, obesity_status, bmi_group 등 논문/데이터셋마다 이름이 다르다)
print(colnames(adipo@meta.data))

# ---- 아래 "obesity_group" 부분을 위에서 확인한 실제 컬럼명으로 바꿔서 실행 ----

# 2) YAP/TAZ 표적 유전자 모듈 점수 계산 (Seurat 표준 함수, 대조유전자 보정 포함)
yap_genes <- c("CCN2", "CCN1", "ANKRD1", "AMOTL2", "BIRC5", "AXL", "F3", "LATS2")
adipo <- AddModuleScore(adipo, features = list(yap_genes), name = "YAP_score", ctrl = 100)

# 3-A) 방법 A: 이미 계산된 점수만 내보내기 → app.py ②탭 "이미 계산된 YAP 점수 CSV 업로드" 모드
score_out <- data.frame(
  sample_id = colnames(adipo),
  group     = adipo$obesity_group,   # <- 실제 컬럼명으로 교체
  yap_score = adipo$YAP_score1
)
write.csv(score_out, "yap_score_by_sample.csv", row.names = FALSE)

# 3-B) 방법 B: 원자료(정규화 발현행렬)를 통째로 내보내기 → app.py ②탭 "원자료 업로드" 모드
#      세포 수가 매우 많으면 파일이 커질 수 있으니, 필요하면 관심 유전자만 골라서 저장해도 된다.
target_genes <- unique(c(yap_genes, "ACTA2", "MYH9", "ROCK1", "RHOA", "FLNA"))
target_genes <- intersect(target_genes, rownames(adipo))

expr <- as.data.frame(GetAssayData(adipo, assay = "RNA", slot = "data")[target_genes, ])
write.csv(expr, "expression_matrix.csv")

meta_out <- data.frame(sample_id = colnames(adipo), group = adipo$obesity_group)  # <- 컬럼명 교체
write.csv(meta_out, "sample_groups.csv", row.names = FALSE)

cat("완료: yap_score_by_sample.csv, expression_matrix.csv, sample_groups.csv 생성됨\n")
