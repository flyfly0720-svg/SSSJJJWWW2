
# -*- coding: utf-8 -*-
"""
비만 → YAP/TAZ 과활성화 → 암 발병 가설 검증 대시보드
====================================================
에너지 과잉 상태 → 지방세포 크기 증가 → 세포골격(액틴) 긴장 상승 → LATS 억제
→ YAP/TAZ 핵 내 축적·과활성화 → (역학적으로) 비만 연관 암종에서 발암 촉진
이라는 가설을, (1) 지방세포 단일세포 데이터와 (2) TCGA 범암 발현 데이터,
두 축으로 나누어 사용자가 직접 검증할 수 있도록 만든 도구.

중요: 이 앱은 업로드된 GEPIA2 boxplot PDF의 이미지(산점도)에서 실제 발현
수치를 추출하지 않았다. 산점도를 텍스트로 변환하면 점의 좌표가 소실되어
신뢰할 수 있는 숫자가 나오지 않기 때문이다. 대신 두 PDF에서 실제로 글자로
읽을 수 있었던 표본 수(num(T), num(N))만 내장했고, 발현값·log2FC는 사용자가
GEPIA2에서 직접 내려받은 자료(또는 아래 템플릿에 채운 값)를 업로드하도록
설계했다. 지방세포 데이터(.rds)도 마찬가지로, Python이 R의 Seurat 객체를
직접 읽을 수 없으므로 R에서 CSV로 내보낸 뒤 업로드하는 구조다.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

from data_reference import (
    CCN1_SIG_DOWN,
    CCN1_SIG_UP,
    CCN2_SIG_DOWN,
    CCN2_SIG_UP,
    CITATIONS,
    DEMO_SOURCE,
    GEPIA2_COHORT,
    TENSION_GENES,
    YAP_TARGET_GENE_ALIASES,
    YAP_TARGET_GENES,
)

st.set_page_config(page_title="비만-YAP/TAZ-암 가설 검증", page_icon="🧬", layout="wide")

# ----------------------------------------------------------------------------
# 공통 유틸
# ----------------------------------------------------------------------------

def cohort_dataframe() -> pd.DataFrame:
    rows = []
    for e in GEPIA2_COHORT.values():
        rows.append(
            dict(
                code=e.code,
                name_kr=e.name_kr,
                name_en=e.name_en,
                n_tumor=e.n_tumor,
                n_normal=e.n_normal,
                obesity_linked=e.obesity_linked,
                rr_per_5bmi=e.rr_per_5bmi,
                evidence=e.evidence,
                note=e.note,
            )
        )
    return pd.DataFrame(rows)


def compute_module_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """expr: index=gene symbol, columns=sample. 유전자별 z-score를 낸 뒤
    목표 유전자 집합에 대해 표본별 평균을 낸다.
    (Seurat AddModuleScore의 대조유전자 보정과는 다른, 단순화된 근사치임을
    README/화면에 명시한다.)"""
    present = [g for g in genes if g in expr.index]
    if not present:
        return None, []
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present


def mannwhitney_table(df: pd.DataFrame, value_col: str, group_col: str) -> pd.DataFrame:
    groups = sorted(df[group_col].dropna().unique().tolist())
    out = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a = df.loc[df[group_col] == groups[i], value_col].dropna()
            b = df.loc[df[group_col] == groups[j], value_col].dropna()
            if len(a) < 2 or len(b) < 2:
                continue
            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            out.append({"비교": f"{groups[i]} vs {groups[j]}", "n1": len(a), "n2": len(b),
                        "U 통계량": round(stat, 2), "p-value": p})
    return pd.DataFrame(out)


if "gene_effect_df" not in st.session_state:
    st.session_state["gene_effect_df"] = None  # Tab3에서 업로드한 암종별 log2FC 저장 → Tab4에서 재사용
if "gepia_demo_df" not in st.session_state:
    st.session_state["gepia_demo_df"] = None

# ----------------------------------------------------------------------------
# 헤더
# ----------------------------------------------------------------------------
st.title("🧬 비만 → YAP/TAZ 과활성화 → 암 발병 가설 검증")
st.caption(
    "에너지 과잉 상태에서 지방세포 크기 증가·세포골격 긴장 상승이 LATS를 억제해 "
    "YAP/TAZ가 과활성화되고, 이것이 비만 연관 암종의 발암을 촉진할 수 있다는 가설을 "
    "단일세포 지방조직 데이터와 TCGA 범암 발현 데이터로 나누어 검증하는 도구입니다."
)

tab0, tab1, tab2, tab3 = st.tabs(
    ["① 가설·데이터 개요", "② 지방세포: 비만 → YAP/TAZ", "③ TCGA 범암: YAP/TAZ → 암", "④ 상관관계 분석"]
)

# ----------------------------------------------------------------------------
# TAB 0: 가설 개요
# ----------------------------------------------------------------------------
with tab0:
    st.subheader("가설의 기전")
    st.markdown(
        """
1. **에너지 과잉** 상태가 지속되면 지방세포가 지질을 저장하며 크기가 커진다(비대, hypertrophy).
2. 세포가 팽창하면서 **세포막·세포골격(F-actin)에 걸리는 기계적 장력이 증가**한다.
3. 높은 장력은 Hippo 경로의 핵심 억제 인산화효소인 **LATS1/2의 활성을 낮춘다**
   (LATS가 YAP/TAZ를 인산화하지 못하면 YAP/TAZ가 분해·세포질 격리에서 벗어난다).
4. 그 결과 **YAP/TAZ가 핵 안에 축적되어 과도하게 활성화**되고, CCN1(CYR61)·CCN2(CTGF) 같은
   표적 유전자의 전사가 늘어난다.
5. YAP/TAZ는 세포 증식·생존 신호를 강하게 유도하는 전사보조인자이므로, 이 경로가 만성적으로
   켜져 있으면 **비만과 역학적으로 연관된 암종(대장암, 자궁내막암, 신장암, 간암 등)**의
   발암·진행에 기여할 수 있다는 것이 이 탐구의 핵심 가설이다.
        """
    )

    st.subheader("두 축으로 나눈 검증 전략")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**축 1 — 비만 → YAP/TAZ (지방세포 단일세포 데이터)**")
        st.markdown(
            "- 데이터: Emont Lab 인체 지방세포 scRNA-seq (`human_adipocytes_lite.rds`)\n"
            "- 확인할 것: 비만군(Obese)이 정상군(Lean)보다 YAP/TAZ 표적 유전자 점수가 높은가?"
        )
    with c2:
        st.markdown("**축 2 — YAP/TAZ → 암 (TCGA 범암 발현 데이터)**")
        st.markdown(
            "- 데이터: GEPIA2 CTGF(CCN2)·CYR61(CCN1) Tumor vs Normal boxplot (TCGA 31개 암종)\n"
            "- 확인할 것: YAP/TAZ 표적 유전자가 종양에서 높게 발현되는 암종일수록, "
            "그 암종이 역학적으로 비만과 강하게 연관되어 있는가?"
        )

    st.info(
        "⚠️ **데이터 접근성에 대한 안내** — 업로드하신 두 GEPIA2 boxplot PDF는 산점도(scatter) "
        "이미지이기 때문에 텍스트로 변환하면 점 하나하나의 좌표(발현값)는 복원되지 않고, "
        "각 암종 아래 적힌 표본 수(num(T), num(N))만 정확히 읽을 수 있었습니다. 이 값은 "
        "②③ 탭 내 템플릿에 이미 반영해 두었습니다. 실제 발현값/log2FC는 GEPIA2 웹사이트에서 "
        "TSV로 내려받아 업로드해 주세요(방법은 ③ 탭 참고). `.rds`(Seurat) 파일도 Python이 R의 "
        "S4 객체를 직접 읽을 수 없어, ② 탭에 있는 R 스크립트로 CSV를 먼저 만들어야 합니다."
    )

    with st.expander("사용한 유전자 목록"):
        alias_str = ", ".join(f"{k}({v})" for k, v in YAP_TARGET_GENE_ALIASES.items())
        st.markdown(f"- **YAP/TAZ 표적 유전자**: {', '.join(YAP_TARGET_GENES)}  \n  (별칭: {alias_str})")
        st.markdown(f"- **세포골격 긴장 관련 유전자**: {', '.join(TENSION_GENES)}")

    with st.expander("참고문헌"):
        for c in CITATIONS:
            st.markdown(f"- {c}")

# ----------------------------------------------------------------------------
# TAB 1: 지방세포 데이터 (비만 → YAP/TAZ)
# ----------------------------------------------------------------------------
with tab1:
    st.subheader("지방세포 단일세포 데이터로 '비만 → YAP/TAZ 활성화' 검증")

    with st.expander("① R에서 CSV 내보내기 (human_adipocytes_lite.rds → CSV)", expanded=False):
        st.markdown("Python은 R의 Seurat(S4) 객체를 직접 읽지 못하므로, R에서 먼저 필요한 값만 CSV로 뽑아야 합니다.")
        st.code(
            """library(Seurat)
adipo <- readRDS("human_adipocytes_lite.rds")

# 메타데이터에서 비만 관련 컬럼 이름을 먼저 확인
colnames(adipo@meta.data)

# YAP/TAZ 표적 유전자 모듈 점수 계산
yap_genes <- c("CCN2","CCN1","ANKRD1","AMOTL2","BIRC5","AXL","F3","LATS2")
adipo <- AddModuleScore(adipo, features = list(yap_genes), name = "YAP_score", ctrl = 100)

# 방법 A: 점수만 내보내기 (앱의 '계산된 점수 업로드' 모드)
out <- data.frame(sample_id = colnames(adipo),
                   group = adipo$obesity_group,   # 실제 컬럼명으로 교체
                   yap_score = adipo$YAP_score1)
write.csv(out, "yap_score_by_sample.csv", row.names = FALSE)

# 방법 B: 원자료(발현행렬)를 통째로 내보내기 (앱의 '원자료 업로드' 모드)
expr <- as.data.frame(GetAssayData(adipo, slot = "data"))
write.csv(expr, "expression_matrix.csv")
meta <- data.frame(sample_id = colnames(adipo), group = adipo$obesity_group)
write.csv(meta, "sample_groups.csv", row.names = FALSE)
""",
            language="r",
        )

    mode = st.radio(
        "데이터 입력 방식",
        ["이미 계산된 YAP 점수 CSV 업로드", "원자료(발현행렬) 업로드 후 앱에서 점수 계산"],
        horizontal=True,
    )

    if mode == "이미 계산된 YAP 점수 CSV 업로드":
        st.caption("필요 컬럼: `sample_id`, `group`(예: Lean/Overweight/Obese), `yap_score`")
        sample_template = pd.DataFrame(
            {"sample_id": ["cell_1", "cell_2", "cell_3"], "group": ["Lean", "Overweight", "Obese"], "yap_score": [np.nan, np.nan, np.nan]}
        )
        st.download_button(
            "📥 입력 양식 CSV 내려받기",
            data=sample_template.to_csv(index=False).encode("utf-8-sig"),
            file_name="yap_score_template.csv",
            mime="text/csv",
        )
        f = st.file_uploader("YAP 점수 CSV 업로드", type=["csv"], key="score_csv")
        if f is not None:
            df = pd.read_csv(f)
            needed = {"group", "yap_score"}
            if not needed.issubset(df.columns):
                st.error(f"다음 컬럼이 필요합니다: {sorted(needed)}")
            else:
                fig = px.box(df, x="group", y="yap_score", points="all", color="group",
                             title="비만 그룹별 YAP/TAZ 표적 유전자 점수")
                st.plotly_chart(fig, use_container_width=True)
                res = mannwhitney_table(df, "yap_score", "group")
                if not res.empty:
                    st.markdown("**Mann-Whitney U 검정 (그룹 간 쌍대 비교)**")
                    st.dataframe(res, use_container_width=True)
                    sig = res[res["p-value"] < 0.05]
                    if len(sig):
                        st.success(
                            f"{len(sig)}개 그룹 쌍에서 p<0.05로 유의한 차이가 관찰됩니다 → "
                            "비만군에서 YAP/TAZ 점수가 유의하게 높다면 가설(축 1)을 지지하는 근거입니다."
                        )
                    else:
                        st.warning("그룹 간 통계적으로 유의한 차이가 확인되지 않았습니다.")

    else:
        st.caption("발현행렬 CSV: 행=유전자 기호(첫 컬럼), 열=샘플. 그룹 매핑 CSV: `sample_id`, `group`.")
        c1, c2 = st.columns(2)
        with c1:
            expr_f = st.file_uploader("발현행렬 CSV", type=["csv"], key="expr_csv")
        with c2:
            grp_f = st.file_uploader("샘플-그룹 매핑 CSV", type=["csv"], key="grp_csv")
        if expr_f is not None and grp_f is not None:
            expr = pd.read_csv(expr_f, index_col=0)
            grp = pd.read_csv(grp_f)
            yap_score, used_yap = compute_module_score(expr, YAP_TARGET_GENES)
            tension_score, used_tension = compute_module_score(expr, TENSION_GENES)
            if yap_score is None:
                st.error("업로드한 발현행렬에서 YAP/TAZ 표적 유전자를 하나도 찾지 못했습니다. 유전자 기호(symbol)를 확인해 주세요.")
            else:
                merged = grp.set_index("sample_id").copy()
                merged["yap_score"] = yap_score
                merged["tension_score"] = tension_score
                merged = merged.reset_index()
                st.caption(f"사용된 YAP/TAZ 유전자: {', '.join(used_yap)} / 세포골격 긴장 유전자: {', '.join(used_tension) if used_tension else '없음'}")

                fig1 = px.box(merged, x="group", y="yap_score", points="all", color="group",
                              title="비만 그룹별 YAP/TAZ 점수(z-score 평균, 단순화된 근사치)")
                st.plotly_chart(fig1, use_container_width=True)

                if used_tension:
                    fig2 = px.scatter(merged, x="tension_score", y="yap_score", color="group",
                                      trendline="ols",
                                      title="세포골격 긴장 점수 vs YAP/TAZ 점수")
                    st.plotly_chart(fig2, use_container_width=True)

                res = mannwhitney_table(merged, "yap_score", "group")
                if not res.empty:
                    st.markdown("**Mann-Whitney U 검정**")
                    st.dataframe(res, use_container_width=True)

# ----------------------------------------------------------------------------
# TAB 2: TCGA 범암 발현 데이터 (YAP/TAZ → 암)
# ----------------------------------------------------------------------------
with tab2:
    st.subheader("TCGA 범암 데이터로 'YAP/TAZ → 암' 검증")

    st.markdown(
        "GEPIA2(http://gepia2.cancer-pku.cn) → **Expression DIY → Box Plots**에서 "
        "`CCN2`(CTGF), `CCN1`(CYR61)을 각각 검색하고 *Match TCGA normal and GTEx data*를 "
        "체크한 뒤, 화면 하단의 **Differential Genes / Download** 버튼으로 암종별 Tumor·Normal "
        "median(또는 log2FC) 값을 내려받아 아래 템플릿에 채워 업로드하세요."
    )

    cohort_df = cohort_dataframe()
    template = cohort_df[["code", "name_kr", "n_tumor", "n_normal"]].copy()
    template["log2FC_CCN2"] = np.nan
    template["log2FC_CCN1"] = np.nan
    template["padj"] = np.nan

    st.download_button(
        "📥 암종별 입력 템플릿 CSV 내려받기 (표본 수는 실제 PDF 값 반영됨)",
        data=template.to_csv(index=False).encode("utf-8-sig"),
        file_name="gepia2_template.csv",
        mime="text/csv",
    )

    with st.expander("GEPIA2 실제 코호트 (두 PDF에서 확인된 표본 수)"):
        st.dataframe(cohort_df[["code", "name_kr", "name_en", "n_tumor", "n_normal"]],
                     use_container_width=True, hide_index=True)

    st.markdown("**아직 자신의 GEPIA2 자료가 없다면?** 실제 발표된 논문의 GEPIA2 분석 결과를 예시로 불러올 수 있습니다.")
    st.caption(
        f"출처: {DEMO_SOURCE}. 이 예시는 연속형 log2FC가 아니라 '유의하게 증가(+1)/감소(-1)'만 "
        "나타내는 방향성 지표이며, 논문에 언급되지 않은 암종은 결측(NaN)으로 둡니다."
    )
    if st.button("🔎 실제 논문 기반 예시 데이터 불러오기"):
        demo = cohort_df[["code", "name_kr", "n_tumor", "n_normal"]].copy()
        demo["log2FC_CCN2"] = demo["code"].apply(
            lambda c: 1.0 if c in CCN2_SIG_UP else (-1.0 if c in CCN2_SIG_DOWN else np.nan)
        )
        demo["log2FC_CCN1"] = demo["code"].apply(
            lambda c: 1.0 if c in CCN1_SIG_UP else (-1.0 if c in CCN1_SIG_DOWN else np.nan)
        )
        st.session_state["gepia_demo_df"] = demo

    f = st.file_uploader("채운 템플릿(또는 동일한 형식의 GEPIA2 결과 CSV) 업로드", type=["csv"], key="gepia_csv")

    df = None
    if f is not None:
        df = pd.read_csv(f)
    elif st.session_state.get("gepia_demo_df") is not None:
        df = st.session_state["gepia_demo_df"]
        st.info("현재 예시(데모) 데이터가 적용되어 있습니다. 직접 업로드하면 그 데이터로 교체됩니다.")

    if df is not None:
        if "code" not in df.columns:
            st.error("`code`(TCGA 암종 약어) 컬럼이 필요합니다.")
        else:
            gene_cols = [c for c in df.columns if c.startswith("log2FC_")]
            if not gene_cols:
                st.error("`log2FC_CCN2`, `log2FC_CYR61` 처럼 `log2FC_`로 시작하는 컬럼이 최소 1개 필요합니다.")
            else:
                merged = df.merge(cohort_df, on="code", how="left", suffixes=("", "_ref"))
                for gcol in gene_cols:
                    gene_label = gcol.replace("log2FC_", "")
                    plot_df = merged.dropna(subset=[gcol]).sort_values(gcol, ascending=False)
                    if plot_df.empty:
                        continue
                    fig = px.bar(
                        plot_df, x="code", y=gcol, color="obesity_linked",
                        color_discrete_map={True: "#d62728", False: "#7f7f7f"},
                        hover_data=["name_kr", "n_tumor", "n_normal", "evidence"],
                        title=f"{gene_label} — 암종별 Tumor vs Normal log2FC (빨강 = IARC 비만 연관 암종)",
                        labels={gcol: "log2FC (Tumor vs Normal)", "code": "TCGA 암종", "obesity_linked": "비만 연관(IARC)"},
                    )
                    st.plotly_chart(fig, use_container_width=True)

                st.session_state["gene_effect_df"] = merged
                st.success("④ 상관관계 분석 탭에서 이 데이터를 바로 사용할 수 있습니다.")

# ----------------------------------------------------------------------------
# TAB 3: 상관관계 분석
# ----------------------------------------------------------------------------
with tab3:
    st.subheader("YAP/TAZ 표적 유전자 발현과 비만-암 역학적 연관성의 상관관계")
    st.caption(
        "x축: BMI 5 kg/m² 증가당 상대위험도(RR, Renehan 2008 외) · "
        "y축: 해당 암종에서 YAP/TAZ 표적 유전자의 Tumor vs Normal log2FC (③ 탭에서 업로드한 값)"
    )

    merged = st.session_state.get("gene_effect_df")
    if merged is None:
        st.warning("먼저 ③ 탭에서 암종별 log2FC 데이터를 업로드해 주세요.")
    else:
        gene_cols = [c for c in merged.columns if c.startswith("log2FC_")]
        gene_choice = st.selectbox("분석할 유전자", gene_cols, format_func=lambda x: x.replace("log2FC_", ""))

        sub = merged.dropna(subset=[gene_choice, "rr_per_5bmi"])
        if len(sub) < 4:
            st.error("상관관계를 계산하려면 RR 값이 있는 암종 중 최소 4개 이상에서 log2FC 값이 필요합니다.")
        else:
            r_pearson, p_pearson = stats.pearsonr(sub["rr_per_5bmi"], sub[gene_choice])
            r_spearman, p_spearman = stats.spearmanr(sub["rr_per_5bmi"], sub[gene_choice])

            c1, c2 = st.columns(2)
            c1.metric("Pearson r", f"{r_pearson:.3f}", f"p = {p_pearson:.4f}")
            c2.metric("Spearman ρ", f"{r_spearman:.3f}", f"p = {p_spearman:.4f}")

            fig = px.scatter(
                sub, x="rr_per_5bmi", y=gene_choice, text="code", trendline="ols",
                hover_data=["name_kr", "evidence"],
                labels={"rr_per_5bmi": "비만-암 상대위험도 (RR, per +5 BMI)", gene_choice: f"{gene_choice.replace('log2FC_', '')} log2FC (T vs N)"},
                title="암종별 비만 연관 위험도 vs YAP/TAZ 표적 유전자 발현 변화",
            )
            fig.update_traces(textposition="top center")
            st.plotly_chart(fig, use_container_width=True)

            st.markdown("**비만 연관 여부(IARC 확정군 vs 미확립군) 간 log2FC 비교**")
            box_df = merged.dropna(subset=[gene_choice])
            fig2 = px.box(box_df, x="obesity_linked", y=gene_choice, points="all", color="obesity_linked",
                          color_discrete_map={True: "#d62728", False: "#7f7f7f"},
                          labels={"obesity_linked": "IARC 비만 연관 암종", gene_choice: f"{gene_choice.replace('log2FC_', '')} log2FC"})
            st.plotly_chart(fig2, use_container_width=True)
            grp_a = box_df.loc[box_df["obesity_linked"], gene_choice].dropna()
            grp_b = box_df.loc[~box_df["obesity_linked"], gene_choice].dropna()
            if len(grp_a) >= 2 and len(grp_b) >= 2:
                u_stat, u_p = stats.mannwhitneyu(grp_a, grp_b, alternative="two-sided")
                st.write(f"Mann-Whitney U 검정: U = {u_stat:.1f}, p = {u_p:.4f} "
                         f"(비만 연관군 n={len(grp_a)}, 미확립군 n={len(grp_b)})")

            st.markdown("### 결과 해석")
            direction = "양의" if r_pearson > 0 else "음의"
            sig_word = "통계적으로 유의한" if p_pearson < 0.05 else "통계적으로 유의하지는 않은"
            st.markdown(
                f"- 비만-암 상대위험도와 {gene_choice.replace('log2FC_', '')} log2FC 사이에는 "
                f"Pearson r = {r_pearson:.3f}({sig_word} {direction} 상관)이 관찰되었습니다.\n"
                "- **주의**: 이 상관관계는 (1) '암종 수준'의 집단 통계(ecological correlation)이지 "
                "개별 환자의 BMI와 유전자 발현을 직접 짝지은 것이 아니며, (2) TCGA의 Tumor vs Normal "
                "비교는 '비만 여부'가 아니라 '종양 조직 여부'를 비교한 것이므로, 이 결과만으로 "
                "'비만이 YAP/TAZ를 통해 암을 유발한다'는 인과관계를 증명할 수는 없습니다. "
                "가설을 뒷받침하는 정황 증거(association) 수준으로 해석해야 합니다."
            )

            st.dataframe(
                sub[["code", "name_kr", "n_tumor", "n_normal", "rr_per_5bmi", "evidence", gene_choice]]
                .sort_values(gene_choice, ascending=False),
                use_container_width=True, hide_index=True,
            )

st.divider()
st.caption("이 앱은 교육·탐구 목적의 도구이며, 임상적 진단·치료 판단에 사용할 수 없습니다.")
