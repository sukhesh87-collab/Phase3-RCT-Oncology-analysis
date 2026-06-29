## ============================================================
## PROJECT 2: QoL and Hemoglobin Secondary Analyses
## NCT00119613 — Darbepoetin Alfa in ES-SCLC
## ============================================================
## Objectives (per SAP Section 6.2 and 6.3):
## 1. Hemoglobin over time by arm (C_HEMAT) -- pharmacological
##    proof of activity: expect higher Hgb in darbepoetin arm
## 2. FACT-F (fatigue) and EQ-5D (utility) by arm and
##    timepoint (A_QOL) -- descriptive only, not powered
## ============================================================

library(haven)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)

setwd("Datasets/SAS files")

ADSL <- readRDS("ADSL.rds")

## ============================================================
## BLOCK 1: LOAD DATASETS AND INSPECT STRUCTURE
## ============================================================

C_HEMAT <- read_sas("C_HEMAT.sas7bdat") %>% clean_names()
A_QOL   <- read_sas("A_QOL.sas7bdat")   %>% clean_names()

cat("C_HEMAT loaded:", nrow(C_HEMAT), "rows\n")
cat("A_QOL loaded:", nrow(A_QOL), "rows\n")

## Always inspect first -- variable names must be confirmed
cat("\n====== C_HEMAT VARIABLES ======\n")
print(names(C_HEMAT))
cat("\n====== C_HEMAT STRUCTURE ======\n")
glimpse(C_HEMAT)

cat("\n====== A_QOL VARIABLES ======\n")
print(names(A_QOL))
cat("\n====== A_QOL STRUCTURE ======\n")
glimpse(A_QOL)

## ============================================================
## BLOCK 2: HEMOGLOBIN ANALYSIS (C_HEMAT)
## ============================================================
## Darbepoetin alfa is an erythropoiesis-stimulating agent --
## its primary pharmacological action is raising hemoglobin.
## We expect to see clearly higher Hgb in the NESP arm over
## time compared to placebo. If present, this confirms the
## drug was biologically active (pharmacodynamic proof),
## even though it didn't improve OS.
##
## Key variables expected in C_HEMAT:
## - txgroup: treatment arm (already in dataset like C_AE)
## - hgb or similar: hemoglobin value (g/dL)
## - week or cycle: visit/timepoint
## ============================================================

## Check if txgroup already in C_HEMAT (likely yes, same structure as C_AE)
if ("txgroup" %in% names(C_HEMAT)) {
  HEMAT <- C_HEMAT %>%
    mutate(ARM = case_when(
      txgroup == "NESP"    ~ "Darbepoetin Alfa",
      txgroup == "PLACEBO" ~ "Placebo",
      TRUE                 ~ as.character(txgroup)
    ))
  cat("\ntxgroup found in C_HEMAT -- using directly\n")
} else {
  HEMAT <- C_HEMAT %>%
    left_join(ADSL %>% select(pt, TRT_LABEL), by = "subjid") %>%
    rename(ARM = TRT_LABEL)
  cat("\nJoined treatment arm from ADSL\n")
}

## Identify hemoglobin variable
## Common names: hgb, hgbval, hemog, hgbrs
hgb_candidates <- names(C_HEMAT)[grepl(
  "^hgb$|hgbval|hgbrs|hemog|haemog|hgb_",
  names(C_HEMAT), ignore.case = TRUE
)]
cat("\nHemoglobin variable candidates:", hgb_candidates, "\n")

## Print summary of each candidate to find the right one
for (v in hgb_candidates) {
  vals <- suppressWarnings(as.numeric(C_HEMAT[[v]]))
  cat(sprintf("  %s: range %.1f - %.1f, mean %.1f\n",
              v, min(vals, na.rm=TRUE),
              max(vals, na.rm=TRUE),
              mean(vals, na.rm=TRUE)))
}
## Hemoglobin in g/dL should range roughly 7-15

## Identify visit/timepoint variable
## Common names: week, cycle, visitwk, visitnum
visit_candidates <- names(C_HEMAT)[grepl(
  "^week$|^cycle$|visitwk|visitnum|^wk$",
  names(C_HEMAT), ignore.case = TRUE
)]
cat("\nVisit/timepoint candidates:", visit_candidates, "\n")
for (v in visit_candidates) {
  cat(sprintf("  %s: unique values: %s\n",
              v, paste(sort(unique(C_HEMAT[[v]])), collapse=", ")))
}

## ---- UPDATE THESE TWO LINES based on output above ----
hgb_var   <- hgb_candidates[1]    ## e.g., "hgb"
visit_var <- visit_candidates[1]  ## e.g., "week"
## -------------------------------------------------------

cat("\nUsing Hgb variable:", hgb_var, "\n")
cat("Using visit variable:", visit_var, "\n")

## Convert Hgb to numeric (may be stored as character)
HEMAT <- HEMAT %>%
  mutate(
    hgb_val   = suppressWarnings(as.numeric(.data[[hgb_var]])),
    visit_num = suppressWarnings(as.numeric(.data[[visit_var]]))
  ) %>%
  filter(!is.na(hgb_val) & !is.na(visit_num) & !is.na(ARM))

## Summarize mean Hgb by arm and visit
hgb_summary <- HEMAT %>%
  group_by(ARM, Visit = visit_num) %>%
  summarize(
    N       = n(),
    Mean    = round(mean(hgb_val, na.rm=TRUE), 2),
    SD      = round(sd(hgb_val, na.rm=TRUE), 2),
    SE      = round(sd(hgb_val, na.rm=TRUE) / sqrt(n()), 3),
    CI_Low  = round(Mean - 1.96 * SE, 2),
    CI_High = round(Mean + 1.96 * SE, 2),
    .groups = "drop"
  ) %>%
  filter(N >= 5) %>%  ## exclude timepoints with very few patients
  arrange(ARM, Visit)

cat("\n====== HEMOGLOBIN BY ARM AND VISIT ======\n")
print(hgb_summary, n = 60)

write.csv(hgb_summary, "Table5_Hemoglobin.csv", row.names = FALSE)
cat("\nTable5_Hemoglobin.csv saved\n")

## ============================================================
## BLOCK 3: FIGURE 3 — MEAN HEMOGLOBIN OVER TIME
## ============================================================
## Expected pattern: Darbepoetin arm rises above placebo
## within first 2-3 cycles and stays elevated throughout
## treatment. This is the pharmacodynamic proof-of-activity
## figure -- the most pharmacologically interpretable result.
## ============================================================

fig3 <- ggplot(
  hgb_summary %>% filter(!is.na(ARM)),
  aes(x     = Visit,
      y     = Mean,
      color = ARM,
      fill  = ARM)
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2.8) +
  geom_ribbon(
    aes(ymin = CI_Low, ymax = CI_High),
    alpha = 0.15, color = NA
  ) +
  scale_color_manual(
    values = c("Darbepoetin Alfa" = "#C0504D",
               "Placebo"          = "#1B3A5C")
  ) +
  scale_fill_manual(
    values = c("Darbepoetin Alfa" = "#C0504D",
               "Placebo"          = "#1B3A5C")
  ) +
  labs(
    title = "Mean Hemoglobin Over Time by Treatment Arm\nSafety Population — NCT00119613",
    x     = "Study Week / Cycle",
    y     = "Mean Hemoglobin (g/dL)",
    color = "",
    fill  = ""
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("Figure3_Hemoglobin.png", fig3,
       width = 10, height = 6, dpi = 180)
cat("\nFigure3_Hemoglobin.png saved\n")

## ============================================================
## BLOCK 4: QOL ANALYSIS (A_QOL)
## ============================================================
## Quality of life instruments:
## FACT-F = fatigue subscale (higher = less fatigue, better)
## FACT-G = general QoL (higher = better overall QoL)
## FACT-A = anemia-specific subscale
## EQ-5D  = health utility (0=death, 1=perfect health)
##          Used in HEOR cost-effectiveness models
##
## Analysis is descriptive only (not powered per SAP).
## We report mean scores at baseline and end of treatment.
##
## Join using 'pt' if A_QOL uses same ID as ADSL.
## If A_QOL uses 'subjid', adjust join key below.

## A_QOL is a WIDE dataset: one row per patient, each QoL
## score at each timepoint is a separate column.
##
## Variable naming convention: {domain}{timepoint}
## Domains:
##   fat   = FACT-F Fatigue total (primary QoL endpoint)
##   fgt   = FACT-G Total overall QoL
##   hsi   = Hemoglobin Symptom Index (ESA-specific)
##   vas   = EQ-5D VAS (0-100 visual analogue scale)
##   eqv   = EQ-5D utility (eqv0=baseline, eqve=end only)
##   fan   = FACT-An Anemia subscale
##   lcs   = Lung Cancer Subscale
## Timepoints:
##   b/0   = baseline
##   w7    = week 7
##   w13   = week 13
##   w24   = week 24
##   w50   = week 50
##   eos/e = end of study
##   chgN  = change from baseline at week N
## ============================================================

ADSL <- readRDS("ADSL.rds")
A_QOL <- read_sas("A_QOL.sas7bdat") %>% clean_names()

cat("A_QOL loaded:", nrow(A_QOL), "patients (wide format)\n")
cat("Columns:", ncol(A_QOL), "\n")

## Add clean treatment arm label (txgroup already in A_QOL)
A_QOL <- A_QOL %>%
  mutate(
    ARM = case_when(
      txgroup == "NESP"    ~ "Darbepoetin Alfa",
      txgroup == "PLACEBO" ~ "Placebo",
      TRUE                 ~ as.character(txgroup)
    )
  )

cat("\nPatients by arm:\n")
print(table(A_QOL$ARM))

## ============================================================
## BLOCK 5: PIVOT FACT-F (FATIGUE) TO LONG FORMAT
## ============================================================
## fatb   = baseline FACT-F score
## fatw7  = week 7 FACT-F score
## fatw13 = week 13 FACT-F score
## fatw24 = week 24 FACT-F score
## fatw50 = week 50 FACT-F score
## fateos = end of study FACT-F score
##
## FACT-F scoring: higher score = LESS fatigue (better)
## Range: 0-52
##
## pivot_longer() converts wide → long:
## Each fat* column becomes a separate row with:
##   - 'Timepoint' identifying which visit
##   - 'FACTF_Score' containing the value
## ============================================================

factf_long <- A_QOL %>%
  select(subjid, ARM,
         fatb, fatw7, fatw13, fatw24, fatw50, fateos) %>%
  pivot_longer(
    cols      = c(fatb, fatw7, fatw13, fatw24, fatw50, fateos),
    names_to  = "Timepoint_Raw",
    values_to = "FACTF_Score"
  ) %>%
  ## Convert raw column name to readable timepoint label
  mutate(
    Timepoint = case_when(
      Timepoint_Raw == "fatb"   ~ "Baseline",
      Timepoint_Raw == "fatw7"  ~ "Week 7",
      Timepoint_Raw == "fatw13" ~ "Week 13",
      Timepoint_Raw == "fatw24" ~ "Week 24",
      Timepoint_Raw == "fatw50" ~ "Week 50",
      Timepoint_Raw == "fateos" ~ "End of Study"
    ),
    ## Numeric order for plotting
    Visit_Num = case_when(
      Timepoint_Raw == "fatb"   ~ 0,
      Timepoint_Raw == "fatw7"  ~ 7,
      Timepoint_Raw == "fatw13" ~ 13,
      Timepoint_Raw == "fatw24" ~ 24,
      Timepoint_Raw == "fatw50" ~ 50,
      Timepoint_Raw == "fateos" ~ 55  ## place after week 50
    ),
    Timepoint = factor(Timepoint,
                       levels = c("Baseline","Week 7","Week 13",
                                  "Week 24","Week 50","End of Study"))
  ) %>%
  filter(!is.na(FACTF_Score))

cat("\nFACT-F long format: ", nrow(factf_long), "records\n")
cat("Records by timepoint:\n")
print(table(factf_long$Timepoint))

## Summary table: mean (SD) by arm and timepoint
factf_summary <- factf_long %>%
  group_by(ARM, Timepoint, Visit_Num) %>%
  summarize(
    N         = n(),
    Mean      = round(mean(FACTF_Score, na.rm=TRUE), 1),
    SD        = round(sd(FACTF_Score, na.rm=TRUE), 1),
    SE        = round(sd(FACTF_Score, na.rm=TRUE)/sqrt(n()), 2),
    CI_Low    = round(Mean - 1.96*SE, 1),
    CI_High   = round(Mean + 1.96*SE, 1),
    .groups   = "drop"
  ) %>%
  arrange(ARM, Visit_Num)

cat("\n====== TABLE 6A: FACT-F FATIGUE SCORE (Higher = Less Fatigue) ======\n")
print(factf_summary %>% select(-Visit_Num, -SE, -CI_Low, -CI_High),
      row.names=FALSE)

## ============================================================
## BLOCK 6: PIVOT FACT-G TOTAL TO LONG FORMAT
## ============================================================
## fgtb, fgtw7, fgtw13, fgtw24, fgtw50, fgteos
## FACT-G Total: overall quality of life composite score
## Higher = better QoL. Range 0-108.
## ============================================================

factg_long <- A_QOL %>%
  select(subjid, ARM,
         fgtb, fgtw7, fgtw13, fgtw24, fgtw50, fgteos) %>%
  pivot_longer(
    cols      = c(fgtb, fgtw7, fgtw13, fgtw24, fgtw50, fgteos),
    names_to  = "Timepoint_Raw",
    values_to = "FACTG_Score"
  ) %>%
  mutate(
    Timepoint = case_when(
      Timepoint_Raw == "fgtb"   ~ "Baseline",
      Timepoint_Raw == "fgtw7"  ~ "Week 7",
      Timepoint_Raw == "fgtw13" ~ "Week 13",
      Timepoint_Raw == "fgtw24" ~ "Week 24",
      Timepoint_Raw == "fgtw50" ~ "Week 50",
      Timepoint_Raw == "fgteos" ~ "End of Study"
    ),
    Visit_Num = case_when(
      Timepoint_Raw == "fgtb"   ~ 0,
      Timepoint_Raw == "fgtw7"  ~ 7,
      Timepoint_Raw == "fgtw13" ~ 13,
      Timepoint_Raw == "fgtw24" ~ 24,
      Timepoint_Raw == "fgtw50" ~ 50,
      Timepoint_Raw == "fgteos" ~ 55
    ),
    Timepoint = factor(Timepoint,
                       levels = c("Baseline","Week 7","Week 13",
                                  "Week 24","Week 50","End of Study"))
  ) %>%
  filter(!is.na(FACTG_Score))

factg_summary <- factg_long %>%
  group_by(ARM, Timepoint, Visit_Num) %>%
  summarize(
    N     = n(),
    Mean  = round(mean(FACTG_Score, na.rm=TRUE), 1),
    SD    = round(sd(FACTG_Score, na.rm=TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(ARM, Visit_Num)

cat("\n====== TABLE 6B: FACT-G TOTAL (Higher = Better QoL) ======\n")
print(factg_summary %>% select(-Visit_Num), row.names=FALSE)

## ============================================================
## BLOCK 7: HEMOGLOBIN SYMPTOM INDEX (HSI)
## ============================================================
## hsib, hsiw7, hsiw13, hsiw24, hsiw50, hsieos
## HSI specifically measures symptoms driven by hemoglobin
## level -- designed specifically for ESA oncology trials.
## Higher = fewer Hgb-related symptoms (better).
## This is arguably the most relevant PRO for this drug.
## ============================================================

hsi_long <- A_QOL %>%
  select(subjid, ARM,
         hsib, hsiw7, hsiw13, hsiw24, hsiw50, hsieos) %>%
  pivot_longer(
    cols      = c(hsib, hsiw7, hsiw13, hsiw24, hsiw50, hsieos),
    names_to  = "Timepoint_Raw",
    values_to = "HSI_Score"
  ) %>%
  mutate(
    Timepoint = case_when(
      Timepoint_Raw == "hsib"   ~ "Baseline",
      Timepoint_Raw == "hsiw7"  ~ "Week 7",
      Timepoint_Raw == "hsiw13" ~ "Week 13",
      Timepoint_Raw == "hsiw24" ~ "Week 24",
      Timepoint_Raw == "hsiw50" ~ "Week 50",
      Timepoint_Raw == "hsieos" ~ "End of Study"
    ),
    Visit_Num = case_when(
      Timepoint_Raw == "hsib"   ~ 0,
      Timepoint_Raw == "hsiw7"  ~ 7,
      Timepoint_Raw == "hsiw13" ~ 13,
      Timepoint_Raw == "hsiw24" ~ 24,
      Timepoint_Raw == "hsiw50" ~ 50,
      Timepoint_Raw == "hsieos" ~ 55
    ),
    Timepoint = factor(Timepoint,
                       levels = c("Baseline","Week 7","Week 13",
                                  "Week 24","Week 50","End of Study"))
  ) %>%
  filter(!is.na(HSI_Score))

hsi_summary <- hsi_long %>%
  group_by(ARM, Timepoint, Visit_Num) %>%
  summarize(
    N    = n(),
    Mean = round(mean(HSI_Score, na.rm=TRUE), 1),
    SD   = round(sd(HSI_Score, na.rm=TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(ARM, Visit_Num)

cat("\n====== TABLE 6C: HEMOGLOBIN SYMPTOM INDEX (Higher = Fewer Symptoms) ======\n")
print(hsi_summary %>% select(-Visit_Num), row.names=FALSE)

## ============================================================
## BLOCK 8: EQ-5D VAS (Visual Analogue Scale)
## ============================================================
## vasb, vasw7, vasw13, vasw24, vasw50, vaseos
## EQ-5D VAS: 0-100 scale, 100=best imaginable health
## This is the component used in HEOR/HTA utilities.
## Note: eqv0/eqve = EQ-5D utility index (only 2 timepoints)
## ============================================================

vas_long <- A_QOL %>%
  select(subjid, ARM,
         vasb, vasw7, vasw13, vasw24, vasw50, vaseos) %>%
  pivot_longer(
    cols      = c(vasb, vasw7, vasw13, vasw24, vasw50, vaseos),
    names_to  = "Timepoint_Raw",
    values_to = "VAS_Score"
  ) %>%
  mutate(
    Timepoint = case_when(
      Timepoint_Raw == "vasb"   ~ "Baseline",
      Timepoint_Raw == "vasw7"  ~ "Week 7",
      Timepoint_Raw == "vasw13" ~ "Week 13",
      Timepoint_Raw == "vasw24" ~ "Week 24",
      Timepoint_Raw == "vasw50" ~ "Week 50",
      Timepoint_Raw == "vaseos" ~ "End of Study"
    ),
    Visit_Num = case_when(
      Timepoint_Raw == "vasb"   ~ 0,
      Timepoint_Raw == "vasw7"  ~ 7,
      Timepoint_Raw == "vasw13" ~ 13,
      Timepoint_Raw == "vasw24" ~ 24,
      Timepoint_Raw == "vasw50" ~ 50,
      Timepoint_Raw == "vaseos" ~ 55
    ),
    Timepoint = factor(Timepoint,
                       levels = c("Baseline","Week 7","Week 13",
                                  "Week 24","Week 50","End of Study"))
  ) %>%
  filter(!is.na(VAS_Score))

vas_summary <- vas_long %>%
  group_by(ARM, Timepoint, Visit_Num) %>%
  summarize(
    N    = n(),
    Mean = round(mean(VAS_Score, na.rm=TRUE), 1),
    SD   = round(sd(VAS_Score, na.rm=TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(ARM, Visit_Num)

cat("\n====== TABLE 6D: EQ-5D VAS (0-100, Higher = Better Health) ======\n")
print(vas_summary %>% select(-Visit_Num), row.names=FALSE)

## EQ-5D Utility Index (only baseline vs end of study)
eq5d_summary <- A_QOL %>%
  select(subjid, ARM, eqv0, eqve) %>%
  filter(!is.na(eqv0) | !is.na(eqve)) %>%
  pivot_longer(
    cols      = c(eqv0, eqve),
    names_to  = "Timepoint_Raw",
    values_to = "EQ5D_Utility"
  ) %>%
  mutate(Timepoint = case_when(
    Timepoint_Raw == "eqv0" ~ "Baseline",
    Timepoint_Raw == "eqve" ~ "End of Study"
  )) %>%
  filter(!is.na(EQ5D_Utility)) %>%
  group_by(ARM, Timepoint) %>%
  summarize(
    N    = n(),
    Mean = round(mean(EQ5D_Utility, na.rm=TRUE), 3),
    SD   = round(sd(EQ5D_Utility, na.rm=TRUE), 3),
    .groups = "drop"
  )

cat("\n====== TABLE 6E: EQ-5D UTILITY INDEX (0=Death, 1=Perfect Health) ======\n")
cat("HEOR NOTE: Used in cost-effectiveness and HTA submissions\n")
print(eq5d_summary, row.names=FALSE)

## ============================================================
## BLOCK 9: CHANGE FROM BASELINE ANALYSIS
## ============================================================
## The chg variables (fatchg7, fatchg13, etc.) give us
## change from baseline directly -- positive = improvement
## for all FACT scores (higher = better).
## This is the standard clinical meaningfulness assessment.
## ============================================================

## FACT-F change from baseline
factf_chg <- A_QOL %>%
  select(subjid, ARM,
         fatchg7, fatchg13, fatchg24, fatchg50) %>%
  pivot_longer(
    cols      = c(fatchg7, fatchg13, fatchg24, fatchg50),
    names_to  = "Timepoint_Raw",
    values_to = "FACTF_Change"
  ) %>%
  mutate(
    Timepoint = case_when(
      Timepoint_Raw == "fatchg7"  ~ "Week 7",
      Timepoint_Raw == "fatchg13" ~ "Week 13",
      Timepoint_Raw == "fatchg24" ~ "Week 24",
      Timepoint_Raw == "fatchg50" ~ "Week 50"
    ),
    Visit_Num = case_when(
      Timepoint_Raw == "fatchg7"  ~ 7,
      Timepoint_Raw == "fatchg13" ~ 13,
      Timepoint_Raw == "fatchg24" ~ 24,
      Timepoint_Raw == "fatchg50" ~ 50
    )
  ) %>%
  filter(!is.na(FACTF_Change)) %>%
  group_by(ARM, Timepoint, Visit_Num) %>%
  summarize(
    N          = n(),
    Mean_Chg   = round(mean(FACTF_Change, na.rm=TRUE), 1),
    SD         = round(sd(FACTF_Change, na.rm=TRUE), 1),
    .groups    = "drop"
  ) %>%
  arrange(ARM, Visit_Num)

cat("\n====== FACT-F CHANGE FROM BASELINE (Positive = Improvement) ======\n")
print(factf_chg %>% select(-Visit_Num), row.names=FALSE)

## ============================================================
## BLOCK 10: FIGURE 4 — FACT-F FATIGUE OVER TIME
## ============================================================

fig4 <- ggplot(
  factf_summary %>% filter(!is.na(ARM)),
  aes(x     = Visit_Num,
      y     = Mean,
      color = ARM,
      fill  = ARM,
      group = ARM)
) +
  geom_ribbon(
    aes(ymin = Mean - SD/sqrt(N),
        ymax = Mean + SD/sqrt(N)),
    alpha = 0.15, color = NA
  ) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Darbepoetin Alfa" = "#C0504D",
               "Placebo"         = "#1B3A5C")
  ) +
  scale_fill_manual(
    values = c("Darbepoetin Alfa" = "#C0504D",
               "Placebo"         = "#1B3A5C")
  ) +
  scale_x_continuous(
    breaks = c(0, 7, 13, 24, 50, 55),
    labels = c("Baseline","Wk 7","Wk 13","Wk 24","Wk 50","EOS")
  ) +
  labs(
    title    = "FACT-F Fatigue Score Over Time by Treatment Arm\nITT Population — NCT00119613",
    subtitle = "Higher score = less fatigue (better). Shaded band = ±1 SE.",
    x        = "Study Timepoint",
    y        = "Mean FACT-F Score",
    color    = "",
    fill     = ""
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("Figure4_FACTF_Fatigue.png", fig4,
       width = 10, height = 6, dpi = 180)
cat("\nFigure4_FACTF_Fatigue.png saved\n")

## ============================================================
## BLOCK 11: FIGURE 5 — EQ-5D VAS OVER TIME
## ============================================================
## EQ-5D VAS is the HEOR-specific utility measure.
## Plotting separately because it has different scale (0-100)
## and is the metric used in HTA cost-effectiveness models.
## ============================================================

fig5 <- ggplot(
  vas_summary %>% filter(!is.na(ARM)),
  aes(x     = Visit_Num,
      y     = Mean,
      color = ARM,
      group = ARM)
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = Mean - SD/sqrt(N),
        ymax = Mean + SD/sqrt(N)),
    width = 1.5, linewidth = 0.8
  ) +
  scale_color_manual(
    values = c("Darbepoetin Alfa" = "#C0504D",
               "Placebo"         = "#1B3A5C")
  ) +
  scale_x_continuous(
    breaks = c(0, 7, 13, 24, 50, 55),
    labels = c("Baseline","Wk 7","Wk 13","Wk 24","Wk 50","EOS")
  ) +
  labs(
    title    = "EQ-5D VAS Score Over Time by Treatment Arm\nITT Population — NCT00119613",
    subtitle = "0 = worst imaginable health, 100 = best imaginable health. Error bars = ±1 SE.",
    x        = "Study Timepoint",
    y        = "Mean EQ-5D VAS Score",
    color    = ""
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("Figure5_EQ5D_VAS.png", fig5,
       width = 10, height = 6, dpi = 180)
cat("\nFigure5_EQ5D_VAS.png saved\n")

## ============================================================
## BLOCK 12: SAVE ALL QOL TABLES
## ============================================================

write.csv(factf_summary %>% select(-Visit_Num, -SE, -CI_Low, -CI_High),
          "Table6a_FACTF.csv",    row.names=FALSE)
write.csv(factg_summary %>% select(-Visit_Num),
          "Table6b_FACTG.csv",    row.names=FALSE)
write.csv(hsi_summary %>% select(-Visit_Num),
          "Table6c_HSI.csv",      row.names=FALSE)
write.csv(vas_summary %>% select(-Visit_Num),
          "Table6d_VAS.csv",      row.names=FALSE)
write.csv(eq5d_summary,
          "Table6e_EQ5D.csv",     row.names=FALSE)
write.csv(factf_chg %>% select(-Visit_Num),
          "Table6f_FACTF_Chg.csv",row.names=FALSE)