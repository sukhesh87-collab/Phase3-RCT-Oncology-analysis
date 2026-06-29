## ============================================================
## PROJECT 2: Phase 3 RCT Analysis — NCT00119613
## Step 2: ADTTE Construction + Primary OS Analysis
## ============================================================
## Pre-requisite: ADSL.rds must exist in your working directory
## (produced by Project2_01_ADSL.R)
## ============================================================

library(haven)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(ggplot2)
library(gtsummary)
library(tableone)
library(janitor)

set.seed(2026)  # reproducibility per SAP

## ============================================================
## BLOCK 1: RELOAD ADSL AND CONFIRM
## ============================================================

setwd("Datasets/SAS files")

ADSL <- readRDS("ADSL.rds")

cat("ADSL reloaded:", nrow(ADSL), "patients,", ncol(ADSL), "variables\n")
cat("ITT population:", sum(ADSL$ITTFL == "Y"), "patients\n")
cat("Treatment arms:\n")
print(table(ADSL$TRT_LABEL))

## ============================================================
## BLOCK 2: BUILD ADTTE (Time-to-Event Analysis Dataset)
## ============================================================
## What is ADTTE?
## ADTTE is the ADaM dataset specifically for time-to-event
## endpoints. It contains one row per patient per endpoint.
## Since we have one primary endpoint (OS), it's one row
## per patient here.
##
## Key ADTTE variables per CDISC ADaM convention:
## USUBJID  = unique subject identifier
## PARAM    = endpoint name ("Overall Survival")
## PARAMCD  = short code ("OS")
## AVAL     = analysis value = time in days (our OS_DAYS)
## CNSR     = censoring indicator
##            0 = event occurred (death)
##            1 = censored (alive at last contact)
##            NOTE: this is opposite of OS_EVENT (1=death)
##            The survival package uses 1=event in Surv()
##            so we use OS_EVENT directly there, not CNSR
## EVNTDESC = description of event or censoring reason
## STARTDT  = start date (randomization)
##
## Why build ADTTE separately from ADSL?
## ADSL has one row per patient with OS variables embedded.
## ADTTE follows the ADaM long format -- one row per patient
## per endpoint. When we have multiple endpoints (OS + PFS
## + TTPR), ADTTE expands to multiple rows per patient.
## Keeping it separate is the CDISC-standard discipline.
## ============================================================

ADTTE <- ADSL %>%
  filter(ITTFL == "Y") %>%   # ITT population only for primary analysis
  mutate(
    PARAM    = "Overall Survival",
    PARAMCD  = "OS",
    AVAL     = OS_DAYS,       # time in days (our primary analysis time variable)
    CNSR     = OS_CNSR,       # 0=event(death), 1=censored -- CDISC convention
    AVALC    = case_when(     # human-readable event description
      OS_EVENT == 1 ~ "Death",
      OS_EVENT == 0 ~ "Censored - Alive at Last Contact"
    )
  ) %>%
  select(subjid, TRT_LABEL, TRT01P, ITTFL, SAFFL,
         PARAM, PARAMCD, AVAL, CNSR, OS_EVENT,
         OS_DAYS, OS_MONTHS, AVALC, age.y, sex.y, b_ecog2.y)

cat("\nADTTE created:", nrow(ADTTE), "rows (one per patient)\n")
cat("Events (deaths):", sum(ADTTE$OS_EVENT == 1, na.rm=TRUE), "\n")
cat("Censored:", sum(ADTTE$OS_EVENT == 0, na.rm=TRUE), "\n")
cat("Overall event rate:", round(100*mean(ADTTE$OS_EVENT, na.rm=TRUE),1), "%\n")

## Quick check by arm -- raw event rates before any analysis
cat("\nEvents by treatment arm:\n")
print(
  ADTTE %>%
    group_by(TRT_LABEL) %>%
    summarize(
      N = n(),
      Deaths = sum(OS_EVENT, na.rm=TRUE),
      Event_Rate_Pct = round(100*mean(OS_EVENT, na.rm=TRUE),1),
      Median_OS_Days = median(OS_DAYS, na.rm=TRUE),
      Median_OS_Months = round(median(OS_DAYS, na.rm=TRUE)/30.4375, 1),
      .groups = "drop"
    )
)

saveRDS(ADTTE, "ADTTE.rds")
cat("\nADTTE saved as ADTTE.rds\n")

## ============================================================
## BLOCK 3: TABLE 1 — BASELINE CHARACTERISTICS (ITT)
## ============================================================
## What is Table 1?
## Table 1 is the standard first table in any clinical trial
## or RWE paper. It shows baseline characteristics of patients
## by treatment arm -- confirms that randomization achieved
## balance (in an RCT, we expect balance by design, unlike
## RWE where we had to CREATE balance via matching).
##
## In an RCT, you generally do NOT show p-values in Table 1
## (unlike observational studies). Why? Because any imbalance
## is purely due to chance -- the trial was randomized --
## and a p-value testing "are these groups different" is
## philosophically wrong in a randomized dataset. This is an
## ICH E9 recommendation. We'll use SMD instead (same as
## Project 1) just for balance assessment.
##
## ============================================================

## Define variables for Table 1
table1_vars <- c("age.y", "sex.y", "b_ecog2.y")

## Build Table 1 using tbl_summary
## by = splits by treatment arm
## statistic = format: continuous as mean(SD), categorical as n(%)
## missing = "no" suppresses missing row (adjust if you have missing)
table1 <- ADSL %>%
  filter(ITTFL == "Y") %>%
  select(TRT_LABEL, age.y, sex.y, b_ecog2.y) %>%
  tbl_summary(
    by = TRT_LABEL,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age.y  ~ "Age, years",
      sex.y  ~ "Sex",
      b_ecog2.y ~ "ECOG Performance Status"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_n() %>%
  bold_labels()

## Print the table in Viewer panel
table1

## Save as CSV for the results report
table1_df <- as.data.frame(table1)
write.csv(table1_df, "Table1_Baseline.csv", row.names = FALSE)
cat("\nTable 1 saved as Table1_Baseline.csv\n")

## ============================================================
## BLOCK 4: KAPLAN-MEIER ANALYSIS
## ============================================================
## What is the Kaplan-Meier estimator?
## KM is a non-parametric method for estimating the survival
## function S(t) = probability of surviving beyond time t.
## "Non-parametric" means it makes no assumptions about the
## shape of the survival distribution -- it just steps down
## at each observed event time, based purely on the data.
##
## How does it work conceptually?
## Start at S(0) = 1.0 (everyone alive at time zero).
## At each event time (death), it calculates:
## S(t) = S(t-1) × (1 - events_at_t / at_risk_at_t)
## This accounts for censored patients appropriately --
## they contribute follow-up time while they're observable,
## then "leave" the risk set when censored, without being
## counted as events.
##
## The number-at-risk table below the KM curve shows how many
## patients remain in follow-up at each time point -- this is
## required in all pharma KM presentations because it lets
## readers judge the reliability of the curve at late time
## points (where few patients remain, curves become unstable).
##
## What is Surv()?
## Surv(time, event) creates a "survival object" -- R's way
## of packaging time + event indicator together for survival
## analysis functions. event must be 1=event, 0=censored.
## ============================================================

## Fit KM curves by treatment arm (ITT population)
km_fit <- survfit(
  Surv(OS_DAYS, OS_EVENT) ~ TRT_LABEL,
  data = ADTTE
)

## Print KM summary -- median OS with 95% CI for each arm
cat("\n====== KAPLAN-MEIER SUMMARY ======\n")
print(km_fit)

## More detailed output: survival at specific time points
cat("\nSurvival probabilities at key time points:\n")
print(summary(km_fit, times = c(90, 180, 270, 365, 547)))
## 90=3mo, 180=6mo, 270=9mo, 365=12mo, 547=18mo

## ============================================================
## BLOCK 5: LOG-RANK TEST
## ============================================================
## What is the log-rank test?
## It's the standard statistical test for comparing survival
## curves between two groups in a clinical trial.
## Null hypothesis: the survival curves are identical in both
## arms (i.e., the drug has no effect on OS).
## It counts observed vs. expected deaths in each group at
## every event time, then aggregates those differences into
## a chi-squared test statistic.
## A p-value < 0.05 means the curves are significantly different.
##
## Why log-rank specifically (not a t-test)?
## Because survival times are right-skewed (many short,
## few very long) and include censored observations --
## neither of which standard t-tests handle correctly.
## Log-rank is specifically designed for time-to-event data
## with censoring.
## ============================================================

logrank_test <- survdiff(
  Surv(OS_DAYS, OS_EVENT) ~ TRT_LABEL,
  data = ADTTE
)

cat("\n====== LOG-RANK TEST ======\n")
print(logrank_test)

## Extract p-value manually for reporting
logrank_pval <- 1 - pchisq(logrank_test$chisq, df = 1)
cat("\nLog-rank p-value:", round(logrank_pval, 4), "\n")

## ============================================================
## BLOCK 6: COX PROPORTIONAL HAZARDS MODEL
## ============================================================
## What does Cox PH add beyond the log-rank test?
## The log-rank test only tells you IF the curves differ
## (p-value). The Cox model tells you HOW MUCH they differ --
## the Hazard Ratio (HR) with 95% CI.
##
## In an RCT, the Cox model is usually simple -- just treatment
## arm as the predictor, because randomization already balanced
## all confounders. We don't need propensity scores (unlike
## RWE sudies where we have observational confounding).
##
## HR interpretation for this trial:
## HR < 1.0 = darbepoetin alfa associated with LOWER death risk
##            (longer survival -- favorable for the drug)
## HR > 1.0 = darbepoetin alfa associated with HIGHER death risk
##            (shorter survival -- potentially harmful)
## HR = 1.0 = no difference between arms
##
## ============================================================

cox_model <- coxph(
  Surv(OS_DAYS, OS_EVENT) ~ TRT_LABEL,
  data = ADTTE
)

cat("\n====== COX PROPORTIONAL HAZARDS MODEL ======\n")
print(summary(cox_model))

## Extract key results
hr  <- exp(coef(cox_model))
ci  <- exp(confint(cox_model))
pval <- summary(cox_model)$coefficients[,"Pr(>|z|)"]

cat("\n=================== PRIMARY OS RESULT ===================\n")
cat("Hazard Ratio (Darbepoetin vs Placebo):", round(hr, 3), "\n")
cat("95% Confidence Interval:", round(ci[1],3), "to", round(ci[2],3), "\n")
cat("Cox model p-value:", round(pval, 4), "\n")
cat("Log-rank p-value:", round(logrank_pval, 4), "\n")
cat("==========================================================\n")

## ============================================================
## BLOCK 7: PROPORTIONAL HAZARDS ASSUMPTION CHECK
## ============================================================
## Schoenfeld residuals test.
## In an RCT context this matters less (the estimate is valid
## even if PH is violated -- we just interpret it as an
## average HR over follow-up rather than a constant HR),
## but standard good practice to always report it.
## ============================================================

ph_test <- cox.zph(cox_model)
cat("\n====== PH ASSUMPTION TEST (Schoenfeld) ======\n")
print(ph_test)

## ============================================================
## BLOCK 8: FIGURE 1 — KAPLAN-MEIER CURVE
## ============================================================
## ggsurvplot() from survminer is the standard function for
## publication-quality KM curves in R.
##
## Key elements of a pharma-standard KM plot:
## - Two curves, one per arm, different colors
## - 95% confidence interval bands (shaded or dashed)
## - Number-at-risk table below the x-axis (REQUIRED in
##   regulatory submissions -- shows how many patients are
##   still being followed at each time point)
## - Median OS marked with horizontal dashed line at S(t)=0.5
## - HR and p-value annotated on the plot
## ============================================================

## Convert days to months for x-axis (standard in oncology)
ADTTE$OS_MONTHS_PLOT <- ADTTE$OS_DAYS / 30.4375

km_fit_months <- survfit(
  Surv(OS_MONTHS_PLOT, OS_EVENT) ~ TRT_LABEL,
  data = ADTTE
)

km_plot <- ggsurvplot(
  km_fit_months,
  data = ADTTE,
  palette = c("#C0504D", "#1B3A5C"),  # red=Darbepoetin, blue=Placebo (alphabetical)
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  risk.table = TRUE,                   # number-at-risk table below plot
  risk.table.height = 0.28,
  risk.table.title = "Number at Risk",
  xlab = "Time from Randomization (Months)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival by Treatment Arm\nITT Population (NCT00119613)",
  legend.labs = c("Darbepoetin Alfa", "Placebo"),
  legend.title = "",
  break.time.by = 3,                   # x-axis ticks every 3 months
  xlim = c(0, 24),
  surv.median.line = "hv",             # horizontal+vertical line at median
  pval = TRUE,                         # show log-rank p-value on plot
  pval.coord = c(1, 0.10),
  ggtheme = theme_classic(base_size = 13)
)

## Save the plot
png("Figure1_KM_Overall_Survival.png",
    width = 2000, height = 1600, res = 180)
print(km_plot)
dev.off()
cat("\nFigure 1 (KM curve) saved\n")

## ============================================================
## BLOCK 9: FIGURE 2 — FOREST PLOT (SUBGROUP ANALYSIS)
## ============================================================
## What is a subgroup forest plot?
## It shows the treatment effect (HR with 95% CI) estimated
## separately in pre-specified patient subgroups -- e.g., by
## ECOG status, age, sex. Each subgroup gets its own row on
## the plot, with the HR as a square and CI as horizontal lines.
## A vertical reference line at HR=1.0 marks "no effect."
##
## Why pre-specify subgroups?
## Post-hoc subgroup hunting is a major source of false
## positives in clinical research -- if you test enough
## subgroups, one will be "significant" by chance alone.
## Pre-specifying the subgroups (in the SAP) and
## treating them as exploratory is the methodologically
## honest approach.
##
## Per our SAP: ECOG 0-1 vs 2+, sex, age <65 vs >=65
## ============================================================

## Create subgroup variables
ADTTE_SG <- ADTTE %>%
  mutate(
    ECOG_GRP = case_when(
      b_ecog2.y <= 1 ~ "ECOG 0-1",
      b_ecog2.y >= 2 ~ "ECOG 2+",
      TRUE      ~ NA_character_
    ),
    AGE_GRP = case_when(
      age.y < 65  ~ "Age < 65",
      age.y >= 65 ~ "Age >= 65",
      TRUE      ~ NA_character_
    ),
    SEX_LBL = case_when(
      tolower(sex.y) %in% c("m","Male")   ~ "Male",
      tolower(sex.y) %in% c("f","Female") ~ "Female",
      TRUE                               ~ as.character(sex.y)
    )
  )

## Function to fit Cox model for one subgroup and extract HR/CI
fit_subgroup <- function(data, subgroup_label) {
  if (nrow(data) < 10 || n_distinct(data$TRT_LABEL) < 2) return(NULL)
  fit <- tryCatch(
    coxph(Surv(OS_DAYS, OS_EVENT) ~ TRT_LABEL, data = data),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  data.frame(
    Subgroup = subgroup_label,
    N        = nrow(data),
    Events   = sum(data$OS_EVENT, na.rm=TRUE),
    HR       = round(exp(coef(fit)), 3),
    Lower    = round(exp(confint(fit))[1], 3),
    Upper    = round(exp(confint(fit))[2], 3)
  )
}

## Run each subgroup
subgroup_results <- bind_rows(
  fit_subgroup(ADTTE_SG,                              "Overall (ITT)"),
  fit_subgroup(filter(ADTTE_SG, ECOG_GRP=="ECOG 0-1"), "ECOG 0-1"),
  fit_subgroup(filter(ADTTE_SG, ECOG_GRP=="ECOG 2+"),  "ECOG 2+"),
  fit_subgroup(filter(ADTTE_SG, AGE_GRP=="Age < 65"),  "Age < 65"),
  fit_subgroup(filter(ADTTE_SG, AGE_GRP=="Age >= 65"), "Age >= 65"),
  fit_subgroup(filter(ADTTE_SG, SEX_LBL=="Male"),      "Male"),
  fit_subgroup(filter(ADTTE_SG, SEX_LBL=="Female"),    "Female")
)

cat("\n====== SUBGROUP ANALYSIS RESULTS ======\n")
print(subgroup_results)

## Build forest plot
subgroup_results$Subgroup <- factor(subgroup_results$Subgroup,
                                    levels = rev(subgroup_results$Subgroup))
forest_plot <- ggplot(subgroup_results, aes(x = HR, y = Subgroup)) +
  geom_point(size = 3, color = "#1B3A5C") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.25, color = "#1B3A5C", linewidth = 0.8) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", HR, Lower, Upper)),
            hjust = -0.1, size = 3.2, color = "gray20") +
  labs(title = "Subgroup Analysis — Overall Survival\nHazard Ratio (Darbepoetin Alfa vs Placebo, ITT)",
       x = "Hazard Ratio (95% CI)", y = "") +
  xlim(0.4, 2.2) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 11))

ggsave("Figure2_Forest_Subgroup.png", forest_plot,
       width = 10, height = 6, dpi = 180)
cat("\nFigure 2 (subgroup forest plot) saved\n")

## ============================================================
## BLOCK 10: TABLE 2 — PRIMARY EFFICACY RESULTS
## ============================================================

## Extract median OS with 95% CI from KM fit
km_summary <- summary(km_fit)$table

## Survival probabilities at key timepoints
surv_6mo  <- summary(km_fit, times = 182.5)
surv_12mo <- summary(km_fit, times = 365)

## Build results table
results_table <- data.frame(
  Endpoint       = c("N randomized", "Deaths, n (%)",
                     "Median OS, days (95% CI)",
                     "Median OS, months",
                     "6-month OS probability (%)",
                     "12-month OS probability (%)",
                     "Hazard Ratio (vs Placebo)",
                     "95% CI",
                     "Log-rank p-value"),
  Darbepoetin    = c(
    "300",
    paste0(km_summary["TRT_LABEL=NESP","events"],
           " (", round(100*km_summary["TRT_LABEL=NESP","events"]/300,1), "%)"),
    paste0(round(km_summary["TRT_LABEL=NESP","median"],0),
           " (", round(km_summary["TRT_LABEL=NESP","0.95LCL"],0),
           "-", round(km_summary["TRT_LABEL=NESP","0.95UCL"],0), ")"),
    as.character(round(km_summary["TRT_LABEL=NESP","median"]/30.4375,1)),
    paste0(round(100*surv_6mo$surv[surv_6mo$strata=="TRT_LABEL=NESP"],1),"%"),
    paste0(round(100*surv_12mo$surv[surv_12mo$strata=="TRT_LABEL=NESP"],1),"%"),
    as.character(round(hr,3)),
    paste0(round(ci[1],3), " to ", round(ci[2],3)),
    as.character(round(logrank_pval,4))
  ),
  Placebo = c(
    "300",
    paste0(km_summary["TRT_LABEL=PLACEBO","events"],
           " (", round(100*km_summary["TRT_LABEL=PLACEBO","events"]/300,1), "%)"),
    paste0(round(km_summary["TRT_LABEL=PLACEBO","median"],0),
           " (", round(km_summary["TRT_LABEL=PLACEBO","0.95LCL"],0),
           "-", round(km_summary["TRT_LABEL=PLACEBO","0.95UCL"],0), ")"),
    as.character(round(km_summary["TRT_LABEL=PLACEBO","median"]/30.4375,1)),
    paste0(round(100*surv_6mo$surv[surv_6mo$strata=="TRT_LABEL=PLACEBO"],1),"%"),
    paste0(round(100*surv_12mo$surv[surv_12mo$strata=="TRT_LABEL=PLACEBO"],1),"%"),
    "Reference",
    "—",
    "—"
  )
)

cat("\n====== TABLE 2: PRIMARY EFFICACY RESULTS ======\n")
print(results_table, row.names = FALSE)

write.csv(results_table, "Table2_OS_Results.csv", row.names = FALSE)
cat("\nTable 2 saved as Table2_OS_Results.csv\n")

## ============================================================
## BLOCK 12: SAVE ALL OUTPUTS
## ============================================================

saveRDS(ADTTE, "ADTTE.rds")
saveRDS(subgroup_results, "subgroup_results.rds")
saveRDS(results_table, "results_table.rds")
write.csv(subgroup_results, "Table3_Subgroup.csv", row.names = FALSE)

cat("\n=================== ANALYSIS COMPLETE ===================\n")
cat("Saved: ADTTE.rds\n")
cat("Saved: Table1_Baseline.csv\n")
cat("Saved: Table2_OS_Results.csv\n")
cat("Saved: Table3_Subgroup.csv\n")
cat("Saved: Figure1_KM_Overall_Survival.png\n")
cat("Saved: Figure2_Forest_Subgroup.png\n")

