## ============================================================
## PROJECT 2: Per-Protocol Flag Derivation + Sensitivity Analysis
## NCT00119613 — Darbepoetin Alfa in ES-SCLC
## ============================================================
## What PP should mean in this trial:
## Per-Protocol = patients who followed the trial protocol
## sufficiently to provide valid efficacy data. This applies
## to BOTH arms -- a placebo patient who completed chemotherapy
## cycles per protocol is equally "per-protocol."
##
## PP derivation approach (pre-specified in SAP):
## Include patients who:
## 1. Received at least 2 cycles of chemotherapy (minimum
##    meaningful treatment per oncology PP convention)
## 2. Had no documented major protocol deviation in C_DISP
##    that would invalidate their efficacy data
## ============================================================

library(haven)
library(dplyr)
library(survival)
library(janitor)

## ============================================================
## BLOCK 1: RELOAD ADSL AND SOURCE DATA
## ============================================================

ADSL  <- readRDS("ADSL.rds")
ADTTE <- readRDS("ADTTE.rds")

## Reload the raw tabulation datasets we need for PP derivation
C_CHEMO <- read_sas("C_CHEMO.sas7bdat") %>% clean_names()
C_DISP  <- read_sas("C_DISP.sas7bdat")  %>% clean_names()

cat("C_CHEMO loaded:", nrow(C_CHEMO), "rows\n")
cat("C_DISP loaded:", nrow(C_DISP), "rows\n")

## First, inspect the variables we need
cat("\nC_CHEMO variables:\n")
print(names(C_CHEMO))

cat("\nC_DISP variables:\n")
print(names(C_DISP))

## ============================================================
## BLOCK 2: EXPLORE CHEMO CYCLE INFORMATION
## ============================================================
## Before deriving anything, understand what's in C_CHEMO:
## - How are cycles recorded? (cycle number variable?)
## - Is there a visit number, cycle count, or date sequence?
## - Are both arms (NESP and placebo) represented?
## ============================================================

## What does the cycle/visit structure look like?
cat("\nC_CHEMO structure summary:\n")
glimpse(C_CHEMO)

## How many records per patient? (tells us the cycle structure)
records_per_patient <- C_CHEMO %>%
  group_by(subjid) %>%
  summarize(n_records = n(), .groups = "drop")

cat("\nRecords per patient in C_CHEMO:\n")
print(summary(records_per_patient$n_records))

## Are both treatment arms in C_CHEMO?
## This is the critical check -- if placebo patients appear
## in C_CHEMO, we can derive PP for both arms
cat("\nTreatment arms represented in C_CHEMO:\n")
chemo_with_arm <- C_CHEMO %>%
  left_join(ADSL %>% select(subjid, TRT_LABEL), by = "subjid")
print(table(chemo_with_arm$TRT_LABEL, useNA = "always"))

## What are the unique values of any cycle/visit variable?
## Adjust 'cycle' to whatever the actual variable is named
## Common names: cycle, cycnum, visit, visitnum, cyc
if ("cycle" %in% names(C_CHEMO)) {
  cat("\nUnique cycle values:\n")
  print(sort(unique(C_CHEMO$cycle)))
} else if ("cycnum" %in% names(C_CHEMO)) {
  cat("\nUnique cycnum values:\n")
  print(sort(unique(C_CHEMO$cycnum)))
} else {
  cat("\nNo 'cycle' or 'cycnum' found -- check variable names above\n")
  cat("Candidates for cycle variable:\n")
  print(names(C_CHEMO)[grepl("cy|vis|num", names(C_CHEMO), ignore.case=TRUE)])
}

## ============================================================
## BLOCK 3: EXPLORE DISPOSITION FOR PROTOCOL DEVIATIONS
## ============================================================
## C_DISP captures why/when patients left the study.
## We want to identify patients with MAJOR protocol deviations
## that would invalidate their efficacy data.
##
## Common coding: a "reason for discontinuation" variable
## with values like "Protocol Violation", "Investigator Decision"
## etc. We'll exclude these from the PP population.
## ============================================================

cat("\nC_DISP structure:\n")
glimpse(C_DISP)

## What are the unique disposition/discontinuation reasons?
## Adjust variable name to match your actual C_DISP output
disp_vars <- names(C_DISP)[grepl("reas|disc|stat|end|compl|dev",
                                   names(C_DISP), ignore.case=TRUE)]
cat("\nPotential disposition reason variables:\n")
print(disp_vars)

## Print unique values of each candidate variable
for (v in disp_vars[1:min(3, length(disp_vars))]) {
  cat("\nUnique values of", v, ":\n")
  print(table(C_DISP[[v]], useNA = "always"))
}

## ============================================================
## BLOCK 4a: DERIVE CYCLE COMPLETION CRITERION
## ============================================================
## Logic: patient completed >= 2 cycles in C_CHEMO
## Why >= 2? One cycle alone is insufficient to assess protocol
## compliance in a chemotherapy trial -- 2 cycles is the
## standard minimum meaningful exposure threshold in oncology PP.
##
## n_distinct(cycle) counts how many DIFFERENT cycle numbers
## a patient has records for -- this is more robust than
## counting rows, since a patient could have multiple rows
## per cycle (e.g., multiple dose administrations per cycle).
## ============================================================

cycles_completed <- C_CHEMO %>%
  group_by(pt) %>%
  summarize(
    n_cycles_distinct = n_distinct(cycle, na.rm = TRUE),
    max_cycle         = max(cycle, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  mutate(
    cycle_criterion = (n_cycles_distinct >= 2)
  )

cat("====== CYCLE COMPLETION SUMMARY ======\n")
cat("Total patients with any chemo record:", nrow(cycles_completed), "\n")
print(table(cycles_completed$n_cycles_distinct))
cat("Patients meeting cycle criterion (>=2 cycles):",
    sum(cycles_completed$cycle_criterion), "\n")
cat("Patients NOT meeting cycle criterion (<2 cycles):",
    sum(!cycles_completed$cycle_criterion), "\n")

## ============================================================
## BLOCK 4b: DERIVE PROTOCOL DEVIATION CRITERION FROM C_DISP
## ============================================================
## We exclude from PP any patient whose primary reason for
## ending treatment was:
##   "Protocol Violation"      (n=3)  -- clear major deviation
##   "Ineligibility Determined" (n=1) -- not eligible, data invalid
##   "Noncompliance"            (n=4) -- didn't follow protocol
##
## Patients with NO disposition record (n = 600 - 479 = 121
## approximately) are assumed to have NO major deviation --
## they are still in follow-up, which is not a deviation.
## replace_na(FALSE) handles these: missing = no deviation.
##
## Note: "Disease Progression" (n=70), "Death" (n=50),
## "Adverse Event" (n=21) are NOT protocol deviations --
## these are valid clinical outcomes, not compliance failures.
## Patients who progressed or died are still per-protocol
## (they received treatment and had the outcome).
## ============================================================

## Define which endtx values constitute major protocol deviations
MAJOR_DEVIATION_VALUES <- c(
  "Protocol Violation",
  "Ineligibility Determined",
  "Noncompliance"
)

major_deviations <- C_DISP %>%
  mutate(
    is_deviation = endtx %in% MAJOR_DEVIATION_VALUES
  ) %>%
  group_by(pt) %>%
  summarize(
    any_deviation = any(is_deviation, na.rm = TRUE),
    deviation_reason = paste(endtx[is_deviation], collapse = "; "),
    .groups = "drop"
  )

cat("\n====== PROTOCOL DEVIATION SUMMARY ======\n")
cat("Patients with any disposition record:", nrow(major_deviations), "\n")
cat("Patients with major protocol deviation:",
    sum(major_deviations$any_deviation), "\n")
cat("\nDeviation reasons:\n")
print(table(major_deviations$deviation_reason[major_deviations$any_deviation]))

## ============================================================
## BLOCK 4c: COMBINE CRITERIA AND DERIVE PPFL
## ============================================================
## PPFL = "Y" if BOTH criteria are met:
##   (1) cycle_criterion == TRUE  (>= 2 cycles completed)
##   (2) any_deviation  == FALSE  (no major protocol deviation)
##
## Patients with no chemo record: cycle_criterion defaults
## to FALSE via replace_na -> excluded from PP (reasonable --
## if we have no dosing record, we can't confirm compliance).
##
## Patients with no disposition record: any_deviation defaults
## to FALSE via replace_na -> included in PP (benefit of doubt --
## no recorded deviation means no known deviation).
## ============================================================

ADSL <- ADSL %>%
  left_join(
    cycles_completed %>% select(pt, n_cycles_distinct, cycle_criterion),
    by = "pt"
  ) %>%
  left_join(
    major_deviations %>% select(pt, any_deviation),
    by = "pt"
  ) %>%
  mutate(
    ## Apply defaults for missing records
    cycle_criterion = replace_na(cycle_criterion, FALSE),
    any_deviation   = replace_na(any_deviation,   FALSE),
    ## Derive PPFL
    PPFL = case_when(
      cycle_criterion & !any_deviation ~ "Y",
      TRUE                             ~ "N"
    )
  )

## ============================================================
## BLOCK 4d: VALIDATE PPFL -- CRITICAL CHECK
## ============================================================
cat("\n====== PPFL VALIDATION ======\n")
cat("Total patients (ITT):", nrow(ADSL), "(expect 600)\n")
cat("PP population (PPFL=Y):", sum(ADSL$PPFL == "Y"), "\n")
cat("Excluded from PP (PPFL=N):", sum(ADSL$PPFL == "N"), "\n")

cat("\nPPFL by treatment arm -- BOTH arms must have Y patients:\n")
print(table(ADSL$TRT_LABEL, ADSL$PPFL, dnn = c("Arm", "PPFL")))
## If both arms show Y patients, proceed to Block 5.
## If only one arm shows Y, something is wrong -- stop and share output.

cat("\nCycle completion by arm:\n")
print(
  ADSL %>%
    group_by(TRT_LABEL) %>%
    summarize(
      N_total    = n(),
      PP_Y       = sum(PPFL == "Y"),
      PP_N       = sum(PPFL == "N"),
      Pct_PP     = round(100 * mean(PPFL == "Y"), 1),
      Mean_cycles = round(mean(n_cycles_distinct, na.rm=TRUE), 1),
      .groups    = "drop"
    )
)

cat("\nReasons for PP exclusion:\n")
cat("  No chemo record or < 2 cycles (cycle_criterion=FALSE):",
    sum(!ADSL$cycle_criterion), "\n")
cat("  Major protocol deviation (any_deviation=TRUE):",
    sum(ADSL$any_deviation), "\n")
cat("  Both reasons:",
    sum(!ADSL$cycle_criterion & ADSL$any_deviation), "\n")

## Save updated ADSL with PPFL
saveRDS(ADSL, "ADSL.rds")
cat("\nADSL updated with PPFL and saved\n")

## ============================================================
## BLOCK 5: PP SENSITIVITY ANALYSIS (Cox PH on PP population)
## ============================================================
## We merge PPFL from ADSL into ADTTE (ADTTE was built before
## PPFL existed, so it doesn't have it yet), filter to PP,
## and rerun the same Cox model as the primary analysis.
##
## Expected result: HR directionally consistent with primary
## ITT analysis, possibly slightly different magnitude since
## PP excludes the least compliant patients (who may have
## worse baseline characteristics or more adverse events).
## ============================================================

## Merge PPFL into ADTTE
ADTTE <- ADTTE %>%
  select(-any_of("PPFL")) %>%  # remove old PPFL column if exists
  left_join(ADSL %>% select(pt, PPFL), by = "pt")

## Confirm both arms present in PP subset
ADTTE_PP <- ADTTE %>% filter(PPFL == "Y")

cat("\n====== PP SENSITIVITY ANALYSIS ======\n")
cat("PP population size:", nrow(ADTTE_PP), "\n")

cat("\nPP by arm:\n")
print(table(ADTTE_PP$TRT_LABEL))

cat("\nEvents by arm (PP):\n")
print(table(ADTTE_PP$TRT_LABEL, ADTTE_PP$OS_EVENT,
            dnn = c("Arm", "Event (1=Death)")))

## Safety check before running model
if (n_distinct(ADTTE_PP$TRT_LABEL) < 2) {
  stop("ERROR: Only one treatment arm in PP population. Check PPFL derivation.")
}
if (sum(ADTTE_PP$OS_EVENT) < 10) {
  warning("Very few events in PP population -- Cox model may be unstable.")
}

## Cox PH model on PP population
cox_pp <- coxph(
  Surv(OS_DAYS, OS_EVENT) ~ TRT_LABEL,
  data = ADTTE_PP
)

hr_pp <- exp(coef(cox_pp))
ci_pp <- exp(confint(cox_pp))
pv_pp <- summary(cox_pp)$coefficients[,"Pr(>|z|)"]

## Log-rank test on PP population
lr_pp <- survdiff(Surv(OS_DAYS, OS_EVENT) ~ TRT_LABEL, data = ADTTE_PP)
lr_pval_pp <- 1 - pchisq(lr_pp$chisq, df = 1)

cat("\n=================== PP SENSITIVITY RESULT ===================\n")
cat("PP HR (Darbepoetin vs Placebo):", round(hr_pp, 3), "\n")
cat("95% CI:", round(ci_pp[1],3), "to", round(ci_pp[2],3), "\n")
cat("Cox p-value:", round(pv_pp, 4), "\n")
cat("Log-rank p-value:", round(lr_pval_pp, 4), "\n")
cat("=============================================================\n")
cat("\nInterpretation:\n")
cat("Compare this HR to your primary ITT HR from Script 02.\n")
cat("If both point in the same direction and are similar in\n")
cat("magnitude, the PP analysis supports robustness of the ITT finding.\n")
cat("If they diverge, this warrants discussion of non-compliance effects.\n")

## Save combined results for Table 2 update
pp_results <- data.frame(
  Analysis  = c("Sensitivity (Per-Protocol)"),
  N         = nrow(ADTTE_PP),
  Events    = sum(ADTTE_PP$OS_EVENT),
  HR        = round(hr_pp, 3),
  CI_Lower  = round(ci_pp[1], 3),
  CI_Upper  = round(ci_pp[2], 3),
  p_value   = round(pv_pp, 4)
)

write.csv(pp_results, "Table_PP_Sensitivity.csv", row.names = FALSE)
saveRDS(ADTTE, "ADTTE.rds")

cat("\nSaved: Table_PP_Sensitivity.csv\n")
cat("Saved: ADTTE.rds (updated with PPFL)\n")
cat("\nPP sensitivity analysis complete.\n")

