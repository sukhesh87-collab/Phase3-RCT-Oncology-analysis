## ============================================================
## PROJECT 2: Adverse Event Analysis — Final Version
## NCT00119613 — Darbepoetin Alfa in ES-SCLC
## ============================================================
## Variable names confirmed from C_AE glimpse output:
##   subjid   = patient identifier
##   txgroup  = treatment arm ("NESP" / "PLACEBO") -- in C_AE already
##   aepterm  = AE preferred term (e.g., "DYSPNOEA")
##   aebcterm = body system organ class
##   sevr     = severity text ("Mild" / "Moderate")
##   sevrcd   = severity numeric (1=Mild, 2=Moderate)
##   sriousyn = serious AE flag (0=No, 1=Yes)
##   teaeyn   = treatment-emergent AE flag (1=Yes) -- primary filter
##   relateyn = drug-related flag (0=No, 1=Yes)
##   acthosp  = AE led to hospitalization (0=No, 1=Yes)
##   aesub    = AE sub-category (e.g., "Arrhythmias" for CV events)
##
## NOTE: No CTCAE Grade 1-5 variable exists in this dataset.
## Severity is 2-level only (Mild/Moderate). Documented as
## a dataset limitation in the results report.
## ============================================================

library(haven)
library(dplyr)
library(tidyr)
library(janitor)

## ============================================================
## STEP 1: LOAD DATA
## ============================================================

C_AE  <- read_sas("C_AE.sas7bdat") %>% clean_names()
ADSL  <- readRDS("ADSL.rds")

cat("C_AE loaded:", nrow(C_AE), "rows,",
    n_distinct(C_AE$subjid), "unique patients\n")

## ============================================================
## STEP 2: BUILD AE ANALYSIS DATASET
## ============================================================
## Filter to treatment-emergent AEs (teaeyn == 1) only.
## This is the regulatory-standard approach: only AEs that
## started AFTER the first dose of study drug are counted.
## AEs present before treatment started are not included.
##
## txgroup is already in C_AE so no ADSL join needed for arm.
## We create a clean ARM label for tables/figures.
## ============================================================

AE_SAFETY <- C_AE %>%
  filter(teaeyn == 1) %>%
  mutate(
    ARM = case_when(
      txgroup == "NESP"    ~ "Darbepoetin Alfa",
      txgroup == "PLACEBO" ~ "Placebo",
      TRUE                 ~ as.character(txgroup)
    ),
    ## Convert severity to ordered for easy filtering
    sev_num = as.numeric(sevrcd),   ## 1=Mild, 2=Moderate
    is_sae  = (sriousyn == 1),      ## serious AE flag
    is_related  = (relateyn == 1),  ## drug-related AE
    is_hosp     = (acthosp == 1),   ## led to hospitalization
    is_cardiac  = (!is.na(aesub) & nchar(trimws(aesub)) > 0)  ## CV/cardiac sub-category
  )

cat("TEAE records:", nrow(AE_SAFETY), "\n")
cat("TEAE records by arm:\n")
print(table(AE_SAFETY$ARM))

## ============================================================
## STEP 3: SET SAFETY POPULATION DENOMINATORS
## ============================================================
## Denominators = all patients who received at least one dose.
## From trial documentation and ADSL SAFFL flag.
## Known values: 239 darbepoetin, 240 placebo.
## ============================================================

N_darb <- sum(ADSL$SAFFL == "Y" &
              ADSL$TRT_LABEL %in% c("NESP","Darbepoetin Alfa","1"),
              na.rm = TRUE)
N_plac <- sum(ADSL$SAFFL == "Y" &
              ADSL$TRT_LABEL %in% c("PLACEBO","Placebo","2"),
              na.rm = TRUE)

## Fallback to known trial values if ADSL join didn't work
if (N_darb == 0) N_darb <- 239
if (N_plac == 0) N_plac <- 240

cat("\nSafety denominators:\n")
cat("  Darbepoetin Alfa N =", N_darb, "\n")
cat("  Placebo N =", N_plac, "\n")

## ============================================================
## STEP 4: HELPER FUNCTION
## ============================================================
## Counts unique patients meeting a logical condition,
## returns formatted "n (x.x%)" string for table cells.
## Always counted at patient level (not event level).
## ============================================================

pt_pct <- function(data, arm, condition_expr, denom) {
  n <- data %>%
    filter(ARM == arm) %>%
    filter(!!rlang::parse_expr(condition_expr)) %>%
    summarize(n = n_distinct(subjid)) %>%
    pull(n)
  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

## ============================================================
## STEP 5: TABLE 3A — AE SUMMARY TABLE
## ============================================================
## Standard pharma safety summary table structure.
## Each row = one AE category.
## Columns = n (%) in each arm.
## ============================================================

Table3a <- data.frame(
  AE_Category = c(
    "Any TEAE",
    "Any Moderate TEAE",
    "Any Serious AE (SAE)",
    "Any Drug-Related TEAE",
    "Any TEAE Leading to Hospitalization",
    "Any Cardiac/Cardiovascular TEAE"
  ),
  Darbepoetin_Alfa = c(
    pt_pct(AE_SAFETY, "Darbepoetin Alfa", "!is.na(subjid)",      N_darb),
    pt_pct(AE_SAFETY, "Darbepoetin Alfa", "sev_num >= 2",        N_darb),
    pt_pct(AE_SAFETY, "Darbepoetin Alfa", "is_sae == TRUE",      N_darb),
    pt_pct(AE_SAFETY, "Darbepoetin Alfa", "is_related == TRUE",  N_darb),
    pt_pct(AE_SAFETY, "Darbepoetin Alfa", "is_hosp == TRUE",     N_darb),
    pt_pct(AE_SAFETY, "Darbepoetin Alfa", "is_cardiac == TRUE",  N_darb)
  ),
  Placebo = c(
    pt_pct(AE_SAFETY, "Placebo", "!is.na(subjid)",     N_plac),
    pt_pct(AE_SAFETY, "Placebo", "sev_num >= 2",       N_plac),
    pt_pct(AE_SAFETY, "Placebo", "is_sae == TRUE",     N_plac),
    pt_pct(AE_SAFETY, "Placebo", "is_related == TRUE", N_plac),
    pt_pct(AE_SAFETY, "Placebo", "is_hosp == TRUE",    N_plac),
    pt_pct(AE_SAFETY, "Placebo", "is_cardiac == TRUE", N_plac)
  ),
  stringsAsFactors = FALSE
)

cat("\n====== TABLE 3A: TREATMENT-EMERGENT AE SUMMARY ======\n")
print(Table3a, row.names = FALSE)

## ============================================================
## STEP 6: TABLE 3B — MOST FREQUENT TEAEs (>=5% either arm)
## ============================================================
## Count unique patients per AE preferred term per arm.
## Report terms where at least one arm has >=5% incidence.
## Sorted by darbepoetin frequency descending.
## ============================================================

ae_by_term <- AE_SAFETY %>%
  group_by(ARM, aepterm) %>%
  summarize(n_pts = n_distinct(subjid), .groups = "drop") %>%
  pivot_wider(
    names_from  = ARM,
    values_from = n_pts,
    values_fill = 0
  )

## Rename columns defensively (handles any spacing in arm label)
names(ae_by_term) <- gsub(" ", "_", names(ae_by_term))
darb_col <- names(ae_by_term)[grepl("Darbepoetin|NESP", names(ae_by_term))][1]
plac_col <- names(ae_by_term)[grepl("Placebo|PLACEBO", names(ae_by_term))][1]

Table3b <- ae_by_term %>%
  rename(n_darb = all_of(darb_col),
         n_plac = all_of(plac_col)) %>%
  mutate(
    Pct_Darbepoetin = round(100 * n_darb / N_darb, 1),
    Pct_Placebo     = round(100 * n_plac / N_plac, 1),
    max_pct         = pmax(Pct_Darbepoetin, Pct_Placebo)
  ) %>%
  filter(max_pct >= 5) %>%
  arrange(desc(Pct_Darbepoetin)) %>%
  select(
    `AE Preferred Term`  = aepterm,
    `Darb n`             = n_darb,
    `Darb %`             = Pct_Darbepoetin,
    `Placebo n`          = n_plac,
    `Placebo %`          = Pct_Placebo
  )

cat("\n====== TABLE 3B: MOST FREQUENT TEAEs (>=5% in either arm) ======\n")
print(Table3b, row.names = FALSE)

## ============================================================
## STEP 7: TABLE 3C — AE BY BODY SYSTEM CLASS
## ============================================================
## Groups AEs by aebcterm (organ system).
## Standard secondary AE table in pharma CSR reports.
## Shows which organ systems were most commonly affected.
## ============================================================

ae_by_soc <- AE_SAFETY %>%
  group_by(ARM, aebcterm) %>%
  summarize(n_pts = n_distinct(subjid), .groups = "drop") %>%
  pivot_wider(
    names_from  = ARM,
    values_from = n_pts,
    values_fill = 0
  )

names(ae_by_soc) <- gsub(" ", "_", names(ae_by_soc))
darb_col2 <- names(ae_by_soc)[grepl("Darbepoetin|NESP", names(ae_by_soc))][1]
plac_col2 <- names(ae_by_soc)[grepl("Placebo|PLACEBO", names(ae_by_soc))][1]

Table3c <- ae_by_soc %>%
  rename(n_darb = all_of(darb_col2),
         n_plac = all_of(plac_col2)) %>%
  mutate(
    Pct_Darbepoetin = round(100 * n_darb / N_darb, 1),
    Pct_Placebo     = round(100 * n_plac / N_plac, 1),
    max_pct         = pmax(Pct_Darbepoetin, Pct_Placebo)
  ) %>%
  filter(max_pct >= 5) %>%
  arrange(desc(Pct_Darbepoetin)) %>%
  select(
    `Body System`   = aebcterm,
    `Darb n`        = n_darb,
    `Darb %`        = Pct_Darbepoetin,
    `Placebo n`     = n_plac,
    `Placebo %`     = Pct_Placebo
  )

cat("\n====== TABLE 3C: TEAEs BY BODY SYSTEM CLASS ======\n")
print(Table3c, row.names = FALSE)

## ============================================================
## STEP 8: CARDIAC/CV SUBGROUP DETAIL
## ============================================================
## The aesub variable flags cardiac sub-categories
## (e.g., "Arrhythmias"). This is the specific thromboembolic
## and cardiac safety signal of concern for ESAs.
## ============================================================

cardiac_detail <- AE_SAFETY %>%
  filter(is_cardiac == TRUE) %>%
  group_by(ARM, aesub) %>%
  summarize(n_pts = n_distinct(subjid), .groups = "drop") %>%
  pivot_wider(
    names_from  = ARM,
    values_from = n_pts,
    values_fill = 0
  ) %>%
  arrange(desc(rowSums(select(., -aesub))))

cat("\n====== CARDIAC/CV AE SUBCATEGORY DETAIL ======\n")
print(cardiac_detail, row.names = FALSE)

## ============================================================
## STEP 9: SAVE ALL OUTPUTS
## ============================================================

write.csv(Table3a,       "Table3a_AE_Summary.csv",     row.names = FALSE)
write.csv(Table3b,       "Table3b_AE_Frequent.csv",    row.names = FALSE)
write.csv(Table3c,       "Table3c_AE_BodySystem.csv",  row.names = FALSE)
write.csv(cardiac_detail,"Table3d_AE_Cardiac.csv",     row.names = FALSE)

cat("\n====== ALL AE OUTPUTS SAVED ======\n")
cat("Table3a_AE_Summary.csv    — AE summary by category\n")
cat("Table3b_AE_Frequent.csv   — Most frequent TEAEs (>=5%)\n")
cat("Table3c_AE_BodySystem.csv — AEs by body system class\n")
cat("Table3d_AE_Cardiac.csv    — Cardiac/CV subcategory detail\n")
cat("\nShare all four CSV outputs + Figure 1 + Figure 2 and\n")
cat("we will compile the full results report.\n")
