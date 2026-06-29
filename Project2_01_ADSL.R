## ============================================================
## PROJECT 2: Phase 3 RCT Analysis
## Darbepoetin Alfa in Extensive-Stage SCLC (NCT00119613)
## Step 1: Package Installation, Data Loading, ADSL Construction
## ============================================================

## ============================================================
## BLOCK 1: INSTALL REQUIRED PACKAGES
## ============================================================
install.packages(c(
  "haven",      # read SAS .sas7bdat files
  "dplyr",      # data manipulation
  "tidyr",      # data reshaping
  "lubridate",  # date arithmetic
  "survival",   # Kaplan-Meier, Cox PH, log-rank
  "survminer",  # KM plots with number-at-risk
  "ggplot2",    # visualization
  "gtsummary",  # Table 1 generation
  "tableone",   # alternative Table 1
  "janitor",    # data cleaning utilities (clean variable names)
  "writexl"     # export tables to Excel for inspection
))

## ============================================================
## BLOCK 2: LOAD PACKAGES INTO CURRENT R SESSION
## ============================================================
library(haven)      # SAS file reader
library(dplyr)      # data manipulation
library(tidyr)      # data reshaping
library(lubridate)  # date arithmetic
library(survival)   # survival analysis
library(survminer)  # KM plots
library(ggplot2)    # visualization
library(gtsummary)  # Table 1
library(tableone)   # Table 1 alternative
library(janitor)    # data cleaning
library(writexl)    # Excel export

## ============================================================
## BLOCK 3: SET YOUR WORKING DIRECTORY
## ============================================================

## EDIT THIS LINE -- put your actual folder path here:
setwd("Datasets/SAS files")

## Confirm it worked -- this shows your current directory:
getwd()

## List all files in that folder to confirm SAS files are there:
list.files(pattern = "\\.sas7bdat$")
## You should see all your C_ and A_ files listed.
## If you see an empty result, your path is wrong or the files
## aren't in that folder -- double-check and re-run.

## ============================================================
## BLOCK 4: LOAD ALL SAS DATASETS INTO R
## ============================================================
## 
##
## Why are we loading ALL datasets first?
## We want to understand the full data landscape before writing any 
## analysis code.
## We'll explore each dataset, then decide which variables to pull
## into ADSL.
## ============================================================

## --- Tabulation datasets (C_ prefix = raw/SDTM-equivalent) ---

## C_KEYVAR: Patient identifiers, treatment arm, randomization
## This is your "spine" -- every patient should appear here.
## Equivalent to SDTM DM (Demographics) domain.
C_KEYVAR <- read_sas("C_KEYVAR.sas7bdat") %>% clean_names()
cat("C_KEYVAR loaded:", nrow(C_KEYVAR), "rows,", ncol(C_KEYVAR), "columns\n")

## C_BCHAR: Baseline characteristics
## Age, sex, ECOG performance status, disease history.
## ECOG is the oncology equivalent of our comorbidity score
## -- measures functional status (0=fully active,
## 4=completely disabled). Standard baseline covariate in all
## oncology trials.
C_BCHAR <- read_sas("C_BCHAR.sas7bdat") %>% clean_names()
cat("C_BCHAR loaded:", nrow(C_BCHAR), "rows,", ncol(C_BCHAR), "columns\n")

## C_DISP: Patient disposition
## Why/when patients left the study (completed, withdrew consent,
## died, lost to follow-up). Drives our censoring logic --
## if a patient didn't die, this tells us when they exited and
## why, which determines their censoring date for OS analysis.
C_DISP <- read_sas("C_DISP.sas7bdat") %>% clean_names()
cat("C_DISP loaded:", nrow(C_DISP), "rows,", ncol(C_DISP), "columns\n")

## C_AE: Adverse events
## One row per adverse event per patient (patients can have
## multiple rows). Contains CTCAE grade, event term, dates,
## seriousness flag. Used for our safety analysis.
C_AE <- read_sas("C_AE.sas7bdat") %>% clean_names()
cat("C_AE loaded:", nrow(C_AE), "rows,", ncol(C_AE), "columns\n")

## C_HEMAT: Hematology
## Hemoglobin and blood count measurements over time.
## Central to this trial -- darbepoetin alfa raises hemoglobin,
## so tracking Hgb over time is a key secondary endpoint.
C_HEMAT <- read_sas("C_HEMAT.sas7bdat") %>% clean_names()
cat("C_HEMAT loaded:", nrow(C_HEMAT), "rows,", ncol(C_HEMAT), "columns\n")

## C_CHEMO: Chemotherapy/study drug exposure
## Dosing records for both the platinum/etoposide chemotherapy
## and the study drug (darbepoetin alfa or placebo).
## Used to derive Safety Population flag (at least one dose).
C_CHEMO <- read_sas("C_CHEMO.sas7bdat") %>% clean_names()
cat("C_CHEMO loaded:", nrow(C_CHEMO), "rows,", ncol(C_CHEMO), "columns\n")

## --- Analysis datasets (A_ prefix = ADaM-equivalent) ---

## A_EENDPT: Efficacy endpoints -- this is our PRIMARY dataset
## Contains the OS variables we need most:
## dth = death indicator (1=died, 0=censored/alive)
## dthdy = days from randomization to death or last contact
## dthwk = same in weeks (convenience variable)
A_EENDPT <- read_sas("A_EENDPT.sas7bdat") %>% clean_names()
cat("A_EENDPT loaded:", nrow(A_EENDPT), "rows,", ncol(A_EENDPT), "columns\n")

## A_EENDES: Efficacy endpoint descriptions
## Companion to A_EENDPT with event description text
A_EENDES <- read_sas("A_EENDES.sas7bdat") %>% clean_names()
cat("A_EENDES loaded:", nrow(A_EENDES), "rows,", ncol(A_EENDES), "columns\n")

## A_EENDFU: Efficacy follow-up
## Follow-up dates and last known alive information
A_EENDFU <- read_sas("A_EENDFU.sas7bdat") %>% clean_names()
cat("A_EENDFU loaded:", nrow(A_EENDFU), "rows,", ncol(A_EENDFU), "columns\n")

## A_SENDPT: Safety endpoints
## Pre-derived safety flags: cardiovascular events,
## thromboembolic events, hemoglobin threshold crossings
A_SENDPT <- read_sas("A_SENDPT.sas7bdat") %>% clean_names()
cat("A_SENDPT loaded:", nrow(A_SENDPT), "rows,", ncol(A_SENDPT), "columns\n")

## A_LAB: Laboratory data (analysis-ready)
A_LAB <- read_sas("A_LAB.sas7bdat") %>% clean_names()
cat("A_LAB loaded:", nrow(A_LAB), "rows,", ncol(A_LAB), "columns\n")

## A_QOL: Quality of life scores
## FACT-F (Fatigue), FACT-G (General), FACT-A (Anemia),
## EQ-5D (health utility). HEOR-relevant PRO endpoints.
A_QOL <- read_sas("A_QOL.sas7bdat") %>% clean_names()
cat("A_QOL loaded:", nrow(A_QOL), "rows,", ncol(A_QOL), "columns\n")

## ============================================================
## BLOCK 5: EXPLORE EACH DATASET STRUCTURE
## ============================================================
## Run each glimpse() separately and study the output carefully
## before moving on. This is where you verify your variable
## names match what the DDT says they should be.
## ============================================================

## Start with C_KEYVAR -- most important, the patient spine:
cat("\n====== C_KEYVAR: Patient Identifiers & Treatment Arm ======\n")
glimpse(C_KEYVAR)

## Look for: patient ID variable (likely 'pt', 'id', or 'usubjid'),
## treatment group variable (likely 'trt', 'grp', or 'txgroup'),
## randomization date if present.
## Note the exact variable names -- you'll need them in Block 6.

cat("\n====== C_BCHAR: Baseline Characteristics ======\n")
glimpse(C_BCHAR)
## Look for: age, sex, ECOG score variable name

cat("\n====== A_EENDPT: OS Endpoint ======\n")
glimpse(A_EENDPT)
## Look for: dth, dthdy, dthwk -- confirm they exist with these names
## Also note the patient ID variable name (should match C_KEYVAR)

cat("\n====== C_DISP: Disposition ======\n")
glimpse(C_DISP)
## Look for: discontinuation reason, dates

cat("\n====== C_AE: Adverse Events ======\n")
glimpse(C_AE)
## Look for: AE term, CTCAE grade variable, seriousness flag

cat("\n====== C_CHEMO: Exposure ======\n")
glimpse(C_CHEMO)
## Look for: dose variable, date of first dose

## ============================================================
## BLOCK 6: CHECK KEY VARIABLES BEFORE BUILDING ADSL
## ============================================================
## Before we merge datasets, we must:
## 1. Confirm the patient ID variable is the same name across
##    all datasets (it's the "key" that links them together)
## 2. Check treatment group coding (what values does it take?
##    "1"/"2", "DARB"/"PBO", "Active"/"Placebo"?)
## 3. Check death indicator coding (is 1=died or is it "Y"/"N"?)
## 4. Confirm row counts make sense (C_KEYVAR should have ~600 rows)
##
## This is the data governance discipline step 
## ============================================================

## How many unique patients in each dataset?
## n_distinct() counts unique values -- use this to confirm
## no unexpected duplicates in our "one row per patient" datasets

cat("\nUnique patients in C_KEYVAR:", n_distinct(C_KEYVAR$subjid), "\n")
cat("Unique patients in C_BCHAR:", n_distinct(C_BCHAR$subjid), "\n")
cat("Unique patients in A_EENDPT:", n_distinct(A_EENDPT$subjid), "\n")

## NOTE: Replace 'subjid' with whatever the actual patient ID
## variable is called in the trial data (from the glimpse() output above)

## What are the treatment group values?
## this confirms the 1:1 randomization

cat("\nTreatment group distribution (C_KEYVAR):\n")
print(table(C_KEYVAR$txgroup))  # replace 'txgroup' with actual variable name in study

## What are the death indicator values?
cat("\nDeath indicator distribution (A_EENDPT):\n")
print(table(A_EENDPT$dth))  # confirm coding: 1=died? 0=alive?

## Distribution of survival times
cat("\nSurvival time (dthdy) summary:\n")
print(summary(A_EENDPT$dthdy))
## Min should be > 0, Max is the longest follow-up in days
## Median gives a rough sense of follow-up before formal analysis

## ============================================================
## BLOCK 7: BUILD ADSL
## ============================================================
## Because our data lives in separate tables (C_KEYVAR has
## treatment arm, C_BCHAR has age/ECOG, A_EENDPT has OS data),
## we need to combine them into one row-per-patient table.
## This is done with a "join" -- linking datasets by a shared
## variable (the patient ID).
##
## We always use left_join starting from C_KEYVAR (the spine)
## so we never accidentally lose patients or create duplicates.
##
## Why are we building ADSL per SAP Section 4?
## ADSL is the foundation of every analysis. By building it
## explicitly -- rather than just filtering on the fly -- we:
## 1. Have one authoritative definition of each population flag
## 2. Can audit every derivation back to source variables
## 3. Follow the CDISC ADaM discipline CRO reviewers expect
## ============================================================

## --- STEP 7a: Start with C_KEYVAR as the spine ---
## Every patient who was randomized must appear exactly once.
## We rename the treatment variable to TRT01P (planned treatment)
## to match CDISC ADaM ADSL naming convention.
## ITTFL: Intent-to-Treat flag -- ALL randomized patients = "Y"
## (In an RCT, by definition, if you were randomized you are ITT)

ADSL <- C_KEYVAR %>%
  rename(TRT01P = txgroup) %>%      # rename to ADaM convention
  mutate(
    ITTFL  = "Y",                   # all randomized = ITT
    TRT01P = as.character(TRT01P)   # ensure text format
  )

cat("ADSL spine created:", nrow(ADSL), "patients\n")

## --- STEP 7b: Attach baseline characteristics from C_BCHAR ---
## We join AGE, SEX, and ECOG score.

ADSL <- ADSL %>%
  left_join(
    C_BCHAR %>% select(subjid, age, sex, b_ecog2), #adjust names
    by = "subjid"
  )

cat("After C_BCHAR join:", ncol(ADSL), "variables\n")

## --- STEP 7c: Derive Safety Population flag from C_CHEMO ---
## Safety Population = received at least one dose of study drug.
## We summarize C_CHEMO to one row per patient (TRUE if any
## dose record exists), then join onto ADSL.

dose_received <- C_CHEMO %>%
  group_by(subjid) %>%
  summarize(
    any_dose = any(!is.na(txdose) & txdose > 0),  # adjust 'dose' to actual variable name
    .groups = "drop"
  ) %>%
  mutate(SAFFL = ifelse(any_dose, "Y", "N"))

ADSL <- ADSL %>%
  left_join(dose_received %>% select(subjid, SAFFL), by = "subjid") %>%
  mutate(SAFFL = replace_na(SAFFL, "N"))  # no dose record = not in safety population

cat("Safety population (SAFFL=Y):", sum(ADSL$SAFFL == "Y"), "patients\n")
cat("Not in safety pop (SAFFL=N):", sum(ADSL$SAFFL == "N"), "patients\n")

## --- STEP 7d: Attach OS endpoint variables from A_EENDPT ---
## This brings in the two variables we need for survival analysis:
## DTHFL = did the patient die? (1=yes/event, 0=no/censored)
## DTHDY = how many days from randomization to death or censoring?
##
## In ADTTE terms:
## AVAL = DTHDY (analysis value = time in days)
## CNSR = inverse of DTHFL (1=censored, 0=event) -- note the
##        flip: CNSR is the OPPOSITE of the event indicator
##        (0 = event occurred, 1 = censored) in survival package

ADSL <- ADSL %>%
  left_join(
    A_EENDPT %>% select(subjid, dth, dthdy, dthwk),
    by = "subjid"
  ) %>%
  rename(
    DTHFL = dth,    # death flag (1=died, 0=alive/censored)
    DTHDY = dthdy,  # days to death or censoring
    DTHWK = dthwk   # weeks to death or censoring
  )

cat("After A_EENDPT join:", ncol(ADSL), "variables\n")
cat("Total deaths (DTHFL=1):", sum(ADSL$DTHFL == 1, na.rm=TRUE), "\n")
cat("Censored (DTHFL=0):", sum(ADSL$DTHFL == 0, na.rm=TRUE), "\n")

## --- STEP 7e: Create CNSR (censoring indicator for survival pkg) ---
## survival package Cox/KM functions expect:
## event indicator where 1=event, 0=censored (opposite of CNSR)
## OR use CNSR convention where 0=event, 1=censored
## We'll create both for clarity:

ADSL <- ADSL %>%
  mutate(
    OS_EVENT = DTHFL,                  # 1=death (event), 0=censored
    OS_CNSR  = ifelse(DTHFL==1, 0, 1), # 0=event, 1=censored (CDISC CNSR convention)
    OS_DAYS  = DTHDY,                  # survival time in days
    OS_MONTHS = round(DTHDY / 30.4375, 1) # convert to months for presentation
  )

## --- STEP 7f: Label treatment arms clearly ---
## The raw treatment variable likely codes arms as numbers or
## short codes. We recode to clear labels for all tables/figures.
## NOTE: Check your actual coding from Block 6 and update
## the values below (e.g., if 1=Darbepoetin and 2=Placebo,
## keep as is; if coded differently, adjust accordingly)

ADSL <- ADSL %>%
  mutate(
    TRT_LABEL = case_when(
      TRT01P == "1" ~ "Darbepoetin Alfa",
      TRT01P == "2" ~ "Placebo",
      TRUE          ~ as.character(TRT01P) # fallback if coding differs
    )
  )

## ============================================================
## BLOCK 8: VALIDATE AND INSPECT THE COMPLETED ADSL
## ============================================================
## Before writing a single analysis line, we quality-check ADSL.
## Confirming row counts, population flags, no unexpected missing
## values, and sensible distributions.
## This step is called "dataset QC" or "data review"
## and is often done by a second analyst independently.
## ============================================================

cat("\n=================== ADSL VALIDATION ===================\n")
cat("Total patients:", nrow(ADSL), "(expect 600)\n")
cat("Variables in ADSL:", ncol(ADSL), "\n")

cat("\nTreatment arm distribution:\n")
print(table(ADSL$TRT_LABEL))  # should be ~300 each

cat("\nAnalysis population flags:\n")
cat("  ITT (ITTFL=Y):", sum(ADSL$ITTFL == "Y"), "\n")
cat("  Safety (SAFFL=Y):", sum(ADSL$SAFFL == "Y"), "\n")

cat("\nOS events by arm:\n")
print(table(ADSL$TRT_LABEL, ADSL$OS_EVENT))
## Rows = treatment arm, Columns = 0 (censored) / 1 (death)

cat("\nSurvival time (OS_DAYS) summary by arm:\n")
print(
  ADSL %>%
    group_by(TRT_LABEL) %>%
    summarize(
      n = n(),
      events = sum(OS_EVENT, na.rm=TRUE),
      event_rate_pct = round(100 * mean(OS_EVENT, na.rm=TRUE), 1),
      median_days = median(OS_DAYS, na.rm=TRUE),
      max_days = max(OS_DAYS, na.rm=TRUE),
      .groups = "drop"
    )
)

cat("\nMissing values check:\n")
print(colSums(is.na(ADSL)))
## Any variable with unexpected NAs should be investigated before
## proceeding to analysis.

cat("\nSample of completed ADSL (first 5 rows):\n")
print(head(ADSL %>% select(subjid, TRT_LABEL, ITTFL, SAFFL, AGE=age.y,
                            SEX=sex.y, ECOG=b_ecog2.y, OS_EVENT, OS_DAYS, OS_MONTHS), 5))

## ============================================================
## BLOCK 9: SAVE THE ADSL DATASET
## ============================================================
## We save ADSL as both an R data file (.rds) and a CSV.
## The locked ADSL would be saved as a .sas7bdat
## and submitted alongside the analysis code. Our .rds serves
## the same purpose within this R-based analysis.
## ============================================================

saveRDS(ADSL, "ADSL.rds")
write.csv(ADSL, "ADSL.csv", row.names = FALSE)

cat("\nADSL saved as ADSL.rds and ADSL.csv\n")
cat("ADSL construction complete. Ready for analysis.\n")
cat("\nNext step: build ADTTE and run primary OS analysis.\n")
