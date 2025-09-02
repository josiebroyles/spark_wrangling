library(data.table)
library(tidyr)
library(dplyr)
 

####################### setting up #######################
setwd("~/Desktop/sebatlab/SPARK")

# clears your current environment!
rm(list = ls())

####################### setting up #######################

area_deprivation_index_2022_12_12.csv = fread("area_deprivation_index_2022-12-12.csv")
asr_2022_12_12.csv = fread("asr_2022-12-12.csv")
background_history_adult_2022_12_12.csv = fread("background_history_adult_2022-12-12.csv")
background_history_child_2022_12_12.csv = fread("background_history_child_2022-12-12.csv")
#background_history_sibling_2022_12_12.csv = fread("background_history_sibling_2022-12-12.csv")
basic_medical_screening_2022_12_12.csv = fread("basic_medical_screening_2022-12-12.csv")
cbcl_1_5_2022_12_12.csv = fread("cbcl_1_5_2022-12-12.csv")
cbcl_6_18_2022_12_12.csv = fread("cbcl_6_18_2022-12-12.csv")
clinical_lab_results_2022_06_03.csv = fread("clinical_lab_results-2022-06-03.csv")
core_descriptive_variables_2022_12_12.csv = fread("core_descriptive_variables_2022-12-12.csv")
dcdq_2022_12_12.csv = fread("dcdq_2022-12-12.csv")
individuals_registration_2022_12_12.csv = fread("individuals_registration_2022-12-12.csv")
iq_2022_12_12.csv = fread("iq_2022-12-12.csv")
predicted_iq_experimental_2022_12_12.csv = fread("predicted_iq_experimental_2022-12-12.csv")
rbsr_2022_12_12.csv = fread("rbsr_2022-12-12.csv")
roles_index_2022_12_12.csv = fread("roles_index_2022-12-12.csv")
scq_2022_12_12.csv = fread("scq_2022-12-12.csv")
srs_2_adult_self_report_2022_12_12.csv = fread("srs-2_adult_self_report_2022-12-12.csv")
vineland_3_2022_12_12.csv = fread("vineland-3_2022-12-12.csv")


####################### functions #######################
clean_percentile <- function(x) {
  x <- gsub("(th|rd|st|nd)", "", x, ignore.case = TRUE)
  x <- gsub(">", "", x)
  x <- gsub("percentile", "", x, ignore.case = TRUE)
  x <- trimws(x)
  as.numeric(x)
}


####################### cleaning data table #######################

za = data.frame(colnames(background_history_adult_2022_12_12.csv))
zc = data.frame(colnames(background_history_child_2022_12_12.csv))


phens = background_history_child_2022_12_12.csv %>%
  select(age_at_eval_months,contains('age_mos'),contains('regress'),
         hand,contains('cog_age'),contains('sped'),contains('intervention'))

zzsnv = clinical_lab_results_2022_06_03.csv %>%
  select(snv_genetic_status,cnv_genetic_status) %>%
  count(cnv_genetic_status)
#count(snv_genetic_status)


zgen = basic_medical_screening_2022_12_12.csv %>%
  select(starts_with("gen_"))


zgen %>% count(gen_test)

cbcl_phens = cbcl_1_5_2022_12_12.csv %>%
  select(subject_sp_id, 9:last_col())



######## Make into one df ##################
df = clinical_lab_results_2022_06_03.csv %>%
  select(subject_sp_id, snv_genetic_status, cnv_genetic_status)


## Phenotypes from basic_medical_screening ##################
df_bms = basic_medical_screening_2022_12_12.csv %>%
  select(subject_sp_id, visaud_blind, visaud_catar, visaud_strab, med_cond_neuro, 
         neuro_sz, birth_def_gastro, birth_def_urogen, birth_def_bone, 
         dev_ld, dev_lang, dev_lang_dis, mood_dep, behav_adhd, mood_ocd,
         mood_anx, sleep_dx, sleep_probs)

# Make vision category using visaud_blind, visaud_catar, visaud_strab
df_bms[, vision := rowSums(.SD, na.rm = TRUE),
       .SDcols = c("visaud_blind", "visaud_catar", "visaud_strab")]
df_bms[, c("visaud_blind", "visaud_catar", "visaud_strab") := NULL]

# Make language_disorder category using dev_lang, dev_lang_dis
df_bms[, language_disorder := rowSums(.SD, na.rm = TRUE), 
        .SDcols = c("dev_lang", "dev_lang_dis")]
df_bms[, c("dev_lang", "dev_lang_dis") := NULL]

# sleepdisturbancescore: 
df_bms[, sleepdisturbancescore := rowSums(.SD, na.rm = TRUE), 
       .SDcols = c("sleep_dx", "sleep_probs")]
df_bms[, c("sleep_dx", "sleep_probs") := NULL] 

#Now make all of those columns binary
df_bms[, c("vision", "language_disorder", "sleepdisturbancescore") :=
         lapply(.SD, function(x) as.integer(x > 0)),
       .SDcols = c("vision", "language_disorder", "sleepdisturbancescore")]

#All NAs are 0
df_bms[is.na(df_bms)] <- 0

#Rename all columns in df_bms
colnames(df_bms)[colnames(df_bms) == "med_cond_neuro"] <- "neuro"
colnames(df_bms)[colnames(df_bms) == "neuro_sz"] <- "neuro_seizure"
colnames(df_bms)[colnames(df_bms) == "birth_def_gastro"] <- "gastro"
colnames(df_bms)[colnames(df_bms) == "birth_def_urogen"] <- "genital"
colnames(df_bms)[colnames(df_bms) == "birth_def_bone"] <- "orthopedic"
colnames(df_bms)[colnames(df_bms) == "dev_ld"] <- "learning_disability"
colnames(df_bms)[colnames(df_bms) == "mood_dep"] <- "depression"
colnames(df_bms)[colnames(df_bms) == "behav_adhd"] <- "attention_deficit_hyperactivity_disorder"
colnames(df_bms)[colnames(df_bms) == "mood_ocd"] <- "obsessive_compulsive_disorder"
colnames(df_bms)[colnames(df_bms) == "mood_anx"] <- "anxiety_disorder"

## Phenotypes from core_descriptive_variables##################
df_cdv = core_descriptive_variables_2022_12_12.csv %>%
  select(subject_sp_id, cognitive_impairment_latest, viq, regress_other_y_n, regress_lang_y_n)

df_cdv[, autistic_regression := rowSums(.SD, na.rm = TRUE), 
       .SDcols = c("regress_other_y_n", "regress_lang_y_n")]
df_cdv[, c("regress_other_y_n", "regress_lang_y_n") := NULL]
df_cdv$autistic_regression[df_cdv$autistic_regression > 1] <- 1
df_cdv$cognitive_impairment_latest <- ifelse(is.na(df_cdv$cognitive_impairment_latest), NA, as.integer(df_cdv$cognitive_impairment_latest))
colnames(df_cdv)[colnames(df_cdv) == "viq"] <- "minimally_verbal"
colnames(df_cdv)[colnames(df_cdv) == "cognitive_impairment_latest"] <- "intellectual_disability"
df_cdv$minimally_verbal <- ifelse(is.na(df_cdv$minimally_verbal), NA,
                                  ifelse(df_cdv$minimally_verbal < 70, 1, 0))





## Phenotypes from vineland-3##################
df_vine = vineland_3_2022_12_12.csv %>%
  select(subject_sp_id, abc_standard, motor_standard)

########### scq
df_scq = scq_2022_12_12.csv %>% 
  select(subject_sp_id, final_score)
colnames(df_scq)[colnames(df_scq) == "final_score"] <- "scq_life_summary_score"


###########sensory 
df_sensory = rbsr_2022_12_12.csv %>% 
  select(subject_sp_id, q06_sensory)
df_sensory$q06_sensory <- ifelse(df_sensory$q06_sensory > 0, 1,
                                 ifelse(df_sensory$q06_sensory == 0, 0, NA))
colnames(df_sensory)[colnames(df_sensory) == "q06_sensory"] <- "sensory_integration_disorder_or_reported_sensory_issues"

##### Merge all data tables
merged_df <- Reduce(function(x, y) merge(x, y, by = "subject_sp_id", all = TRUE),
                    list(df, df_bms, df_cdv, df_scq, df_sensory, df_vine))

overlap_rows <- df[!is.na(snv_genetic_status) & !is.na(cnv_genetic_status), ]
nrow(overlap_rows)  # how many overlaps

colnames(merged_df) 
colnames(searchlight_df)


# Find the common columns in the order of searchlight_df
common_cols <- intersect(colnames(searchlight_df), colnames(merged_df))

# Remove the ID and the two genetic status columns from the list (we'll put them first)
common_cols <- setdiff(common_cols, c("subject_sp_id", "snv_genetic_status", "cnv_genetic_status"))

# Reorder merged_df: ID + genetic columns first, then matching columns in searchlight order
merged_df <- merged_df[, c("subject_sp_id", "snv_genetic_status", "cnv_genetic_status", common_cols), with = FALSE]


df_for_lda = merged_df %>%
  mutate(
    snv_genetic_status = na_if(snv_genetic_status, ""),
    cnv_genetic_status = na_if(cnv_genetic_status, "")
  ) %>%
  mutate(genetic_status = coalesce(snv_genetic_status, cnv_genetic_status)) %>% 
  select(-snv_genetic_status, -cnv_genetic_status) %>% 
  rename("sfari_id" = "subject_sp_id") %>% 
  select(sfari_id, genetic_status, everything())


###### Save it for use in LDA 
fwrite(df_for_lda, "~/Desktop/sebatlab/spark_wrangling/df_spark_phenotypes.tsv", sep = "\t")
