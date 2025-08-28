###################################################################################################
# Load necessary libraries
###################################################################################################
library(caret)      # For confusion matrix
library(conflicted)
library(data.table) # For high-performance data manipulation
#library(lmerTest)
library(MASS)       # For LDA and other statistical methods
library(scatterplot3d)
library(tidyverse)  # This loads dplyr, ggplot2, and other core tidyverse packages

###################################################################################################
# Set up
###################################################################################################
# Specify preferences for conflicting functions
# I use dplyr functions like these ones. When you load other libraries, sometimes these function names are overridden by some other library's function with the same names.
# So I use this to ensure that the dplyr functions are the ones that are called (otherwise you can just write e.g. dplyr::select).
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("first", "dplyr")
# Remove all variables in current environment
#rm(list = ls())
#setwd("/Users/jamesguevara/sebatlab Dropbox/James Guevara/searchlight_analysis")

setwd("/Users/josie/Desktop/sebatlab/searchlight_multimodal/data/")
searchlight_df = fread("df_phenotypes.tsv") %>% 
  mutate(genetic_status = as.factor(genetic_status)) %>%
  select(-all_of(cols_to_remove))

setwd("/Users/josie/Desktop/sebatlab/spark_wrangling")
df_phenotypes = fread("df_spark_phenotypes.tsv")%>%
  mutate(genetic_status = as.factor(genetic_status)) %>%
  select(all_of(cols)) %>%
  mutate(prop_na = rowMeans(is.na(.))) %>% # Proportion of NAs per row
  filter(prop_na < 0.40) %>% # Remove samples with missingness greater than 0.40
  filter(!is.na(genetic_status) & genetic_status != "") %>% # Remove samples that don't have genetic_status
  select(-prop_na) %>%
  droplevels()

###################################################################################################
# Functions
###################################################################################################
# Function to get total number of values (NA or not) for each column
is_binary <- function(column) {
  unique_vals <- unique(na.omit(column))  # Get unique values, ignoring NAs
  length(unique_vals) == 2  # Return TRUE if there are exactly two unique values
}

# Function to calculate mode
get_mode <- function(column) {
  uniq_vals <- na.omit(column)
  uniq_vals[which.max(tabulate(match(column, uniq_vals)))]
}

impute_column <- function(column) {
  if (is_binary(column)) {
    # Impute with mode for binary columns
    mode_value <- get_mode(column)
    column[is.na(column)] <- mode_value
  } else {
    # Impute with mean for non-binary columns
    mean_value <- mean(column, na.rm = TRUE)
    column[is.na(column)] <- mean_value
  }
  return(column)
}

######### Vectors that may come in use later #################

cols_to_remove <- c(
  "best_latest_walking",
  "best_latest_talking",
  "infectious_disease",
  "surgery",
  "respiratory",
  "heart",
  "kidney",
  "endocrinologic",
  "structural_birth_defects",
  "joint_problems",
  "dermatologic",
  "allergy",
  "autoimmune",
  "immunodeficiency",
  "no_reported_diagnoses",
  "asd",
  "unable_to_walk",
  "developmental_regression_due_to_seizures",
  "urinary_incontinence",
  "developmental_delay_not_id",
  "lang_calc",
  "ketogenic_diet",
  "seizures_in_sleep",
  "ever_needed_emergency_intervention",
  "count_current_medications",
  "childrens_sleep_habits_totalscore",
  "qi_disability_total_score"
)


cols <- c(
  "genetic_status",
  "sfari_id",
  "vision",
  "neuro",
  "neuro_seizure",
  "gastro",
  "genital",
  "orthopedic",
  "learning_disability",
  "language_disorder",
  "intellectual_disability",
  "depression",
  "minimally_verbal",
  "attention_deficit_hyperactivity_disorder",
  "obsessive_compulsive_disorder",
  "anxiety_disorder",
  "sensory_integration_disorder_or_reported_sensory_issues",
  "autistic_regression",
  "scq_life_summary_score",
  "abc_standard",
  "motor_standard",
  "sleepdisturbancescore"
)

########### Summary on SPARK data ############## 
#### Shows that there is very limited data for minimally_verbal and motor_standard #######

# Create summary table (before imputation)
summary_table = df_phenotypes %>%
  summarize(
    num_phenotypes = ncol(df_phenotypes) - 2,
    unique_genetic_status = n_distinct(genetic_status),
    num_samples = n()
  )
print(summary_table)

# Number of samples per genetic_status
samples_per_status <- df_phenotypes %>%
  group_by(genetic_status) %>%
  summarize(samples_per_genetic_status = n())

# Coverage per phenotype 
df_phenotype_coverage = df_phenotypes %>%
  group_by(genetic_status) %>%
  summarise(across(everything(), ~ mean(!is.na(.)), .names = "{col}" ))
# Reshape the dataframe from wide to long format
df_long_original <- df_phenotype_coverage %>%
  pivot_longer(cols = -genetic_status,  # Keep genetic_status as is, reshape all other columns
               names_to = "Variable",   # New column name for former column names
               values_to = "Value")     # New column name for values
df_long_original$Variable <- factor(df_long_original$Variable, levels = unique(df_long_original$Variable))
# Plot the heatmap
ggplot(df_long_original, aes(x = Variable, y = genetic_status, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  labs(x = "Variables", y = "Genetic Status", fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))




# Heatmap of phenotype measure files per genetic_status
 carrier_count_by_measure = fread("/Users/josie/Desktop/sebatlab/searchlight_multimodal/data/carrier_count_by_measure.csv") %>% distinct()
# Melt the data
melted_data = melt(carrier_count_by_measure, id.vars = "measure")
searchlight_core_descriptive_variables_gene_counts = carrier_count_by_measure %>%
  filter(measure == "searchlight_core_descriptive_variables") %>%
  pivot_longer(cols = -1, names_to = "gene", values_to = "count") %>%
  select(gene, count) %>%
  arrange(desc(count))
gene_order <- searchlight_core_descriptive_variables_gene_counts$gene
# Find the position where total count falls below 100
threshold_position_100 <- which(searchlight_core_descriptive_variables_gene_counts$count < 100)[1]
threshold_position_50  <- which(searchlight_core_descriptive_variables_gene_counts$count < 50)[1]
threshold_position_25  <- which(searchlight_core_descriptive_variables_gene_counts$count < 25)[1]
# Calculate total counts for each measure (phenotype) and order them
measure_order <- melted_data %>%
  group_by(measure) %>%
  summarize(total = sum(value, na.rm = TRUE)) %>%
  arrange(total) %>%
  pull(measure)
# Convert variable and measure to factors with levels ordered by total counts
melted_data$variable <- factor(melted_data$variable, levels = gene_order)
melted_data$measure  <- factor(melted_data$measure, levels = measure_order)
# Create the heatmap
ggplot(melted_data, aes(x = variable, y = measure, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_vline(xintercept = threshold_position_100, color = "blue") +
  geom_vline(xintercept = threshold_position_50, color = "green") +
  geom_vline(xintercept = threshold_position_25, color = "orange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes (ordered by count)",
       y = "Measures (ordered by count)",
       fill = "Count",
       title = "Heatmap of Gene-Measure Counts",
       subtitle = "Blue line indicates where gene count falls below 100\nPurple line indicates where gene count falls below 50\nOrange line indicates where gene count falls below 25"
  )

###################################################################################################
# Basic imputation
###################################################################################################
# Mean (or mode) imputation
df_imputed = df_phenotypes %>%
  mutate(across(
    .cols = -c(genetic_status, sfari_id),  # Exclude first 2 columns and specific columns
    .fns = impute_column
  ))
tmp_is_factors = data.frame(sapply(df_imputed, is.factor)) # Ensure genetic_status is a factor

searchlight_imputed <- searchlight_df %>% 
  mutate(across(-c(genetic_status, sfari_id), impute_column)) %>%
  filter(!is.na(genetic_status), genetic_status != "") %>%
  mutate(genetic_status = factor(genetic_status),         # make it a factor
         genetic_status = droplevels(genetic_status))     # drop unused levels


###################################################################################################
# LDA (Trained on Spark)
###################################################################################################

# Run PCA on phenotypes ("phenotypic PCA")
pca_result = prcomp(df_imputed %>% select(-genetic_status, -sfari_id), scale. = TRUE)
df_imputed_with_pcs = cbind(df_imputed, pca_result$x)

# Run LDA on phenotypes
lda_model = lda(genetic_status ~ ., data = df_imputed %>% select(-sfari_id) ) 
lda_pred = predict(lda_model) # Predict (on the same data that was trained on...)
# Create a confusion matrix to evaluate the predictions
confusion_matrix = confusionMatrix(lda_pred$class, df_imputed$genetic_status)
# Visualize the confusion matrix
print(confusion_matrix)

# Create a combined data.table with the genetic_status, sfari_id, phenotypes, PCs, and LDs, and demographics...
df_imputed_with_pcs_lds = cbind(df_imputed_with_pcs, lda_pred$x)
is.data.table(df_imputed_with_pcs_lds)

###################################################################################################
# Spark trained LDA tested on Spark data | Accuracy : 0.4706   
###################################################################################################
# 1. Predict on df_imputed
lda_pred_test <- predict(lda_model, newdata = df_imputed %>% select(-sfari_id, -genetic_status))

# 2. Align factor levels between prediction and reference
common_levels <- intersect(levels(lda_pred_test$class), levels(df_imputed$genetic_status))
pred_aligned <- factor(lda_pred_test$class, levels = common_levels)
ref_aligned  <- factor(df_imputed$genetic_status, levels = common_levels)

# 3. Compute confusion matrix
confusion_matrix <- confusionMatrix(pred_aligned, ref_aligned)
print(confusion_matrix)

###################################################################################################
# Spark trained LDA tested on Searchlight data | Accuracy : 0
###################################################################################################

# Define predictors (excluding ID and outcome)
predictor_cols <- setdiff(colnames(searchlight_imputed), c("sfari_id", "genetic_status"))

# 1. Predict on searchlight_imputed
lda_pred_test <- predict(lda_model, newdata = searchlight_imputed %>% select(-sfari_id, -genetic_status))

# 2. Align factor levels between prediction and reference
common_levels <- intersect(levels(lda_pred_test$class), levels(searchlight_imputed$genetic_status))
pred_aligned <- factor(lda_pred_test$class, levels = common_levels)
ref_aligned  <- factor(searchlight_imputed$genetic_status, levels = common_levels)

# 3. Compute confusion matrix
confusion_matrix <- confusionMatrix(pred_aligned, ref_aligned)
print(confusion_matrix)

###################################################################################################
# LDA (Trained on Searchlight)
###################################################################################################

# Run PCA on phenotypes ("phenotypic PCA")
pca_result = prcomp(searchlight_imputed %>% select(-genetic_status, -sfari_id), scale. = TRUE)
df_imputed_with_pcs = cbind(searchlight_imputed, pca_result$x)

# Run LDA on phenotypes
lda_model = lda(genetic_status ~ ., data = searchlight_imputed %>% select(-sfari_id) ) 
lda_pred = predict(lda_model) # Predict (on the same data that was trained on...)
# Create a confusion matrix to evaluate the predictions
confusion_matrix = confusionMatrix(lda_pred$class, searchlight_imputed$genetic_status)
# Visualize the confusion matrix
print(confusion_matrix)

# Create a combined data.table with the genetic_status, sfari_id, phenotypes, PCs, and LDs, and demographics...
df_imputed_with_pcs_lds = cbind(df_imputed_with_pcs, lda_pred$x)
is.data.table(df_imputed_with_pcs_lds)

###################################################################################################
# Searchlight trained LDA tested on Spark data | Accuracy : 0.1818
###################################################################################################

# Define predictors (excluding ID and outcome)
predictor_cols <- setdiff(colnames(df_imputed), c("sfari_id", "genetic_status"))

# 1. Predict on searchlight_imputed
lda_pred_test <- predict(lda_model, newdata = df_imputed %>% select(-sfari_id, -genetic_status))

# 2. Align factor levels between prediction and reference
common_levels <- intersect(levels(lda_pred_test$class), levels(df_imputed$genetic_status))
pred_aligned <- factor(lda_pred_test$class, levels = common_levels)
ref_aligned  <- factor(df_imputed$genetic_status, levels = common_levels)

# 3. Compute confusion matrix
confusion_matrix <- confusionMatrix(pred_aligned, ref_aligned)
print(confusion_matrix)

