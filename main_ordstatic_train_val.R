library(readr)
library(MASS)
library(pROC)
library(ggplot2)
library(foreach)
library(doParallel)
source("functions_with_median.R")

# Declare features and units

features <- c(
  "qnet",	"inal_auc",	"ical_auc",	"apd90",	"apd50",	"apd_tri",	"cad90",	"cad50",
  "cad_tri",	"dvmdt_repol",	"dvmdt_peak",	"vm_peak",	"vm_dia",	"ca_peak",	"ca_dia"
)

units <-c ("", "", "", "", "", "", "", "",
           "", "", "", "", "", "", "")

# Declare paths for the file input
filepath_training <- "data/ordstatic_zscore_train.csv"
filepath_testing <- "data/ordstatic_zscore_val.csv"

# Set feature dimension (how many feature pair)
dimension <-1  # can be changed to 2 (double-biomarker); or 3; 4; 5; 6

# Create pairsdf with all unique combinations
pairsdf <- pairsdfinitfun(features = features, units = units, dimension = dimension)

# Create the results folder
results_folder <- "results/single_biom"

# Choose whether data needs to be normalized.
is_normalized <- FALSE  #

# Choose how many attempts are required before skipping the fitting process (when divergence is found)
max_attempts <- 10000 # 1000

# Choose how many tests are required for evaluating model performance
num_tests <- 10000 #1, 1000, 10000

# Register parallel backend
numCores <- 6

# Check if the folder exists
if (!dir.exists(results_folder)) {
  # The folder does not exist, so create it
  dir.create(results_folder, recursive = TRUE)
  cat("Folder created:", results_folder, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = results_folder, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}

# Execute the tasks in parallel
cl <- makeCluster(numCores)
registerDoParallel(cl)

is_loocv <- FALSE
summarydf <- foreach( pair_id = 1:nrow(pairsdf),
                      .combine = 'rbind',
                      .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
                        # Select the row corresponding to pair_id
                        pair_row <- pairsdf[pair_id, ]
                        
                        # Extract vectors of features and units
                        features_vector <- pair_row[grep("feature_", names(pair_row))]
                        units_vector <- pair_row[grep("unit_", names(pair_row))]
                        
                        # Print the current pair_id and features being processed
                        print(paste(c(pair_id, features_vector), collapse = "_"))
                        
                        # Read in the training and testing datasets
                        training <- read_csv(filepath_training, show_col_types = FALSE)
                        testing <- read_csv(filepath_testing, show_col_types = FALSE)
                        
                        result <- run_all_testing (results_folder = results_folder,
                                                    training = training,
                                                    testing = testing,
                                                    features_vector = features_vector,
                                                    units_vector = units_vector,
                                                    is_normalized = is_normalized,
                                                    max_attempts = max_attempts,
                                                    num_tests = num_tests)
                        
                        # Return the result
                        result
                      }
# -------------------------------------------------------
# Split wide summarydf into two output CSV files:
#   summary.csv        — worst-case CI metrics (original behaviour)
#   summary_median.csv — median (50th percentile) metrics
# Shared metadata columns (logLik, Alpha_i, Beta_i) appear in both files.
# -------------------------------------------------------
med_cols  <- grep("_med$", names(summarydf), value = TRUE)
base_cols <- setdiff(names(summarydf), med_cols)   # all non-median columns

# -------------------------------------------------------
# Desired column order for BOTH output files.
# "feature_col_name" is used as a placeholder for the first column
# (actual name depends on dimension, e.g. "feature_1").
# -------------------------------------------------------
desired_order <- c(
  "FEATURE_COL_PLACEHOLDER",
  "logLik", "Alpha_1", "Alpha_2", "Beta_1",
  "AUC1", "AUC2", "pairwise",
  "LR_pos_th1", "LR_pos_th2", "LR_neg_th1", "LR_neg_th2",
  "Classification_error",
  "Sensitivity1", "Sensitivity2",
  "Specificity1", "Specificity2",
  "Accuracy_th1",    "Accuracy_th2",
  "F1score_th1",     "F1score_th2",
  "Rank_score",
  "AUC1_training", "AUC2_training",
  "LR_pos_th1_training", "LR_pos_th2_training",
  "LR_neg_th1_training", "LR_neg_th2_training",
  "Classification_error_training",
  "Sensitivity1_training", "Sensitivity2_training",
  "Specificity1_training", "Specificity2_training",
  "F1score_th1_training",     "F1score_th2_training",
  "Rank_score_training"
)

# Helper: reorder a dataframe by desired_order, appending any leftover cols
reorder_cols <- function(df, desired, feature_col) {
  ordered  <- sub("FEATURE_COL_PLACEHOLDER", feature_col, desired)
  present  <- intersect(ordered, names(df))
  leftover <- setdiff(names(df), present)
  df[, c(present, leftover), drop = FALSE]
}

# --- summary.csv ---
feature_col_name <- names(summarydf[, base_cols])[1]   # e.g. "feature_1"
summarydf_out    <- reorder_cols(summarydf[, base_cols], desired_order, feature_col_name)
write.csv(summarydf_out,
          file.path(results_folder, "summary.csv"),
          row.names = FALSE)

# --- summary_median.csv ---
meta_cols        <- intersect(
  c("logLik", grep("^Alpha_|^Beta_", names(summarydf), value = TRUE)),
  names(summarydf)
)
median_sel_cols  <- intersect(c(feature_col_name, meta_cols, med_cols), names(summarydf))
summarydf_median <- summarydf[, median_sel_cols, drop = FALSE]
names(summarydf_median) <- sub("_med$", "", names(summarydf_median))

# Rename Classification_error_median -> Classification_error to match summary.csv
names(summarydf_median)[names(summarydf_median) == "Classification_error_median"] <-
  "Classification_error"
names(summarydf_median)[names(summarydf_median) == "Classification_error_median_training"] <-
  "Classification_error_training"

summarydf_median <- reorder_cols(summarydf_median, desired_order, feature_col_name)
write.csv(summarydf_median,
          file.path(results_folder, "summary_median.csv"),
          row.names = FALSE)

stopCluster(cl)
