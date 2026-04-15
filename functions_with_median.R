# Function to calculate sigmoid
sigmoid <- function(x){
  return(1.0/(1.0 + exp(-x)))
}

# Function to calculate variable B
calculate_B1 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(sigmoid(z1) - 0.5 * sigmoid(z2))
}

calculate_B2 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(1.0 + sigmoid(z1) - 2.0 * sigmoid(z2))
}

run_all_testing <- function(results_folder = 'results',
                            training = 'trainingdf',
                            testing = 'testingdf',
                            features_vector = 'features',
                            units_vector = 'units',
                            is_normalized = FALSE,
                            max_attempts = 1000,
                            num_tests = 10000){
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features_vector) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename,"training","log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  training <- training[, training_cols]
  testing <- testing[, testing_cols]
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  training$label <- as.integer(factor(training$label, levels = levels, labels = values))
  testing$label <- as.integer(factor(testing$label, levels = levels, labels = values))
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  # Remove errors
  training <- na.omit(training)
  testing <- na.omit(testing)
  
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first m columns of ordinaldf
    means <- sapply(training[1:dimension], mean)
    sds <- sapply(training[1:dimension], sd)
    tempdf <- training
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, training[1:dimension], means, sds)
    training <- tempdf
    tempdf <- testing
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, testing[1:dimension], means, sds)
    testing <- tempdf
  }
  
  # =======================================================
  # Testing performance stochastically depending on "num_tests"
  # =======================================================
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)
  
  # Specify the number of attempts
  fit_results <- fitandgetTMS(training_data = training,
                              testing_data = testing,
                              formula = formula,
                              dimension = dimension,
                              max_attempts = max_attempts,
                              logfile = logfile)
  converged <- fit_results$converged
  all_or_loocv <- FALSE  # defined here so it is available in both converged and non-converged branches
  if (converged) {
    training <- fit_results$training_data
    testing <- fit_results$testing_data
    mod <- fit_results$mod
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    training <- cbind(training,predict(mod, newdata = training, type = "probs"))
    training["pred"] <- predict(mod, newdata = training, type = "class")
    testing <- cbind(testing,predict(mod, newdata = testing, type = "probs"))
    testing["pred"] <- predict(mod, newdata = testing, type = "class")
    if (dimension == 2) {
      scatterplotfun(data = training,
                     mod = mod,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = TRUE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
      scatterplotfun(data = testing,
                     mod = mod,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = FALSE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
    }
    # Constants of the TMS
    # if (!is_single) {
    th1 <- (alphas[2] - alphas[1] ) / 2.0 + log(1.0 - 2.0 * exp(alphas[1] - alphas[2]))
    th2 <- -th1
    
    # Plot TMS for training dataset
    label_colors <- c("low" = "green4", "intermediate" = "orange", "high" = "red")
    training$risk <- ifelse(training$label == 1, "low",
                            ifelse(training$label == 2, "intermediate", "high"))
    
    # Plot TMS for testing dataset
    testing$risk <- ifelse(testing$label == 1, "low",
                           ifelse(testing$label == 2, "intermediate", "high"))
    
    # if (!is_single) {
    tms_name <- "TMS"
    filename_training <- paste(results_folder,"/",paste(log_filename, "training_development_tms.jpg", sep="_"), sep = "")
    filename_testing <- paste(results_folder,"/",paste(log_filename, "testing_development_tms.jpg", sep="_"), sep = "")
    
    # Plot TMS for training dataset
    tmsplotfun(data = training,
               th1 = th1,
               th2 = th2,
               label_colors = label_colors,
               title = "Training dataset",
               file_name = filename_training,
               tms_name = "TMS")
    tmsplotfun(data = testing,
               th1 = th1,
               th2 = th2,
               label_colors = label_colors,
               title = "Testing dataset",
               file_name = filename_testing,
               tms_name = "TMS")
  } else {
    if (dimension == 2) {
      scatterplotfun(data = training,
                     mod = NA,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = TRUE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
      scatterplotfun(data = testing,
                     mod = NA,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = FALSE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
    }
    temp_return <- returnNA(log_filename,dimension)
    temp_return["th1"] <- NA
    temp_return["th2"] <- NA
    temp_return["Rank_score_training"] <- NA
    temp_return["Accuracy_th1_training"] <- NA
    temp_return["Accuracy_th2_training"] <- NA
    temp_return["AUC1_training"] <- NA
    temp_return["AUC2_training"] <- NA
    temp_return["Sensitivity1_training"] <- NA
    temp_return["Sensitivity2_training"] <- NA
    temp_return["Specificity1_training"] <- NA
    temp_return["Specificity2_training"] <- NA
    temp_return["LR_pos_th1_training"] <- NA
    temp_return["LR_pos_th2_training"] <- NA
    temp_return["LR_neg_th1_training"] <- NA
    temp_return["LR_neg_th2_training"] <- NA
    temp_return["F1score_th1_training"] <- NA
    temp_return["F1score_th2_training"] <- NA
    temp_return["Classification_error_training"] <- NA
    
    # Median columns — all NA when model did not converge
    temp_return["Rank_score_med"]                            <- NA
    temp_return["Accuracy_th1_med"]                       <- NA
    temp_return["Accuracy_th2_med"]                       <- NA
    temp_return["AUC1_med"]                            <- NA
    temp_return["AUC2_med"]                            <- NA
    temp_return["Sensitivity1_med"]                    <- NA
    temp_return["Sensitivity2_med"]                    <- NA
    temp_return["Specificity1_med"]                    <- NA
    temp_return["Specificity2_med"]                    <- NA
    temp_return["LR_pos_th1_med"]                    <- NA
    temp_return["LR_pos_th2_med"]                    <- NA
    temp_return["LR_neg_th1_med"]                    <- NA
    temp_return["LR_neg_th2_med"]                    <- NA
    temp_return["F1score_th1_med"]                        <- NA
    temp_return["F1score_th2_med"]                        <- NA
    temp_return["Classification_error_median_med"]           <- NA
    temp_return["Classification_error_mean_med"]             <- NA
    temp_return["pairwise_med"]      <- NA
    temp_return["Rank_score_training_med"]                   <- NA
    temp_return["Accuracy_th1_training_med"]              <- NA
    temp_return["Accuracy_th2_training_med"]              <- NA
    temp_return["AUC1_training_med"]                   <- NA
    temp_return["AUC2_training_med"]                   <- NA
    temp_return["Sensitivity1_training_med"]           <- NA
    temp_return["Sensitivity2_training_med"]           <- NA
    temp_return["Specificity1_training_med"]           <- NA
    temp_return["Specificity2_training_med"]           <- NA
    temp_return["LR_pos_th1_training_med"]           <- NA
    temp_return["LR_pos_th2_training_med"]           <- NA
    temp_return["LR_neg_th1_training_med"]           <- NA
    temp_return["LR_neg_th2_training_med"]           <- NA
    temp_return["F1score_th1_training_med"]               <- NA
    temp_return["F1score_th2_training_med"]               <- NA
    temp_return["Classification_error_median_training_med"]  <- NA
    temp_return["Classification_error_mean_training_med"]    <- NA
    
    return(temp_return)
  }
  
  # Preallocate the pmeasures dataframe for testing drugs
  pmeasures_testing <- vector("list", num_tests)
  
  # Combine all data
  all_data <- rbind(training,testing)
  
  # Calculate the classification error from the whole dataset
  predicted_labels <- predict(mod, newdata = testing, type = "class")
  pred_err_testing <- (abs(as.integer(predicted_labels) - as.integer(testing$label)))
  
  # Calculate the pmeasures for the whole dataset
  training_drugs <- unique(training$drug_name)
  testing_drugs <- unique(testing$drug_name)
  training$is_training <- 1
  testing$is_training <- 0
  # Iterate through the num_tests tests
  set.seed(1)
  for (i in 1:num_tests) {
    # Sample rows for testing data
    testing_data <- do.call(rbind, lapply(unique(testing$drug_name), function(drug) {
      drug_data <- testing[testing$drug_name == drug, ]
      drug_data[sample(nrow(drug_data), 1), ]
    }))
    
    # Sample rows for training data
    training_data <- do.call(rbind, lapply(unique(training$drug_name), function(drug) {
      drug_data <- training[training$drug_name == drug, ]
      drug_data[sample(nrow(drug_data), 1), ]
    }))
    
    temp_data <- rbind(training_data, testing_data)
    
    # Calculate all pmeasures
    pmeasures_testing[[i]] <- pmeasuresfun_loocv(data = temp_data,
                                                 label_values = values,
                                                 all_or_loocv = FALSE,
                                                 pairwise_test = TRUE)
  } # i-th test
  
  # Combine all pmeasures results at once
  pmeasures_testing <- do.call(rbind, pmeasures_testing)
  
  writepmeasuresfun(pmeasures = pmeasures_testing,
                    pred_error = pred_err_testing,
                    logfile = logfile,
                    all_or_loocv = all_or_loocv,
                    is_loocv = is_loocv,
                    th1 = th1,
                    th2 = th2)
  
  # Compute pmeasures for whole training data
  pmeasures_training <- pmeasuresfun_loocv(data = training,
                                           label_values = values,
                                           all_or_loocv = TRUE,
                                           pairwise_test = FALSE)
  
  # Calculate the classification error from the whole training dataset
  predicted_labels_training <- predict(mod, newdata = training, type = "class")
  pred_err_training <- (abs(as.integer(predicted_labels_training) - as.integer(training$label)))
  
  writepmeasuresfun(pmeasures = pmeasures_training,
                    pred_error = pred_err_training,
                    logfile = logfile,
                    all_or_loocv = TRUE,
                    is_loocv = is_loocv,
                    th1 = th1,
                    th2 = th2)
  
  # Close the connection to the text file
  close(logfile)
  
  # Rank score for the OLR model
  rank_score_testing <- rankscorefun(pmeasures = pmeasures_testing,
                                     pred_error = pred_err_testing,
                                     is_normalized = TRUE)
  rank_score_training <- rankscorefun_training(pmeasures = pmeasures_training,
                                               pred_error = pred_err_training,
                                               is_normalized = TRUE)
  
  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename
  
  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_th1 = quantile(pmeasures_testing$Accuracy_th1, 0.025),
    Accuracy_th2 = quantile(pmeasures_testing$Accuracy_th2, 0.025),
    AUC1 = quantile(pmeasures_testing$AUC1, 0.025),
    AUC2 = quantile(pmeasures_testing$AUC2, 0.025),
    Sensitivity1 = quantile(pmeasures_testing$Sensitivity1, 0.025),
    Sensitivity2 = quantile(pmeasures_testing$Sensitivity2, 0.025),
    Specificity1 = quantile(pmeasures_testing$Specificity1, 0.025),
    Specificity2 = quantile(pmeasures_testing$Specificity2, 0.025),
    LR_pos_th1 = quantile(pmeasures_testing$LR_pos_th1, 0.025),
    LR_pos_th2 = quantile(pmeasures_testing$LR_pos_th2, 0.025),
    LR_neg_th1 = quantile(pmeasures_testing$LR_neg_th1, 0.975),
    LR_neg_th2 = quantile(pmeasures_testing$LR_neg_th2, 0.975),
    F1score_th1 = quantile(pmeasures_testing$F1score_th1, 0.025),
    F1score_th2 = quantile(pmeasures_testing$F1score_th2, 0.025),
    Classification_error = mean(pred_err_testing) + 1.96 * sd(pred_err_testing) / sqrt(length(pred_err_testing)),
    pairwise = quantile(pmeasures_testing$Pairwise, 0.025),
    Rank_score = rank_score_testing
  )
  # Add normalized log likelihood value
  summarydf['logLik'] <- logLik(mod)/nrow(training)
  # Add TMS thresholds
  summarydf['th1'] <- th1
  summarydf['th2'] <- th2
  
  # Add Alphas
  for (i in 1:length(alphas)) {
    summarydf[[paste0("Alpha_", i)]] <- alphas[i]
  }
  
  # Add Betas
  for (i in 1:length(betas)) {
    summarydf[[paste0("Beta_", i)]] <- betas[i]
  }
  
  # Add training scores
  summarydf["Rank_score_training"] <- rank_score_training
  summarydf["Accuracy_th1_training"] <- quantile(pmeasures_training$Accuracy_th1, 0.025)
  summarydf["Accuracy_th2_training"] <- quantile(pmeasures_training$Accuracy_th2, 0.025)
  summarydf["AUC1_training"] <- quantile(pmeasures_training$AUC1, 0.025)
  summarydf["AUC2_training"] <- quantile(pmeasures_training$AUC2, 0.025)
  summarydf["Sensitivity1_training"] <- quantile(pmeasures_training$Sensitivity1, 0.025)
  summarydf["Sensitivity2_training"] <- quantile(pmeasures_training$Sensitivity2, 0.025)
  summarydf["Specificity1_training"] <- quantile(pmeasures_training$Specificity1, 0.025)
  summarydf["Specificity2_training"] <- quantile(pmeasures_training$Specificity2, 0.025)
  summarydf["LR_pos_th1_training"] <- quantile(pmeasures_training$LR_pos_th1, 0.025)
  summarydf["LR_pos_th2_training"] <- quantile(pmeasures_training$LR_pos_th2, 0.025)
  summarydf["LR_neg_th1_training"] <- quantile(pmeasures_training$LR_neg_th1, 0.975)
  summarydf["LR_neg_th2_training"] <- quantile(pmeasures_training$LR_neg_th2, 0.975)
  summarydf["F1score_th1_training"] <- quantile(pmeasures_training$F1score_th1, 0.025)
  summarydf["F1score_th2_training"] <- quantile(pmeasures_training$F1score_th2, 0.025)
  summarydf["Classification_error_training"] <- mean(pred_err_training) + 1.96 * sd(pred_err_training) / sqrt(length(pred_err_training))
  
  # ================================================
  # Median-based metrics (testing & training)
  # ------------------------------------------------
  # These columns carry a "_med" suffix so training_and_validation.R
  # can cleanly split the wide row into summary.csv and summary_median.csv.
  # ------------------------------------------------
  rank_score_testing_med  <- rankscorefun_median(pmeasures  = pmeasures_testing,
                                                 pred_error = pred_err_testing,
                                                 is_normalized = TRUE)
  rank_score_training_med <- rankscorefun_training_median(pmeasures  = pmeasures_training,
                                                          pred_error = pred_err_training,
                                                          is_normalized = TRUE)
  
  # Testing — median metrics
  summarydf["Rank_score_med"]                            <- rank_score_testing_med
  summarydf["Accuracy_th1_med"]                       <- quantile(pmeasures_testing$Accuracy_th1,   0.5)
  summarydf["Accuracy_th2_med"]                       <- quantile(pmeasures_testing$Accuracy_th2,   0.5)
  summarydf["AUC1_med"]                            <- quantile(pmeasures_testing$AUC1,        0.5)
  summarydf["AUC2_med"]                            <- quantile(pmeasures_testing$AUC2,        0.5)
  summarydf["Sensitivity1_med"]                    <- quantile(pmeasures_testing$Sensitivity1, 0.5)
  summarydf["Sensitivity2_med"]                    <- quantile(pmeasures_testing$Sensitivity2, 0.5)
  summarydf["Specificity1_med"]                    <- quantile(pmeasures_testing$Specificity1, 0.5)
  summarydf["Specificity2_med"]                    <- quantile(pmeasures_testing$Specificity2, 0.5)
  summarydf["LR_pos_th1_med"]                    <- quantile(pmeasures_testing$LR_pos_th1, 0.5)
  summarydf["LR_pos_th2_med"]                    <- quantile(pmeasures_testing$LR_pos_th2, 0.5)
  summarydf["LR_neg_th1_med"]                    <- quantile(pmeasures_testing$LR_neg_th1, 0.5)
  summarydf["LR_neg_th2_med"]                    <- quantile(pmeasures_testing$LR_neg_th2, 0.5)
  summarydf["F1score_th1_med"]                        <- quantile(pmeasures_testing$F1score_th1,    0.5)
  summarydf["F1score_th2_med"]                        <- quantile(pmeasures_testing$F1score_th2,    0.5)
  summarydf["Classification_error_median_med"]           <- median(pred_err_testing)
  summarydf["Classification_error_mean_med"]             <- mean(pred_err_testing)
  summarydf["pairwise_med"]      <- quantile(pmeasures_testing$Pairwise,          0.5)
  
  # Training — median metrics
  summarydf["Rank_score_training_med"]                   <- rank_score_training_med
  summarydf["Accuracy_th1_training_med"]              <- quantile(pmeasures_training$Accuracy_th1,   0.5)
  summarydf["Accuracy_th2_training_med"]              <- quantile(pmeasures_training$Accuracy_th2,   0.5)
  summarydf["AUC1_training_med"]                   <- quantile(pmeasures_training$AUC1,        0.5)
  summarydf["AUC2_training_med"]                   <- quantile(pmeasures_training$AUC2,        0.5)
  summarydf["Sensitivity1_training_med"]           <- quantile(pmeasures_training$Sensitivity1, 0.5)
  summarydf["Sensitivity2_training_med"]           <- quantile(pmeasures_training$Sensitivity2, 0.5)
  summarydf["Specificity1_training_med"]           <- quantile(pmeasures_training$Specificity1, 0.5)
  summarydf["Specificity2_training_med"]           <- quantile(pmeasures_training$Specificity2, 0.5)
  summarydf["LR_pos_th1_training_med"]           <- quantile(pmeasures_training$LR_pos_th1, 0.5)
  summarydf["LR_pos_th2_training_med"]           <- quantile(pmeasures_training$LR_pos_th2, 0.5)
  summarydf["LR_neg_th1_training_med"]           <- quantile(pmeasures_training$LR_neg_th1, 0.5)
  summarydf["LR_neg_th2_training_med"]           <- quantile(pmeasures_training$LR_neg_th2, 0.5)
  summarydf["F1score_th1_training_med"]               <- quantile(pmeasures_training$F1score_th1,    0.5)
  summarydf["F1score_th2_training_med"]               <- quantile(pmeasures_training$F1score_th2,    0.5)
  summarydf["Classification_error_median_training_med"]  <- median(pred_err_training)
  summarydf["Classification_error_mean_training_med"]    <- mean(pred_err_training)
  
  return(summarydf)
}

run_all_training <- function(results_folder = 'results',
                             training = 'trainingdf',
                             testing = 'testingdf',
                             features_vector = 'features',
                             units_vector = 'units',
                             is_normalized = FALSE,
                             max_attempts = 1000){
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features_vector) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename,"training","log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  training <- training[, training_cols]
  testing <- testing[, testing_cols]
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  training$label <- as.integer(factor(training$label, levels = levels, labels = values))
  testing$label <- as.integer(factor(testing$label, levels = levels, labels = values))
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  # Remove errors
  training <- na.omit(training)
  testing <- na.omit(testing)
  
  # Combine all data
  all_data <- rbind(training,testing)
  
  # Normalize the data
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first 'dimension' columns of all_data
    means <- sapply(all_data[1:dimension], mean)
    sds <- sapply(all_data[1:dimension], sd)
    tempdf <- all_data
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, all_data[1:dimension], means, sds)
    all_data <- tempdf
  }
  
  # ================================
  # Training performance on all data
  # ================================
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)
  
  # Specify the number of attempts
  # max_attempts <- 10000
  fit_results <- fitandgetTMS(training_data = all_data,
                              testing_data = all_data,
                              formula = formula,
                              dimension = dimension,
                              max_attempts = max_attempts,
                              logfile = logfile)
  converged <- fit_results$converged
  if (converged) {
    all_data <- fit_results$training_data
    mod <- fit_results$mod
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    all_data <- cbind(all_data,predict(mod, newdata = all_data, type = "probs"))
    all_data["pred"] <- predict(mod, newdata = all_data, type = "class")
    all_or_loocv <- TRUE
    if (dimension == 2) {
      scatterplotfun(data = all_data,
                     mod = mod,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = TRUE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
    }
    # Constants of the TMS
    # if (!is_single) {
    th1 <- (alphas[2] - alphas[1] ) / 2.0 + log(1.0 - 2.0 * exp(alphas[1] - alphas[2]))
    th2 <- -th1
    
    # Plot TMS for training dataset
    label_colors <- c("low" = "green4", "intermediate" = "orange", "high" = "red")
    all_data$risk <- ifelse(all_data$label == 1, "low",
                            ifelse(all_data$label == 2, "intermediate", "high"))
    
    # if (!is_single) {
    tms_name <- "TMS"
    filename_training <- paste(results_folder,"/",paste(log_filename, "dataset_tms.jpg", sep="_"), sep = "")
    
    # Plot TMS for training dataset
    tmsplotfun(data = all_data,
               th1 = th1,
               th2 = th2,
               label_colors = label_colors,
               title = "Training all dataset",
               file_name = filename_training,
               tms_name = "TMS")
  } else {
    if (dimension == 2) {
      scatterplotfun(data = all_data,
                     mod = NA,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = TRUE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
    }
    return(returnNA(log_filename,dimension))
  }
  
  # Preallocate the pmeasures dataframe
  pmeasures <- data.frame()
  
  # Calculate the classification error from the whole training dataset
  predicted_labels <- predict(mod, newdata = all_data, type = "class")
  pred_err <- (abs(as.integer(predicted_labels) - as.integer(all_data$label)))
  
  # Calculate the pmeasures for the whole dataset
  training_drugs <- unique(training$drug_name)
  testing_drugs <- unique(testing$drug_name)
  all_data$is_training <- ifelse(all_data$drug_name %in% training_drugs, 1, 0)
  for (sample_id in training$Sample_ID[training$drug_name == training_drugs[1]]) {
    temp_data <- subset(all_data, all_data$Sample_ID == sample_id)
    temppmeasure <- pmeasuresfun_loocv(data = temp_data,
                                       label_values = values,
                                       all_or_loocv = TRUE,
                                       pairwise_test = TRUE)
    pmeasures <- rbind(pmeasures,temppmeasure)
  }
  
  writepmeasuresfun(pmeasures = pmeasures,
                    pred_error = pred_err,
                    logfile = logfile,
                    all_or_loocv = all_or_loocv,
                    is_loocv = is_loocv,
                    th1 = th1,
                    th2 = th2)
  
  # Close the connection to the text file
  close(logfile)
  
  # Rank score for the OLR model
  rank_score <- rankscorefun(pmeasures = pmeasures,
                             pred_error = pred_err,
                             is_normalized = TRUE)
  
  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename
  
  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_th1 = quantile(pmeasures$Accuracy_th1, 0.025),
    Accuracy_th2 = quantile(pmeasures$Accuracy_th2, 0.025),
    AUC1 = quantile(pmeasures$AUC1, 0.025),
    AUC2 = quantile(pmeasures$AUC2, 0.025),
    Sensitivity1 = quantile(pmeasures$Sensitivity1, 0.025),
    Sensitivity2 = quantile(pmeasures$Sensitivity2, 0.025),
    Specificity1 = quantile(pmeasures$Specificity1, 0.025),
    Specificity2 = quantile(pmeasures$Specificity2, 0.025),
    LR_pos_th1 = quantile(pmeasures$LR_pos_th1, 0.025),
    LR_pos_th2 = quantile(pmeasures$LR_pos_th2, 0.025),
    LR_neg_th1 = quantile(pmeasures$LR_neg_th1, 0.975),
    LR_neg_th2 = quantile(pmeasures$LR_neg_th2, 0.975),
    F1score_th1 = quantile(pmeasures$F1score_th1, 0.025),
    F1score_th2 = quantile(pmeasures$F1score_th2, 0.025),
    Classification_error = mean(pred_err) + 1.96 * sd(pred_err) / sqrt(length(pred_err)),
    pairwise = quantile(pmeasures$Pairwise, 0.025),
    Rank_score = rank_score
  )
  # Add normalized log likelihood value
  summarydf['logLik'] <- logLik(mod)/nrow(all_data)
  
  # Add Alphas
  for (i in 1:length(alphas)) {
    summarydf[[paste0("Alpha_", i)]] <- alphas[i]
  }
  
  # Add Betas
  for (i in 1:length(betas)) {
    summarydf[[paste0("Beta_", i)]] <- betas[i]
  }
  
  
  return(summarydf)
}

run_all_loocv <- function(results_folder = 'results',
                          training = 'trainingdf',
                          testing = 'testingdf',
                          features_vector = 'features',
                          units_vector = 'units',
                          is_normalized = FALSE,
                          max_attempts = 1000) {
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features_vector) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename,"loocv","log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  training <- training[, training_cols]
  testing <- testing[, testing_cols]
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  training$label <- as.integer(factor(training$label, levels = levels, labels = values))
  testing$label <- as.integer(factor(testing$label, levels = levels, labels = values))
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  # Remove errors
  training <- na.omit(training)
  testing <- na.omit(testing)
  
  # Combine all data
  all_data <- rbind(training,testing)
  
  # Normalize the data
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first 'dimension' columns of all_data
    means <- sapply(all_data[1:dimension], mean)
    sds <- sapply(all_data[1:dimension], sd)
    tempdf <- all_data
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, all_data[1:dimension], means, sds)
    all_data <- tempdf
  }
  
  # ================================
  # LOOCV performance on all data
  # ================================
  # All drugs
  drugs <- unique(all_data$drug_name)
  
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)
  
  # Specify the number of attempts
  # max_attempts <- 10000
  
  # Preallocated loocv data
  loocv_data <- data.frame()
  for (drug in drugs) {
    converged <- FALSE
    fit_results <- fitandgetTMS(training_data = subset(all_data,all_data$drug_name != drug),
                                testing_data = subset(all_data,all_data$drug_name == drug),
                                formula = formula,
                                dimension = dimension,
                                max_attempts = max_attempts,
                                logfile = logfile)
    converged <- fit_results$converged
    if (converged){
      testing_data <- fit_results$testing_data
      mod <- fit_results$mod
      probs <- predict(mod, newdata = testing_data, type = "probs")
      testing_data["1"] <- probs[1]
      testing_data["2"] <- probs[2]
      testing_data["3"] <- probs[3]
      testing_data["pred"] <- predict(mod, newdata = testing_data, type = "class")
      loocv_data <- rbind(loocv_data,testing_data)
    } else {
      write(sprintf("Cross validating with %s failed! Skipping...",drug),logfile)
    }
  }
  
  # Check whether the loocv_data contains three TdP risk samples
  if (length(unique(loocv_data$label))!=3 | length(unique(loocv_data$pred))!=3) {
    return(returnNA(log_filename))
  }
  
  # Preallocate the pmeasures dataframe
  pmeasures <- data.frame()
  
  # Calculate the classification error from the whole loocv data
  pred_err <- (abs(as.integer(loocv_data$pred) - as.integer(loocv_data$label)))
  
  # Calculate the pmeasures for the whole dataset
  training_drugs <- unique(training$drug_name)
  testing_drugs <- unique(testing$drug_name)
  loocv_data$is_training <- ifelse(loocv_data$drug_name %in% training_drugs, 1, 0)
  
  for (sample_id in training$Sample_ID[training$drug_name == training_drugs[1]]) {
    temp_data <- subset(loocv_data, loocv_data$Sample_ID == sample_id)
    temppmeasure <- pmeasuresfun_loocv(data = temp_data,
                                       label_values = values,
                                       all_or_loocv = TRUE,
                                       pairwise_test = TRUE)
    pmeasures <- rbind(pmeasures,temppmeasure)
  }
  all_or_loocv <- TRUE
  writepmeasuresfun(pmeasures = pmeasures,
                    pred_error = pred_err,
                    logfile = logfile,
                    all_or_loocv = all_or_loocv,
                    is_loocv = is_loocv,
                    th1 = NA,
                    th2 = NA)
  
  # Close the connection to the text file
  close(logfile)
  
  # Rank score for the OLR model
  rank_score <- rankscorefun(pmeasures = pmeasures,
                             pred_error = pred_err,
                             is_normalized = TRUE)
  
  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename
  
  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_th1 = quantile(pmeasures$Accuracy_th1, 0.025),
    Accuracy_th2 = quantile(pmeasures$Accuracy_th2, 0.025),
    AUC1 = quantile(pmeasures$AUC1, 0.025),
    AUC2 = quantile(pmeasures$AUC2, 0.025),
    Sensitivity1 = quantile(pmeasures$Sensitivity1, 0.025),
    Sensitivity2 = quantile(pmeasures$Sensitivity2, 0.025),
    Specificity1 = quantile(pmeasures$Specificity1, 0.025),
    Specificity2 = quantile(pmeasures$Specificity2, 0.025),
    LR_pos_th1 = quantile(pmeasures$LR_pos_th1, 0.025),
    LR_pos_th2 = quantile(pmeasures$LR_pos_th2, 0.025),
    LR_neg_th1 = quantile(pmeasures$LR_neg_th1, 0.975),
    LR_neg_th2 = quantile(pmeasures$LR_neg_th2, 0.975),
    F1score_th1 = quantile(pmeasures$F1score_th1, 0.025),
    F1score_th2 = quantile(pmeasures$F1score_th2, 0.025),
    Classification_error = mean(pred_err) + 1.96 * sd(pred_err) / sqrt(length(pred_err)),
    pairwise = quantile(pmeasures$Pairwise, 0.025),
    Rank_score = rank_score
  )
  
  return(summarydf)
}

TMS <- function(alphas, betas, features_row){
  z1 <- alphas[1] - sum(betas * features_row)
  z2 <- alphas[2] - sum(betas * features_row)
  
  return((z1 + z2) / 2.0)
}

pairwisefun<- function(fulltable){
  cmb <- combn(seq_len(nrow(fulltable)), 2)
  mergedtable<-cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
  validpairidx <- (mergedtable[,4]!=mergedtable[,10])&(!mergedtable[,6]|!mergedtable[,12])
  correctidx1 <- ((mergedtable[,4]>mergedtable[,10])&(mergedtable[,5]>mergedtable[,11]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,5]<mergedtable[,11])) #when predicted class are different
  correctidx2 <- (mergedtable[,5]==3)&(mergedtable[,11]==3)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both high
  correctidx3 <- (mergedtable[,5]==1)&(mergedtable[,11]==1)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both low
  correctidx4 <- (mergedtable[,5]==2)&(mergedtable[,11]==2)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both intermediate
  correctidx <- correctidx1|correctidx2|correctidx3|correctidx4
  sum(validpairidx&correctidx)/sum(validpairidx)
}

aucrocfun_loocv <- function(data, label_values){
  auc_scores <- c()
  for (class_label in label_values) {
    if (class_label == 1) {
      actual <- as.integer(data$label != class_label)
    } else {
      actual <- as.integer(data$label == class_label)
    }
    predicted_prob <- data[,c("1","2","3")]
    if (class_label == 1) {
      predicted_prob <- predicted_prob[, 2] + predicted_prob[, 3]
    } else {
      predicted_prob <- predicted_prob[, class_label]
    }
    roc_obj <- roc(actual,
                   unlist(predicted_prob),
                   direction = "<",
                   quiet = TRUE)
    auc_score <- auc(roc_obj)
    auc_scores <- c(auc_scores, auc_score)
  }
  return(auc_scores)
}

pmeasuresfun_loocv <- function(data, label_values, tms_name = "TMS", all_or_loocv = TRUE, pairwise_test = TRUE){
  # Calculate performance measures
  if (all_or_loocv) {
    testing_data <- data
  } else {
    testing_data <- subset(data, data$is_training == 0)
  }
  # Create confusion matrix
  confusion_matrix <- table(testing_data$pred, testing_data$label)
  
  # Calculate the accuracy for each class
  tp_cipa_1 = confusion_matrix[2,2] + confusion_matrix[2,3] + confusion_matrix[3,2] + confusion_matrix[3,3]
  tn_cipa_1 = confusion_matrix[1,1]
  fp_cipa_1 = confusion_matrix[2,1] + confusion_matrix[3,1]
  fn_cipa_1 = confusion_matrix[1,2] + confusion_matrix[1,3]
  
  tp_cipa_2 = confusion_matrix[3,3]
  tn_cipa_2 = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[2,1] + confusion_matrix[2,2]
  fp_cipa_2 = confusion_matrix[3,1] + confusion_matrix[3,2]
  fn_cipa_2 = confusion_matrix[1,3] + confusion_matrix[2,3]
  
  F1score_th1 <- 2.0 * tp_cipa_1 / (2.0 * tp_cipa_1 + fp_cipa_1 + fn_cipa_1)
  F1score_th2 <- 2.0 * tp_cipa_2 / (2.0 * tp_cipa_2 + fp_cipa_2 + fn_cipa_2)
  
  Accuracy_th1 <- (tp_cipa_1 + tn_cipa_1) / sum(confusion_matrix)
  Accuracy_th2 <- (tp_cipa_2 + tn_cipa_2) / sum(confusion_matrix)
  
  Sensitivity1 <- tp_cipa_1 / (tp_cipa_1 + fn_cipa_1)
  Sensitivity2 <- tp_cipa_2 / (tp_cipa_2 + fn_cipa_2)
  
  Specificity1 <- tn_cipa_1 / (tn_cipa_1 + fp_cipa_1)
  Specificity2 <- tn_cipa_2 / (tn_cipa_2 + fp_cipa_2)
  
  # Add random number to LR+ and LR- to prevent zero devision
  u <- 1e-6
  sd <- 1e-12
  
  lr_positive_cipa_1 <- (Sensitivity1 + rnorm(1, mean = u, sd = sd)) / (1 - Specificity1 + rnorm(1, mean = u, sd = sd))
  lr_positive_cipa_2 <- (Sensitivity2 + rnorm(1, mean = u, sd = sd)) / (1 - Specificity2 + rnorm(1, mean = u, sd = sd))
  
  lr_negative_cipa_1 <- (1 - Sensitivity1 + rnorm(1, mean = u, sd = sd)) / (Specificity1 + rnorm(1, mean = u, sd = sd))
  lr_negative_cipa_2 <- (1 - Sensitivity2 + rnorm(1, mean = u, sd = sd)) / (Specificity2 + rnorm(1, mean = u, sd = sd))
  
  auc_scores <- aucrocfun_loocv(testing_data, label_values)
  
  # Calculate the pairwise classification error
  if (pairwise_test) {
    data["risk_label"] <- as.integer(data$label)
    data["risk_pred"] <- as.integer(data$pred)
    data <- data[,c("Sample_ID","drug_name",tms_name,"risk_label","risk_pred","is_training")]
    pairwise <- pairwisefun(data)
  } else {
    pairwise <- NA
  }
  
  # Fill the pmeasures row by row
  pmeasures <- data.frame(
    TP_cipa_1 = tp_cipa_1,
    TP_cipa_2 = tp_cipa_2,
    TN_cipa_1 = tn_cipa_1,
    TN_cipa_2 = tn_cipa_2,
    FP_cipa_1 = fp_cipa_1,
    FP_cipa_2 = fp_cipa_2,
    FN_cipa_1 = fn_cipa_1,
    FN_cipa_2 = fn_cipa_2,
    F1score_th1 = F1score_th1,
    F1score_th2 = F1score_th2,
    Accuracy_th1 = Accuracy_th1,
    Accuracy_th2 = Accuracy_th2,
    AUC1 = auc_scores[1],
    AUC2 = auc_scores[3],
    Sensitivity1 = Sensitivity1,
    Sensitivity2 = Sensitivity2,
    Specificity1 = Specificity1,
    Specificity2 = Specificity2,
    LR_pos_th1 = lr_positive_cipa_1,
    LR_pos_th2 = lr_positive_cipa_2,
    LR_neg_th1 = lr_negative_cipa_1,
    LR_neg_th2 = lr_negative_cipa_2,
    Pairwise = pairwise
  )
  
  return(pmeasures)
}

tmsplotfun <- function(data, th1, th2, label_colors, title, file_name, tms_name, tms_unit){
  # Known TdP risk classification (CiPA reference list)
  high_risk_drugs <- c("sotalol","quinidine","dofetilide","bepridil",
                       "vandetanib","ibutilide","disopyramide","azimilide")
  low_risk_drugs  <- c("verapamil","ranolazine","mexiletine","diltiazem",
                       "tamoxifen","nitrendipine","nifedipine","metoprolol","loratadine")

  # Order drugs by label (low -> intermediate -> high) for y-axis
  data$drug_name <- factor(data$drug_name,
                           levels = unique(data$drug_name[order(data$label)]))

  # Build y-axis label color vector matching factor level order
  drug_levels    <- levels(data$drug_name)
  yaxis_colors   <- ifelse(drug_levels %in% high_risk_drugs, "red",
                    ifelse(drug_levels %in% low_risk_drugs,  "green4", "orange"))

  tms  <- tms_name
  plot <- ggplot(data, aes_string(x = tms_name, y = "drug_name", fill = "risk")) +
    geom_boxplot(color = "black", width = 0.5, size = 0.2,
                 outlier.size = 0.5, outlier.shape = NA) +
    labs(title = title, x = tms, y = "") +
    geom_vline(xintercept = th1, linetype = "dashed", color = "blue", size = 1) +
    geom_vline(xintercept = th2, linetype = "dashed", color = "red",  size = 1) +
    scale_fill_manual(values = label_colors) +
    theme(
      plot.title       = element_text(size = 16),
      axis.title.x     = element_text(size = 14),
      axis.title.y     = element_text(size = 14),
      axis.text.x      = element_text(size = 12),
      axis.text.y      = element_text(size = 12, colour = yaxis_colors, face = "bold"),
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 10)
    )
  ggsave(file_name, plot, width = 6, height = 4, dpi = 900)
}

scatterplotfun <- function(data, 
                           mod = NA, 
                           feature_1, 
                           feature_2, 
                           unit_1, 
                           unit_2, 
                           is_converged, 
                           is_training, 
                           is_legend, 
                           results_folder,
                           is_normalized, 
                           all_or_loocv = TRUE){
  # Check the column index
  data <- data.frame(data)
  idx_model_1 <- as.integer(which(colnames(data) == feature_1))
  idx_model_2 <- as.integer(which(colnames(data) == feature_2))
  idx_label <- as.integer(which(colnames(data) == "label"))
  
  if (is_converged) {
    # Variables from Ordinal Logistic Regression model
    alpha1 <- as.numeric(mod$zeta[1])
    alpha2 <- as.numeric(mod$zeta[2])
    beta1 <- as.numeric(mod$coefficients[1])
    beta2 <- as.numeric(mod$coefficients[2])
    
    # Some descriptions of decision boundaries
    m <- - beta1 / beta2
    c1 <- - 1.0 / beta2 * (- alpha1 + log(1.0 - 2.0 * exp(alpha1 - alpha2)))
    c2 <- 1.0 / beta2 * (alpha1 + log(exp(alpha2 - alpha1) - 2.0))
  }
  
  # Create a meshgrid for contour plotting testing dataset
  x <- seq(min(data[,idx_model_1]), max(data[,idx_model_1]), length.out = 100)
  y <- seq(min(data[,idx_model_2]), max(data[,idx_model_2]), length.out = 100)
  if (is_converged) {
    z1 <- outer(x, y, Vectorize(function(x, y) calculate_B1(alpha1, alpha2, beta1, beta2, x, y)))
    z2 <- outer(x, y, Vectorize(function(x, y) calculate_B2(alpha1, alpha2, beta1, beta2, x, y)))
  }
  if (is_training) {
    title <- "(A) Training dataset"
    file_name <- "training"
  } else {
    title <- "(B) Testing dataset"
    file_name <- "testing"
  }
  if (is_normalized) {
    unit_1 <- ""
    unit_2 <- ""
  }
  if (!all_or_loocv) {
    file_name <- paste(file_name,"development",sep = "_")
  }
  
  jpeg(paste(results_folder,"/",paste(feature_1,feature_2,file_name,"dataset.jpg",sep = "_"), sep = ""),quality = 100, units = "in", width = 5, height = 5, res = 900)
  plot(data[,idx_model_1], 
       data[,idx_model_2], 
       xlab = paste(feature_1, unit_1, sep = " "), 
       ylab = paste(feature_2, unit_2, sep = " "),
       main = title,
       cex.axis = 1.5, 
       cex.lab = 1.5, 
       cex.main = 1.5, 
       cex = 0.5)
  if (is_legend) {
    legend("bottomright", legend = c("Low", "Intermediate", "High"), fill = c("green4", "orange", "red"))
  }
  points(data[,idx_model_1][data$label == "1"], data[,idx_model_2][data$label == "1"], col = "green4", cex = 0.5)
  points(data[,idx_model_1][data$label == "2"], data[,idx_model_2][data$label == "2"], col = "orange", cex = 0.5)
  points(data[,idx_model_1][data$label == "3"], data[,idx_model_2][data$label == "3"], col = "red", cex = 0.5)
  if (is_converged) {
    abline(a = c1, b = m, col = "blue", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
    abline(a = c2, b = m, col = "red", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
  }
  dev.off()
}

pairsdfinitfun <- function(features, units, dimension) {
  if (dimension > length(features)) {
    stop("Dimension cannot be greater than the number of features.")
  }
  
  # Calculate all possible combinations of features and units
  feature_combinations <- combn(features, dimension, simplify = FALSE)
  unit_combinations <- combn(units, dimension, simplify = FALSE)
  
  # Initialize an empty list to store data frames for each combination
  pairsdf_list <- list()
  
  # Loop through each combination and create a data frame
  for (i in seq_along(feature_combinations)) {
    feature_combination <- feature_combinations[[i]]
    unit_combination <- unit_combinations[[i]]
    
    # Create a data frame with the feature names and units for this combination
    feature_names <- paste0("feature_", seq_len(dimension))
    unit_names <- paste0("unit_", seq_len(dimension))
    tempdf <- setNames(data.frame(t(feature_combination)), feature_names)
    unitdf <- setNames(data.frame(t(unit_combination)), unit_names)
    
    # Combine features and units into one data frame
    combinedf <- cbind(tempdf, unitdf)
    
    # Add the combined data frame to the list
    pairsdf_list[[i]] <- combinedf
  }
  
  # Combine all data frames in the list into a single data frame
  pairsdf <- do.call("rbind", pairsdf_list)
  
  return(pairsdf)
}

solvepair <- function(results_folder,
                      pair_id,
                      filepath_training,
                      filepath_testing,
                      features_vector,
                      units_vector,
                      is_normalized,
                      max_attempts = 1000,
                      is_loocv){
  # Print the current pair_id and features being processed
  print(paste(c(pair_id, features_vector), collapse = "_"))
  
  # Read in the training and testing datasets
  training <- read_csv(filepath_training, show_col_types = FALSE)
  testing <- read_csv(filepath_testing, show_col_types = FALSE)
  
  # Perform all
  if (is_loocv) {
    resultdf <- run_all_loocv(results_folder = results_folder,
                              training = training,
                              testing = testing,
                              features_vector = features_vector,
                              units_vector = units_vector,
                              is_normalized = is_normalized,
                              max_attempts = max_attempts)
  } else {
    resultdf <- run_all_training(results_folder = results_folder,
                                 training = training,
                                 testing = testing,
                                 features_vector = features_vector,
                                 units_vector = units_vector,
                                 is_normalized = is_normalized,
                                 max_attempts = max_attempts)
  }
  
  return(resultdf)
}

writepmeasuresfun <- function(pmeasures, 
                              pred_error,
                              logfile,
                              all_or_loocv = TRUE,
                              is_loocv = TRUE,
                              th1 = NA,
                              th2 = NA){
  
  # Print the Thresholds into logfile
  write("=======================",logfile)
  write(sprintf("TMS thresholds"),logfile)
  write("=======================",logfile)
  write(paste(sprintf('Threshold_1: %.4f ', th1)), logfile)
  write(paste(sprintf('Threshold_2: %.4f ', th2)), logfile)
  
  # Print the pmeasures dataframe into logfile
  if (is_loocv) {
    write("==============================",logfile)
    write("Performance measures for loocv data",logfile)
    write("==============================",logfile)
  } else if (all_or_loocv ) {
    write("==============================",logfile)
    write("Performance measures for whole training data",logfile)
    write("==============================",logfile)
  } else {
    write("==============================",logfile)
    write("Performance measures for testing data based on num_tests",logfile)
    write("==============================",logfile)
  }
  write(sprintf('AUC1: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC1),
                quantile(pmeasures$AUC1, 0.025),
                quantile(pmeasures$AUC1, 0.975)),
        logfile)
  write(sprintf('AUC2: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC2),
                quantile(pmeasures$AUC2, 0.025),
                quantile(pmeasures$AUC2, 0.975)),
        logfile)
  write(sprintf('Accuracy_th1: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_th1),
                quantile(pmeasures$Accuracy_th1, 0.025),
                quantile(pmeasures$Accuracy_th1, 0.975)),
        logfile)
  write(sprintf('Accuracy_th2: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_th2),
                quantile(pmeasures$Accuracy_th2, 0.025),
                quantile(pmeasures$Accuracy_th2, 0.975)),
        logfile)
  write(sprintf('Sensitivity1: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity1),
                quantile(pmeasures$Sensitivity1, 0.025),
                quantile(pmeasures$Sensitivity1, 0.975)),
        logfile)
  write(sprintf('Sensitivity2: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity2),
                quantile(pmeasures$Sensitivity2, 0.025),
                quantile(pmeasures$Sensitivity2, 0.975)),
        logfile)
  write(sprintf('Specificity1: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity1),
                quantile(pmeasures$Specificity1, 0.025),
                quantile(pmeasures$Specificity1, 0.975)),
        logfile)
  write(sprintf('Specificity2: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity2),
                quantile(pmeasures$Specificity2, 0.025),
                quantile(pmeasures$Specificity2, 0.975)),
        logfile)
  write(sprintf('LR_pos_th1: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_pos_th1),
                quantile(pmeasures$LR_pos_th1, 0.025),
                quantile(pmeasures$LR_pos_th1, 0.975)),
        logfile)
  write(sprintf('LR_pos_th2: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_pos_th2),
                quantile(pmeasures$LR_pos_th2, 0.025),
                quantile(pmeasures$LR_pos_th2, 0.975)),
        logfile)
  write(sprintf('LR_neg_th1: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_neg_th1),
                quantile(pmeasures$LR_neg_th1, 0.025),
                quantile(pmeasures$LR_neg_th1, 0.975)),
        logfile)
  write(sprintf('LR_neg_th2: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_neg_th2),
                quantile(pmeasures$LR_neg_th2, 0.025),
                quantile(pmeasures$LR_neg_th2, 0.975)),
        logfile)
  write(sprintf('F1score_th1: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_th1),
                quantile(pmeasures$F1score_th1, 0.025),
                quantile(pmeasures$F1score_th1, 0.975)),
        logfile)
  write(sprintf('F1score_th2: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_th2),
                quantile(pmeasures$F1score_th2, 0.025),
                quantile(pmeasures$F1score_th2, 0.975)),
        logfile)
  write(sprintf('Classification_error: %.4f (%.4f, %.4f)',
                mean(pred_error),
                mean(pred_error) - 1.96 * sd(pred_error) / sqrt(length(pred_error)),
                mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))),
        logfile)
  if (all_or_loocv) {
    write(sprintf('Pairwise: NA'),
          logfile)
  } else {
    write(sprintf('Pairwise: %.4f (%.4f, %.4f)',
                  median(pmeasures$Pairwise),
                  quantile(pmeasures$Pairwise, 0.025),
                  quantile(pmeasures$Pairwise, 0.975)),
          logfile)
  }
  
}

rankscorefun <- function(pmeasures,
                         pred_error,
                         is_normalized = TRUE){
  Performance_measures <- c("AUC1", "AUC2",
                            "LR_pos_th1", "LR_pos_th2",
                            "LR_neg_th1", "LR_neg_th2",
                            "Pairwise","Classification_error")
  Performance_levels <- c("Excellent_performance", 
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  # A dataframe for the weights of performance measures
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0)
  )
  
  # A dataframe for the weights of performance levels
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA # Not acceptable performance is removed
  
  # Initialize the performance dataframe for the model
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA, NA)
  )
  
  # Check the performance level for each performance measure 
  # by looking at the 95% confidence interval that match the CiPA's criteria
  
  # AUC1
  score_to_check <- quantile(pmeasures$AUC1, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC2
  score_to_check <- quantile(pmeasures$AUC2, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th1
  score_to_check <- quantile(pmeasures$LR_pos_th1, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th2
  score_to_check <- quantile(pmeasures$LR_pos_th2, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th1
  score_to_check <- quantile(pmeasures$LR_neg_th1, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th2
  score_to_check <- quantile(pmeasures$LR_neg_th2, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Pairwise
  score_to_check <- quantile(pmeasures$Pairwise, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error
  score_to_check <- mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Calculate rank score
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  
  return(as.numeric(rank_score))
}

rankscorefun_median <- function(pmeasures,
                                pred_error,
                                is_normalized = TRUE){
  # Identical to rankscorefun but uses median (50th percentile) instead of
  # worst-case CI (2.5th / 97.5th) for every performance measure.
  Performance_measures <- c("AUC1", "AUC2",
                            "LR_pos_th1", "LR_pos_th2",
                            "LR_neg_th1", "LR_neg_th2",
                            "Pairwise","Classification_error")
  Performance_levels <- c("Excellent_performance",
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0)
  )
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA
  
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA, NA)
  )
  
  # AUC1
  score_to_check <- quantile(pmeasures$AUC1, 0.5)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC2
  score_to_check <- quantile(pmeasures$AUC2, 0.5)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th1
  score_to_check <- quantile(pmeasures$LR_pos_th1, 0.5)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th2
  score_to_check <- quantile(pmeasures$LR_pos_th2, 0.5)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th1 (lower is better → median, same direction check)
  score_to_check <- quantile(pmeasures$LR_neg_th1, 0.5)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th2
  score_to_check <- quantile(pmeasures$LR_neg_th2, 0.5)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Pairwise
  score_to_check <- quantile(pmeasures$Pairwise, 0.5)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error — use median(pred_error) as the point estimate
  score_to_check <- median(pred_error)
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  return(as.numeric(rank_score))
}

rankscorefun_training <- function(pmeasures,
                                  pred_error,
                                  is_normalized = TRUE){
  Performance_measures <- c("AUC1", "AUC2",
                            "LR_pos_th1", "LR_pos_th2",
                            "LR_neg_th1", "LR_neg_th2",
                            "Classification_error")
  Performance_levels <- c("Excellent_performance", 
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  # A dataframe for the weights of performance measures
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0)
  )
  
  # A dataframe for the weights of performance levels
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA # Not acceptable performance is removed
  
  # Initialize the performance dataframe for the model
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA)
  )
  
  # Check the performance level for each performance measure 
  # by looking at the 95% confidence interval that match the CiPA's criteria
  
  # AUC1
  score_to_check <- quantile(pmeasures$AUC1, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC2
  score_to_check <- quantile(pmeasures$AUC2, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th1
  score_to_check <- quantile(pmeasures$LR_pos_th1, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th2
  score_to_check <- quantile(pmeasures$LR_pos_th2, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th1
  score_to_check <- quantile(pmeasures$LR_neg_th1, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th2
  score_to_check <- quantile(pmeasures$LR_neg_th2, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error
  score_to_check <- mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Calculate rank score
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  
  return(as.numeric(rank_score))
}

rankscorefun_training_median <- function(pmeasures,
                                         pred_error,
                                         is_normalized = TRUE){
  # Identical to rankscorefun_training but uses median (50th percentile)
  # instead of worst-case CI for every performance measure.
  Performance_measures <- c("AUC1", "AUC2",
                            "LR_pos_th1", "LR_pos_th2",
                            "LR_neg_th1", "LR_neg_th2",
                            "Classification_error")
  Performance_levels <- c("Excellent_performance",
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0)
  )
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA
  
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA)
  )
  
  # AUC1
  score_to_check <- quantile(pmeasures$AUC1, 0.5)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC2
  score_to_check <- quantile(pmeasures$AUC2, 0.5)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th1
  score_to_check <- quantile(pmeasures$LR_pos_th1, 0.5)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th2
  score_to_check <- quantile(pmeasures$LR_pos_th2, 0.5)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th1
  score_to_check <- quantile(pmeasures$LR_neg_th1, 0.5)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th2
  score_to_check <- quantile(pmeasures$LR_neg_th2, 0.5)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error — use median(pred_error) as the point estimate
  score_to_check <- median(pred_error)
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  return(as.numeric(rank_score))
}


fitandgetTMS <- function(training_data,testing_data,formula,dimension,max_attempts,logfile){
  converged <- FALSE
  
  # First trial without start option
  tryCatch({
    # Attempt to fit the model without specifying start
    mod <- polr(formula, data = training_data, Hess = TRUE)
    converged <- TRUE
  }, error = function(e) {
    write(sprintf("First trial without 'start' failed. Retrying..."),logfile)
  })
  
  # If the first trial fails, attempt to fit with random starting values
  attempts <- as.integer(0)
  set.seed(1)
  if (!converged) {
    for (attemp_idx in 1:max_attempts) {
      # Generate random starting values
      start_values <- rnorm(dimension + 2, mean = 0.0, sd = 100.0)
      
      # Attempt to fit the model with random starting values
      tryCatch({
        mod <- polr(formula, data = training_data, start = start_values, Hess = TRUE)
        converged <- TRUE
        break
      }, error = function(e) {
        attempts <- attempts + 1
      })
    }
  }
  
  # If all attempts fail, print a message
  if (!converged) {
    write(sprintf("Could not fit the model after %d attempts",max_attempts),logfile)
    return(
      list(
        training_data = training_data,
        testing_data = testing_data,
        mod = NA,
        converged = converged
      )
    )
    
  } else {
    # Make predictions for training dataset
    predicted_labels <- predict(mod, newdata = training_data, type = "class")
    
    # Variables from Ordinal Logistic Regression model
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    
    # if (!is_single) {
    # Apply TMS function row-wise efficiently
    training_data$TMS <- apply(training_data[, 1:dimension], 1, function(row) TMS(alphas, betas, row))
    testing_data$TMS <- apply(testing_data[, 1:dimension], 1, function(row) TMS(alphas, betas, row))
    # }
    return(
      list(
        training_data = training_data,
        testing_data = testing_data,
        mod = mod,
        converged = converged
      )
    )
  }
}

returnNA <- function(log_filename,dimension = NA){
  summarydf <- data.frame(
    Feature_Pair = log_filename,
    Accuracy_th1 = NA,
    Accuracy_th2 = NA,
    AUC1 = NA,
    AUC2 = NA,
    Sensitivity1 = NA,
    Sensitivity2 = NA,
    Specificity1 = NA,
    Specificity2 = NA,
    LR_pos_th1 = NA,
    LR_pos_th2 = NA,
    LR_neg_th1 = NA,
    LR_neg_th2 = NA,
    F1score_th1 = NA,
    F1score_th2 = NA,
    Classification_error = NA,
    pairwise = NA,
    Rank_score = NA,
    th1 = NA,
    th2 = NA
  )
  
  if (!is.na(dimension)) {
    # Add log likelihood value
    summarydf['logLik'] <- NA 
    
    # Add Alphas
    for (i in 1:2) {
      summarydf[[paste0("Alpha_", i)]] <- NA
    }
    
    # Add Betas
    for (i in 1:dimension) {
      summarydf[[paste0("Beta_", i)]] <- NA
    }
  }
  
  return (summarydf)
}

try_polr<-function(datadf, predname, penalty=0, tol=NA, maxit=1e6){
  if(length(penalty)>1 || penalty!=0)
    print("Ignoring penalty parameter for polr function!")
  if(!is.na(tol))
    print("Ignoring tolerance parameter for polr function!")
  
  # get starting guess (polr's default often fails)
  u <- as.integer(table(datadf$class))
  u <- (cumsum(u)/sum(u))[1:(length(u)-1)]
  zetas <- qlogis(u)
  s0 <- c(0, zetas[1], log(diff(zetas)))
  
  try({
    polr(as.formula(paste0("class~",predname)), data=datadf, start=s0, Hess=TRUE, control=list(maxit=maxit)) # note: large maxit is very slow but seems necessary for convergence
  })
}

try_glm<-function(datadf, predname, penalty=0, tol=NA, maxit=1e6){
  if(length(penalty)>1 || penalty!=0)
    print("Ignoring penalty parameter for polr function!")
  if(!is.na(tol))
    print("Ignoring tolerance parameter for polr function!")
  
  try({
    glm(as.formula(paste0("class~",predname)), data=datadf, family=binomial(link="logit"), control=list(maxit=maxit))
  })
}

try_repolr<-function(datadf,predname,penalty=0,tol=1e-10,maxit=1e6){
  
  try({
    repolr(as.formula(paste0("class~",predname)), subjects="drug",data=datadf, categories=3,times=1:maxsamp,corr.mod="independence",alpha=0.9,diffmeth="analytic", fit.opt=c(cmaxit = maxit, omaxit = maxit, ctol = tol, otol = tol, h = 0.01))
  })
  
}#repolr

try_ordLORgee<-function(datadf, predname, penalty=0,tol=1e-10,maxit=1e6){
  options(warn=1)
  try({
    ordLORgee(as.formula(paste0("class~",predname)), id=datadf$drug,data=datadf, link="logit", repeated=datadf$sample,LORstr="independence",control=LORgee.control(verbose=FALSE))
  })
  
}#ordLORgee

try_lrm<-function(datadf, predname, penalty, tol=1e-10, maxit=1e6){
  if(class(penalty)=="list" || length(penalty)==1){ # specified penalty
    fout<-try({
      lrm(as.formula(paste0("class~",predname)), data=datadf, penalty=penalty, x=TRUE, y=TRUE, tol=tol, maxit=maxit)
    })
    return(fout)
  }
  
  pen_vec<-penalty # candidate penalties to try if penalty=0 doesn't work
  fout<-try({
    lrm(as.formula(paste0("class~",predname)), data=datadf, penalty=0, x=TRUE, y=TRUE, tol=tol, maxit=maxit)
  })
  if(inherits(fout, "try-error")){ # try all the penalties
    print("Unpenalized regression didn't work, using first penalty that works!")
    i<-1
    while(inherits(fout, "try-error") && i<=length(pen_vec)){
      pen<-pen_vec[i]
      i<-i+1
      fout<-try({
        lrm(as.formula(paste0("class~",predname)), data=datadf, penalty=pen, x=TRUE, y=TRUE, tol=tol, maxit=maxit)
      })
    }
    #print(sprintf("Using penalty = %g",pen))
  }
  if(inherits(fout, "try-error"))
    return(fout)
  
  p <- pentrace(fout, penalty=pen_vec, tol=tol, maxit=maxit, noaddzero=TRUE)
  if(inherits(p, "try-error"))
    return(p)
  #print(sprintf("Using penalty = %g",p$penalty))
  update(fout, penalty=p$penalty)
}


combine_drug <- function(num_of_drugs, drug_candidates) {
  # Number of elements to choose (excluding L1 and L2)
  # rand_elements <- num_of_drugs-2 # 8 drugs = 2 fixed and 6 random
  rand_elements <- num_of_drugs # All random
  
  # Generate combinations with repetitions for the remaining elements
  elements <- drug_candidates$drug
  combinations <- combinations(n = length(elements), r = rand_elements, v = elements, repeats.allowed = FALSE)
  
  # # Filter combinations to ensure no more than 3 L's in each combination
  # #filtered_combinations <- subset(combinations, rowSums(combinations == "L1" | combinations == "L2") <= 3)
  check_counts <- function(combination) {
    # 8 drugs = L<2; 2<I<3; 2<H<3 
    max_drugs <- num_of_drugs / 3
    # sum(combination %in% c("tamoxifen", "loratadine")) >= floor(max_drugs)-2 &&
    # sum(combination %in% c("tamoxifen", "loratadine")) < ceiling(max_drugs)-1 &&
    sum(combination %in% drug_candidates$drug[drug_candidates$risk==0]) >= floor(max_drugs) &&
      sum(combination %in% drug_candidates$drug[drug_candidates$risk==0]) <= ceiling(max_drugs) &&
      sum(combination %in% drug_candidates$drug[drug_candidates$risk==1]) >= floor(max_drugs) &&
      sum(combination %in% drug_candidates$drug[drug_candidates$risk==1]) <= ceiling(max_drugs) &&
      sum(combination %in% drug_candidates$drug[drug_candidates$risk==2]) >= floor(max_drugs) &&
      sum(combination %in% drug_candidates$drug[drug_candidates$risk==2]) <= ceiling(max_drugs)
  }
  # Filter combinations to ensure no more than 3 H in each combination
  filtered_combinations <- subset(combinations, apply(combinations, 1, check_counts))
  # # Add L1 and L2 to the filtered combinations
  # final_combinations <- cbind("verapamil","ranolazine",filtered_combinations)
  if (num_of_drugs == 28) {
    final_combinations <- combinations
  } else {
    final_combinations <- cbind(filtered_combinations)
  }
  return(final_combinations)
}

open_csv_file <- function(filename) {
  
  # Check if the file exists
  if (!file.exists(filename)) {
    cat(sprintf("File '%s' does not exist. Skipping...\n", filename))
    next
  }
  
  # Try to read the file and handle potential errors
  tryCatch({
    data <- read.csv(filename)
    # cat(sprintf("Successfully opened file: '%s'\n", filename))
    return(data) # Return the data if successful
  }, error = function(e) {
    cat(sprintf("File '%s' is corrupted or unreadable. Skipping...\n", filename))
  })
  
  cat("No valid files were found.\n")
  return(NULL) # Return NULL if no valid file was found
}

#--- function to try ordinal logistic regression
try_lrm_loocv <- function(datadf, tol = 1e-10, maxit = 1e6) {
  # Create and assign datadist object
  dd <- datadist(datadf)
  options(datadist = "dd")
  print("Check try lrm")
  try({ lrm(CiPA ~ get(metric), data = datadf, penalty = 0, x = TRUE, y = TRUE, tol = tol, maxit = maxit) })
}

evaluate_drug_combination <- function(drug_candidates, drug_pair, drug_index, simdir, metrics_file){
  # print(drug_candidates)
  # print(drug_pair)
  # print(drug_index)
  # print(simdir)
  # print(metrics_file)
  # Define measures
  measurevec <- c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error",
                  "AUC1_LOOCV","AUC2_LOOCV","Pairwise_LOOCV","LR1plus_LOOCV","LR1minus_LOOCV","LR2plus_LOOCV","LR2minus_LOOCV","Mean_error_LOOCV",
                  "AUC1_training","AUC2_training","Pairwise_training","LR1plus_training","LR1minus_training","LR2plus_training","LR2minus_training","Mean_error_training",
                  "Threshold1", "Threshold2", "Normalized_logLik")
  columns <- c(to_vec(for(i in 1:drug_size) paste0("drug",i)), measurevec) 
  result_df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(result_df) <- columns
  
  drugtable <- data.frame(drug = c(drug_pair))
  for (row in 1:nrow(drugtable)) {
    drugtable$CiPA[row] = drug_candidates$risk[drug_candidates$drug==drugtable$drug[row]]
  }
  drugnames <- as.character(drugtable$drug)
  
  outdir_vec <- simdir
  
  for (outdir in outdir_vec) {
    # read in dataset
    infile <- metrics_file
    df <- read.csv(infile)
    df <- df[df$drug != "control" & df$dose != 0, ]
    df <- df[!is.na(df$max_dv), ] # depolarization failures
    df <- df[(df$dose) %in% c("1", "2", "3", "4"), ]                  #very important cause some drugs may have additional doses in metrics.rds!!!!
    if (input == "uncertainty") {
      qNettable <- aggregate(df[, metric], list(df$drug, df$sample), mean)  #name is qNettable, but it could be any metric
      colnames(qNettable) <- c("drug", "sample", metric)
    } else {
      qNettable <- aggregate(df[, metric], list(df$drug), mean)
      colnames(qNettable) <- c("drug", "qNet")
      qNettable$sample <- NA
    }
    # print(head(qNettable))
    # print(class(qNettable$qNet))
    if (metric == "qNet")
      qNettable$qNet <- qNettable$qNet / 1000
    df <- merge(qNettable, drugtable[, c("drug", "CiPA")], by.x = "drug", by.y = "drug")
    # }
    
    #organize df
    df$drug <- factor(df$drug, levels = drugnames)
    df <- df[order(df$drug, df$sample), ]
    df$CiPA <- ordered(df$CiPA)
    
    # perform logistic regression and leave-one-out cross validation for each drug
    errdf <- data.frame()
    probdf <- data.frame()
    allprobdf <- data.frame()   #corresponding to cvdf but without LOOCV
    cverrdf <- data.frame()
    cvprobdf <- data.frame()
    cvdf <- data.frame()                              #used to store individual sample's probs
    is_training_all_error <- FALSE
    
    for (drug in c(NA, drugnames)) {
      
      # fit logistic regression model
      if (is.na(drug)) {
        print("training on all drugs")
        datadf <- df[, c("drug", "CiPA", "sample", metric)]
        testdf <- datadf
        traindf <- datadf
      } else {
        print(sprintf("cross validating with %s", drug))
        testdf <- datadf[datadf$drug == drug, ]
        traindf <- datadf[datadf$drug != drug, ]
        if (nrow(testdf) == 0) {
          print("no data to test! skipping...")
          next
        }
      }
      
      # lmod <- try_lrm_loocv(traindf)
      dd <- datadist(traindf)
      options(datadist = "dd")
      lmod <- try({ lrm(CiPA ~ get(metric), data = traindf, penalty = 0, x = TRUE, y = TRUE, tol = 1e-10, maxit = 1e6) })
      # print(lmod)
      
      if (inherits(lmod, "try-error")) {
        if (is.na(drug)) {
          print("fitting all drugs failed? skipping...")
          is_training_all_error <- TRUE
          break
        } else {
          print(sprintf("cross validating with %s failed! skipping...", drug))
          next
        }
      }
      
      # save coefficients
      print(sprintf("Convergence failure: %s", lmod$fail))
      cf <- coefficients(lmod)
      cfvec <- c()
      for (kint in 1:(length(cf) - 1))
        cfvec[[paste0("intercept", kint)]] <- cf[[kint]]
      cfvec[["slope"]] <- cf[[length(cf)]]
      
      #use formal math for threshold 1 and 2
      ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
      t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
      t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]
      
      # get training/prediction error
      y0 <- testdf$CiPA
      if (is.na(drug)) {
        probs <- predict(lmod, type = "fitted.ind")
      } else {
        probs <- matrix(predict(lmod, newdata = testdf, type = "fitted.ind"), ncol = length(levels(y0)))
      }
      
      #sometimes APD90 is NA, giving rise to NA in probs
      idx <- apply(probs, 1, function(x) any(is.na(x)))
      probs[idx, ] <- matrix(rep(c(0, 0, 1), sum(idx)), nrow = sum(idx), byrow = T) #this is to assume all NAs are due to repolarization failure, and thus high risk
      
      yPred <- apply(probs, 1, function(x) which.max(x))
      pred_err <- mean(abs(yPred - as.integer(y0)))
      
      # Loglikelihood value of the model
      logLikVal <- numeric(length(y0))
      
      for (i in 1:length(y0)) {
        logLikVal[i] <- log(probs[i,y0[i]])
      }
      
      # Normalized loglikelihood value
      nrows_trainingdf <- nrow(traindf)
      normalized_logLik <- logLikVal / nrows_trainingdf
      
      # append to data frame
      if (is.na(drug)) {
        newrow <- data.frame(t(cfvec), error = pred_err)
        errdf <- rbind(errdf, newrow)
        
        # allprobdf <- cbind(traindf[, 1:3], probs)
        # colnames(allprobdf) <- c("drug", "CiPA", "sample", "low_prob", "inter_prob", "high_prob")
        allprobdf <- cbind(traindf[, 1:3], probs, abs(yPred - as.integer(y0)), normalized_logLik)
        colnames(allprobdf) <- c("drug", "CiPA", "sample", "low_prob", "inter_prob", "high_prob", "pred_err", "normalized_logLik")
        
      } else {
        newrow <- data.frame(drug = drug, t(cfvec), error = pred_err)
        cverrdf <- rbind(cverrdf, newrow)
        if (outputI) {                            #make detailed prob table here
          testdf <- cbind(testdf, high_prob = probs[, 3], inter_prob = probs[, 2], low_prob = probs[, 1], pred_err = abs(yPred - as.integer(y0)))
          cvdf <- rbind(cvdf, testdf)
        } #if outputI
      }
      
      # detailed results
      newcols <- sapply(levels(y0), function(s) paste0("predict_", s), USE.NAMES = F)
      predictions <- factor(newcols[yPred], levels = newcols)
      tmpdf <- as.data.frame(tapply(1:nrow(testdf), list(drug = as.character(testdf$drug), predict = predictions), FUN = length))
      tmpdf[is.na(tmpdf)] <- 0
      tmpdf$drug <- factor(rownames(tmpdf), levels = drugnames)
      
      tmpdf <- merge(drugtable[, c("drug", "CiPA")], tmpdf)
      tmpdf <- tmpdf[order(tmpdf$drug), ]
      rownames(tmpdf) <- NULL
      if (is.na(drug)) {
        probdf <- rbind(probdf, tmpdf)
        
      } else {
        cvprobdf <- rbind(cvprobdf, tmpdf)
      }
      
    } # for drug
    
    if (!is_training_all_error) {
      outfile <- paste0(outdir, "/",metric, "_", drug_index, "_training_probs.csv")   #here training means "trained on all drugs"
      write.csv(probdf, outfile, row.names = F, quote = F)
      # print(head(probdf))
      
      outfile <- paste0(outdir, "/",metric,  "_", drug_index, "_training_errors.csv")
      list_columns <- sapply(errdf, is.list) # Identify columns in errdf that are of type list
      if (any(list_columns)) { # Convert list columns to character (or other suitable type) if any
        errdf[list_columns] <- lapply(errdf[list_columns], as.numeric)
      }
      write.csv(errdf, outfile, row.names = F, quote = F)
      # print(head(errdf))
      
      outfile <- paste0(outdir, "/",metric,  "_", drug_index, "_LOOCV_probs.csv")
      write.csv(cvprobdf, outfile, row.names = F, quote = F)
      # print(head(cvprobdf))
      
      outfile <- paste0(outdir, "/",metric,  "_", drug_index, "_LOOCV_errors.csv")
      list_columns <- sapply(cverrdf, is.list) # Identify columns in errdf that are of type list
      if (any(list_columns)) { # Convert list columns to character (or other suitable type) if any
        cverrdf[list_columns] <- lapply(cverrdf[list_columns], as.numeric)
      }
      write.csv(cverrdf, outfile, row.names = F, quote = F)
      # print(head(cverrdf))
      
      if (outputI) {
        outfile <- paste0(outdir,"/", metric,  "_", drug_index, "_LOOCV_allprobes.csv")
        write.csv(cvdf, outfile, row.names = F, quote = F)
        # print(head(cvdf))
        outfile <- paste0(outdir,"/", metric,  "_", drug_index, "_training_allprobes.csv")
        write.csv(allprobdf, outfile, row.names = F, quote = F)
      }
    }
  } # for outdir
  
  if (!is_training_all_error) {
    
    # Initialize data frames
    all28df <- data.frame(metric=c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error",
                                   "AUC1_LOOCV","AUC2_LOOCV","Pairwise_LOOCV","LR1plus_LOOCV","LR1minus_LOOCV","LR2plus_LOOCV","LR2minus_LOOCV","Mean_error_LOOCV",
                                   "AUC1_training","AUC2_training","Pairwise_training","LR1plus_training","LR1minus_training","LR2plus_training","LR2minus_training","Mean_error_training",
                                   "Threshold1", "Threshold2", "Normalized_logLik"))
    colordf <- data.frame(
      drug = c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
               "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
               "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine"),
      classidx=c(rep(2,8),rep(1,11),rep(0,9)),
      coloridx=c(rep(1,8),rep(3,11),rep(2,9)),
      isTraining=c(1,1,1,1,0,0,0,0,
                   1,1,1,1,0,0,0,0,0,0,0,
                   1,1,1,1,0,0,0,0,0)
    )
    
    # Loop over datasets and metrics
    for(dataset in datasetvec){
      for(metric in metricv){
        altstr <- paste0(dataset,"_",metric)
        altvec <- 0
        valdf <- read.delim(paste0(simdir,"/",metric,"_", drug_index, "_LOOCV_allprobes.csv"),sep=",",as.is=T)
        trainingdf <- read.delim(paste0(simdir,"/",metric,"_", drug_index, "_training_allprobes.csv"),sep=",",as.is=T)
        file.remove(paste0(simdir,"/",metric,"_", drug_index, "_LOOCV_allprobes.csv"))
        file.remove(paste0(simdir,"/",metric,"_", drug_index, "_training_allprobes.csv"))
        
        # Classify data
        valdf$classidx <- valdf$CiPA
        valdf$class <- 2 - valdf$classidx 
        
        trainingdf$classidx <- trainingdf$CiPA
        trainingdf$class <- 2 - trainingdf$classidx 
        
        # Prepare data for ROC analysis
        forROC1 <- valdf[,c("inter_prob","classidx","sample", "drug")]
        forROC1$inter_prob <- rowSums(valdf[,c("inter_prob","high_prob")])
        colnames(forROC1)[1] <- "positive_prob"
        forROC1$class <- forROC1$classidx > 0  # now high and intermediate are both "TRUE" or "positive" 
        
        forROC1_training <- trainingdf[,c("inter_prob","classidx","sample", "drug")]
        forROC1_training$inter_prob <- rowSums(trainingdf[,c("inter_prob","high_prob")])
        colnames(forROC1_training)[1] <- "positive_prob"
        forROC1_training$class <- forROC1_training$classidx > 0  # now high and intermediate are both "TRUE" or "positive" 
        
        # Generate predictions and labels
        forROC1predictions <- do.call(cbind,by(forROC1[,"positive_prob"],forROC1$sample,function(x) x))
        forROC1labels <- do.call(cbind,by(forROC1$class,forROC1$sample,function(x) x))
        
        forROC1predictions_training <- do.call(cbind,by(forROC1_training[,"positive_prob"],forROC1_training$sample,function(x) x))
        forROC1labels_training <- do.call(cbind,by(forROC1_training$class,forROC1_training$sample,function(x) x))
        
        # Perform ROC analysis
        ROC1obj <- prediction(forROC1predictions, forROC1labels)        
        AUC1 <- performance(ROC1obj,"auc")  #the value is in y.values
        AUC1dist <- formatC(quantile(unlist(AUC1@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(AUC1dist[2]," (",AUC1dist[1]," - ",AUC1dist[3],")"))
        
        ROC1obj_training <- prediction(forROC1predictions_training, forROC1labels_training)        
        AUC1_training <- performance(ROC1obj_training,"auc")  #the value is in y.values
        AUC1dist_training <- formatC(quantile(unlist(AUC1_training@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(AUC1dist_training[2]," (",AUC1dist_training[1]," - ",AUC1dist_training[3],")"))
        
        # Prepare data for second ROC analysis
        forROC2 <- valdf[,c("high_prob","classidx","sample","drug")]
        colnames(forROC2)[1] <- "positive_prob"
        forROC2$class <- forROC2$classidx > 1                         
        
        forROC2_training <- trainingdf[,c("high_prob","classidx","sample","drug")]
        colnames(forROC2_training)[1] <- "positive_prob"
        forROC2_training$class <- forROC2_training$classidx > 1
        
        # Generate predictions and labels
        forROC2predictions <- do.call(cbind,by(forROC2[,"positive_prob"],forROC2$sample,function(x) x))  
        forROC2labels <- do.call(cbind,by(forROC2$class,forROC2$sample,function(x) x))
        
        forROC2predictions_training <- do.call(cbind,by(forROC2_training[,"positive_prob"],forROC2_training$sample,function(x) x))  
        forROC2labels_training <- do.call(cbind,by(forROC2_training$class,forROC2_training$sample,function(x) x))
        
        # Perform second ROC analysis
        ROC2obj <- prediction(forROC2predictions, forROC2labels)        
        AUC2 <- performance(ROC2obj,"auc")  #the value is in y.values 
        AUC2dist <- formatC(quantile(unlist(AUC2@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(AUC2dist[2]," (",AUC2dist[1]," - ",AUC2dist[3],")"))
        
        ROC2obj_training <- prediction(forROC2predictions_training, forROC2labels_training)        
        AUC2_training <- performance(ROC2obj_training,"auc")  #the value is in y.values 
        AUC2dist_training <- formatC(quantile(unlist(AUC2_training@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(AUC2dist_training[2]," (",AUC2dist_training[1]," - ",AUC2dist_training[3],")"))
        
        # Pairwise ranking
        pairwisefun <- function(fulltable){
          cmb <- combn(seq_len(nrow(fulltable)), 2)
          mergedtable <- cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
          validpairidx <- (mergedtable[,7]!=mergedtable[,16])&(!mergedtable[,9]|!mergedtable[,18])
          correctidx1 <- ((mergedtable[,7]>mergedtable[,16])&(mergedtable[,6]<mergedtable[,15]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,6]>mergedtable[,15])) #when predicted class are different
          correctidx2 <- (mergedtable[,6]==1)&(mergedtable[,15]==1)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,3]>mergedtable[,12]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,3]<mergedtable[,12]))) #when predicted class are both high
          correctidx3 <- (mergedtable[,6]==3)&(mergedtable[,15]==3)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both low
          correctidx4 <- (mergedtable[,6]==2)&(mergedtable[,15]==2)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both intermediate
          correctidx <- correctidx1|correctidx2|correctidx3|correctidx4
          sum(validpairidx&correctidx)/sum(validpairidx)
        }
        tempdf <- valdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
        tempdf$pred <- apply(valdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
        tempdf <- merge(tempdf,colordf,by="drug")
        distpairwise <- by(tempdf, valdf$sample, pairwisefun)
        pairwisedist <- formatC(quantile(unlist(distpairwise),prob=c(0.025,0.5,0.975),na.rm=TRUE),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(pairwisedist[2]," (",pairwisedist[1]," - ",pairwisedist[3],")"))
        
        tempdf_training <- trainingdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
        tempdf_training$pred <- apply(trainingdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
        tempdf_training <- merge(tempdf_training,colordf,by="drug")
        distpairwise_training <- by(tempdf_training, trainingdf$sample, pairwisefun)
        pairwisedist_training <- formatC(quantile(unlist(distpairwise_training),prob=c(0.025,0.5,0.975),na.rm=TRUE),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(pairwisedist_training[2]," (",pairwisedist_training[1]," - ",pairwisedist_training[3],")"))
        
        # Calculate sensitivity and specificity
        sens1 <- with(valdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
        spec1 <-  with(valdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
        LRplus1 <- sens1/(1-spec1)
        LRminus1 <- (1-sens1)/spec1
        
        sens2 <- with(valdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
        spec2 <-  with(valdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
        LRplus2 <- sens2/(1-spec2)
        LRminus2 <- (1-sens2)/spec2
        
        sens1_training <- with(trainingdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
        spec1_training <-  with(trainingdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
        LRplus1_training <- sens1_training/(1-spec1_training)
        LRminus1_training <- (1-sens1_training)/spec1_training
        
        sens2_training <- with(trainingdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
        spec2_training <-  with(trainingdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
        LRplus2_training <- sens2_training/(1-spec2_training)
        LRminus2_training <- (1-sens2_training)/spec2_training
        
        # Calculate LR using random sampling
        sens1UQ <- by(valdf, valdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
        spec1UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
        LRplus1UQ <- (sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(1-spec1UQ+rnorm(length(sens1UQ),1e-6,1e-12))
        LRplus1dist <- formatC(quantile(unlist(LRplus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRplus1dist[2]," (",LRplus1dist[1]," - ",LRplus1dist[3],")"))
        LRminus1UQ <- (1-sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(spec1UQ+rnorm(length(sens1UQ),1e-12,1e-12))
        LRminus1dist <- formatC(quantile(unlist(LRminus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRminus1dist[2]," (",LRminus1dist[1]," - ",LRminus1dist[3],")"))
        
        sens2UQ <- by(valdf, valdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
        spec2UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
        LRplus2UQ <- (sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(1-spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
        LRplus2dist <- formatC(quantile(unlist(LRplus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRplus2dist[2]," (",LRplus2dist[1]," - ",LRplus2dist[3],")"))
        LRminus2UQ <- (1-sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
        LRminus2dist <- formatC(quantile(unlist(LRminus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRminus2dist[2]," (",LRminus2dist[1]," - ",LRminus2dist[3],")"))
        
        sens1UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
        spec1UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
        LRplus1UQ_training <- (sens1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))/(1-spec1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))
        LRplus1dist_training <- formatC(quantile(unlist(LRplus1UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRplus1dist_training[2]," (",LRplus1dist_training[1]," - ",LRplus1dist_training[3],")"))
        LRminus1UQ_training <- (1-sens1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))/(spec1UQ_training+rnorm(length(sens1UQ_training),1e-12,1e-12))
        LRminus1dist_training <- formatC(quantile(unlist(LRminus1UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRminus1dist_training[2]," (",LRminus1dist_training[1]," - ",LRminus1dist_training[3],")"))
        
        sens2UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
        spec2UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
        LRplus2UQ_training <- (sens2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))/(1-spec2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))
        LRplus2dist_training <- formatC(quantile(unlist(LRplus2UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRplus2dist_training[2]," (",LRplus2dist_training[1]," - ",LRplus2dist_training[3],")"))
        LRminus2UQ_training <- (1-sens2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))/(spec2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))
        LRminus2dist_training <- formatC(quantile(unlist(LRminus2UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRminus2dist_training[2]," (",LRminus2dist_training[1]," - ",LRminus2dist_training[3],")"))
        
        # Calculate mean error
        meanerror <- mean(valdf$pred_err)
        lowererror <- meanerror - 1.96* sd(valdf$pred_err)/sqrt(dim(valdf)[1]); uppererror <- meanerror+1.96*sd(valdf$pred_err)/sqrt(dim(valdf)[1])
        errordist <- formatC(c(meanerror,lowererror,uppererror),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(errordist[1]," (",errordist[2]," - ",errordist[3],")"))
        
        meanerror_training <- mean(trainingdf$pred_err)
        lowererror_training <- meanerror_training - 1.96* sd(trainingdf$pred_err)/sqrt(dim(trainingdf)[1]); uppererror_training <- meanerror_training+1.96*sd(trainingdf$pred_err)/sqrt(dim(trainingdf)[1])
        errordist_training <- formatC(c(meanerror_training,lowererror_training,uppererror_training),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(errordist_training[1]," (",errordist_training[2]," - ",errordist_training[3],")"))
        
        # LOOCV results for ranking
        altvec <- c(altvec,paste0(AUC1dist[2]))
        altvec <- c(altvec,paste0(AUC2dist[2]))
        altvec <- c(altvec,paste0(pairwisedist[2]))
        altvec <- c(altvec,paste0(LRplus1dist[2]))
        altvec <- c(altvec,paste0(LRminus1dist[2]))
        altvec <- c(altvec,paste0(LRplus2dist[2]))
        altvec <- c(altvec,paste0(LRminus2dist[2]))
        altvec <- c(altvec,paste0(errordist[1]))
        
        # altvec <- c(altvec,paste0(AUC1dist[1]))
        # altvec <- c(altvec,paste0(AUC2dist[1]))
        # altvec <- c(altvec,paste0(pairwisedist[1]))
        # altvec <- c(altvec,paste0(LRplus1dist[1]))
        # altvec <- c(altvec,paste0(LRminus1dist[3]))
        # altvec <- c(altvec,paste0(LRplus2dist[1]))
        # altvec <- c(altvec,paste0(LRminus2dist[3]))
        # altvec <- c(altvec,paste0(errordist[3]))
        
        
        # Training results for ranking
        altvec <- c(altvec,paste0( AUC1dist_training[2]))
        altvec <- c(altvec,paste0( AUC2dist_training[2]))
        altvec <- c(altvec,paste0(pairwisedist_training[2]))
        altvec <- c(altvec,paste0(LRplus1dist_training[2]))
        altvec <- c(altvec,paste0(LRminus1dist_training[2]))
        altvec <- c(altvec,paste0(LRplus2dist_training[2]))
        altvec <- c(altvec,paste0(LRminus2dist_training[2]))
        altvec <- c(altvec,paste0(errordist_training[1]))
        
        # altvec <- c(altvec,paste0( AUC1dist_training[1]))
        # altvec <- c(altvec,paste0( AUC2dist_training[1]))
        # altvec <- c(altvec,paste0(pairwisedist_training[1]))
        # altvec <- c(altvec,paste0(LRplus1dist_training[1]))
        # altvec <- c(altvec,paste0(LRminus1dist_training[3]))
        # altvec <- c(altvec,paste0(LRplus2dist_training[1]))
        # altvec <- c(altvec,paste0(LRminus2dist_training[3]))
        # altvec <- c(altvec,paste0(errordist_training[3]))
        
        # Thresholds
        altvec <- c(altvec,paste0(t1))
        altvec <- c(altvec,paste0(t2))
        
        # The normalized loglikelihood value
        altvec <- c(altvec,paste0(mean(normalized_logLik)))
        
        # Update altvec and all28df
        altvec <- altvec[-1]
        all28df[,altstr] <- altvec
      }
    }
    
    # Write results to file
    
    table7 <- data.frame(metric=rep(measurevec,each=length(datasetvec)), 
                         dataset=rep(datasetvec, length(measurevec)),
                         qNet=NA)
    table7[,-1:-2] <- matrix(unlist(t(all28df[,-1])),nrow= length(measurevec)*length(datasetvec),byrow=T)
    outfile <- paste0(simdir,"/",drug_index,"_Table7.txt")
    write.table(table7, outfile, col.names=T,row.names=F,sep="\t",quote=F)
    table7 <- subset(table7, select=-c(dataset))
    #table7 <- subset(table7, select=-c(1))
    names(table7) <- NULL
    table7 <- setNames(data.frame(t(table7[,-1])), table7[,1])
    drugnames <- t(drugnames)
    colnames(drugnames) <- c(to_vec(for(i in 1:drug_size) paste0("drug_",i)))
    
    final_result <- cbind(drugnames, table7[1,])
    #print(final_result)
    result_df <- rbind(result_df, final_result)
    #print(result_df)
    outfile <- paste0(outdir, "/",metric, "metrics.csv")
    write.csv(result_df, outfile ,row.names=F)
  } #if is_training_all_error
}

rankscorefun_loocv <- function(pmeasures,
                               is_normalized = TRUE){
  Performance_measures <- c("AUC1", "AUC2",
                            "LR_pos_th1", "LR_pos_th2",
                            "LR_neg_th1", "LR_neg_th2",
                            "Pairwise","Classification_error")
  Performance_levels <- c("Excellent_performance", 
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  # A dataframe for the weights of performance measures
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0)
  )
  
  # A dataframe for the weights of performance levels
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA # Not acceptable performance is removed
  
  # Initialize the performance dataframe for the model
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA, NA)
  )
  
  # Check the performance level for each performance measure 
  # by looking at the 95% confidence interval that match the CiPA's criteria
  
  # AUC1
  score_to_check <- pmeasures$AUC1_LOOCV
  # print("AUC1_LOOCV")
  # print(score_to_check)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC2
  score_to_check <- pmeasures$AUC2_LOOCV
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th1
  score_to_check <- pmeasures$LR1plus_LOOCV
  # print("LR1plus_LOOCV")
  # print(score_to_check)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_pos_th2
  score_to_check <- pmeasures$LR2plus_LOOCV
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_pos_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th1
  score_to_check <- pmeasures$LR1minus_LOOCV
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_neg_th2
  score_to_check <- pmeasures$LR2minus_LOOCV
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_neg_th2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Pairwise
  score_to_check <- pmeasures$Pairwise_LOOCV
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error
  score_to_check <- pmeasures$Mean_error_LOOCV
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Calculate rank score
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  
  return(as.numeric(rank_score))
}
