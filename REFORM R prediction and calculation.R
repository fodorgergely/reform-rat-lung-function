# ===============================================================================
# REFERENCE EQUATIONS FOR RESPIRATORY MECHANICS AND END-EXPIRATORY LUNG VOLUME IN RATS (REFORM)
# Lung mechanics and volume prediction sheet
#
# Please cite the original article if using this calculator.
#
# DISCLAIMER:
# To prevent extrapolation beyond observed data, we do not recommend not using the 
# calculator beyond mass ranges actually observed for each strain-sex combination:
#   Sprague-Dawley males (160-750g)
#   Sprague-Dawley females (160-380g)
#   Wistar males (160-540g)
#   Wistar females (160-370g)
#
# STATISTICAL APPROACH:
# - Raw, G, H use log-normal distributions: log(parameter) ~ Normal(μ, σ²)
# - EELV uses normal distribution: parameter ~ Normal(μ, σ²) 
# - Z-scores calculated as: z = (observed - predicted_μ) / predicted_σ
# - Percentiles based on normal distribution assumptions
#
# ===============================================================================


# ===============================================================================
# 0. Load the required packages
# built under R version 4.2.1, gamlss version 5.4-22
# ===============================================================================

library(gamlss)
library(dplyr)

# ===============================================================================
# 1. Load the trained models
# Trained models can be examined with summary(model_name), 
# where model_name should be replaced with appropriate model
# ===============================================================================

model_Raw <- readRDS("REFORM_Raw.rds")
model_G <- readRDS("REFORM_G.rds")
model_H <- readRDS("REFORM_H.rds")
model_EELV <- readRDS("REFORM_EELV.rds")

#summary(model_Raw)

# ===============================================================================
# 2. First we create the structure of the data frame for predictions
#    - Categorical values are recoded to numeric values for ease of use
#      - Strain: Sprague Dawley as 0, Wistar as 1
#      - Sex: Female as 0, Male as 1
#      - PEEP: in cmH2O, 0-6, with the exception of 5 cmH2O (no training data)
# ===============================================================================

newdata_struct <- 
  data.frame(
    rat = character(),                                       # optional
    strain = factor(levels = c("0", "1")),                   # SPRD as 0, Wistar as 1
    sex = factor(levels = c("0", "1")),                      # Female as 0, Male as 1
    peep = factor(levels = c("0", "1", "2", "3", "4", "6")), # level of PEEP in cmH2O (0, 1, 2, 3, 4, 6 cmH2O)
    mass = double(),                                         # body mass in grams
    Raw = double(),                                          # Measured airway resistance in cmH2O.s/l, optional, used only for z-scores
    G = double(),                                            # Measured tissue damping in cmH2O/l, optional, used only for z-scores
    H = double(),                                            # Measured tissue elastance in cmH2O/l, optional, used only for z-scores
    EELV = double(),                                         # Measured end-expiratory lung volume in ml, optional, used only for z-scores
    stringsAsFactors = FALSE
  )

# ===============================================================================
# 3. Adding new data to the data frame
#    - When adding new data manually, you can use a c() operator for each variable. 
#        The order of elements in these c() operators does matter, 
#        i.e., the first element in each c() will correspond to first rat,
#        the second element will correspond to second rat, etc.
#    - In case measured data for one variable is missing, please use NA, 
#        then predictions will be calculated, but not z-scores.
#    - Sample data presented here is the same as the two examples in the Online Data Supplement
# ===============================================================================

# Defining new data, coding and use of proper factor levels is a must.
new_rows <- data.frame(
  rat    = c("example_rat_001","example_rat_002"),
  strain = factor(c("0","1"), levels = c("0","1")),
  sex    = factor(c("0","1"), levels = c("0","1")),
  peep   = factor(c("2","0"), levels = c("0","1","2","3","4","6")),
  mass   = c(230,170),
  Raw    = c(41.82, NA),
  G      = c(NA, NA),
  H      = c(NA, NA),
  EELV   = c(NA, 2)
)

newdata <- bind_rows(newdata_struct, new_rows) # creating the actual prediction dataset

# ===============================================================================
# GAMLSS PREDICTION SCRIPT WITH PERCENTILES
# For standalone prediction of REFORM models without original training data
# All lines below need to be run at least once per session to define the functions
# Includes 5th and 95th percentiles as LLN and ULN
# Instructions for use are to be found below
# ===============================================================================

# Core prediction function
gamlss_predict <- function(model_path, newdata, what = "mu", type = "response") {
  # Predict GAMLSS model parameters
  
  model <- readRDS(model_path)
  coef_name <- paste0(what, ".coefficients")
  coefs <- model[[coef_name]]
  
  if (is.null(coefs)) {
    stop(paste("Parameter", what, "not found in model"))
  }
  
  results <- numeric(nrow(newdata))
  
  for (i in 1:nrow(newdata)) {
    row <- newdata[i, ]
    prediction <- coefs["(Intercept)"]
    
    
    if (what == "mu" && "sqrt(mass)" %in% names(coefs)) {
      # For mu: use sqrt(mass)
      prediction <- prediction + coefs["sqrt(mass)"] * sqrt(row$mass)
    } else if (what == "sigma" && "mass" %in% names(coefs)) {
      # For sigma: use mass (no sqrt)
      prediction <- prediction + coefs["mass"] * row$mass
    }
    
    # Add peep terms
    peep_level <- as.character(row$peep)
    peep_coef_name <- paste0("peep", peep_level)
    if (peep_coef_name %in% names(coefs)) {
      prediction <- prediction + coefs[peep_coef_name]
    }
    
    # Add sex term
    sex_level <- as.character(row$sex)
    sex_coef_name <- paste0("sex", sex_level)
    if (sex_coef_name %in% names(coefs)) {
      prediction <- prediction + coefs[sex_coef_name]
    }
    
    # Add strain term
    strain_level <- as.character(row$strain)
    strain_coef_name <- paste0("strain", strain_level)
    if (strain_coef_name %in% names(coefs)) {
      prediction <- prediction + coefs[strain_coef_name]
    }
    
    # Apply link function
    if (type == "response" && what == "sigma") {
      prediction <- exp(prediction)  # Sigma uses log link
    }
    # Mu uses identity link (no transformation)
    
    results[i] <- prediction
  }
  
  return(results)
}

# Function to calculate percentiles
calculate_percentiles <- function(mu_link, sigma_response, percentiles = c(0.05, 0.95), is_log_transformed = FALSE) {
  # Calculate percentiles for normal distribution
  # percentiles: vector of percentiles (e.g., c(0.05, 0.95) for 5th and 95th)
  # is_log_transformed: TRUE for Raw, G, H; FALSE for EELV
  
  # Convert percentiles to z-scores
  z_scores <- qnorm(percentiles)
  
  # Calculate percentiles
  if (is_log_transformed) {
    # For log-transformed: calculate on log scale then transform back
    log_percentiles <- mu_link + z_scores * sigma_response
    result_percentiles <- exp(log_percentiles)
  } else {
    # For direct: calculate directly
    result_percentiles <- mu_link + z_scores * sigma_response
  }
  
  return(result_percentiles)
}

# Main function to predict all REFORM parameters with percentiles
predict_reform <- function(newdata, 
                           model_dir = ".", 
                           model_prefix = "REFORM_") {
  # Predict all REFORM parameters with z-scores and percentiles
  # Raw, G, H use log transformation; EELV uses direct
  
  parameters <- c("Raw", "G", "H", "EELV")
  log_params <- c("Raw", "G", "H")
  
  # Build model paths
  model_paths <- list()
  for (param in parameters) {
    model_paths[[param]] <- file.path(model_dir, paste0(model_prefix, param, ".rds"))
  }
  
  results <- newdata
  
  cat("Processing REFORM models...\n")
  
  for (param in parameters) {
    model_path <- model_paths[[param]]
    
    if (!file.exists(model_path)) {
      warning(paste("Model file not found:", model_path))
      next
    }
    
    # Get predictions
    mu_pred <- gamlss_predict(model_path, newdata, "mu", "link")
    sigma_pred <- gamlss_predict(model_path, newdata, "sigma", "response")
    
    # Store basic results
    results[[paste0(param, "_mu")]] <- mu_pred
    results[[paste0(param, "_sigma")]] <- sigma_pred
    
    # Calculate response scale prediction
    if (param %in% log_params) {
      results[[paste0(param, "_predicted")]] <- exp(mu_pred)
    } else {
      results[[paste0(param, "_predicted")]] <- mu_pred
    }
    
    # Calculate 5th and 95th percentiles (LLN and ULN)
    is_log <- param %in% log_params
    percentiles_5_95 <- matrix(nrow = nrow(newdata), ncol = 2)
    
    for (i in 1:nrow(newdata)) {
      percs <- calculate_percentiles(mu_pred[i], sigma_pred[i], c(0.05, 0.95), is_log)
      percentiles_5_95[i, ] <- percs
    }
    
    results[[paste0(param, "_LLN")]] <- percentiles_5_95[, 1]  # 5th percentile
    results[[paste0(param, "_ULN")]] <- percentiles_5_95[, 2]  # 95th percentile
    
    # Calculate z-scores if measured data available
    if (param %in% names(newdata)) {
      measured_values <- newdata[[param]]
      z_scores <- rep(NA, length(measured_values))
      non_missing <- !is.na(measured_values)
      
      if (any(non_missing)) {
        if (param %in% log_params) {
          # For log models: z = (log(measured) - mu_link) / sigma_response
          log_measured <- log(measured_values[non_missing])
          z_scores[non_missing] <- (log_measured - mu_pred[non_missing]) / sigma_pred[non_missing]
        } else {
          # For direct models: z = (measured - mu_link) / sigma_response
          z_scores[non_missing] <- (measured_values[non_missing] - mu_pred[non_missing]) / sigma_pred[non_missing]
        }
      }
      
      results[[paste0(param, "_z_score")]] <- z_scores
    }
    
    cat("  ✓", param, "processed\n")
  }
  
  return(results)
}

# Create summary report with percentiles
create_summary <- function(results, subject_id_col = "rat") {
  # Create readable summary of predictions, percentiles, and z-scores
  
  cat("=== REFORM PREDICTION SUMMARY ===\n\n")
  
  parameters <- c("Raw", "G", "H", "EELV")
  param_info <- list(
    Raw = list(name = "Airway Resistance", units = "cmH2O.s/L"),
    G = list(name = "Tissue Damping", units = "cmH2O/L"),
    H = list(name = "Tissue Elastance", units = "cmH2O/L"),
    EELV = list(name = "End-Expiratory Lung Volume", units = "mL")
  )
  
  for (i in 1:nrow(results)) {
    # Subject header
    if (subject_id_col %in% names(results)) {
      cat("SUBJECT:", results[[subject_id_col]][i], "\n")
    } else {
      cat("SUBJECT", i, "\n")
    }
    
    # Demographics with correct factor handling
    strain_val <- as.character(results$strain[i])
    sex_val <- as.character(results$sex[i])
    peep_val <- as.character(results$peep[i])
    
    strain_name <- ifelse(strain_val == "0", "Sprague-Dawley", "Wistar")
    sex_name <- ifelse(sex_val == "0", "Female", "Male")
    
    cat("Demographics: ", strain_name, ", ", sex_name, ", PEEP ", peep_val, " cmH2O, Mass ", results$mass[i], " g\n\n")
    
    # Predictions with percentiles
    cat("Predictions (Mean ± SD, [5th-95th percentiles]):\n")
    for (param in parameters) {
      pred_col <- paste0(param, "_predicted")
      sigma_col <- paste0(param, "_sigma")
      lln_col <- paste0(param, "_LLN")
      uln_col <- paste0(param, "_ULN")
      
      if (all(c(pred_col, sigma_col, lln_col, uln_col) %in% names(results))) {
        pred_val <- results[[pred_col]][i]
        sigma_val <- results[[sigma_col]][i]
        lln_val <- results[[lln_col]][i]
        uln_val <- results[[uln_col]][i]
        info <- param_info[[param]]
        
        cat("  ", info$name, ": ", round(pred_val, 2), " ± ", round(sigma_val, 3), 
            " [", round(lln_val, 2), "-", round(uln_val, 2), "] ", info$units, "\n", sep = "")
      }
    }
    
    # Z-scores and interpretation
    z_cols <- paste0(parameters, "_z_score")
    available_z <- intersect(z_cols, names(results))
    
    if (length(available_z) > 0) {
      cat("\nMeasured vs Predicted:\n")
      for (z_col in available_z) {
        param <- gsub("_z_score", "", z_col)
        z_val <- results[[z_col]][i]
        measured_val <- results[[param]][i]
        lln_val <- results[[paste0(param, "_LLN")]][i]
        uln_val <- results[[paste0(param, "_ULN")]][i]
        
        if (!is.na(z_val)) {
          # Interpretation based on percentiles
          if (measured_val >= lln_val && measured_val <= uln_val) {
            status <- "Normal (5th-95th percentile)"
            flag <- "✓"
          } else if (measured_val < lln_val) {
            status <- "Below LLN (< 5th percentile)"
            flag <- "⬇"
          } else {
            status <- "Above ULN (> 95th percentile)"
            flag <- "⬆"
          }
          
          info <- param_info[[param]]
          cat("  ", info$name, ": ", round(measured_val, 2), " ", info$units, 
              " (z = ", round(z_val, 2), ", ", status, ") ", flag, "\n", sep = "")
        }
      }
    }
    
    cat("\n")
  }
  
  cat("Legend:\n")
  cat("  ✓ Normal: within 5th-95th percentile range\n")
  cat("  ⬇ Low: below 5th percentile (LLN)\n") 
  cat("  ⬆ High: above 95th percentile (ULN)\n")
  cat("  LLN = Lower Limit of Normal, ULN = Upper Limit of Normal\n")
}

# Export results to CSV
export_results <- function(results, filename = "reform_predictions.csv") {
  # Export predictions to CSV file
  write.csv(results, filename, row.names = FALSE)
  cat("Results exported to:", filename, "\n")
  
  # Show column summary
  param_cols <- grep("_(predicted|LLN|ULN|z_score)$", names(results), value = TRUE)
  cat("Exported columns include:\n")
  cat("  - Predictions: *_predicted\n")
  cat("  - Lower limits: *_LLN (5th percentile)\n")
  cat("  - Upper limits: *_ULN (95th percentile)\n") 
  cat("  - Z-scores: *_z_score\n")
}

# Check available models
check_models <- function(model_dir = ".", model_prefix = "REFORM_") {
  # Check which model files are available
  
  parameters <- c("Raw", "G", "H", "EELV")
  available <- c()
  
  cat("Checking REFORM models:\n")
  
  for (param in parameters) {
    filepath <- file.path(model_dir, paste0(model_prefix, param, ".rds"))
    if (file.exists(filepath)) {
      available <- c(available, param)
      cat("  ✓", param, "\n")
    } else {
      cat("  ✗", param, "(missing)\n")
    }
  }
  
  cat("\nReady to predict:", length(available), "of 4 parameters\n")
  return(available)
}

# Function to get just the percentiles table
get_percentiles_table <- function(results) {
  # Extract just the percentiles for easy viewing
  
  parameters <- c("Raw", "G", "H", "EELV")
  param_info <- list(
    Raw = list(name = "Airway Resistance", units = "cmH2O.s/L"),
    G = list(name = "Tissue Damping", units = "cmH2O/L"),
    H = list(name = "Tissue Elastance", units = "cmH2O/L"),
    EELV = list(name = "End-Expiratory Lung Volume", units = "mL")
  )
  
  cat("=== PERCENTILES TABLE ===\n\n")
  
  for (i in 1:nrow(results)) {
    if ("rat" %in% names(results)) {
      cat("Subject:", results$rat[i], "\n")
    } else {
      cat("Subject", i, "\n")
    }
    
    cat("Parameter                    Predicted   LLN(5th)   ULN(95th)  Units\n")
    cat("--------------------------------------------------------\n")
    
    for (param in parameters) {
      pred_col <- paste0(param, "_predicted")
      lln_col <- paste0(param, "_LLN")
      uln_col <- paste0(param, "_ULN")
      
      if (all(c(pred_col, lln_col, uln_col) %in% names(results))) {
        pred_val <- results[[pred_col]][i]
        lln_val <- results[[lln_col]][i]
        uln_val <- results[[uln_col]][i]
        info <- param_info[[param]]
        
        cat(sprintf("%-25s %9.2f %10.2f %10.2f  %s\n", 
                    info$name, pred_val, lln_val, uln_val, info$units))
      }
    }
    cat("\n")
  }
}

# ===============================================================================
# USAGE EXAMPLES

# REFORM PREDICTION SCRIPT WITH PERCENTILES LOADED
# Key features:
# - mu uses sqrt(mass), sigma uses mass
# - Raw, G, H use log transformation
# - EELV uses direct modeling
# - Automatic z-score calculation
# - 5th and 95th percentiles (LLN/ULN)

# Main functions:
# - check_models()                    Check available model files
# - predict_reform(newdata)           Get all predictions, percentiles, z-scores
# - create_summary(results)           Create readable report
# - get_percentiles_table(results)    Show percentiles table
# - export_results(results)           Export to CSV

# Quick start:
# results <- predict_reform(newdata)
# create_summary(results)
# get_percentiles_table(results)
# export_results(results, "predictions_REFORM.csv")
# ===============================================================================

# Get predictions with percentiles
results <- predict_reform(newdata)

# Full summary with percentiles
create_summary(results)

# Just the percentiles table
get_percentiles_table(results)

# Export results
export_results(results, "predictions_REFORM.csv")
