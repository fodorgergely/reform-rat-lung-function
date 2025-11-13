# REFERENCE EQUATIONS FOR RESPIRATORY MECHANICS AND END-EXPIRATORY LUNG VOLUME IN RATS (REFORM)
## Overview
Predictive equations for respiratory mechanics parameters in laboratory rats, enabling calculation of predicted values, z-scores, and reference percentiles.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17055086.svg)](https://doi.org/10.5281/zenodo.17055086)


**Status**: Manuscript under review.

## Citation
If you use this calculator, please cite:
> Fodor GH et al. Standardized lung function reference values in rats for translational respiratory research. *(under review)*. 2025.

## Grant support
This work was supported by the National Research, Development and Innovation Office of Hungary in the form of two grants: OTKA-NKFIH FK134274 and K138032.

## Quick Start
```r
# Load required packages
library(gamlss)
library(dplyr)

# Run the calculator
source("REFORM_calculator.R")
results <- predict_reform(newdata)
create_summary(results)
```

## Requirements
- R
- gamlss package
- dplyr package

## Parameters Predicted
- **Raw**: Airway resistance (cmH2O.s/L)
- **G**: Tissue damping (cmH2O/L)  
- **H**: Tissue elastance (cmH2O/L)
- **EELV**: End-expiratory lung volume (mL)

## Input Variables
- **strain**: Sprague-Dawley (0) or Wistar (1)
- **sex**: Female (0) or Male (1)
- **peep**: PEEP level in cmH2O (0, 1, 2, 3, 4, 6)
- **mass**: Body mass in grams

## DISCLAIMER
Do not use beyond observed mass ranges:
- Sprague-Dawley males: 160-750g
- Sprague-Dawley females: 160-380g
- Wistar males: 160-540g  
- Wistar females: 160-370g

## Statistical Approach
- Raw, G, H: log-normal distributions
- EELV: normal distribution
- Z-scores: (observed - predicted_μ) / predicted_σ
- Percentiles: 5th and 95th for reference limits

## Files
- `REFORM R prediction and calculation.R` - Main prediction script implemented in R
- `REFORM_*.rds` - Trained GAMLSS models
- `REFORM Excel calculator.xlsx` - Prediction implemented in Excel
- Example data included

## Contact
[fodor.gergely@med.u-szeged.hu] for questions or issues

## License
MIT License - see LICENSE file
