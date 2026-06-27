# 12-CpG Methylation SNaPshot Multiplex Assay for Forensic Age Estimation from Semen

[![DOI](https://pfst.cf2.poecdn.net/base/image/2b99cfd2da8ae0e5e1a03c9cd25c3f18db1a62ad6d3c63a3c087908beb0b798e?pmaid=632892806)](https://doi.org/)
[![License: MIT](https://pfst.cf2.poecdn.net/base/image/0f22612a092cd1947cd99e9a61a24fc5c955c4a2986d89bc0991dfcd0a8d42f2?pmaid=632892807)](LICENSE)

## Overview

This repository provides trained age estimation models for estimating chronological age from semen/sperm DNA methylation data using a 12-CpG SNaPshot multiplex panel. The assay is designed for direct deployment on capillary electrophoresis (CE) instruments commonly available in forensic DNA laboratories.

**Key performance:**

- **SVR model:** MAE = 2.82 years (independent test set, *n* = 76); CV-MAE = 2.86 ± 0.44 years (100×10-fold nested CV, *n* = 260)
- **MLR model:** MAE = 3.41 years (test set); CV-MAE = 3.24 ± 0.47 years (nested CV); transparent linear equation
- **MQR model:** MAE = 3.27 years (test set); CV-MAE = 3.22 ± 0.47 years (nested CV); provides 80% estimation intervals

**Population:** Han Chinese males (*n* = 260, aged 22–67 years)  

**Tissue:** Semen (sperm genomic DNA after differential lysis)  

**Platform:** SNaPshot single-base extension on ABI 3130 Genetic Analyzer

## Repository Structure

```text
semen-age-snapshot-12cpg/
├── README.md
├── LICENSE
├── R/
│   └── estimate.R              # Estimate function & model info printer
├── models/
│   ├── mlr_full.rds            # Full-sample MLR model (stats::lm)
│   ├── svr_full.rds            # Full-sample SVR model (e1071::svm)
│   ├── mqr_full.rds            # Full-sample MQR models (quantreg::rq, tau=0.10/0.50/0.90)
│   └── model_metadata.rds      # Marker info, hyperparameters, performance metrics
├── example/
│   ├── example_input.csv       # Example methylation input (5 samples)
│   └── run_example.R           # Demonstration script
```

## Quick Start

### Requirements

- R ≥ 4.1.0
- Required packages: `e1071` (for SVR), `quantreg` (for MQR)

```r
install.packages(c("e1071", "quantreg"))
```

### Usage

```r
# Source the estimation function
source("R/estimate_age.R")

# View model information
print_model_info(model_dir = "models/")

# Prepare your data: data.frame with 12 methylation ratio columns
my_data <- read.csv("my_data.csv", row.names = 1)

# Estimate age (SVR recommended)
results <- estimate_semen_age(my_data, model_dir = "models/", method = "SVR")
print(results)

# Get 80% estimation intervals (MQR only)
results_pi <- estimate_semen_age(my_data, model_dir = "models/", method = "MQR",
                                 estimation_interval = TRUE)
print(results_pi)
```

## Input Data Format

Your input must be a `data.frame` with the following **12 columns** (exact names required):

| Column Name      | Marker ID      | Direction | Detection |
| ---------------- | -------------- | --------- | --------- |
| `cg19998819`     | cg19998819     | Positive  | C/(C+T)   |
| `cg11262154`     | cg11262154     | Positive  | C/(C+T)   |
| `cg12277678`     | cg12277678     | Positive  | C/(C+T)   |
| `cg04123357`     | cg04123357     | Positive  | C/(C+T)   |
| `cg18037145`     | cg18037145     | Positive  | G/(G+A)   |
| `cg25187042`     | cg25187042     | Positive  | C/(C+T)   |
| `cg20602007`     | cg20602007     | Positive  | C/(C+T)   |
| `cg21843517`     | cg21843517     | Positive  | C/(C+T)   |
| `cg13872326`     | cg13872326     | Negative  | C/(C+T)   |
| `chr1.19339447`  | chr1:19339447  | Negative  | G/(G+A)   |
| `chr13.93039685` | chr13:93039685 | Positive  | C/(C+T)   |
| `chr19.18610678` | chr19:18610678 | Negative  | C/(C+T)   |

- All values must be **methylation ratios** between 0 and 1
- "Direction" indicates the age-correlation direction (positive = hypermethylation with age)
- "Detection" indicates the formula used to calculate the methylation ratio from SNaPshot peak heights

## Model Details

### SVR (Recommended)

- **Algorithm:** Support Vector Regression with radial basis function (RBF) kernel
- **Package:** `e1071::svm(type = "eps-regression")`
- **Hyperparameters:** cost = 128, gamma = 0.005, epsilon = 1.0
- **Internal scaling:** The `svm()` function applies automatic feature scaling (`scale = TRUE` by default); scaling parameters are stored within the model object
- **Performance:** CV-MAE = 2.86 ± 0.44 years, CV-RMSE = 3.75 ± 0.57, CV-R2 = 0.850 ± 0.057

### MLR (Interpretable reference)

- **Algorithm:** Ordinary least squares multiple linear regression

- **Package:** `stats::lm()`

- **Equation:**

  ```apache
  Age = −17.11 + 28.43×cg19998819 + 32.53×cg11262154 − 14.49×cg12277678
        + 17.81×cg04123357 + 28.60×cg18037145 − 18.07×cg25187042
        + 7.37×cg20602007 + 14.34×cg21843517 − 23.68×cg13872326
        − 13.93×chr1:19339447 + 4.67×chr13:93039685 + 33.87×chr19:18610678
  ```

- **Performance:** CV-MAE = 3.24 ± 0.47 years, CV-RMSE = 4.11 ± 0.55, CV-R2 = 0.822 ± 0.057

### MQR (Robust + uncertainty quantification)

- **Algorithm:** Multiple quantile regression at τ = 0.10, 0.50, 0.90
- **Package:** `quantreg::rq()`
- **Point estimate:** Median regression (τ = 0.50)
- **80% prediction interval:** [τ = 0.10, τ = 0.90]
- **Performance:** CV-MAE = 3.22 ± 0.47 years (median), CV-RMSE = 4.11 ± 0.57 (median), CV-R2 = 0.823 ± 0.053 (median)
- **Note:** Empirical PI under-coverage reflects estimation variance in moderately-sized training data; interpret intervals as approximate guidance

## Important Limitations

1. **Population specificity:** Models were trained exclusively on Han Chinese males. Application to other populations requires independent validation.
2. **Platform specificity:** Models are calibrated for the SNaPshot (CE-based) methylation quantification platform. Direct application of methylation values from other platforms (pyrosequencing, MPS, microarray) will introduce systematic bias. Cross-platform use requires platform-specific recalibration (see paper for details).
3. **Age range:** Training data spans 22–67 years, with limited representation above 50 years (10% of samples). Extrapolation beyond this range is not recommended.
4. **DNA input:** Reliable results require ≥ 10 ng genomic DNA input for bisulfite conversion (inter-replicate SD ≤ 1.55 years).
5. **Sample type:** Models were developed for ejaculated semen processed by differential lysis to enrich for sperm DNA. Performance on mixed body fluid samples (e.g., semen–vaginal secretion mixtures) has not been evaluated.
6. **Age-dependent bias:** Systematic under-estimation occurs for older donors (> 50 years) due to regression toward the mean in the skewed training distribution.
7. **Confounders:** Effects of smoking, BMI, chronic diseases, and fertility status were not assessed.

## Sensitivity and Stability

- **Minimum DNA input:** ≥ 10 ng for single-measurement casework application (preliminary, *n* = 5 donors)
- **Stain stability:** Preliminary evidence (*n* = 5, 118-day ambient storage) suggests estimation deviations ≤ 2 years from fresh semen

## Citation

If you use these models in your research, please cite:

> [Authors] (2026). Development and validation of a 12-CpG methylation SNaPshot multiplex assay for forensic age estimation from semen. *[Journal]*, *[Volume]*, [Pages]. https://doi.org/[DOI]

## Related Publications

- Xiao C, Li Y, Chen M, Yi S, Huang D (2023). Improved age estimation from semen using sperm-specific age-related CpG markers. *Forensic Sci Int Genet* 67:102941.
- Li Y, Liu X, Chen M, Yi S, He X, Xiao C, Huang D (2025). DNA methylation-based age estimation from semen: genome-wide marker identification and model development. *Forensic Sci Int Genet* 76:103215.

## License

This project is licensed under the MIT License. See LICENSE for details.

## Contact

For questions regarding the models or assay implementation, please open an issue or contact xiaochao369@hust.edu.cn.

## LICENSE

MIT License

Copyright (c) 2026 Huazhong University of Science and Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
