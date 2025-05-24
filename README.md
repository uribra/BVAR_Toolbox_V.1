# Bayesian Vector Autoregression (BVAR) Toolbox 

# BVAR_Toolbox_V.1

**BVAR_Toolbox_V.1** is a MATLAB toolbox for estimating **Bayesian Vector Autoregression (BVAR)** models with multiple prior options. It includes tools for **forecast error variance decomposition (FEVD)** and **structural impulse response functions (SIRFs)**, widely used in macroeconomic analysis.

---

## 🔍 Features

- Estimate Standard OLS VAR and BVAR models
- Support for multiple Bayesian priors
- Compute:
  - **Forecast Error Variance Decomposition (FEVD)**
  - **Structural Impulse Response Functions (SIRFs)**
- MATLAB-native, well-commented, and easy to modify

---

## 📦 Priors Supported

| Code | Description |
|------|-------------|
| `00` | Standard OLS VAR |
| `01` | Diffuse (Jeffrey’s) Prior |
| `02` | *(Planned)* Original Minnesota Prior |
| `03` | Independent Normal-Diffuse Prior (Minnesota moments) |
| `04` | Normal-Inverted-Wishart Prior (Kadiyala & Karlsson, 1997) |
| `05` | Dummy Observations Prior (Banbura et al., 2010) |

---

## 📘 Example: Stock & Watson (2001)

This toolbox includes a working example replicating part of the analysis from:

> **Stock, J. H., & Watson, M. W. (2001)**  
> *Vector Autoregressions*, Journal of Economic Perspectives.

The example uses macroeconomic U.S. time series data (e.g., GDP, inflation, interest rates) to demonstrate:

- Estimation of BVAR model with a Normal-Inverted-Wishart prior
- Structural impulse response functions (SIRFs) to monetary policy shocks
- Forecast error variance decomposition (FEVD) of output and inflation


[IRF_3.pdf](https://github.com/user-attachments/files/20421876/IRF_3.pdf)

---

## ⚙️ Installation

```bash
git clone https://github.com/your-username/BVAR_Toolbox_V.1.git
addpath(genpath('path_to/BVAR_Toolbox_V.1'));

