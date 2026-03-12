# Bacteroides thetaiotaomicron Regulates Host Homeostasis

This repository contains the complete analytical pipeline, R scripts, and processed data matrices for our study on how *Bacteroides thetaiotaomicron* (Bt) functions as a metabolic architect to anchor host systemic homeostasis.

## 📌 Project Overview
Our study utilizes a **Cross-scale Mechanistic Validation** approach:
1. **In Vivo Discovery (ST001683):** Multi-organ metabolomics (cecum, feces, serum, urine) to reveal the "Metabolic Anchor" effect and entropy reduction.
2. **In Vitro Validation (ST001688):** Large-scale pure culture data ($N=1049$ samples) to mechanistically prove the "Substrate Sink" effect (Tryptophan consumption and IPA/5-HI production).

## 📂 Repository Structure
* `/Scripts/`
  * `01_In_Vivo_Systemic_Analysis.R` - Pipeline for ST001683 (PCA, CV, AUC, MCI, FGSEA).
  * `02_In_Vitro_Mechanistic_Validation.R` - Robust computational parser and Log2FC analysis for ST001688.
* `/Data/Cleaned_Matrices/`
  * *Note: Raw files from the ST001688 repository contained severe RTF formatting artifacts. We have provided the sanitized, ready-to-analyze CSV matrices here to ensure total reproducibility.*
  * `ST001688_HILIC_POS_clean.csv` (Tryptophan extraction)
  * `ST001688_RP_POS_clean.csv` (IPA extraction)
  * `ST001688_RP_NEG_clean.csv` (5-HIAA extraction)

## 📊 Key Output
The pipeline automatically generates all main figures and supplementary tables used in the manuscript, including the highly deterministic dose-response coupling and the cross-scale Log2FC validation plot (Figure 12).

## 🛠 Usage
All scripts are written in R (v4.4.1+). Ensure that the required packages (`dplyr`, `ggplot2`, `ggpubr`, `tidyr`, `pROC`, `fgsea`, `pheatmap`) are installed before running the pipeline.
