# ST001683-Metabolic-Anchoring-

# Metabolic Anchoring and Competitive Substrate Diversion by *Bacteroides thetaiotaomicron* (Project ST001683)

This repository contains the R scripts and processed data for the study:  
**"Bacteroides thetaiota
omicron regulates host homeostasis through metabolic anchoring and competitive substrate diversion"**

## 🔬 Overview
This project characterizes a novel metabolic phenomenon termed **"Metabolic Anchoring"** using a high-resolution multi-organ metabolomics approach. We utilize the ST001683 dataset (Bt-colonized vs. Germ-Free mice) to demonstrate how *B. thetaiotaomicron* (Bt) dictates host systemic metabolism.

### Key Findings:
- **Metabolic Anchoring:** Bt colonization compresses global metabolic variance (CV) by up to 73.0%.
- **Mathematical Invariant:** Discovery of a constant suppression effect on the IDO pathway (Metabolic Control Index, **MCI ≈ -1.30**).
- **Bioenergetic Pivot:** Identification of a systemic shift from glycolysis to the TCA cycle (flux ratio = 6.15).

## 📂 Repository Structure
- `/Scripts`: R scripts for data preprocessing, PCA, MCI calculation, and FGSEA analysis.
- `/Data`: Link to the processed abundance matrices (Cecal, Feces, Serum, and Urine).
- `/Figures`: High-resolution versions of Figures 1-11.

## 🛠 Prerequisites
To reproduce the analysis, you will need **R (v4.5.2 or later)** and the following packages:
```R
install.packages(c("ggplot2", "dplyr", "pROC", "pheatmap", "ggpubr"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fgsea")

🚀 Quick Start
Clone the repository:
git clone https://github.com/bingmoon/ST001683-Metabolic-Anchoring-/tree/main
Run the main analysis script:
source("Scripts/analysis.R")
📊 Data Source
Raw data were acquired from the Metabolomics Workbench under Project ID: ST001683.

✉️ Contact
Corresponding Author: Hongbo He (moxizhen@gmail.com)
Institution: West China Hospital, Sichuan University
⚖️ License
This project is licensed under the MIT License - see the LICENSE file for details.
