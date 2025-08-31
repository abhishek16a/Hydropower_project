# Hydropower Decision Support System

This repository contains the implementation of a **probabilistic decision support framework** for hydropower load management in Nepal.  
The project integrates **hydrological data, load forecasting, and text-derived signals** into a unified system that helps operators manage risks such as **spill, shortfall, and high load**.

---

## Project Overview
Hydropower in Nepal is highly dependent on **rainfall, river discharge, and reservoir levels**, while operations are further influenced by **maintenance notices, weather bulletins, and policy updates**.  
This project demonstrates how to combine these diverse signals using:

- **Gaussian Process Regression (GPR):** for daily load forecasting with uncertainty intervals.  
- **Bayesian Networks (BN):** for risk reasoning and decision inference.  
- **Topic Modelling (LDA + keyword flags):** to convert unstructured notices into structured features.  

---

##  Repository Structure
Hydropower_project/
│
├── 02_clean_features.ipynb # Feature engineering from hydrology + text
├── 03_text_topics.ipynb # Text preprocessing & topic modelling
├── 04_exploratory_analysis.ipynb # Exploratory plots & data checks
├── 05_gp_model.ipynb # Gaussian Process regression model
├── 06_bn_model.ipynb # Bayesian Network structure & inference
├── 06_decision_support.ipynb # Translate risks into operator advisories
├── 07_bayesian_network.ipynb # Extended BN modelling
├── 08_evaluation_reporting.ipynb # Model metrics & evaluation plots
├── 09_reporting_figures_and_tables.ipynb # Final publication-ready outputs
│
├── data/ # Input datasets (hydrology, rainfall, text)
├── figures/ # Exported figures used in paper/report
├── appendix/ # Additional experimental figures
│
└── README.md # Project description (this file)


---
# Results(summary)
Forecasting Accuracy: GPR achieved RMSE ≈ 129 MW and MAPE < 6% on test data.

Risk Alerts: BN achieved high recall (> 0.95) for shortfall risks, though precision was modest.

Decision Support: The system generated daily advisories such as “Spill advisory,” “Monitor inflows,” and “Shortfall mitigation.”

##  Installation & Requirements
Clone the repository:
```bash
git clone https://github.com/abhishek16a/Hydropower_project.git
cd Hydropower_project
