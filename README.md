# The specific role of parental humiliation in the guilt–suicidal ideation pathway: a Moderated Network Analysis

This repository contains the R code for the analyses of the experimental bachelor's thesis *"The specific role of parental humiliation in the guilt–suicidal ideation pathway: A Moderated Network Analysis"* (University of Rome Tor Vergata, supervisor Prof. Rodolfo Rossi).

## Overview

The study investigates the network structure of depressive symptoms in a community sample of Italian late adolescents and tests whether early family adversity moderates specific symptom–symptom associations.

- **Symptoms:** the nine items of the Patient Health Questionnaire-9 (PHQ-9; `phq01`–`phq09`), used as network nodes.
- **Moderators:** the total score and the single items of a reformulated 5-item version of the Risky Families Questionnaire (RFQ; `rfq_tot`, `rfq02`, `rfq05`, `rfq07`, `rfq12`, `rfq13`).
- **Sample:** N = 969 adolescents with complete PHQ-9 data for the Gaussian Graphical Model; the Moderated Network Models use complete cases for each moderator.

**Main finding:** the MNM with the RFQ total score selected a single three-way interaction (RFQ × guilt × suicidal ideation). When the items were tested separately, only parental verbal abuse (insults, humiliation, threats; `rfq02`) showed this interaction.

## Repository structure

```
.
├── README.md
└── code/
    └── NA_dt_ACE_depression.R   # complete analysis pipeline
```

## Analysis pipeline

The script `code/NA_dt_ACE_depression.R` runs the following steps:

1. **Data preparation** – selection of PHQ-9 items and RFQ moderators; listwise deletion of incomplete cases.
2. **Exploratory Graph Analysis (EGA)** – dimensionality of the PHQ-9 (graphical LASSO, Walktrap algorithm) and bootstrap EGA (parametric).
3. **Gaussian Graphical Model (GGM)** – network estimation with EBICglasso (`bootnet`) and visualization (`qgraph`).
4. **Accuracy and stability** – nonparametric bootstrap of edge weights, case-dropping bootstrap and correlation-stability (CS) coefficients for Strength and Expected Influence, centrality plots.
5. **Sample size simulation** – `netSimulator` with N = 200–1,000.
6. **Moderated Network Models (MNM)** – `mgm` with each RFQ variable as moderator (EBIC, γ = 0.5, AND rule, scaled variables), inspection of three-way interactions and factor graphs (standard and nodewise).
7. **Conditional networks** – networks conditioned on three levels of `rfq02` (never, sometimes, often).
8. **Bootstrap of the MNM** – `resample` for the models with `rfq_tot` and `rfq02`.

## Requirements

R packages: `bootnet`, `qgraph`, `networktools`, `mgm`, `tidyverse`, `EGAnet`, `igraph` (`parallel` is part of base R).

```r
install.packages(c("bootnet", "qgraph", "networktools", "mgm",
                   "tidyverse", "EGAnet", "igraph"))
```


## Authors

Paolo Amicone, Tommaso B. Jannini, Paolo Ossola, Martino Belvederi Murri, Grazia Terrone, Cinzia Niolu, Adriano Schimmenti, Patrizia Velotti, Rodolfo Rossi
