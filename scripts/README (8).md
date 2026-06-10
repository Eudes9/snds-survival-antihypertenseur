# Étude pharmaco-épidémiologique — Abandon du traitement antihypertenseur

[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Reproducible](https://img.shields.io/badge/Reproducible-Yes-brightgreen)](data/00_simulate_snds.R)

## Contexte

Ce projet simule une étude pharmaco-épidémiologique de type **SNDS** (Système National des Données de Santé) portant sur l'**abandon du traitement antihypertenseur** dans les 12 mois suivant l'initiation.

Il reproduit les standards méthodologiques appliqués dans les études du **GIS EPI-PHARE** (ANSM/CNAM) sur les données médico-administratives françaises.

> ⚠️ Les données sont entièrement **simulées** à des fins pédagogiques. Aucune donnée réelle de patient n'est utilisée.

---

## Question de recherche

> **Quel est le délai d'abandon du traitement antihypertenseur dans les 12 mois suivant l'initiation, et existe-t-il une différence significative entre les inhibiteurs de l'enzyme de conversion (IEC) et les antagonistes des récepteurs de l'angiotensine II (ARA2) ?**

---

## Structure du projet

```
snds-survival-antihypertenseur/
│
├── data/
│   ├── 00_simulate_snds.R          # Simulation base de type SNDS
│   └── *.rds                       # Données simulées (générées localement)
│
├── scripts/
│   ├── 01_data_quality.R           # Contrôle qualité
│   ├── 02_cohort_selection.R       # Construction de la cohorte + flow chart
│   ├── 03_propensity_score.R       # Score de propension + appariement
│   ├── 04_survival_analysis.R      # Kaplan-Meier + modèle de Cox
│   └── 05_sensitivity_analysis.R  # Analyses de sensibilité
│
├── report/
│   └── rapport_final.Rmd           # Rapport reproductible complet
│
├── shiny/
│   └── app.R                       # Dashboard interactif R Shiny
│
└── outputs/
    ├── figures/                    # Graphiques (KM, forest plots, love plot...)
    └── tables/                     # Tableaux (Tableau 1, HR, Cox...)
```

---

## Méthodologie

### Population
- Patients de 40 à 80 ans, primo-initiateurs d'un traitement antihypertenseur
- Fenêtre d'étude : 2020–2023 | Suivi maximum : 365 jours

### Données simulées — Structure SNDS
| Table simulée | Équivalent SNDS réel | Contenu |
|---|---|---|
| `ir_ben_r` | IR_BEN_R | Caractéristiques bénéficiaires |
| `er_pha_f` | ER_PHA_F (DCIR) | Délivrances médicaments |
| `mco_b`    | MCO_B (PMSI) | Hospitalisations |
| `ald`      | Table ALD | Affections longue durée |

### Pipeline analytique
1. **Contrôle qualité** — complétude, cohérence temporelle, valeurs aberrantes, doublons
2. **Sélection de la cohorte** — critères d'inclusion/exclusion, diagramme de flux
3. **Score de propension** — modèle logistique, appariement 1:1 nearest neighbor (caliper = 0.2 SD), love plot
4. **Analyse de survie** — Kaplan-Meier, test log-rank, modèle de Cox multivarié, test de Schoenfeld
5. **Analyses de sensibilité** — définitions alternatives, sous-groupes, risques compétitifs, E-value

---

## Packages R utilisés

```r
# Analyse de survie
library(survival)     # Cox, Kaplan-Meier
library(survminer)    # Visualisation KM

# Score de propension
library(MatchIt)      # Appariement
library(tableone)     # Tableau 1

# Manipulation & visualisation
library(dplyr)
library(ggplot2)
library(broom)        # Tidying des modèles
library(lubridate)

# Rapport
library(rmarkdown)
library(knitr)
library(kableExtra)
```

---

## Reproductibilité

Tous les scripts utilisent `set.seed(20260101)` pour garantir la reproductibilité totale des résultats.

Pour reproduire l'analyse :

```r
# 1. Simuler les données
source("data/00_simulate_snds.R")

# 2. Exécuter le pipeline dans l'ordre
source("scripts/01_data_quality.R")
source("scripts/02_cohort_selection.R")
source("scripts/03_propensity_score.R")
source("scripts/04_survival_analysis.R")
source("scripts/05_sensitivity_analysis.R")

# 3. Générer le rapport
rmarkdown::render("report/rapport_final.Rmd")
```

---

## Résultats principaux

*(générés après exécution du pipeline)*

- **Taux d'abandon global** : ~40% à 12 mois
- **Classe avec meilleure observance** : IEC/ARA2 (vs BB, diurétiques)
- **Facteurs associés à l'abandon** : âge jeune, absence de comorbidité, faible nombre de délivrances antérieures

---

## Références méthodologiques

- Austin PC (2011). An Introduction to Propensity Score Methods. *Multivariate Behavioral Research*.
- Cox DR (1972). Regression Models and Life-Tables. *JRSS-B*.
- Hernán MA, Robins JM (2016). Using Big Data to Emulate a Target Trial. *AJE*.
- VanderWeele TJ, Ding P (2017). Sensitivity Analysis: Introducing the E-Value. *Ann Intern Med*.
- Fine JP, Gray RJ (1999). A Proportional Hazards Model for Competing Risks. *JASA*.

---

## Auteur

**Jean-Eude GBADA**
Master Économétrie & Data Science — Aix-Marseille School of Economics
Mastère IA & Big Data — ESGI Paris

📧 jeaneudes.gbada@gmail.com | [LinkedIn](#) | [GitHub](#)
