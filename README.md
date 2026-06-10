\# Étude pharmaco-épidémiologique — Abandon du traitement antihypertenseur



\## Contexte

Étude de survie sur données simulées de type SNDS (Système National des

Données de Santé), portant sur l'abandon du traitement antihypertenseur

dans les 12 mois suivant l'initiation.



Reproduit les standards méthodologiques du GIS EPI-PHARE (ANSM/CNAM).



> Les données sont entièrement simulées à des fins pédagogiques.



\## Structure du projet

snds-survival-antihypertenseur/

├── data/

│   └── 00\_simulate\_snds.R

├── scripts/

│   ├── 01\_data\_quality.R

│   ├── 02\_cohort\_selection.R

│   ├── 03\_propensity\_score.R

│   ├── 04\_survival\_analysis.R

│   └── 05\_sensitivity\_analysis.R

├── report/

│   └── rapport\_final.Rmd

└── shiny/

└── app.R



\## Méthodologie

\- Score de propension (MatchIt) — appariement 1:1

\- Courbes de Kaplan-Meier + test log-rank

\- Modèle de Cox multivarié + test de Schoenfeld

\- 5 analyses de sensibilité (E-value, risques compétitifs, sous-groupes)



\## Reproductibilité

```r

source("data/00\_simulate\_snds.R")

source("scripts/01\_data\_quality.R")

source("scripts/02\_cohort\_selection.R")

source("scripts/03\_propensity\_score.R")

source("scripts/04\_survival\_analysis.R")

source("scripts/05\_sensitivity\_analysis.R")

rmarkdown::render("report/rapport\_final.Rmd")

shiny::runApp("shiny/app.R")

```



\## Auteur

\*\*Jean-Eude GBADA\*\*

Master Économétrie \& Data Science — Aix-Marseille School of Economics

Mastère IA \& Big Data — ESGI Paris

