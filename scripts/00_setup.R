
# =============================================================================
# FICHIER  : 00_setup.R
# kableExtra doit être chargé EN DERNIER car il importe stringr
# =============================================================================

packages <- c(
  "stringr",      # chargé en premier
  "dplyr",
  "lubridate",
  "ggplot2",
  "tidyr",
  "survival",
  "survminer",
  "MatchIt",
  "tableone",
  "broom",
  "purrr",
  "rmarkdown",
  "knitr",
  "kableExtra"    # chargé en dernier
)

packages_manquants <- packages[!packages %in% installed.packages()[,"Package"]]
if (length(packages_manquants) > 0) {
  install.packages(packages_manquants)
}

suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

cat("Tous les packages sont chargés ✓\n")