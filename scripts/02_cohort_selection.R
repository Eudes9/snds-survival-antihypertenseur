# =============================================================================
# FICHIER  : 02_cohort_selection.R
# OBJECTIF : Construction de la cohorte analytique
#            Reproduit le processus de sélection appliqué dans les études
#            pharmaco-épidémiologiques sur le SNDS (EPI-PHARE, CNAM)
#
# CRITÈRES D'INCLUSION :
#   - Primo-initiation d'un antihypertenseur entre 01/01/2020 et 31/12/2022
#   - Âge entre 40 et 80 ans à la date d'inclusion
#   - Au moins 6 mois de données disponibles avant l'inclusion (wash-out)
#
# CRITÈRES D'EXCLUSION :
#   - Traitement antihypertenseur dans les 6 mois précédant l'inclusion
#   - Insuffisance rénale terminale (contre-indication à certaines classes)
#   - Données manquantes sur les variables essentielles
#
# OUTPUT :
#   - Cohorte finale analysable
#   - Diagramme de flux (flow chart) — standard des publications
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tibble)

cat("=== CONSTRUCTION DE LA COHORTE ===\n\n")

snds_brut <- readRDS("data/snds_analytique.rds")
cat(sprintf("Base brute : %d patients\n\n", nrow(snds_brut)))

# Suivi du flux de patients à chaque étape
flux <- tibble(
  etape       = character(),
  n_patients  = integer(),
  n_exclus    = integer(),
  raison      = character()
)

enregistrer_flux <- function(flux, etape, n, n_avant, raison = "") {
  bind_rows(flux, tibble(
    etape      = etape,
    n_patients = n,
    n_exclus   = n_avant - n,
    raison     = raison
  ))
}

# --- ÉTAPE 0 : Population initiale -------------------------------------------
n_initial <- nrow(snds_brut)
flux <- enregistrer_flux(flux, "Population initiale", n_initial, n_initial,
                         "Ensemble des patients avec ≥1 délivrance d'antihypertenseur")

# --- ÉTAPE 1 : Critère d'âge -------------------------------------------------
cohorte <- snds_brut %>%
  filter(between(age_inclusion, 40, 80))

flux <- enregistrer_flux(flux, "Critère d'âge", nrow(cohorte), n_initial,
                         "Exclusion : âge < 40 ans ou > 80 ans")

cat(sprintf("Après critère d'âge (40-80 ans)    : %d patients (--%d)\n",
            nrow(cohorte), n_initial - nrow(cohorte)))

# --- ÉTAPE 2 : Exclusion insuffisance rénale terminale -----------------------
# Dans le vrai SNDS : repérage via codes CIM-10 N18.5 ou actes de dialyse
# Ici : on utilise le proxy ald_insuf_renale avec score élevé
cohorte_2 <- cohorte %>%
  filter(!(ald_insuf_renale == 1 & score_charlson >= 4))

flux <- enregistrer_flux(flux, "Exclusion IRC terminale", nrow(cohorte_2),
                         nrow(cohorte),
                         "Insuffisance rénale sévère (score Charlson ≥4 + IRC)")

cat(sprintf("Après exclusion IRC terminale       : %d patients (--%d)\n",
            nrow(cohorte_2), nrow(cohorte) - nrow(cohorte_2)))

# --- ÉTAPE 3 : Données complètes sur variables essentielles ------------------
vars_essentielles <- c("age_inclusion", "ben_sex_cod", "classe_atc",
                       "delai_observe", "abandon", "score_charlson")

cohorte_3 <- cohorte_2 %>%
  filter(complete.cases(select(., all_of(vars_essentielles))))

flux <- enregistrer_flux(flux, "Données complètes", nrow(cohorte_3),
                         nrow(cohorte_2),
                         "Exclusion : données manquantes sur variables essentielles")

cat(sprintf("Après exclusion données manquantes  : %d patients (--%d)\n",
            nrow(cohorte_3), nrow(cohorte_2) - nrow(cohorte_3)))

# --- ÉTAPE 4 : Délai de suivi minimum ----------------------------------------
# On exige au moins 30 jours de suivi pour inclure un patient
SUIVI_MIN <- 28

cohorte_finale <- cohorte_3 %>%
  filter(delai_observe >= SUIVI_MIN)

flux <- enregistrer_flux(flux, "Cohorte finale", nrow(cohorte_finale),
                         nrow(cohorte_3),
                         sprintf("Exclusion : suivi < %d jours", SUIVI_MIN))

cat(sprintf("Cohorte finale analysable           : %d patients (--%d)\n\n",
            nrow(cohorte_finale), nrow(cohorte_3) - nrow(cohorte_finale)))


# =============================================================================
# DIAGRAMME DE FLUX (FLOW CHART)
# Standard indispensable dans toute publication pharmaco-épidémiologique
# =============================================================================

cat("--- Génération du diagramme de flux ---\n")

flux_plot <- flux %>%
  mutate(etape_num = row_number(),
         label_boite = sprintf("%s\nn = %s",
                                etape,
                                format(n_patients, big.mark = " ")),
         label_exclus = if_else(
           n_exclus > 0,
           sprintf("Exclus : %s\n(%s)", raison,
                   format(n_exclus, big.mark = " ")),
           NA_character_
         ))

# Visualisation du flow chart avec ggplot2
p_flux <- ggplot() +
  # Boîtes principales
  geom_rect(data = flux_plot,
            aes(xmin = 0.1, xmax = 0.9,
                ymin = -etape_num - 0.35,
                ymax = -etape_num + 0.35),
            fill = "#EBF5FB", color = "#2C7BB6", linewidth = 0.8) +
  # Texte dans les boîtes
  geom_text(data = flux_plot,
            aes(x = 0.5, y = -etape_num, label = label_boite),
            size = 3.2, fontface = "bold", hjust = 0.5) +
  # Flèches entre les boîtes
  geom_segment(data = flux_plot %>% filter(etape_num < max(etape_num)),
               aes(x = 0.5, xend = 0.5,
                   y = -etape_num - 0.35,
                   yend = -etape_num - 0.65),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "#2C7BB6", linewidth = 0.7) +
  # Annotations d'exclusion
  geom_text(data = flux_plot %>% filter(!is.na(label_exclus)),
            aes(x = 1.05, y = -etape_num, label = label_exclus),
            size = 2.8, hjust = 0, color = "#D7191C") +
  xlim(0, 1.8) +
  theme_void() +
  labs(
    title   = "Diagramme de flux — Construction de la cohorte",
    caption = "Données simulées de type SNDS | Projet pharmaco-épidémiologie"
  ) +
  theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.caption = element_text(size = 8, color = "grey50"))

ggsave("outputs/figures/flow_chart_cohorte.png",
       plot = p_flux, width = 9, height = 7, dpi = 150)

cat("   Flow chart exporté : outputs/figures/flow_chart_cohorte.png\n\n")


# =============================================================================
# DESCRIPTION DE LA COHORTE FINALE
# =============================================================================

cat("=== DESCRIPTION DE LA COHORTE FINALE ===\n")
cat(sprintf("  N total                  : %d\n", nrow(cohorte_finale)))
cat(sprintf("  Événements (abandons)    : %d (%.1f%%)\n",
            sum(cohorte_finale$abandon),
            100 * mean(cohorte_finale$abandon)))
cat(sprintf("  Âge moyen (±SD)          : %.1f ± %.1f ans\n",
            mean(cohorte_finale$age_inclusion),
            sd(cohorte_finale$age_inclusion)))
cat(sprintf("  Proportion femmes        : %.1f%%\n",
            100 * mean(cohorte_finale$ben_sex_cod == 2)))
cat(sprintf("  Délai médian suivi       : %.0f jours\n",
            median(cohorte_finale$delai_observe)))
cat("=========================================\n\n")

# Sauvegarde
saveRDS(cohorte_finale, "data/cohorte_finale.rds")
saveRDS(flux, "outputs/flux_patients.rds")

cat(">>> Cohorte finale sauvegardée : data/cohorte_finale.rds\n")
cat(">>> Sélection de la cohorte terminée.\n")
