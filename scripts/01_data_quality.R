# =============================================================================
# FICHIER  : 01_data_quality.R
# OBJECTIF : Contrôle qualité rigoureux de la base analytique
#            Reproduit les vérifications standards appliquées au SNDS
#
# CHECKS RÉALISÉS :
#   1. Complétude des variables clés
#   2. Cohérence temporelle (dates)
#   3. Valeurs aberrantes (âge, délais)
#   4. Doublons
#   5. Distribution des variables — détection d'anomalies
#   6. Rapport qualité exporté
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

cat("=== CONTRÔLE QUALITÉ — BASE SNDS ANALYTIQUE ===\n\n")

# --- Chargement --------------------------------------------------------------
snds <- readRDS("data/snds_analytique.rds")
cat(sprintf("Base chargée : %d observations, %d variables\n\n",
            nrow(snds), ncol(snds)))

rapport_qualite <- list()


# =============================================================================
# CHECK 1 : COMPLÉTUDE DES VARIABLES CLÉS
# =============================================================================

cat("--- CHECK 1 : Complétude des variables ---\n")

vars_cles <- c("id_patient", "age_inclusion", "ben_sex_cod", "classe_atc",
               "delai_observe", "abandon", "date_index",
               "score_charlson", "ald_diabete")

completude <- snds %>%
  select(all_of(vars_cles)) %>%
  summarise(across(everything(),
                   list(
                     n_manquant  = ~sum(is.na(.)),
                     pct_manquant = ~round(100 * mean(is.na(.)), 2)
                   ))) %>%
  pivot_longer(everything(),
               names_to  = c("variable", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value)

print(completude)

# Alerte si > 5% manquants sur variable clé
vars_problematiques <- completude %>%
  filter(manquant > 5) %>%
  pull(variable)

if (length(vars_problematiques) > 0) {
  warning(sprintf("ALERTE : Variables avec >5%% manquants : %s",
                  paste(vars_problematiques, collapse = ", ")))
} else {
  cat("   OK : Aucune variable clé avec >5% de valeurs manquantes\n\n")
}

rapport_qualite[["completude"]] <- completude


# =============================================================================
# CHECK 2 : COHÉRENCE TEMPORELLE
# =============================================================================

cat("--- CHECK 2 : Cohérence temporelle ---\n")

DATE_DEBUT <- as.Date("2020-01-01")
DATE_FIN   <- as.Date("2023-12-31")

incoherences_dates <- snds %>%
  summarise(
    n_date_avant_debut  = sum(date_index < DATE_DEBUT, na.rm = TRUE),
    n_date_apres_fin    = sum(date_index > DATE_FIN, na.rm = TRUE),
    n_delai_negatif     = sum(delai_observe <= 0, na.rm = TRUE),
    n_delai_sup_suivi   = sum(delai_observe > 366, na.rm = TRUE)
  )

print(incoherences_dates)

if (any(incoherences_dates > 0)) {
  warning("ALERTE : Incohérences temporelles détectées — vérifier la simulation")
} else {
  cat("   OK : Aucune incohérence temporelle\n\n")
}

rapport_qualite[["coherence_dates"]] <- incoherences_dates


# =============================================================================
# CHECK 3 : VALEURS ABERRANTES
# =============================================================================

cat("--- CHECK 3 : Valeurs aberrantes ---\n")

# Âge
age_stats <- snds %>%
  summarise(
    age_min    = min(age_inclusion),
    age_max    = max(age_inclusion),
    age_mean   = round(mean(age_inclusion), 1),
    age_median = median(age_inclusion),
    n_age_hors_cible = sum(age_inclusion < 40 | age_inclusion > 80)
  )

cat("Statistiques d'âge :\n")
print(age_stats)

# Délai de survie
delai_stats <- snds %>%
  summarise(
    delai_min    = min(delai_observe),
    delai_max    = max(delai_observe),
    delai_mean   = round(mean(delai_observe), 1),
    delai_median = median(delai_observe),
    n_delai_1j   = sum(delai_observe <= 1)  # Délais très courts suspects
  )

cat("\nStatistiques de délai (jours) :\n")
print(delai_stats)

rapport_qualite[["aberrants"]] <- list(age = age_stats, delai = delai_stats)


# =============================================================================
# CHECK 4 : DOUBLONS
# =============================================================================

cat("\n--- CHECK 4 : Doublons ---\n")

n_doublons <- snds %>%
  group_by(id_patient) %>%
  filter(n() > 1) %>%
  nrow()

cat(sprintf("   Nombre de doublons sur id_patient : %d\n", n_doublons))

if (n_doublons > 0) {
  warning("ALERTE : Doublons détectés — dédoublonnage nécessaire")
} else {
  cat("   OK : Aucun doublon\n\n")
}

rapport_qualite[["doublons"]] <- n_doublons


# =============================================================================
# CHECK 5 : DISTRIBUTIONS ET DÉTECTION D'ANOMALIES
# =============================================================================

cat("--- CHECK 5 : Distributions ---\n")

# Distribution des classes thérapeutiques
dist_classes <- snds %>%
  count(classe_atc) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))

cat("\nDistribution classes thérapeutiques :\n")
print(dist_classes)

# Distribution de l'abandon par classe
abandon_classe <- snds %>%
  group_by(classe_atc) %>%
  summarise(
    n            = n(),
    n_abandons   = sum(abandon),
    taux_abandon = round(100 * mean(abandon), 1),
    delai_median = median(delai_observe),
    .groups = "drop"
  ) %>%
  arrange(desc(taux_abandon))

cat("\nTaux d'abandon par classe :\n")
print(abandon_classe)

# Alerte si une classe a un taux d'abandon > 3x la moyenne
taux_moyen <- mean(snds$abandon)
classes_anormales <- abandon_classe %>%
  filter(taux_abandon / 100 > 3 * taux_moyen) %>%
  pull(classe_atc)

if (length(classes_anormales) > 0) {
  warning(sprintf("ALERTE : Taux d'abandon anormalement élevé pour : %s",
                  paste(classes_anormales, collapse = ", ")))
}

rapport_qualite[["distributions"]] <- list(
  classes  = dist_classes,
  abandon  = abandon_classe
)


# =============================================================================
# CHECK 6 : GRAPHIQUES DE VALIDATION
# =============================================================================

cat("\n--- CHECK 6 : Export des graphiques de validation ---\n")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# Distribution de l'âge
p_age <- ggplot(snds, aes(x = age_inclusion)) +
  geom_histogram(bins = 20, fill = "#2C7BB6", color = "white", alpha = 0.85) +
  geom_vline(xintercept = c(40, 80), linetype = "dashed",
             color = "red", linewidth = 0.7) +
  labs(
    title    = "Distribution de l'âge à l'inclusion",
    subtitle = sprintf("N = %d patients | Lignes rouges = limites d'inclusion (40-80 ans)",
                       nrow(snds)),
    x        = "Âge (années)",
    y        = "Nombre de patients",
    caption  = "Données simulées — Projet pharmaco-épidémiologie SNDS"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/qc_distribution_age.png",
       plot = p_age, width = 8, height = 5, dpi = 150)

# Distribution des délais
p_delai <- ggplot(snds, aes(x = delai_observe,
                             fill = factor(abandon,
                                           labels = c("Censuré", "Abandon")))) +
  geom_histogram(bins = 30, alpha = 0.75, position = "stack") +
  scale_fill_manual(values = c("#2C7BB6", "#D7191C")) +
  labs(
    title    = "Distribution des délais de suivi",
    subtitle = "Bleu = censuré (fin de suivi) | Rouge = abandon (événement)",
    x        = "Délai (jours)",
    y        = "Nombre de patients",
    fill     = "Statut",
    caption  = "Données simulées"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave("outputs/figures/qc_distribution_delais.png",
       plot = p_delai, width = 8, height = 5, dpi = 150)

cat("   Graphiques exportés dans outputs/figures/\n")


# =============================================================================
# RAPPORT DE SYNTHÈSE QUALITÉ
# =============================================================================

cat("\n=== SYNTHÈSE DU CONTRÔLE QUALITÉ ===\n")
cat(sprintf("  Observations analysées    : %d\n", nrow(snds)))
cat(sprintf("  Variables clés complètes  : %s\n",
            ifelse(length(vars_problematiques) == 0, "OUI", "NON")))
cat(sprintf("  Incohérences temporelles  : %d\n",
            sum(incoherences_dates)))
cat(sprintf("  Doublons                  : %d\n", n_doublons))
cat(sprintf("  Taux d'abandon global     : %.1f%%\n",
            100 * mean(snds$abandon)))
cat(sprintf("  Délai médian de suivi     : %.0f jours\n",
            median(snds$delai_observe)))
cat("=====================================\n\n")

saveRDS(rapport_qualite, "outputs/rapport_qualite.rds")
cat(">>> Rapport qualité sauvegardé : outputs/rapport_qualite.rds\n")
cat(">>> Contrôle qualité terminé.\n")
