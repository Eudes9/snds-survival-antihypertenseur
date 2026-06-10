# =============================================================================
# FICHIER  : 03_propensity_score.R
# OBJECTIF : Calcul du score de propension et appariement
#            Méthode centrale en pharmaco-épidémiologie sur données
#            observationnelles pour contrôler le biais d'indication
#
# CONTEXTE :
#   Dans le SNDS, le choix de la classe thérapeutique par le médecin
#   n'est pas aléatoire. Il dépend des caractéristiques du patient
#   (âge, comorbidités, antécédents). C'est le "biais d'indication".
#   Le score de propension permet de rééquilibrer les groupes comme
#   dans un essai randomisé — c'est le principe de l'"essai émulé"
#   (Hernán & Robins, 2016).
#
# STRATÉGIE :
#   Comparaison IEC vs ARA2 (question clinique principale)
#   Méthode : appariement 1:1 par nearest neighbor sur le logit du PS
#   Évaluation : SMD avant/après appariement (seuil standard : <0.10)
#
# RÉFÉRENCES :
#   - Austin PC (2011). An Introduction to Propensity Score Methods.
#     Multivariate Behavioral Research, 46(3), 399-424.
#   - Hernán MA, Robins JM (2016). Using Big Data to Emulate a Target Trial.
#     American Journal of Epidemiology, 183(8), 758-764.
# =============================================================================

library(dplyr)
library(MatchIt)
library(tableone)
library(ggplot2)
library(tidyr)

cat("=== SCORE DE PROPENSION ET APPARIEMENT ===\n\n")

cohorte <- readRDS("data/cohorte_finale.rds")


# =============================================================================
# SOUS-COHORTE : IEC vs ARA2
# Comparaison principale — deux classes aux mécanismes différents
# =============================================================================

sous_cohorte <- cohorte %>%
  filter(classe_atc %in% c("IEC", "ARA2")) %>%
  mutate(
    # Variable de traitement binaire : 1 = IEC (exposé), 0 = ARA2 (référence)
    traitement_iec = as.integer(classe_atc == "IEC"),
    # Recodage pour MatchIt
    age_std        = scale(age_inclusion)[, 1],  # Standardisation de l'âge
    sexe_bin       = as.integer(ben_sex_cod == 1)  # 1 = Homme
  )

cat(sprintf("Sous-cohorte IEC vs ARA2 : %d patients\n", nrow(sous_cohorte)))
cat(sprintf("  IEC  : %d patients\n", sum(sous_cohorte$traitement_iec == 1)))
cat(sprintf("  ARA2 : %d patients\n\n", sum(sous_cohorte$traitement_iec == 0)))


# =============================================================================
# ÉTAPE 1 : TABLEAU 1 AVANT APPARIEMENT
# Standard de toute publication pharmaco-épidémiologique
# =============================================================================

cat("--- Tableau 1 : Caractéristiques avant appariement ---\n")

vars_tableau1 <- c("age_inclusion", "sexe_bin", "ald_diabete",
                   "ald_insuf_card", "ald_insuf_renale",
                   "score_charlson", "nb_hospitalisations", "ben_cmu_top")

vars_categoriques <- c("sexe_bin", "ald_diabete", "ald_insuf_card",
                       "ald_insuf_renale", "ben_cmu_top")

tab1_avant <- CreateTableOne(
  vars        = vars_tableau1,
  strata      = "classe_atc",
  data        = sous_cohorte,
  factorVars  = vars_categoriques,
  test        = TRUE,
  smd         = TRUE
)

print(tab1_avant, smd = TRUE, showAllLevels = FALSE,
      quote = FALSE, noSpaces = TRUE)


# =============================================================================
# ÉTAPE 2 : MODÈLE DE SCORE DE PROPENSION
# Régression logistique : P(IEC | covariables)
# =============================================================================

cat("\n--- Modèle logistique pour le score de propension ---\n")

# Formule incluant toutes les variables de confusion pré-spécifiées
formule_ps <- traitement_iec ~ age_inclusion + sexe_bin +
              ald_diabete + ald_insuf_card + ald_insuf_renale +
              score_charlson + nb_hospitalisations + ben_cmu_top

modele_ps <- glm(formule_ps,
                 data   = sous_cohorte,
                 family = binomial(link = "logit"))

cat("\nRésumé du modèle de score de propension :\n")
print(summary(modele_ps)$coefficients)

# Ajout du score de propension à la base
sous_cohorte <- sous_cohorte %>%
  mutate(
    ps_score = predict(modele_ps, type = "response"),
    ps_logit = qlogis(ps_score)  # Logit du PS (utilisé pour l'appariement)
  )

cat(sprintf("\nScore de propension — Distribution :\n"))
cat(sprintf("  IEC  : médiane = %.3f (IQR: %.3f - %.3f)\n",
            median(sous_cohorte$ps_score[sous_cohorte$traitement_iec == 1]),
            quantile(sous_cohorte$ps_score[sous_cohorte$traitement_iec == 1], 0.25),
            quantile(sous_cohorte$ps_score[sous_cohorte$traitement_iec == 1], 0.75)))
cat(sprintf("  ARA2 : médiane = %.3f (IQR: %.3f - %.3f)\n",
            median(sous_cohorte$ps_score[sous_cohorte$traitement_iec == 0]),
            quantile(sous_cohorte$ps_score[sous_cohorte$traitement_iec == 0], 0.25),
            quantile(sous_cohorte$ps_score[sous_cohorte$traitement_iec == 0], 0.75)))


# =============================================================================
# ÉTAPE 3 : VÉRIFICATION DU SUPPORT COMMUN (OVERLAP)
# Condition nécessaire pour la validité de l'appariement
# =============================================================================

p_overlap <- ggplot(sous_cohorte,
                    aes(x = ps_score,
                        fill = factor(traitement_iec,
                                      labels = c("ARA2 (référence)",
                                                 "IEC (exposé)")))) +
  geom_density(alpha = 0.55, color = "white") +
  scale_fill_manual(values = c("#2C7BB6", "#D7191C")) +
  geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  labs(
    title    = "Distribution du score de propension avant appariement",
    subtitle = "Un chevauchement satisfaisant est nécessaire pour la validité de l'analyse",
    x        = "Score de propension P(IEC | covariables)",
    y        = "Densité",
    fill     = "Groupe",
    caption  = "Lignes pointillées = zones de faible support commun (PS < 0.05 ou > 0.95)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("outputs/figures/ps_overlap_avant.png",
       plot = p_overlap, width = 8, height = 5, dpi = 150)


# =============================================================================
# ÉTAPE 4 : APPARIEMENT 1:1 PAR NEAREST NEIGHBOR
# =============================================================================

cat("\n--- Appariement 1:1 (nearest neighbor, caliper = 0.2 SD du logit PS) ---\n")

# Le caliper de 0.2 SD est le standard recommandé (Austin, 2011)
set.seed(20260101)

match_result <- matchit(
  formula  = formule_ps,
  data     = sous_cohorte,
  method   = "nearest",
  distance = "glm",
  link     = "logit",
  ratio    = 1,
  caliper  = 0.2,        # Caliper standard : 0.2 SD du logit PS
  std.caliper = TRUE,    # Caliper exprimé en SD
  replace  = FALSE       # Sans remise — évite la surreprésentation
)

cat("\nRésumé de l'appariement :\n")
print(summary(match_result, standardize = TRUE))

# Extraction de la base appariée
donnees_appariees <- match.data(match_result)

cat(sprintf("\nBase appariée : %d paires (N total = %d)\n",
            nrow(donnees_appariees) / 2,
            nrow(donnees_appariees)))


# =============================================================================
# ÉTAPE 5 : TABLEAU 1 APRÈS APPARIEMENT + VÉRIFICATION BALANCE
# Critère de succès : SMD < 0.10 pour toutes les covariables
# =============================================================================

cat("\n--- Tableau 1 : Caractéristiques après appariement ---\n")

tab1_apres <- CreateTableOne(
  vars        = vars_tableau1,
  strata      = "classe_atc",
  data        = donnees_appariees,
  factorVars  = vars_categoriques,
  test        = FALSE,   # Pas de test p après appariement (Love, 2002)
  smd         = TRUE
)

print(tab1_apres, smd = TRUE, showAllLevels = FALSE,
      quote = FALSE, noSpaces = TRUE)

# Extraction des SMD pour le love plot
smd_avant <- data.frame(
  variable  = vars_tableau1,
  smd       = as.numeric(ExtractSmd(tab1_avant)),
  periode   = "Avant appariement"
)

smd_apres <- data.frame(
  variable  = vars_tableau1,
  smd       = as.numeric(ExtractSmd(tab1_apres)),
  periode   = "Après appariement"
)

smd_comparison <- bind_rows(smd_avant, smd_apres) %>%
  filter(!is.na(smd))


# =============================================================================
# LOVE PLOT — Visualisation standard de la balance des covariables
# =============================================================================

p_love <- ggplot(smd_comparison,
                 aes(x = smd, y = reorder(variable, smd),
                     color = periode, shape = periode)) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             color = "#D7191C", linewidth = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = c("Avant appariement" = "#D7191C",
                                "Après appariement" = "#2C7BB6")) +
  scale_shape_manual(values = c("Avant appariement" = 16,
                                "Après appariement" = 17)) +
  labs(
    title    = "Love plot — Balance des covariables",
    subtitle = "Ligne rouge = seuil de déséquilibre résiduel (SMD = 0.10)",
    x        = "Différence standardisée moyenne (SMD)",
    y        = "Covariable",
    color    = NULL,
    shape    = NULL,
    caption  = "SMD < 0.10 indique un bon équilibre entre les groupes (Austin, 2011)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("outputs/figures/love_plot_balance.png",
       plot = p_love, width = 8, height = 6, dpi = 150)

cat("\n   Love plot exporté : outputs/figures/love_plot_balance.png\n")

# Vérification automatique du seuil SMD
smd_problematiques <- smd_apres %>%
  filter(smd > 0.10)

if (nrow(smd_problematiques) > 0) {
  warning(sprintf("ALERTE : Déséquilibre résiduel (SMD > 0.10) pour : %s",
                  paste(smd_problematiques$variable, collapse = ", ")))
  cat("   ATTENTION : Déséquilibre résiduel détecté — envisager un ajustement résiduel\n")
} else {
  cat("   OK : Toutes les covariables ont un SMD < 0.10 après appariement\n")
}


# =============================================================================
# SAUVEGARDE
# =============================================================================

saveRDS(donnees_appariees, "data/donnees_appariees.rds")
saveRDS(match_result,      "outputs/match_result.rds")
saveRDS(smd_comparison,    "outputs/smd_comparison.rds")

cat("\n>>> Données appariées sauvegardées : data/donnees_appariees.rds\n")
cat(">>> Score de propension terminé.\n")
