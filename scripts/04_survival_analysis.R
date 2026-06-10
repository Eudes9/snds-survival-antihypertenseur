# =============================================================================
# FICHIER  : 04_survival_analysis.R
# OBJECTIF : Analyse de survie complète
#            - Courbes de Kaplan-Meier
#            - Modèle de Cox (univarié et multivarié)
#            - Vérification de l'hypothèse des risques proportionnels
#            - Hazard ratios avec intervalles de confiance à 95%
#
# ÉVÉNEMENT : Abandon du traitement antihypertenseur
# CENSURE   : Fin de suivi à 365 jours ou fin de la période d'étude
#
# RÉFÉRENCES :
#   - Cox DR (1972). Regression Models and Life-Tables. JRSS-B, 34(2), 187-220.
#   - Kaplan EL, Meier P (1958). Nonparametric Estimation from Incomplete
#     Observations. JASA, 53(282), 457-481.
#   - Grambsch PM, Therneau TM (1994). Proportional Hazards Tests and
#     Diagnostics Based on Weighted Residuals. Biometrika, 81(3), 515-526.
# =============================================================================

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(broom)
library(tidyr)

cat("=== ANALYSE DE SURVIE ===\n\n")

cohorte     <- readRDS("data/cohorte_finale.rds")
donnees_app <- readRDS("data/donnees_appariees.rds")

dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# PARTIE A : KAPLAN-MEIER — COHORTE COMPLÈTE (6 classes)
# =============================================================================

cat("--- Kaplan-Meier : Toutes les classes thérapeutiques ---\n")

surv_obj_global <- Surv(time  = cohorte$delai_observe,
                        event = cohorte$abandon)

km_global <- survfit(surv_obj_global ~ classe_atc,
                     data      = cohorte,
                     conf.type = "log-log")

medianes <- surv_median(km_global)
cat("\nMédianes de survie par classe (jours) :\n")
print(medianes)

p_km_global <- ggsurvplot(
  km_global,
  data              = cohorte,
  fun               = "event",
  pval              = TRUE,
  conf.int          = TRUE,
  conf.int.alpha    = 0.12,
  risk.table        = TRUE,
  risk.table.height = 0.30,
  palette           = c("#2C7BB6","#D7191C","#1A9641",
                        "#FDAE61","#762A83","#F4A582"),
  xlab              = "Temps depuis l'initiation (jours)",
  ylab              = "Probabilité cumulée d'abandon",
  title             = "Courbes de Kaplan-Meier — Abandon du traitement antihypertenseur",
  xlim              = c(0, 60),
  break.time.by     = 20,
  risk.table.col    = "strata",
  ggtheme           = theme_minimal(base_size = 11)
)

png("outputs/figures/km_toutes_classes.png",
    width = 1200, height = 900, res = 120)
print(p_km_global)
dev.off()  # ← Fermeture obligatoire

cat("   KM global exporté : outputs/figures/km_toutes_classes.png\n")


# =============================================================================
# PARTIE B : KAPLAN-MEIER — IEC vs ARA2 (DONNÉES APPARIÉES)
# =============================================================================

cat("\n--- Kaplan-Meier : IEC vs ARA2 (données appariées) ---\n")

surv_app <- Surv(time  = donnees_app$delai_observe,
                 event = donnees_app$abandon)

km_app <- survfit(surv_app ~ classe_atc,
                  data      = donnees_app,
                  conf.type = "log-log")

logrank_test <- survdiff(surv_app ~ classe_atc, data = donnees_app)
p_logrank    <- 1 - pchisq(logrank_test$chisq,
                           df = length(logrank_test$n) - 1)

cat(sprintf("Test log-rank IEC vs ARA2 (données appariées) : p = %.4f\n", p_logrank))
cat(sprintf("Interprétation : %s\n",
            ifelse(p_logrank < 0.05,
                   "Différence significative de survie entre les deux classes",
                   "Pas de différence significative de survie")))

p_km_compare <- ggsurvplot(
  km_app,
  data              = donnees_app,
  fun               = "event",
  pval              = TRUE,
  conf.int          = TRUE,
  conf.int.alpha    = 0.15,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#D7191C", "#2C7BB6"),
  legend.labs       = c("ARA2", "IEC"),
  xlab              = "Temps depuis l'initiation (jours)",
  ylab              = "Probabilité cumulée d'abandon",
  title             = "IEC vs ARA2 — Abandon thérapeutique (données appariées)",
  subtitle          = sprintf("N = %d paires | Score de propension 1:1",
                              nrow(donnees_app) / 2),
  xlim              = c(0, 60),
  break.time.by     = 20,
  risk.table.col    = "strata",
  ggtheme           = theme_minimal(base_size = 11)
)

png("outputs/figures/km_iec_ara2_apparie.png",
    width = 900, height = 700, res = 120)
print(p_km_compare)
dev.off()  # ← Fermeture obligatoire

cat("   KM IEC vs ARA2 exporté : outputs/figures/km_iec_ara2_apparie.png\n")


# =============================================================================
# PARTIE C : MODÈLE DE COX UNIVARIÉ
# =============================================================================

cat("\n--- Modèle de Cox univarié ---\n")

vars_cox <- c("classe_atc", "age_inclusion", "ben_sex_cod",
              "ald_diabete", "ald_insuf_card", "score_charlson",
              "nb_hospitalisations", "ben_cmu_top")

cox_univarie <- lapply(vars_cox, function(var) {
  formule <- as.formula(paste("Surv(delai_observe, abandon) ~", var))
  modele  <- coxph(formule, data = cohorte)
  tidy(modele, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(variable = var)
}) %>%
  bind_rows() %>%
  select(variable, term, estimate, conf.low, conf.high, p.value) %>%
  rename(HR = estimate, IC_bas = conf.low, IC_haut = conf.high) %>%
  mutate(
    HR           = round(HR, 3),
    IC_bas       = round(IC_bas, 3),
    IC_haut      = round(IC_haut, 3),
    p.value      = round(p.value, 4),
    significatif = ifelse(p.value < 0.05, "***", "")
  )

cat("\nRésultats Cox univarié :\n")
print(cox_univarie)

write.csv(cox_univarie, "outputs/tables/cox_univarie.csv",
          row.names = FALSE)


# =============================================================================
# PARTIE D : MODÈLE DE COX MULTIVARIÉ
# =============================================================================

cat("\n--- Modèle de Cox multivarié ---\n")

cox_multi <- coxph(
  Surv(delai_observe, abandon) ~
    classe_atc + age_inclusion + factor(ben_sex_cod) +
    ald_diabete + ald_insuf_card + score_charlson +
    nb_hospitalisations + ben_cmu_top,
  data = cohorte,
  ties = "efron"
)

cat("\nRésumé du modèle de Cox multivarié :\n")
print(summary(cox_multi))

hr_table <- tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    Covariable = term,
    HR         = estimate,
    IC_2.5     = conf.low,
    IC_97.5    = conf.high,
    p_valeur   = p.value
  ) %>%
  mutate(
    `HR (IC 95%)` = sprintf("%.2f (%.2f–%.2f)", HR, IC_2.5, IC_97.5),
    p_valeur      = ifelse(p_valeur < 0.001, "<0.001",
                           as.character(round(p_valeur, 3)))
  )

cat("\nTableau publication (HR avec IC 95%) :\n")
print(hr_table %>% select(Covariable, `HR (IC 95%)`, p_valeur))

write.csv(hr_table, "outputs/tables/cox_multivariate_hr.csv",
          row.names = FALSE)


# =============================================================================
# PARTIE E : VÉRIFICATION HYPOTHÈSE DES RISQUES PROPORTIONNELS
# =============================================================================

cat("\n--- Test de l'hypothèse des risques proportionnels (Schoenfeld) ---\n")

test_schoenfeld <- cox.zph(cox_multi)
cat("\nTest de Schoenfeld :\n")
print(test_schoenfeld)

n_violations <- sum(test_schoenfeld$table[, "p"] < 0.05, na.rm = TRUE)

if (n_violations == 0) {
  cat("\n   OK : Hypothèse des risques proportionnels vérifiée\n")
} else {
  cat(sprintf("\n   ATTENTION : %d covariable(s) violent(s) l'hypothèse PH\n",
              n_violations))
}

png("outputs/figures/residus_schoenfeld.png",
    width = 1200, height = 800, res = 120)
plot(test_schoenfeld)
dev.off()  # ← Fermeture obligatoire
cat("   Résidus de Schoenfeld exportés\n")


# =============================================================================
# PARTIE F : FOREST PLOT DES HAZARD RATIOS
# =============================================================================

p_forest <- hr_table %>%
  filter(!grepl("Intercept", Covariable)) %>%
  mutate(
    Covariable = gsub("factor\\(ben_sex_cod\\)", "Sexe : ", Covariable),
    Covariable = gsub("classe_atc", "Classe : ", Covariable)
  ) %>%
  ggplot(aes(x = HR, y = reorder(Covariable, HR),
             xmin = IC_2.5, xmax = IC_97.5)) +
  geom_pointrange(color = "#2C7BB6", linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#D7191C", linewidth = 0.7) +
  scale_x_log10() +
  labs(
    title    = "Forest plot — Hazard ratios du modèle de Cox multivarié",
    subtitle = "Événement : abandon du traitement antihypertenseur",
    x        = "Hazard Ratio (échelle logarithmique) | Ligne rouge = HR nul",
    y        = NULL,
    caption  = "IC 95% | Méthode d'Efron | Données simulées de type SNDS"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/forest_plot_cox.png",
       plot = p_forest, width = 9, height = 7, dpi = 150)

cat("\n   Forest plot exporté : outputs/figures/forest_plot_cox.png\n")

cat("\n>>> Analyse de survie terminée.\n")
cat(">>> Résultats dans outputs/tables/ et outputs/figures/\n")