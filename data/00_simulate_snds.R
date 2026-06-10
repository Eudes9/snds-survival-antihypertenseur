# =============================================================================
# PROJET : Étude de survie pharmaco-épidémiologique
#          Abandon du traitement antihypertenseur — Analyse de type SNDS
# =============================================================================
# FICHIER  : 00_simulate_snds.R
# AUTEUR   : Jean-Eude GBADA
# DATE     : 2026
# OBJECTIF : Simuler une base médico-administrative de type SNDS
#            reproduisant la structure des tables DCIR et PMSI
#            (remboursements médicaments + hospitalisations)
#
# STRUCTURE SNDS SIMULÉE :
#   - Table bénéficiaires  (IR_BEN_R)   : caractéristiques patients
#   - Table délivrances    (ER_PHA_F)   : remboursements médicaments
#   - Table hospitalisations (MCO_B)    : séjours hospitaliers
#   - Table ALD                         : affections longue durée
#
# NOTE : Ces données sont entièrement simulées à des fins pédagogiques.
#        Elles ne contiennent aucune donnée réelle de patient.
# =============================================================================

# --- 0. Environnement --------------------------------------------------------

set.seed(20260101)  # Reproductibilité totale

library(dplyr)
library(lubridate)
library(stringr)

N_PATIENTS <- 5000          # Taille de cohorte réaliste pour une étude pilote
DATE_DEBUT  <- as.Date("2020-01-01")   # Début de la fenêtre d'étude
DATE_FIN    <- as.Date("2023-12-31")   # Fin de la fenêtre d'étude
SUIVI_MAX   <- 365                     # Suivi maximum en jours (1 an)


# =============================================================================
# TABLE 1 : IR_BEN_R — Référentiel bénéficiaires
# Contient les caractéristiques socio-démographiques des assurés
# =============================================================================

cat(">>> Simulation table bénéficiaires (IR_BEN_R)...\n")

ir_ben_r <- tibble(
  # Identifiant pseudonymisé (dans le vrai SNDS : NIR pseudonymisé)
  id_patient = str_pad(1:N_PATIENTS, width = 10, pad = "0"),

  # Caractéristiques démographiques
  ben_nai_ann  = sample(1940:1983, N_PATIENTS, replace = TRUE),  # Année de naissance
  ben_sex_cod  = sample(c(1, 2), N_PATIENTS, replace = TRUE,     # 1=Homme, 2=Femme
                        prob = c(0.48, 0.52)),

  # Région de résidence (codes INSEE simplifiés)
  ben_reg_cod  = sample(c("11","24","27","28","32","44","52","53",
                           "75","76","84","93","94"), N_PATIENTS, replace = TRUE),

  # Régime d'assurance maladie
  org_cle_new  = sample(c("01","02","03","09"), N_PATIENTS, replace = TRUE,
                         prob = c(0.75, 0.12, 0.08, 0.05)),

  # Indicateur CMU-C / complémentaire santé solidaire
  ben_cmu_top  = rbinom(N_PATIENTS, 1, prob = 0.07),

  # Date d'entrée dans le référentiel
  ben_dcd_dte  = as.Date(NA)  # Date de décès (NA = vivant)
) %>%
  mutate(
    # Âge à la date index (début de suivi)
    age_inclusion = 2020 - ben_nai_ann,

    # Décès simulés (5% de mortalité sur la période — réaliste HTA)
    ben_dcd_dte = if_else(
      rbinom(N_PATIENTS, 1, prob = 0.05) == 1,
      DATE_DEBUT + sample(1:1460, N_PATIENTS, replace = TRUE),
      as.Date(NA)
    ),

    # On restreint la cohorte aux 40-80 ans (cible HTA)
    eligible_age = between(age_inclusion, 40, 80)
  )

cat(sprintf("   %d bénéficiaires simulés dont %d éligibles (40-80 ans)\n",
            N_PATIENTS, sum(ir_ben_r$eligible_age)))


# =============================================================================
# TABLE 2 : ALD — Affections Longue Durée
# Dans le SNDS, l'ALD est un proxy important des comorbidités chroniques
# =============================================================================

cat(">>> Simulation table ALD...\n")

# Codes ALD pertinents pour l'HTA et ses comorbidités
# ALD 1  = Accident vasculaire cérébral
# ALD 3  = Artériopathies chroniques
# ALD 5  = Insuffisance cardiaque
# ALD 8  = Diabète (très fréquent chez hypertendu)
# ALD 12 = HTA sévère (ALD supprimée en 2011 mais proxy utile)
# ALD 19 = Insuffisance rénale chronique
# ALD 30 = Autres maladies graves

ald <- ir_ben_r %>%
  select(id_patient) %>%
  mutate(
    ald_diabete      = rbinom(N_PATIENTS, 1, prob = 0.22),  # 22% — prévalence réaliste
    ald_insuf_card   = rbinom(N_PATIENTS, 1, prob = 0.08),
    ald_insuf_renale = rbinom(N_PATIENTS, 1, prob = 0.06),
    ald_avc          = rbinom(N_PATIENTS, 1, prob = 0.04),
    ald_arterio      = rbinom(N_PATIENTS, 1, prob = 0.05),

    # Score de comorbidité agrégé (inspiré Charlson simplifié)
    score_charlson = ald_diabete + ald_insuf_card * 2 +
                     ald_insuf_renale * 2 + ald_avc * 2 + ald_arterio,
    charlson_cat   = case_when(
      score_charlson == 0 ~ "0",
      score_charlson == 1 ~ "1",
      score_charlson >= 2 ~ "2+"
    )
  )

cat(sprintf("   Distribution score Charlson : 0=%d | 1=%d | 2+=%d\n",
            sum(ald$score_charlson == 0),
            sum(ald$score_charlson == 1),
            sum(ald$score_charlson >= 2)))


# =============================================================================
# TABLE 3 : ER_PHA_F — Délivrances de médicaments (DCIR)
# Table centrale du SNDS : un enregistrement par délivrance en pharmacie
# Variables clés : code CIP13, date délivrance, nombre boîtes, montant remboursé
# =============================================================================

cat(">>> Simulation table délivrances médicaments (ER_PHA_F)...\n")

# Classes d'antihypertenseurs selon classification ATC
# C02 = Antihypertenseurs
# C03 = Diurétiques
# C07 = Bêtabloquants
# C08 = Inhibiteurs calciques
# C09 = Agents agissant sur le système rénine-angiotensine (IEC/ARA2)

classes_atc <- c("IEC", "ARA2", "BB", "ICC", "DIURETIQUE", "ASSOCIATION")

# Probabilités d'exposition par classe (réalistes France 2020)
prob_classes <- c(0.28, 0.25, 0.18, 0.15, 0.09, 0.05)

# Codes CIP13 simulés par classe (dans le vrai SNDS : vrais codes CIP)
cip_par_classe <- list(
  IEC         = c("3400936741234", "3400937852341", "3400938963452"),
  ARA2        = c("3400939074563", "3400940185674", "3400941296785"),
  BB          = c("3400942307896", "3400943418907", "3400944529018"),
  ICC         = c("3400945630129", "3400946741230", "3400947852341"),
  DIURETIQUE  = c("3400948963452", "3400949074563"),
  ASSOCIATION = c("3400950185674", "3400951296785")
)

# Génération des ordonnances initiales (première délivrance = inclusion)
er_pha_f_init <- ir_ben_r %>%
  filter(eligible_age) %>%
  mutate(
    # Classe thérapeutique assignée à l'inclusion
    classe_atc     = sample(classes_atc, sum(eligible_age),
                            replace = TRUE, prob = prob_classes),
    # Date de première délivrance (date index)
    date_index     = DATE_DEBUT + sample(0:1095, sum(eligible_age), replace = TRUE),
    # Conditionnement : 28 ou 30 comprimés (le plus fréquent)
    nb_boites      = 1,
    duree_boite    = sample(c(28, 30), sum(eligible_age), replace = TRUE),
    # Montant remboursé simulé (en euros)
    montant_rembourse = round(runif(sum(eligible_age), 2.5, 18.5), 2),
    # Taux de remboursement
    taux_remb      = sample(c(65, 100), sum(eligible_age),
                            replace = TRUE, prob = c(0.70, 0.30))
  )

# Simulation des renouvellements (jusqu'à l'abandon ou fin de suivi)
# Logique : chaque patient renouvelle jusqu'à ce qu'il arrête
# L'écart entre deux délivrances > 1.5x la durée de la boîte = abandon

generer_renouvellements <- function(df) {
  renouvellements <- list()
  
  for (i in seq_len(nrow(df))) {
    patient       <- df[i, ]
    date_courante <- patient$date_index
    historique    <- list(patient)
    
    # Probabilité d'abandon à chaque renouvellement selon la classe
    prob_abandon_base <- case_when(
      patient$classe_atc == "IEC"         ~ 0.06,
      patient$classe_atc == "ARA2"        ~ 0.06,
      patient$classe_atc == "BB"          ~ 0.10,
      patient$classe_atc == "ICC"         ~ 0.08,
      patient$classe_atc == "DIURETIQUE"  ~ 0.12,
      patient$classe_atc == "ASSOCIATION" ~ 0.05,
      TRUE                                ~ 0.08
    )
    
    for (j in 1:12) {
      # Abandon plus probable au 1er renouvellement
      prob_j <- if (j == 1) prob_abandon_base * 1.5 else prob_abandon_base
      
      # Tirage : abandon ou renouvellement ?
      if (runif(1) < prob_j) break
      
      # Délai jusqu'au prochain renouvellement
      # Variabilité réaliste : boîte + délai administratif (0-15j)
      delai_renouvellement <- patient$duree_boite + sample(0:15, 1)
      date_courante        <- date_courante + delai_renouvellement
      
      # Censure si hors fenêtre
      if (date_courante > patient$date_index + SUIVI_MAX) break
      if (date_courante > DATE_FIN) break
      
      nouveau            <- patient
      nouveau$date_index <- date_courante
      historique         <- c(historique, list(nouveau))
    }
    
    renouvellements[[i]] <- bind_rows(historique)
  }
  
  bind_rows(renouvellements)
}

cat("   Génération des renouvellements (patience ~30s)...\n")
er_pha_f <- generer_renouvellements(er_pha_f_init)

cat(sprintf("   %d délivrances générées pour %d patients\n",
            nrow(er_pha_f), nrow(er_pha_f_init)))


# =============================================================================
# TABLE 4 : MCO_B — Hospitalisations (PMSI MCO)
# Résumés de sortie anonymisés : diagnostics CIM-10, durée séjour, mode sortie
# =============================================================================

cat(">>> Simulation table hospitalisations (MCO_B)...\n")

# Codes CIM-10 pertinents pour l'HTA et ses complications
# I10   = HTA essentielle
# I11   = HTA avec atteinte cardiaque
# I20   = Angine de poitrine
# I21   = Infarctus aigu du myocarde
# I50   = Insuffisance cardiaque
# I63   = Infarctus cérébral
# N18   = Insuffisance rénale chronique
# E11   = Diabète type 2

# 20% des patients ont au moins une hospitalisation sur la période
patients_hospit <- sample(er_pha_f_init$id_patient,
                          size = round(nrow(er_pha_f_init) * 0.20))

mco_b <- tibble(
  id_patient     = patients_hospit,
  date_entree    = DATE_DEBUT + sample(1:1460, length(patients_hospit), replace = TRUE),
  duree_sejour   = sample(1:21, length(patients_hospit), replace = TRUE,
                           prob = c(rep(0.08, 3), rep(0.06, 4),
                                    rep(0.04, 7), rep(0.02, 7))),
  dp_cim10       = sample(c("I10","I11","I20","I21","I50","I63","N18","E11"),
                           length(patients_hospit), replace = TRUE,
                           prob = c(0.25,0.15,0.12,0.10,0.12,0.08,0.10,0.08)),
  mode_sortie    = sample(c("domicile","transfert","deces","HAD"),
                           length(patients_hospit), replace = TRUE,
                           prob = c(0.82, 0.10, 0.04, 0.04)),
  ghm            = paste0("0", sample(5:9, length(patients_hospit), replace = TRUE),
                           "M", sample(10:99, length(patients_hospit), replace = TRUE))
) %>%
  mutate(date_sortie = date_entree + duree_sejour)

cat(sprintf("   %d hospitalisations simulées pour %d patients distincts\n",
            nrow(mco_b), length(unique(mco_b$id_patient))))


# =============================================================================
# CONSTRUCTION DE LA TABLE ANALYTIQUE FINALE
# Jointure de toutes les tables — 1 ligne par patient
# =============================================================================

cat(">>> Construction de la table analytique finale...\n")

# Calcul du délai jusqu'à l'abandon pour chaque patient
survie_par_patient <- er_pha_f %>%
  group_by(id_patient) %>%
  summarise(
    date_index          = min(date_index),
    derniere_delivrance = max(date_index),
    nb_delivrances      = n(),
    classe_atc          = first(classe_atc),
    duree_boite         = first(duree_boite),
    .groups = "drop"
  ) %>%
  mutate(
    # Délai = écart entre première et dernière délivrance + durée boîte
    # Si 1 seule délivrance : délai = durée boîte seulement
    # Si plusieurs : délai couvre toute la période de traitement
    delai_observe = case_when(
      nb_delivrances == 1 ~ as.numeric(duree_boite),
      TRUE ~ as.numeric(derniere_delivrance - date_index) + duree_boite
    ),
    delai_observe = pmin(delai_observe, SUIVI_MAX),
    # Abandon = moins de 3 délivrances ET pas arrivé en fin de suivi
    abandon = as.integer(nb_delivrances <= 2 & delai_observe < 300)
  )# Calcul correct : on recalcule depuis er_pha_f directement
dates_par_patient <- er_pha_f %>%
  group_by(id_patient) %>%
  summarise(
    date_premiere   = min(date_index),
    date_derniere   = max(date_index),
    nb_delivrances  = n(),
    classe_atc      = first(classe_atc),
    duree_boite     = first(duree_boite),
    .groups = "drop"
  )

survie_par_patient <- dates_par_patient %>%
  mutate(
    delai_observe = as.integer(date_derniere - date_premiere) + duree_boite,
    delai_observe = pmin(delai_observe, SUIVI_MAX),
    abandon       = as.integer(nb_delivrances <= 2 & delai_observe < 300)
  ) %>%
  rename(date_index = date_premiere)

# Table analytique finale
snds_analytique <- survie_par_patient %>%
  left_join(ir_ben_r %>% select(id_patient, age_inclusion, ben_sex_cod,
                                 ben_reg_cod, ben_cmu_top, ben_dcd_dte),
            by = "id_patient") %>%
  left_join(ald %>% select(id_patient, ald_diabete, ald_insuf_card,
                            ald_insuf_renale, score_charlson, charlson_cat),
            by = "id_patient") %>%
  left_join(
    mco_b %>%
      group_by(id_patient) %>%
      summarise(nb_hospitalisations = n(),
                hospit_cardiovasc   = as.integer(any(dp_cim10 %in%
                                       c("I20","I21","I50","I63"))),
                .groups = "drop"),
    by = "id_patient"
  ) %>%
  mutate(
    # Variables dérivées utiles pour l'analyse
    sexe_label     = if_else(ben_sex_cod == 1, "Homme", "Femme"),
    age_cat        = cut(age_inclusion,
                         breaks = c(39, 49, 59, 69, 80),
                         labels = c("40-49", "50-59", "60-69", "70-80")),
    cmu_label      = if_else(ben_cmu_top == 1, "Oui", "Non"),
    nb_hospit_cat  = case_when(
      is.na(nb_hospitalisations) ~ "0",
      nb_hospitalisations == 1   ~ "1",
      nb_hospitalisations >= 2   ~ "2+"
    ),
    nb_hospitalisations   = replace_na(nb_hospitalisations, 0),
    hospit_cardiovasc     = replace_na(hospit_cardiovasc, 0),

    # Troncature du délai à la date de décès si applicable
    delai_observe  = if_else(
      !is.na(ben_dcd_dte),
      pmin(delai_observe,
           as.numeric(ben_dcd_dte - date_index)),
      delai_observe
    ),
    # Événement composite : abandon OU décès
    evenement_composite = as.integer(abandon == 1 |
                                     (!is.na(ben_dcd_dte) &
                                      as.numeric(ben_dcd_dte - date_index) <= SUIVI_MAX))
  ) %>%
  # Exclusions finales
  filter(
    delai_observe > 0,          # Exclure délais nuls (erreurs de codage)
    age_inclusion >= 40,         # Critère d'inclusion : 40 ans minimum
    age_inclusion <= 80          # Critère d'inclusion : 80 ans maximum
  )

cat(sprintf("   Table analytique finale : %d patients, %d variables\n",
            nrow(snds_analytique), ncol(snds_analytique)))


# =============================================================================
# STATISTIQUES DE VALIDATION DE LA SIMULATION
# =============================================================================

cat("\n=== RAPPORT DE VALIDATION DE LA SIMULATION ===\n")
cat(sprintf("Patients inclus            : %d\n", nrow(snds_analytique)))
cat(sprintf("Événements (abandons)      : %d (%.1f%%)\n",
            sum(snds_analytique$abandon),
            100 * mean(snds_analytique$abandon)))
cat(sprintf("Délai médian (jours)       : %.0f\n",
            median(snds_analytique$delai_observe)))
cat(sprintf("Proportion femmes          : %.1f%%\n",
            100 * mean(snds_analytique$ben_sex_cod == 2)))
cat(sprintf("Âge moyen                  : %.1f ans\n",
            mean(snds_analytique$age_inclusion)))
cat(sprintf("Proportion diabétiques     : %.1f%%\n",
            100 * mean(snds_analytique$ald_diabete)))
cat(sprintf("Proportion hospitalisés    : %.1f%%\n",
            100 * mean(snds_analytique$nb_hospitalisations > 0)))
cat("\nDistribution des classes thérapeutiques :\n")
print(table(snds_analytique$classe_atc))
cat("==============================================\n\n")


# =============================================================================
# SAUVEGARDE
# =============================================================================

saveRDS(snds_analytique, file = "data/snds_analytique.rds")
saveRDS(er_pha_f,        file = "data/er_pha_f.rds")
saveRDS(mco_b,           file = "data/mco_b.rds")
saveRDS(ir_ben_r,        file = "data/ir_ben_r.rds")
saveRDS(ald,             file = "data/ald.rds")

cat(">>> Fichiers sauvegardés dans data/\n")
cat(">>> Simulation terminée avec succès.\n")
