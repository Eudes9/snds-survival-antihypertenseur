# =============================================================================
# FICHIER  : shiny/app.R
# OBJECTIF : Dashboard interactif — Étude pharmaco-épidémiologique SNDS
#            Visualisation des résultats pour équipes scientifiques
#            et décideurs non techniques
#
# STRUCTURE :
#   - Onglet 1 : Vue d'ensemble de la cohorte
#   - Onglet 2 : Courbes de Kaplan-Meier interactives
#   - Onglet 3 : Modèle de Cox — Hazard ratios
#   - Onglet 4 : Analyses de sensibilité
#
# LANCEMENT :
#   setwd("C:/Users/jeane/Documents/Alternance/snds-survival-antihypertenseur")
#   shiny::runApp("shiny/app.R")
# =============================================================================

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tableone)
library(broom)
library(tidyr)
library(DT)

# =============================================================================
# CHARGEMENT DES DONNÉES
# =============================================================================

# Définition du chemin racine du projet
# Shiny cherche depuis shiny/ — on remonte d'un niveau
ROOT <- dirname(getwd())

# Si lancé depuis la racine du projet
if (file.exists("data/cohorte_finale.rds")) {
  ROOT <- getwd()
}

# Chargement des données
cohorte        <- readRDS(file.path(ROOT, "data/cohorte_finale.rds"))
donnees_app    <- readRDS(file.path(ROOT, "data/donnees_appariees.rds"))
smd_comparison <- readRDS(file.path(ROOT, "outputs/smd_comparison.rds"))
res_sensi      <- readRDS(file.path(ROOT, "outputs/resultats_sensibilite.rds"))
# Palette de couleurs cohérente avec le rapport
PALETTE <- c(
  "ARA2"        = "#2C7BB6",
  "ASSOCIATION" = "#762A83",
  "BB"          = "#D7191C",
  "DIURETIQUE"  = "#FDAE61",
  "ICC"         = "#1A9641",
  "IEC"         = "#F4A582"
)

# =============================================================================
# UI — INTERFACE UTILISATEUR
# =============================================================================

ui <- dashboardPage(
  skin = "blue",

  # En-tête
  dashboardHeader(
    title = "EPI-PHARE | Antihypertenseurs",
    titleWidth = 320
  ),

  # Sidebar
  dashboardSidebar(
    width = 320,
    sidebarMenu(
      menuItem("Vue d'ensemble",
               tabName = "overview",
               icon    = icon("users")),
      menuItem("Kaplan-Meier",
               tabName = "km",
               icon    = icon("chart-line")),
      menuItem("Modèle de Cox",
               tabName = "cox",
               icon    = icon("table")),
      menuItem("Sensibilités",
               tabName = "sensi",
               icon    = icon("flask"))
    ),

    hr(),

    # Filtres globaux
    h4("Filtres globaux", style = "padding-left:15px; color:#ECF0F1"),

    selectInput(
      "filtre_classe",
      "Classes thérapeutiques",
      choices  = c("Toutes", sort(unique(cohorte$classe_atc))),
      selected = "Toutes"
    ),

    sliderInput(
      "filtre_age",
      "Tranche d'âge",
      min   = 40,
      max   = 80,
      value = c(40, 80),
      step  = 5
    ),

    radioButtons(
      "filtre_sexe",
      "Sexe",
      choices  = c("Tous" = "tous", "Hommes" = "1", "Femmes" = "2"),
      selected = "tous"
    ),

    hr(),
    p("Données simulées de type SNDS",
      style = "padding-left:15px; color:#95A5A6; font-size:11px"),
    p("Jean-Eude GBADA | 2026",
      style = "padding-left:15px; color:#95A5A6; font-size:11px")
  ),

  # Corps principal
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F8F9FA; }
      .box { border-radius: 6px; }
      .info-box { border-radius: 6px; }
      .small-box { border-radius: 6px; }
    "))),

    tabItems(

      # =========================================================
      # ONGLET 1 : VUE D'ENSEMBLE
      # =========================================================
      tabItem(
        tabName = "overview",

        # Indicateurs clés
        fluidRow(
          valueBoxOutput("vbox_n_total",    width = 3),
          valueBoxOutput("vbox_n_abandons", width = 3),
          valueBoxOutput("vbox_taux",       width = 3),
          valueBoxOutput("vbox_delai",      width = 3)
        ),

        fluidRow(
          # Distribution des classes
          box(
            title  = "Distribution des classes thérapeutiques",
            status = "primary",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_classes", height = "320px")
          ),

          # Distribution de l'âge
          box(
            title  = "Distribution de l'âge à l'inclusion",
            status = "primary",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_age", height = "320px")
          )
        ),

        fluidRow(
          # Taux d'abandon par classe
          box(
            title  = "Taux d'abandon à 12 mois par classe thérapeutique",
            status = "warning",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_taux_abandon", height = "320px")
          ),

          # Tableau descriptif
          box(
            title  = "Caractéristiques de la cohorte filtrée",
            status = "info",
            solidHeader = TRUE,
            width  = 6,
            DTOutput("table_descriptive")
          )
        )
      ),

      # =========================================================
      # ONGLET 2 : KAPLAN-MEIER
      # =========================================================
      tabItem(
        tabName = "km",

        fluidRow(
          box(
            title  = "Paramètres",
            status = "primary",
            solidHeader = TRUE,
            width  = 3,

            checkboxGroupInput(
              "km_classes",
              "Classes à afficher",
              choices  = sort(unique(cohorte$classe_atc)),
              selected = sort(unique(cohorte$classe_atc))
            ),

            hr(),

            radioButtons(
              "km_fun",
              "Type de courbe",
              choices = c(
                "Probabilité d'abandon"  = "event",
                "Probabilité de maintien" = "surv"
              ),
              selected = "event"
            ),

            hr(),

            checkboxInput("km_ci",    "Intervalles de confiance", value = TRUE),
            checkboxInput("km_table", "Tableau des patients à risque", value = TRUE),

            hr(),

            sliderInput(
              "km_xlim",
              "Période d'observation (jours)",
              min   = 30,
              max   = 365,
              value = 365,
              step  = 30
            )
          ),

          box(
            title  = "Courbes de Kaplan-Meier",
            status = "primary",
            solidHeader = TRUE,
            width  = 9,
            plotOutput("plot_km", height = "520px")
          )
        ),

        fluidRow(
          box(
            title  = "Médianes de survie par classe",
            status = "info",
            solidHeader = TRUE,
            width  = 12,
            DTOutput("table_medianes")
          )
        )
      ),

      # =========================================================
      # ONGLET 3 : MODÈLE DE COX
      # =========================================================
      tabItem(
        tabName = "cox",

        fluidRow(
          box(
            title  = "Forest plot — Hazard ratios ajustés (IC 95%)",
            status = "primary",
            solidHeader = TRUE,
            width  = 7,
            plotOutput("plot_forest", height = "480px")
          ),

          box(
            title  = "Tableau des hazard ratios",
            status = "info",
            solidHeader = TRUE,
            width  = 5,
            DTOutput("table_hr")
          )
        ),

        fluidRow(
          box(
            title  = "Balance des covariables — Love plot",
            status = "warning",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_love", height = "380px")
          ),

          box(
            title  = "Interprétation des résultats",
            status = "success",
            solidHeader = TRUE,
            width  = 6,
            h4("Résultat principal"),
            p("Aucune différence significative d'abandon entre IEC et ARA2
              après appariement par score de propension
              (HR = 0,89 ; IC 95% : 0,72–1,10 ; p = 0,282)."),
            hr(),
            h4("Classes associées à un risque élevé d'abandon"),
            tags$ul(
              tags$li("Diurétiques : HR = 1,96 (IC 95% : 1,53–2,50) ***"),
              tags$li("Bêtabloquants : HR = 1,62 (IC 95% : 1,31–1,99) ***"),
              tags$li("ICC : HR = 1,30 (IC 95% : 1,03–1,63) *")
            ),
            hr(),
            h4("Validité du modèle"),
            p("Hypothèse des risques proportionnels vérifiée
              (test de Schoenfeld global : p = 0,92)."),
            p("Concordance = 0,595.")
          )
        )
      ),

      # =========================================================
      # ONGLET 4 : ANALYSES DE SENSIBILITÉ
      # =========================================================
      tabItem(
        tabName = "sensi",

        fluidRow(
          # S1 — Définitions alternatives
          box(
            title  = "S1 — Définitions alternatives de l'abandon",
            status = "primary",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_s1", height = "320px")
          ),

          # S2 — Sans comorbidité
          box(
            title  = "S2 — Patients sans comorbidité (Charlson = 0)",
            status = "primary",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_s2", height = "320px")
          )
        ),

        fluidRow(
          # S4 — Risques compétitifs
          box(
            title  = "S4 — Modèle de risques compétitifs",
            status = "warning",
            solidHeader = TRUE,
            width  = 6,
            plotOutput("plot_s4", height = "320px")
          ),

          # S5 — E-value
          box(
            title  = "S5 — E-value : robustesse aux confondants non mesurés",
            status = "success",
            solidHeader = TRUE,
            width  = 6,
            br(),
            valueBoxOutput("vbox_evalue",    width = 12),
            valueBoxOutput("vbox_evalue_ic", width = 12),
            br(),
            p("Un facteur de confusion non mesuré devrait avoir une association
              d'au moins cette valeur avec l'exposition ET avec l'événement
              pour annuler l'association observée.",
              style = "padding: 10px; color: #555")
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER — LOGIQUE SERVEUR
# =============================================================================

server <- function(input, output, session) {

  # --- Données filtrées réactives -------------------------------------------
  cohorte_filtre <- reactive({
    df <- cohorte

    if (input$filtre_classe != "Toutes") {
      df <- df %>% filter(classe_atc == input$filtre_classe)
    }

    df <- df %>%
      filter(
        age_inclusion >= input$filtre_age[1],
        age_inclusion <= input$filtre_age[2]
      )

    if (input$filtre_sexe != "tous") {
      df <- df %>% filter(ben_sex_cod == as.integer(input$filtre_sexe))
    }

    df
  })

  # =========================================================
  # ONGLET 1 : VALUE BOXES
  # =========================================================

  output$vbox_n_total <- renderValueBox({
    valueBox(
      value    = format(nrow(cohorte_filtre()), big.mark = " "),
      subtitle = "Patients inclus",
      icon     = icon("users"),
      color    = "blue"
    )
  })

  output$vbox_n_abandons <- renderValueBox({
    valueBox(
      value    = format(sum(cohorte_filtre()$abandon), big.mark = " "),
      subtitle = "Abandons thérapeutiques",
      icon     = icon("user-times"),
      color    = "red"
    )
  })

  output$vbox_taux <- renderValueBox({
    taux <- round(100 * mean(cohorte_filtre()$abandon), 1)
    valueBox(
      value    = paste0(taux, "%"),
      subtitle = "Taux d'abandon à 12 mois",
      icon     = icon("percent"),
      color    = "orange"
    )
  })

  output$vbox_delai <- renderValueBox({
    valueBox(
      value    = paste0(median(cohorte_filtre()$delai_observe), "j"),
      subtitle = "Délai médian de suivi",
      icon     = icon("calendar"),
      color    = "green"
    )
  })

  # =========================================================
  # ONGLET 1 : GRAPHIQUES
  # =========================================================

  output$plot_classes <- renderPlot({
    cohorte_filtre() %>%
      count(classe_atc) %>%
      mutate(pct = round(100 * n / sum(n), 1)) %>%
      ggplot(aes(x = reorder(classe_atc, n), y = n,
                 fill = classe_atc)) +
      geom_col(alpha = 0.85, width = 0.7) +
      geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
                hjust = -0.1, size = 3.5) +
      scale_fill_manual(values = PALETTE) +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(x = NULL, y = "Nombre de patients") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            panel.grid.major.y = element_blank())
  })

  output$plot_age <- renderPlot({
    cohorte_filtre() %>%
      ggplot(aes(x = age_inclusion,
                 fill = factor(abandon,
                               labels = c("Censuré", "Abandon")))) +
      geom_histogram(bins = 20, alpha = 0.8, position = "stack") +
      scale_fill_manual(values = c("#2C7BB6", "#D7191C")) +
      labs(x = "Âge (années)", y = "Nombre de patients", fill = "Statut") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
  })

  output$plot_taux_abandon <- renderPlot({
    cohorte_filtre() %>%
      group_by(classe_atc) %>%
      summarise(
        taux = round(100 * mean(abandon), 1),
        .groups = "drop"
      ) %>%
      ggplot(aes(x = reorder(classe_atc, taux), y = taux,
                 fill = classe_atc)) +
      geom_col(alpha = 0.85, width = 0.7) +
      geom_text(aes(label = paste0(taux, "%")),
                hjust = -0.2, size = 4, fontface = "bold") +
      scale_fill_manual(values = PALETTE) +
      geom_hline(yintercept = 100 * mean(cohorte_filtre()$abandon),
                 linetype = "dashed", color = "grey40") +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(x = NULL, y = "Taux d'abandon (%)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            panel.grid.major.y = element_blank())
  })

  output$table_descriptive <- renderDT({
    cohorte_filtre() %>%
      group_by(classe_atc) %>%
      summarise(
        N              = n(),
        `Âge moyen`    = round(mean(age_inclusion), 1),
        `% Femmes`     = round(100 * mean(ben_sex_cod == 2), 1),
        `% Diabète`    = round(100 * mean(ald_diabete), 1),
        `% Abandon`    = round(100 * mean(abandon), 1),
        `Délai médian` = median(delai_observe),
        .groups = "drop"
      ) %>%
      rename(`Classe ATC` = classe_atc) %>%
      datatable(
        options = list(pageLength = 6, dom = "t"),
        rownames = FALSE
      )
  })

  # =========================================================
  # ONGLET 2 : KAPLAN-MEIER
  # =========================================================

  output$plot_km <- renderPlot({
    req(length(input$km_classes) >= 2)

    df_km <- cohorte %>%
      filter(classe_atc %in% input$km_classes)

    km_fit <- survfit(
      Surv(delai_observe, abandon) ~ classe_atc,
      data      = df_km,
      conf.type = "log-log"
    )

    palette_filtre <- PALETTE[names(PALETTE) %in% input$km_classes]

    p <- ggsurvplot(
      km_fit,
      data           = df_km,
      fun            = input$km_fun,
      conf.int       = input$km_ci,
      conf.int.alpha = 0.12,
      risk.table     = input$km_table,
      risk.table.height = 0.28,
      palette        = palette_filtre,
      xlim           = c(0, input$km_xlim),
      break.time.by  = max(30, round(input$km_xlim / 6 / 30) * 30),
      pval           = TRUE,
      xlab           = "Temps depuis l'initiation (jours)",
      ylab           = ifelse(input$km_fun == "event",
                              "Probabilité cumulée d'abandon",
                              "Probabilité de maintien du traitement"),
      ggtheme        = theme_minimal(base_size = 12)
    )

    print(p)
  })

  output$table_medianes <- renderDT({
    req(length(input$km_classes) >= 2)
    df_km <- cohorte %>%
      filter(classe_atc %in% input$km_classes)

    km_fit <- survfit(
      Surv(delai_observe, abandon) ~ classe_atc,
      data      = df_km,
      conf.type = "log-log"
    )
    surv_median(km_fit) %>%
      mutate(strata = gsub("classe_atc=", "", strata)) %>%
      rename(
        `Classe ATC`     = strata,
        `Médiane (j)`    = median,
        `IC bas (j)`     = lower,
        `IC haut (j)`    = upper
      ) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  # =========================================================
  # ONGLET 3 : COX
  # =========================================================

  output$plot_forest <- renderPlot({
    cox_multi <- coxph(
      Surv(delai_observe, abandon) ~
        classe_atc + age_inclusion + factor(ben_sex_cod) +
        ald_diabete + ald_insuf_card + score_charlson +
        nb_hospitalisations + ben_cmu_top,
      data = cohorte,
      ties = "efron"
    )

    tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(
        term = gsub("factor\\(ben_sex_cod\\)", "Sexe : ", term),
        term = gsub("classe_atc", "Classe : ", term),
        significatif = case_when(
          p.value < 0.001 ~ "p < 0,001",
          p.value < 0.05  ~ "p < 0,05",
          TRUE            ~ "NS"
        )
      ) %>%
      ggplot(aes(x = estimate,
                 y = reorder(term, estimate),
                 xmin = conf.low,
                 xmax = conf.high,
                 color = significatif)) +
      geom_pointrange(size = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 color = "#D7191C", linewidth = 0.7) +
      scale_x_log10() +
      scale_color_manual(values = c(
        "p < 0,001" = "#D7191C",
        "p < 0,05"  = "#FDAE61",
        "NS"        = "#2C7BB6"
      )) +
      labs(
        x     = "Hazard Ratio (IC 95%) — échelle log",
        y     = NULL,
        color = "Significativité"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
  })

  output$table_hr <- renderDT({
    cox_multi <- coxph(
      Surv(delai_observe, abandon) ~
        classe_atc + age_inclusion + factor(ben_sex_cod) +
        ald_diabete + ald_insuf_card + score_charlson +
        nb_hospitalisations + ben_cmu_top,
      data = cohorte,
      ties = "efron"
    )

    tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(
        `HR (IC 95%)` = sprintf("%.2f (%.2f–%.2f)",
                                estimate, conf.low, conf.high),
        p             = ifelse(p.value < 0.001, "<0,001",
                               as.character(round(p.value, 3))),
        Covariable    = gsub("classe_atc", "", term),
        Covariable    = gsub("factor\\(ben_sex_cod\\)2", "Femmes", Covariable)
      ) %>%
      select(Covariable, `HR (IC 95%)`, p) %>%
      datatable(
        options  = list(pageLength = 12, dom = "t"),
        rownames = FALSE
      ) %>%
      formatStyle(
        "p",
        backgroundColor = styleEqual(
          c("<0,001", "0.027", "0.028"),
          c("#FDEDEC", "#FEF9E7", "#FEF9E7")
        )
      )
  })

  output$plot_love <- renderPlot({
    smd_comparison %>%
      ggplot(aes(x = smd,
                 y = reorder(variable, smd),
                 color = periode,
                 shape = periode)) +
      geom_point(size = 3.5) +
      geom_vline(xintercept = 0.10, linetype = "dashed",
                 color = "#D7191C", linewidth = 0.7) +
      geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
      scale_color_manual(values = c(
        "Avant appariement" = "#D7191C",
        "Après appariement" = "#2C7BB6"
      )) +
      labs(
        x     = "Différence standardisée moyenne (SMD)",
        y     = NULL,
        color = NULL,
        shape = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
  })

  # =========================================================
  # ONGLET 4 : SENSIBILITÉS
  # =========================================================

  output$plot_s1 <- renderPlot({
    res_sensi$s1_definition %>%
      mutate(
        term = gsub("classe_atc", "", term),
        `HR` = estimate
      ) %>%
      ggplot(aes(x = HR, y = reorder(term, HR),
                 xmin = conf.low, xmax = conf.high,
                 color = analyse)) +
      geom_pointrange(position = position_dodge(width = 0.4), size = 0.6) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 color = "grey40") +
      scale_x_log10() +
      scale_color_manual(values = c("#2C7BB6", "#D7191C")) +
      labs(x = "HR (IC 95%)", y = NULL, color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  })

  output$plot_s2 <- renderPlot({
    res_sensi$s2_sans_comorbidite %>%
      mutate(term = gsub("classe_atc", "", term)) %>%
      ggplot(aes(x = estimate, y = reorder(term, estimate),
                 xmin = conf.low, xmax = conf.high)) +
      geom_pointrange(color = "#2C7BB6", size = 0.6) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 color = "#D7191C") +
      scale_x_log10() +
      labs(
        x        = "HR (IC 95%)",
        y        = NULL,
        subtitle = "Patients sans comorbidité (Charlson = 0)"
      ) +
      theme_minimal(base_size = 11)
  })

  output$plot_s4 <- renderPlot({
    res_sensi$s4_risques_competitifs %>%
      filter(cause == "Abandon (cause-specific)") %>%
      mutate(term = gsub("classe_atc", "", term)) %>%
      ggplot(aes(x = estimate, y = reorder(term, estimate),
                 xmin = conf.low, xmax = conf.high)) +
      geom_pointrange(color = "#762A83", size = 0.6) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 color = "#D7191C") +
      scale_x_log10() +
      labs(
        x        = "HR cause-specific (IC 95%)",
        y        = NULL,
        subtitle = "Événement : abandon (décès = risque compétiteur)"
      ) +
      theme_minimal(base_size = 11)
  })

  output$vbox_evalue <- renderValueBox({
    valueBox(
      value    = round(res_sensi$s5_evalue$e_value_hr, 2),
      subtitle = "E-value (HR ponctuel IEC vs ARA2)",
      icon     = icon("shield-alt"),
      color    = "green"
    )
  })

  output$vbox_evalue_ic <- renderValueBox({
    valueBox(
      value    = round(res_sensi$s5_evalue$e_value_ic, 2),
      subtitle = "E-value (borne inférieure IC 95%)",
      icon     = icon("info-circle"),
      color    = "blue"
    )
  })
}

# =============================================================================
# LANCEMENT
# =============================================================================

shinyApp(ui = ui, server = server)
