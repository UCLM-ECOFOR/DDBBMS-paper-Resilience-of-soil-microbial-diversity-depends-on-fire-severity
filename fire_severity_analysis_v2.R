# =============================================================================
# Fire severity effects on Pinus halepensis regeneration and soil fertility
# Outdoor mesocosm (ecotron) experiment — RStudio analysis script
# Target journal: Forest Ecology and Management / JEMA (Q1)
# British English throughout
# =============================================================================
# Data files required (place in working directory):
#   Base_datos_Pinos_R.xlsx
#   Fisico-quimicos_R.xlsx
#   Datalogger_medias_final_R.xlsx
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
pkgs <- c(
  "readxl", "dplyr", "tidyr", "tibble", "stringr", "purrr",
  "ggplot2", "patchwork", "ggbeeswarm", "ggpubr",
  "car", "FSA", "rstatix", "effectsize", "outliers",
  "lme4", "lmerTest", "emmeans", "performance",
  "FactoMineR", "factoextra",
  "corrplot", "RColorBrewer"
)
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

options(scipen = 6, digits = 4)
set.seed(2024)

# ── 1. DATA IMPORT ────────────────────────────────────────────────────────────

## 1a. Pine growth data --------------------------------------------------------
pine_raw <- Base_datos_Pinos_R

# Rename by position — avoids encoding issues with special characters
names(pine_raw) <- c("id", "treatment", "aerial_dw", "n_buds", "root_dw")

pine <- pine_raw %>%
  filter(!is.na(treatment), treatment %in% c("C", "B", "A")) %>%
  mutate(
    treatment = factor(treatment,
                       levels = c("C", "B", "A"),
                       labels = c("Control", "Low severity", "High severity")),
    # Derived allocation metrics
    RSR       = root_dw  / aerial_dw,          # Root:shoot ratio (key metric)
    total_dw  = aerial_dw + root_dw,           # Total plant dry weight
    log_aerial = log(aerial_dw + 0.001),       # Log-transform for normality check
    log_root   = log(root_dw   + 0.001)
  )

cat("Pine data — sample sizes per treatment:\n")
print(table(pine$treatment))

## 1b. Soil physico-chemical data ----------------------------------------------
soil_raw <- Fisico_quimicos_R

# Rename all columns by position
names(soil_raw) <- c(
  "date", "treatment", "time",
  "sand", "silt", "clay",
  "pH", "EC",    # EC = electrical conductivity (CE1:5)
  "OM", "N",     # organic matter, total nitrogen
  "CN",          # C:N ratio
  "NO3",         # mineral nitrogen (NaN at 1y)
  "P", "K", "Na", "Ca", "Mg",
  "K_Mg", "Ca_Mg",
  "C", "CEC"     # carbon, cation exchange capacity
)

soil <- soil_raw %>%
  mutate(
    treatment = str_trim(treatment),
    time      = str_trim(time),
    treatment = factor(treatment,
                       levels = c("UB", "LOW", "HIGH"),
                       labels = c("Control", "Low severity", "High severity")),
    time = factor(time,
                  levels = c("PRE", "POST", "1 y after"),
                  labels = c("Pre-fire", "Post-fire", "One year after"))
  ) %>%
  # Assign replicate ID (rep 1 or 2) within each treatment × time combination
  # This tracks the same lysimeter across time points
  group_by(treatment, time) %>%
  mutate(rep = row_number(),
         unit_id = paste(treatment, rep, sep = "_rep")) %>%
  ungroup()

cat("\nSoil data — observations per time × treatment:\n")
print(table(soil$treatment, soil$time))
cat("NOTE: n = 2 replicates per cell — LMM power is limited; results are indicative.\n")

## 1c. Datalogger (environmental homogeneity) ----------------------------------
dl_raw <- Datalogger_medias_final_R
names(dl_raw) <- c("treatment", "date", "pot_wt_g",
                   "air_temp_C", "rel_hum_pct", "PAR", "soil_temp_C")

datalog <- dl_raw %>%
  mutate(
    treatment = str_trim(treatment),
    treatment = factor(treatment,
                       levels = c("CONTROL", "BAJA", "ALTA"),
                       labels = c("Control", "Low severity", "High severity")),
    # Replace "NAN" strings that came in as character
    air_temp_C  = suppressWarnings(as.numeric(as.character(air_temp_C))),
    rel_hum_pct = suppressWarnings(as.numeric(as.character(rel_hum_pct)))
  )

# ── 2. COLOUR PALETTE (consistent across all figures) ────────────────────────
pal <- c("Control"       = "#2E86AB",
         "Low severity"  = "#E9724C",
         "High severity" = "#C5283D")

pal_time <- c("Pre-fire"       = "#5C6B73",
              "Post-fire"      = "#E07A5F",
              "One year after" = "#3D405B")

# ── 3. PINE DATA — DESCRIPTIVE STATISTICS ────────────────────────────────────
pine_vars <- c("aerial_dw", "root_dw", "RSR", "n_buds", "total_dw")

desc_pine <- pine %>%
  group_by(treatment) %>%
  summarise(
    n = n(),
    across(all_of(pine_vars),
           list(mean   = ~mean(., na.rm = TRUE),
                sd     = ~sd(.,   na.rm = TRUE),
                se     = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),
                median = ~median(., na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

write.csv(desc_pine, "Table_S1_pine_descriptive_stats.csv", row.names = FALSE)
cat("\nPine descriptive statistics exported.\n")

# ── 4. PINE DATA — OUTLIER DETECTION ─────────────────────────────────────────
# Grubbs test on aerial DW within High severity (known extreme values)
cat("\n--- Grubbs test: aerial dry weight in High severity ---\n")
grb <- grubbs.test(pine$aerial_dw[pine$treatment == "High severity"])
print(grb)

# Flag |z| > 3 within treatment
pine <- pine %>%
  group_by(treatment) %>%
  mutate(
    z_aerial = as.numeric(scale(aerial_dw)),
    z_root   = as.numeric(scale(root_dw)),
    outlier  = abs(z_aerial) > 3 | abs(z_root) > 3
  ) %>%
  ungroup()

cat("\nFlagged potential outliers (|z| > 3):\n")
print(pine %>% filter(outlier) %>%
        select(id, treatment, aerial_dw, root_dw, z_aerial, z_root))

# Primary analysis retains ALL observations (biologically valid).
# Sensitivity analysis without outliers is run in Section 8.

# ── 5. PINE DATA — NORMALITY & HOMOSCEDASTICITY ───────────────────────────────
cat("\n====== ASSUMPTION CHECKS — PINE DATA ======\n")

norm_pine <- map_dfr(pine_vars, function(v) {
  map_dfr(levels(pine$treatment), function(trt) {
    x <- pine[[v]][pine$treatment == trt]
    x <- x[!is.na(x)]
    sw <- shapiro.test(x)
    tibble(variable = v, treatment = trt,
           W = round(sw$statistic, 4),
           p_SW = round(sw$p.value, 4),
           normal = sw$p.value > 0.05)
  })
})

# Which variables are normal across ALL groups?
normal_pine_vars <- norm_pine %>%
  group_by(variable) %>%
  summarise(all_normal = all(normal), .groups = "drop") %>%
  filter(all_normal) %>% pull(variable)

cat("Variables with normality in all groups:", paste(normal_pine_vars, collapse = ", "), "\n")

# Levene's test
lev_pine <- map_dfr(pine_vars, function(v) {
  d <- pine[!is.na(pine[[v]]), ]
  lev <- leveneTest(as.formula(paste(v, "~ treatment")), data = d)
  tibble(variable = v,
         F_Levene = round(lev$`F value`[1], 3),
         p_Levene = round(lev$`Pr(>F)`[1],  4),
         homogen  = lev$`Pr(>F)`[1] > 0.05)
})

write.csv(norm_pine, "Table_S2_normality_tests.csv", row.names = FALSE)
write.csv(lev_pine,  "Table_S3_levene_tests.csv",    row.names = FALSE)

# ── 6. PINE DATA — BETWEEN-TREATMENT COMPARISONS ─────────────────────────────
# Decision: normality + homoscedasticity → ANOVA + Tukey
#           otherwise                    → Kruskal-Wallis + Dunn (Bonferroni)

cat("\n====== BETWEEN-TREATMENT TESTS — PINE ======\n")

pine_stat_summary <- map_dfr(pine_vars, function(v) {
  d      <- pine[!is.na(pine[[v]]), ]
  is_norm  <- v %in% normal_pine_vars
  is_homo  <- lev_pine$homogen[lev_pine$variable == v]
  method   <- ifelse(is_norm & is_homo, "ANOVA + Tukey", "Kruskal-Wallis + Dunn")

  if (is_norm & is_homo) {
    mod <- aov(as.formula(paste(v, "~ treatment")), data = d)
    sm  <- summary(mod)
    p_overall <- sm[[1]]$`Pr(>F)`[1]
    F_chi     <- sm[[1]]$`F value`[1]
    es        <- eta_squared(mod, partial = FALSE)$Eta2
    cat(sprintf("\n%s [%s]: F = %.3f, p = %.4f, eta² = %.3f\n",
                v, method, F_chi, p_overall, es))
    print(TukeyHSD(mod)$treatment)
  } else {
    kw <- kruskal.test(as.formula(paste(v, "~ treatment")), data = d)
    p_overall <- kw$p.value
    F_chi     <- kw$statistic
    es        <- kruskal_effsize(d, as.formula(paste(v, "~ treatment")))$effsize
    cat(sprintf("\n%s [%s]: chi2 = %.3f, df = %d, p = %.4f, epsilon2 = %.3f\n",
                v, method, F_chi, kw$parameter, p_overall, es))
    dn <- dunnTest(as.formula(paste(v, "~ treatment")), data = d, method = "bonferroni")
    print(dn$res)
  }
  tibble(variable = v, method = method,
         statistic = round(F_chi, 3), p_value = round(p_overall, 4),
         effect_size = round(es, 3))
})

write.csv(pine_stat_summary, "Table_S4_pine_statistical_tests.csv", row.names = FALSE)

# ── 7. PINE FIGURES ───────────────────────────────────────────────────────────

## Helper: compact letter display from Dunn test
add_letters_dunn <- function(df, var) {
  dn <- dunnTest(as.formula(paste(var, "~ treatment")),
                 data = df[!is.na(df[[var]]),], method = "bonferroni")
  # Build compact letter display manually
  pmat <- dn$res
  # Use rstatix for CLD
  stat_df <- df[!is.na(df[[var]]),] %>%
    kruskal_test(as.formula(paste(var, "~ treatment"))) %>%
    { df[!is.na(df[[var]]),] } %>%
    dunn_test(as.formula(paste(var, "~ treatment")), p.adjust.method = "bonferroni")
  cld <- stat_df %>% add_xy_position(x = "treatment")
  return(cld)
}

## Helper: boxplot + beeswarm with significance brackets
make_box <- function(df, var, ylab, log_scale = FALSE) {
  p <- ggplot(df[!is.na(df[[var]]),],
              aes(x = treatment, y = .data[[var]], fill = treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
                 colour = "black", linewidth = 0.4) +
    geom_beeswarm(aes(colour = treatment), cex = 1.8, size = 1.6,
                  alpha = 0.75, method = "swarm") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = "Fire severity treatment", y = ylab) +
    theme_bw(base_size = 11) +
    theme(axis.title   = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
  if (log_scale) p <- p + scale_y_log10()
  p
}

# Significance bracket helper using rstatix
add_sig <- function(p_obj, df, var) {
  stat_res <- df[!is.na(df[[var]]),] %>%
    dunn_test(as.formula(paste(var, "~ treatment")),
              p.adjust.method = "bonferroni") %>%
    filter(p.adj < 0.05) %>%
    add_xy_position(x = "treatment")
  if (nrow(stat_res) > 0) {
    p_obj <- p_obj + stat_pvalue_manual(stat_res, label = "p.adj.signif",
                                         tip.length = 0.01, hide.ns = TRUE,
                                         bracket.nudge.y = 0.02)
  }
  p_obj
}

# --- Fig 1: Above-ground and root dry weight (main biomass result) -----------
p_aerial <- make_box(pine, "aerial_dw", "Above-ground dry weight (g)", log_scale = FALSE)
p_root   <- make_box(pine, "root_dw",   "Root dry weight (g)", log_scale = FALSE)

# DESPUÉS de crear p_aerial y p_root (y después de add_sig())
p_aerial <- add_sig(p_aerial, pine, "aerial_dw") +
  coord_cartesian(ylim = c(0, 20))          # recorta vista; outliers siguen en cálculos

p_root <- add_sig(p_root, pine, "root_dw") +
  coord_cartesian(ylim = c(0, 8))           # root max medio = 2.51 g, outliers ~12 g

fig1 <- (p_aerial | p_root) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Fig. 1. Above-ground (A) and root (B) dry weight (g, log10 scale) of",
      "Pinus halepensis seedlings one year after fire treatments.",
      "Boxes show the interquartile range; central line = median.",
      "Individual observations are overlaid as points.",
      "Significant pairwise differences (Dunn's test, Bonferroni-adjusted) are",
      "indicated by brackets (* p < 0.05; ** p < 0.01; *** p < 0.001)."
    )
  ) &
  theme(plot.caption = element_text(size = 8, hjust = 0))
print(fig1)
ggsave("Fig1_Biomass_aerial_root.pdf", fig1, width = 9, height = 5, dpi = 300)
ggsave("Fig1_Biomass_aerial_root.png", fig1, width = 9, height = 5, dpi = 300)

# --- Fig 2: Root:shoot ratio and number of buds (allocation & vigour) ---------
p_rsr   <- make_box(pine, "RSR",    "Root:shoot ratio (dry weight)")
p_buds  <- make_box(pine, "n_buds", "Number of regenerative buds")
p_rsr   <- add_sig(p_rsr,  pine, "RSR")
p_buds  <- add_sig(p_buds, pine, "n_buds")

fig2 <- (p_rsr | p_buds) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Fig. 2. Root:shoot ratio (A) and number of regenerative buds (B) of",
      "Pinus halepensis seedlings one year after fire treatments.",
      "Notation as in Fig. 1."
    )
  ) &
  theme(plot.caption = element_text(size = 8, hjust = 0))

ggsave("Fig2_RSR_buds.pdf", fig2, width = 9, height = 5, dpi = 300)
ggsave("Fig2_RSR_buds.png", fig2, width = 9, height = 5, dpi = 300)

# --- Fig 3: Stacked biomass partitioning -------------------------------------


# 1. Pre-calcular estadísticas y posiciones acumuladas
library(dplyr)
alloc_summary <- alloc %>%
  group_by(treatment, component) %>%
  summarise(
    mean_dw = mean(dw, na.rm = TRUE),
    se_dw   = sd(dw, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  # Asegurar el orden de los factores para que Root quede abajo y Above-ground arriba
  mutate(component = factor(component, levels = c("Root", "Above-ground"))) %>%
  arrange(treatment, component) %>%
  group_by(treatment) %>%
  mutate(
    # Calcular la posición Y exacta para la barra de error (tope de cada segmento)
    ypos = cumsum(mean_dw)
  )

# 2. Generar el gráfico con los datos pre-calculados
fig3 <- ggplot(alloc_summary, aes(x = treatment, y = mean_dw, fill = component)) +
  geom_col(colour = "black", width = 0.55, linewidth = 0.35) +
  geom_errorbar(aes(ymin = ypos - se_dw, ymax = ypos + se_dw),
                width = 0.2, linewidth = 0.5, colour = "black") +
  scale_fill_manual(values = c("Above-ground" = "#4B9CD3",
                               "Root"         = "#1F5C87"),
                    name = "Plant component") +
  labs(x = "Fire severity treatment",
       y = "Mean dry weight (g)") +
  theme_bw(base_size = 11) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position  = "top")

# Guardar sin errores
ggsave("Fig3_Biomass_partitioning.png", fig3, width = 7, height = 5, dpi = 300)
print(fig3)

# ── 8. PINE — SENSITIVITY ANALYSIS (excluding |z| > 3 outliers) ──────────────
cat("\n====== SENSITIVITY ANALYSIS — PINE (outliers excluded) ======\n")
pine_no_out <- filter(pine, !outlier)
cat(sprintf("Removed %d observation(s); remaining n = %d\n",
            sum(pine$outlier), nrow(pine_no_out)))

for (v in c("aerial_dw", "root_dw", "RSR")) {
  kw <- kruskal.test(as.formula(paste(v, "~ treatment")), data = pine_no_out)
  cat(sprintf("  %s: chi2 = %.3f, p = %.4f\n", v, kw$statistic, kw$p.value))
}

# ── 9. SOIL DATA — DESCRIPTIVE STATISTICS ─────────────────────────────────────
soil_vars <- c("pH", "EC", "OM", "N", "CN", "P", "K", "Ca", "Mg", "CEC", "C")
# N-NO3 excluded (missing at 1y after)

desc_soil <- soil %>%
  group_by(treatment, time) %>%
  summarise(across(all_of(soil_vars),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd   = ~sd(.,   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

write.csv(desc_soil, "Table_S5_soil_descriptive_stats.csv", row.names = FALSE)

# ── 10. SOIL DATA — LINEAR MIXED MODELS ───────────────────────────────────────
# Design: treatment (3) × time (3), random intercept for unit_id (lysimeter)
# NOTE: n = 2 per cell (6 lysimeters total). Models are estimable but
#       degrees of freedom are very small — interpret with caution.
# lmerTest provides Satterthwaite df approximation for p-values.

cat("\n====== LINEAR MIXED MODELS — SOIL DATA ======\n")
cat("Model structure: y ~ treatment * time + (1 | unit_id)\n")
cat("Satterthwaite df approximation; emmeans pairwise comparisons.\n\n")

lmm_results <- list()
lmm_anova   <- list()

for (v in soil_vars) {
  cat(sprintf("--- %s ---\n", v))
  d <- soil[!is.na(soil[[v]]), ]
  if (nrow(d) < 6) { cat("  Insufficient data — skipped.\n"); next }

  mod <- tryCatch(
    lmer(as.formula(paste(v, "~ treatment * time + (1|unit_id)")),
         data = d, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { cat("  Model failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(mod)) next

  # Singular fit warning
  if (isSingular(mod)) cat("  Warning: singular fit (random effect near zero boundary).\n")

  # Type III ANOVA table (Satterthwaite)
  av <- anova(mod, type = "III")
  lmm_anova[[v]] <- av
  print(av)

  # Marginal means and pairwise comparisons by time point
  em_trt <- emmeans(mod, pairwise ~ treatment | time, adjust = "bonferroni")
  em_time <- emmeans(mod, pairwise ~ time | treatment, adjust = "bonferroni")
  lmm_results[[v]] <- list(model       = mod,
                            anova       = av,
                            emm_trt     = em_trt,
                            emm_time    = em_time)
  cat("  Pairwise: treatment | time\n")
  print(em_trt$contrasts)
  cat("  Pairwise: time | treatment\n")
  print(em_time$contrasts)
  cat("\n")
}

# Compile ANOVA summary table
lmm_summary_tbl <- map_dfr(names(lmm_anova), function(v) {
  av <- lmm_anova[[v]]
  tibble(
    variable         = v,
    F_treatment      = round(av["treatment",      "F value"], 3),
    p_treatment      = round(av["treatment",      "Pr(>F)"],  4),
    F_time           = round(av["time",           "F value"], 3),
    p_time           = round(av["time",           "Pr(>F)"],  4),
    F_interaction    = round(av["treatment:time", "F value"], 3),
    p_interaction    = round(av["treatment:time", "Pr(>F)"],  4)
  )
})

write.csv(lmm_summary_tbl, "Table1_LMM_ANOVA_soil.csv", row.names = FALSE)
cat("\nLMM ANOVA summary exported to Table1_LMM_ANOVA_soil.csv\n")

# ── 11. SOIL FIGURES ──────────────────────────────────────────────────────────

## Key variables for "fertilising effect" narrative
# EC, OM, N, P, K, pH — primary axis of manuscript argument

key_soil <- c("EC", "OM", "N", "P", "K", "pH")

key_labels <- c(
  EC   = "Electrical conductivity (dS m\u207b\u00b9)",
  OM   = "Organic matter (%)",
  N    = "Total nitrogen (%)",
  P    = expression(Available~P~(mg~kg^-1)),
  K    = expression(Exchangeable~K~(cmol[c]~kg^-1)),
  pH   = "pH (1:5 H\u2082O)"
)

## Fig 4: Line plot — key soil variables over time by treatment ----------------
soil_long_key <- soil %>%
  select(treatment, time, unit_id, all_of(key_soil)) %>%
  pivot_longer(all_of(key_soil), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = key_soil))

soil_sum_key <- soil_long_key %>%
  group_by(treatment, time, variable) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            se   = sd(value,   na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

fig4 <- ggplot(soil_sum_key,
               aes(x = time, y = mean, colour = treatment,
                   group = treatment, shape = treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3,
             labeller = as_labeller(c(
               EC  = "Electrical conductivity (dS/m)",
               OM  = "Organic matter (%)",
               N   = "Total nitrogen (%)",
               P   = "Available P (mg/kg)",
               K   = "Exchangeable K (cmol/kg)",
               pH  = "pH"
             ))) +
  scale_colour_manual(values = pal, name = "Fire severity") +
  scale_shape_manual(values = c(16, 17, 15), name = "Fire severity") +
  labs(x = "Sampling occasion",
       y = "Mean value (\u00b1 SE)",
       caption = paste(
         "Fig. 4. Temporal dynamics of key soil properties across fire severity",
         "treatments in the ecotron experiment (n = 2 per treatment per occasion).",
         "Points represent treatment means \u00b1 SE. EC: electrical conductivity;",
         "OM: organic matter; N: total nitrogen; P: available phosphorus;",
         "K: exchangeable potassium."
       )) +
  theme_bw(base_size = 10) +
  theme(axis.title       = element_text(face = "bold"),
        legend.position  = "bottom",
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(size = 7.5, hjust = 0))

ggsave("Fig4_Soil_temporal_dynamics.pdf", fig4, width = 10, height = 7, dpi = 300)
ggsave("Fig4_Soil_temporal_dynamics.png", fig4, width = 10, height = 7, dpi = 300)

## Fig 5: Box + jitter of soil vars at 1 year after (compare with pine data) --
soil_1y <- filter(soil, time == "One year after")

make_soil_box <- function(v, ylab) {
  ggplot(soil_1y[!is.na(soil_1y[[v]]),],
         aes(x = treatment, y = .data[[v]], fill = treatment)) +
    geom_col(data = soil_1y %>%
               group_by(treatment) %>%
               summarise(m = mean(.data[[v]], na.rm = TRUE), .groups = "drop"),
             aes(x = treatment, y = m), alpha = 0.6,
             inherit.aes = FALSE, fill = c("#2E86AB","#E9724C","#C5283D"),
             colour = "black", linewidth = 0.35, width = 0.55) +
    geom_point(size = 3, shape = 21, colour = "black",
               fill = "white", stroke = 0.8) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 10) +
    theme(axis.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
}

p_EC  <- make_soil_box("EC",  "Electrical conductivity (dS/m)")
p_OM  <- make_soil_box("OM",  "Organic matter (%)")
p_N   <- make_soil_box("N",   "Total nitrogen (%)")
p_P   <- make_soil_box("P",   "Available P (mg/kg)")
p_K   <- make_soil_box("K",   "Exchangeable K (cmol/kg)")
p_pH  <- make_soil_box("pH",  "pH")

fig5 <- (p_EC | p_OM | p_N) / (p_P | p_K | p_pH) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Fig. 5. Key soil physico-chemical properties one year after fire",
      "treatments in the ecotron experiment (n = 2 per treatment).",
      "Bars represent treatment means; individual observations shown as points."
    )
  ) &
  theme(plot.caption = element_text(size = 7.5, hjust = 0))

ggsave("Fig5_Soil_1year.pdf", fig5, width = 11, height = 7, dpi = 300)
ggsave("Fig5_Soil_1year.png", fig5, width = 11, height = 7, dpi = 300)

# ── 12. SOIL ← FIRE CHANGE: delta from pre-fire ───────────────────────────────
# Calculate the absolute change (POST - PRE) and (1y - PRE) per replicate

soil_wide <- soil %>%
  select(treatment, time, unit_id, all_of(soil_vars)) %>%
  pivot_longer(all_of(soil_vars), names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = time, values_from = value)

names(soil_wide)[names(soil_wide) == "Pre-fire"]       <- "PRE"
names(soil_wide)[names(soil_wide) == "Post-fire"]      <- "POST"
names(soil_wide)[names(soil_wide) == "One year after"] <- "ONEYEAR"

soil_delta <- soil_wide %>%
  mutate(
    delta_post   = POST    - PRE,   # immediate fire effect
    delta_1y     = ONEYEAR - PRE    # one-year net change
  )

delta_summary <- soil_delta %>%
  group_by(treatment, variable) %>%
  summarise(
    mean_delta_post = mean(delta_post, na.rm = TRUE),
    mean_delta_1y   = mean(delta_1y,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(variable %in% key_soil)

write.csv(soil_delta,   "Table_S6_soil_delta_values.csv",   row.names = FALSE)
write.csv(delta_summary,"Table_S7_soil_delta_summary.csv",  row.names = FALSE)

## Fig 6: Heatmap of delta (1y - PRE) across treatments ----------------------
fig6_data <- delta_summary %>%
  select(treatment, variable, mean_delta_1y) %>%
  mutate(variable = factor(variable, levels = rev(key_soil)))

fig6 <- ggplot(fig6_data,
               aes(x = treatment, y = variable, fill = mean_delta_1y)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = round(mean_delta_1y, 3)), size = 3.2) +
  scale_fill_gradient2(
    low      = "#2E86AB", mid = "white", high = "#C5283D",
    midpoint = 0, name = "\u0394 (1y \u2212 PRE)"
  ) +
  labs(x = "Fire severity treatment",
       y = "Soil variable",
       caption = paste(
         "Fig. 6. Change in key soil properties from pre-fire baseline to",
         "one year after burning (\u0394 = 1y after \u2212 pre-fire value),",
         "averaged across replicates (n = 2). Red = increase; blue = decrease."
       )) +
  theme_bw(base_size = 11) +
  theme(axis.title   = element_text(face = "bold"),
        plot.caption = element_text(size = 7.5, hjust = 0))

ggsave("Fig6_Soil_delta_heatmap.pdf", fig6, width = 7, height = 5, dpi = 300)
ggsave("Fig6_Soil_delta_heatmap.png", fig6, width = 7, height = 5, dpi = 300)

# ── 13. CORRELATION: SOIL (1y) × PINE GROWTH ─────────────────────────────────
# Aggregate pine to treatment means, join with soil means at 1y after

pine_means <- pine %>%
  group_by(treatment) %>%
  summarise(across(c(aerial_dw, root_dw, RSR, n_buds),
                   mean, na.rm = TRUE), .groups = "drop")

soil_1y_means <- soil %>%
  filter(time == "One year after") %>%
  group_by(treatment) %>%
  summarise(across(all_of(soil_vars), mean, na.rm = TRUE), .groups = "drop")

corr_data <- left_join(pine_means, soil_1y_means, by = "treatment")

cat("\n====== CORRELATION: SOIL (1y after) × PINE GROWTH ======\n")
cat("NOTE: n = 3 treatment-level means — Spearman r reported descriptively.\n\n")

# Spearman r between each soil var and each pine var
cor_mat <- map_dfr(soil_vars, function(sv) {
  map_dfr(c("aerial_dw", "root_dw", "RSR", "n_buds"), function(pv) {
    x <- corr_data[[sv]]
    y <- corr_data[[pv]]
    if (sum(!is.na(x) & !is.na(y)) < 3) return(NULL)
    sp <- cor.test(x, y, method = "spearman", exact = FALSE)
    tibble(soil_var  = sv, pine_var  = pv,
           rho       = round(sp$estimate, 3),
           p_value   = round(sp$p.value,  3))
  })
})

print(cor_mat)
write.csv(cor_mat, "Table_S8_soil_pine_correlations.csv", row.names = FALSE)

## Fig 7: Scatter plots — EC, N, P vs aerial DW (narrative focus) -------------
library(ggplot2)
library(patchwork)

make_scatter <- function(sv, pv, xlab, ylab, xlim_max = NULL) {
  # Si no le pasamos un límite, calculamos uno que deje espacio a la derecha
  if (is.null(xlim_max)) {
    lims <- range(corr_data[[sv]], na.rm = TRUE)
    xlim_max <- lims[2] + (lims[2] - lims[1]) * 0.25 # Añade un 25% de espacio extra a la derecha
  }
  
  ggplot(corr_data, aes(x = .data[[sv]], y = .data[[pv]], colour = treatment, label = treatment)) +
    geom_point(size = 4, shape = 19) +
    geom_text(hjust = -0.2, vjust = 0.5, size = 3.2, show.legend = FALSE) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_continuous(limits = c(NA, xlim_max)) + # Expande el eje X para que quepa la etiqueta
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 10) +
    theme(axis.title = element_text(face = "bold"))
}

# Al llamarlas, el gráfico ya dejará aire a la derecha para que el punto rojo muestre su etiqueta
p_ec_aerial <- make_scatter("EC", "aerial_dw", "Electrical conductivity (dS/m)", "Above-ground DW (g)")
p_n_aerial  <- make_scatter("N", "aerial_dw", "Total nitrogen (%)", "Above-ground DW (g)")
p_p_aerial  <- make_scatter("P", "aerial_dw", "Available P (mg/kg)", "Above-ground DW (g)")
p_ec_rsr    <- make_scatter("EC", "RSR", "Electrical conductivity (dS/m)", "Root:shoot ratio")

# Ensamblaje igual que antes...

fig7 <- (p_ec_aerial | p_n_aerial) / (p_p_aerial | p_ec_rsr) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Fig. 7. Relationships between key soil properties (one year after burning)",
      "and Pinus halepensis seedling growth metrics across fire severity treatments.",
      "Each point represents a treatment mean (n = 3 treatment levels).",
      "Spearman correlation coefficients are reported in the text."
    )
  ) &
  theme(plot.caption = element_text(size = 7.5, hjust = 0))
print(fig7)
ggsave("Fig7_Soil_pine_correlations.pdf", fig7, width = 9, height = 7, dpi = 300)
ggsave("Fig7_Soil_pine_correlations.png", fig7, width = 9, height = 7, dpi = 300)

# ── 14. SOIL PCA ──────────────────────────────────────────────────────────────
pca_vars_soil <- c("pH","EC","OM","N","P","K","Ca","Mg","CEC","C")

pca_soil_data <- soil %>%
  select(treatment, time, all_of(pca_vars_soil)) %>%
  drop_na()

pca_res <- PCA(pca_soil_data[, pca_vars_soil], scale.unit = TRUE, graph = FALSE)

cat("\n====== SOIL PCA ======\n")
cat("Eigenvalues (components with > 1 variance):\n")
print(pca_res$eig[pca_res$eig[,1] > 1, ])

# Group labels for biplot
groups <- interaction(pca_soil_data$treatment, pca_soil_data$time, sep = " | ")

fig8 <- fviz_pca_biplot(
  pca_res,
  habillage    = groups,
  addEllipses  = FALSE,
  label        = "var",
  col.var      = "black",
  repel        = TRUE,
  palette      = c("#2E86AB","#5BA4CF","#A8D5E2",
                   "#E9724C","#F0A07B","#F9CFBA",
                   "#C5283D","#E07A90","#F3BACB"),
  title        = ""
) +
  labs(
    x       = paste0("PC1 (", round(pca_res$eig[1,2],1), "%)"),
    y       = paste0("PC2 (", round(pca_res$eig[2,2],1), "%)"),
    colour  = "Treatment | Time",
    shape   = "Treatment | Time",
    caption = paste(
      "Fig. 8. Principal component analysis (PCA) biplot of soil physico-chemical",
      "properties across fire severity treatments and sampling occasions.",
      "Variable loadings (arrows) indicate the direction of maximum variation.",
      "Variance explained by each axis is shown in parentheses."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(legend.text  = element_text(size = 7),
        plot.caption = element_text(size = 7.5, hjust = 0))

ggsave("Fig8_Soil_PCA.pdf", fig8, width = 9, height = 7, dpi = 300)
ggsave("Fig8_Soil_PCA.png", fig8, width = 9, height = 7, dpi = 300)

# ── 15. ENVIRONMENTAL HOMOGENEITY (DATALOGGER) ────────────────────────────────
cat("\n====== ENVIRONMENTAL CONDITIONS (Datalogger) ======\n")

env_summary <- datalog %>%
  group_by(treatment) %>%
  summarise(
    n               = n(),
    air_temp_mean   = round(mean(air_temp_C,   na.rm = TRUE), 2),
    air_temp_sd     = round(sd(air_temp_C,     na.rm = TRUE), 2),
    rel_hum_mean    = round(mean(rel_hum_pct,  na.rm = TRUE), 2),
    rel_hum_sd      = round(sd(rel_hum_pct,    na.rm = TRUE), 2),
    soil_temp_mean  = round(mean(soil_temp_C,  na.rm = TRUE), 2),
    soil_temp_sd    = round(sd(soil_temp_C,    na.rm = TRUE), 2),
    pot_wt_mean_g   = round(mean(pot_wt_g,     na.rm = TRUE), 1),
    .groups = "drop"
  )
print(env_summary)
write.csv(env_summary, "Table_S9_environmental_conditions.csv", row.names = FALSE)

# Kruskal-Wallis to confirm homogeneity
for (ev in c("air_temp_C", "rel_hum_pct", "soil_temp_C")) {
  d <- datalog[!is.na(datalog[[ev]]) & is.finite(datalog[[ev]]), ]
  if (nrow(d) < 4 || length(unique(d$treatment)) < 2) next
  kw <- kruskal.test(as.formula(paste(ev, "~ treatment")), data = d)
  cat(sprintf("  %s: chi2 = %.3f, p = %.4f  %s\n",
              ev, kw$statistic, kw$p.value,
              ifelse(kw$p.value > 0.05,
                     "(NS — conditions homogeneous across treatments)",
                     "(SIGNIFICANT — investigate)")))
}

# ── 16. MANUSCRIPT SUMMARY TABLES ─────────────────────────────────────────────

## Table 1 (Main text): Pine biomass — mean ± SE per treatment + test results --
tbl1 <- pine %>%
  group_by(treatment) %>%
  summarise(
    n                   = n(),
    `Above-ground DW (g)` = sprintf("%.3f \u00b1 %.3f",
                                     mean(aerial_dw, na.rm=T),
                                     sd(aerial_dw,   na.rm=T)/sqrt(n())),
    `Root DW (g)`       = sprintf("%.3f \u00b1 %.3f",
                                   mean(root_dw,  na.rm=T),
                                   sd(root_dw,    na.rm=T)/sqrt(n())),
    `Root:shoot ratio`  = sprintf("%.3f \u00b1 %.3f",
                                   mean(RSR,  na.rm=T),
                                   sd(RSR,    na.rm=T)/sqrt(n())),
    `Number of buds`    = sprintf("%.1f \u00b1 %.1f",
                                   mean(n_buds, na.rm=T),
                                   sd(n_buds,   na.rm=T)/sqrt(n())),
    `Total DW (g)`      = sprintf("%.3f \u00b1 %.3f",
                                   mean(total_dw, na.rm=T),
                                   sd(total_dw,   na.rm=T)/sqrt(n())),
    .groups = "drop"
  )

write.csv(tbl1, "Table1_Pine_biomass_summary.csv", row.names = FALSE)

## Table 2 (Main text): Soil properties by treatment at each time point -------
tbl2 <- soil %>%
  select(treatment, time, all_of(key_soil)) %>%
  pivot_longer(all_of(key_soil), names_to = "variable", values_to = "value") %>%
  group_by(variable, treatment, time) %>%
  summarise(
    mean_se = sprintf("%.3f \u00b1 %.3f",
                      mean(value, na.rm=TRUE),
                      sd(value, na.rm=TRUE)/sqrt(n())),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = c(treatment, time), values_from = mean_se)

write.csv(tbl2, "Table2_Soil_properties_summary.csv", row.names = FALSE)

# ── 17. FINAL SUMMARY ─────────────────────────────────────────────────────────
cat("\n\n====== ANALYSIS COMPLETE ======\n")
cat("Output files generated:\n")
outputs <- c(
  "FIGURES:",
  "  Fig1_Biomass_aerial_root.pdf/png",
  "  Fig2_RSR_buds.pdf/png",
  "  Fig3_Biomass_partitioning.pdf/png",
  "  Fig4_Soil_temporal_dynamics.pdf/png",
  "  Fig5_Soil_1year.pdf/png",
  "  Fig6_Soil_delta_heatmap.pdf/png",
  "  Fig7_Soil_pine_correlations.pdf/png",
  "  Fig8_Soil_PCA.pdf/png",
  "",
  "MAIN TEXT TABLES:",
  "  Table1_Pine_biomass_summary.csv",
  "  Table1_LMM_ANOVA_soil.csv",
  "  Table2_Soil_properties_summary.csv",
  "",
  "SUPPLEMENTARY:",
  "  Table_S1_pine_descriptive_stats.csv",
  "  Table_S2_normality_tests.csv",
  "  Table_S3_levene_tests.csv",
  "  Table_S4_pine_statistical_tests.csv",
  "  Table_S5_soil_descriptive_stats.csv",
  "  Table_S6_soil_delta_values.csv",
  "  Table_S7_soil_delta_summary.csv",
  "  Table_S8_soil_pine_correlations.csv",
  "  Table_S9_environmental_conditions.csv"
)
cat(paste(outputs, collapse = "\n"), "\n")
