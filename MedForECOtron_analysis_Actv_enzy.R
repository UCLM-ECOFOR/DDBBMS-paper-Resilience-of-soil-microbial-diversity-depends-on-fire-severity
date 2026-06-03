>   # ── Full saturated model (quadratic × treatment) ──────────────────────────
  >   m_full <- tryCatch(
    +     lmer(reformulate(c("poly1 * BurnSev", "poly2 * BurnSev", "(1|Monolith)"),
                           +                      response = v),
               +          data = dat, REML = FALSE,
               +          control = lmerControl(optimizer = "bobyqa")),
    +     error = function(e) { cat("  Full model failed:", e$message, "\n"); NULL }
    +   )
  # ── Candidate models for AICc-based selection ─────────────────────────────
  >   m1 <- lmer(reformulate(c("poly1 + BurnSev", "(1|Monolith)"), response = v),
                 +              data = dat, REML = FALSE,
                 +              control = lmerControl(optimizer = "bobyqa"))
  m2 <- lmer(reformulate(c("poly1 + poly2 + BurnSev", "(1|Monolith)"), response = v),
             +              data = dat, REML = FALSE,
             +              control = lmerControl(optimizer = "bobyqa"))
  >   m3 <- lmer(reformulate(c("poly1", "(1|Monolith)"), response = v),
                 +              data = dat, REML = FALSE,
                 +              control = lmerControl(optimizer = "bobyqa"))
  m4 <- lmer(reformulate(c("poly1 + poly2", "(1|Monolith)"), response = v),
             +              data = dat, REML = FALSE,
             +              control = lmerControl(optimizer = "bobyqa"))
  m5 <- lmer(reformulate(c("BurnSev", "(1|Monolith)"), response = v),
             +              data = dat, REML = FALSE,
             +              control = lmerControl(optimizer = "bobyqa"))
  m0 <- lmer(reformulate("(1|Monolith)", response = v),
             +              data = dat, REML = FALSE,
             +              control = lmerControl(optimizer = "bobyqa"))
  print(lmm_anova)
  $GLU
  Type III Analysis of Variance Table with Satterthwaite's method
         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
poly1   0.28599 0.28599     1    40  2.2794 0.1389592    
poly2   1.61143 1.61143     1    40 12.8436 0.0009095 ***
BurnSev 1.38160 0.69080     2     3  5.5059 0.0990692 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$UR
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
  poly1   1.7178 1.71780     1    44  4.8375 0.03315 *
    BurnSev 1.9737 0.98686     2    44  2.7791 0.07302 .
  ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  $PHP
  Type III Analysis of Variance Table with Satterthwaite's method
        Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
poly1   0.4367 0.43671     1    43  1.9195 0.1730592    
poly2   1.7776 1.77757     1    43  7.8129 0.0077211 ** 
BurnSev 4.2918 2.14589     2    43  9.4318 0.0004015 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$DHA
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
  poly1   2111.48 2111.48     1    44 93.6278 1.826e-12 ***
    BurnSev  421.25  210.63     2    44  9.3397 0.0004161 ***
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  $BSR
  Type III Analysis of Variance Table with Satterthwaite's method
      Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
poly1 14.299  14.299     1    39   8.838  0.005036 ** 
poly2 33.084  33.084     1    39  20.449 5.583e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$GRSP
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)    
  poly1  875858  875858     1    40   6.831  0.01257 *  
    poly2 4701572 4701572     1    40  36.668 3.95e-07 ***
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  $Cmic
  Type III Analysis of Variance Table with Satterthwaite's method
     Sum Sq Mean Sq NumDF DenDF F value Pr(>F)

>   calc_AICc <- function(model) {
+     k   <- length(fixef(model)) + 1
+     n   <- nobs(model)
+     AIC(model) + (2 * k * (k + 1)) / (n - k - 1)
+   }
>   ctrl_ml   <- lmerControl(optimizer = "bobyqa",
+                            check.conv.singular = .makeCC("ignore", tol = 1e-4))
>   ctrl_reml <- lmerControl(optimizer = "bobyqa",
+                            check.conv.singular = .makeCC("ignore", tol = 1e-4))
>   lmm_results   <- list()
>   lmm_summaries <- list()
>   lmm_anova     <- list()
>   lmm_posthoc   <- list()
>   for (v in enzyme_vars) {
+     
+     cat("\n══════════════════════════\n")
+     cat("Variable:", v, "\n")
+     cat("══════════════════════════\n")
+     
+     dat <- df_post %>%
+       filter(!is.na(.data[[v]])) %>%
+       mutate(
+         poly1 = poly(TSF_months, 2)[, 1],
+         poly2 = poly(TSF_months, 2)[, 2]
+       )
+     
+     f1 <- as.formula(paste(v, "~ poly1 + BurnSev + (1|Monolith)"))
+     f2 <- as.formula(paste(v, "~ poly1 + poly2 + BurnSev + (1|Monolith)"))
+     f3 <- as.formula(paste(v, "~ poly1 + (1|Monolith)"))
+     f4 <- as.formula(paste(v, "~ poly1 + poly2 + (1|Monolith)"))
+     f5 <- as.formula(paste(v, "~ BurnSev + (1|Monolith)"))
+     f0 <- as.formula(paste(v, "~ 1 + (1|Monolith)"))
+     
+     formulas <- list(
+       "Linear_TSF + BurnSev"    = f1,
+       "Quadratic_TSF + BurnSev" = f2,
+       "Linear_TSF only"         = f3,
+       "Quadratic_TSF only"      = f4,
+       "BurnSev only"            = f5,
+       "Null (random only)"      = f0
+     )
+     
+     models_ml <- lapply(formulas, function(f) {
+       tryCatch(
+         lmer(f, data = dat, REML = FALSE, control = ctrl_ml),
+         error = function(e) NULL
+       )
+     })
+     models_ml <- Filter(Negate(is.null), models_ml)
+     
+     aic_table <- data.frame(
+       Model = names(models_ml),
+       AIC   = sapply(models_ml, AIC),
+       AICc  = sapply(models_ml, calc_AICc)
+     ) %>%
+       mutate(DAICc = AICc - min(AICc)) %>%
+       arrange(AICc)
+     
+     cat("  Model selection (AICc):\n"); print(aic_table)
+     
+     best_name    <- aic_table$Model[1]
+     best_formula <- formulas[[best_name]]
+     
+     best_model_reml <- tryCatch(
+       lmer(best_formula, data = dat, REML = TRUE, control = ctrl_reml),
+       error = function(e) { cat("  REML refit failed:", e$message, "\n"); NULL }
+     )
+     
+     if (is.null(best_model_reml)) {
+       cat("  Skipping", v, "— model could not be fitted.\n")
+       next
+     }
+     
+     lmm_results[[v]]   <- best_model_reml
+     lmm_summaries[[v]] <- summary(best_model_reml)
+     lmm_anova[[v]]     <- anova(best_model_reml, type = "III")
+     
+     cat("\n  Best model:", best_name, "\n")
+     cat("\n  ANOVA (type III, Satterthwaite):\n")
+     print(lmm_anova[[v]])
+     
+     r2_vals <- performance::r2(best_model_reml)
+     cat(sprintf("  Marginal R² = %.3f | Conditional R² = %.3f\n",
+                 r2_vals$R2_marginal, r2_vals$R2_conditional))
+     
+     if (grepl("BurnSev", best_name)) {
+       emm <- emmeans(best_model_reml, ~ BurnSev)
+       ph  <- pairs(emm, adjust = "tukey")
+       lmm_posthoc[[v]] <- ph
+       cat("\n  Post-hoc (Tukey):\n"); print(ph)
+     }
+     
+   }  # ← único cierre del bucle for

══════════════════════════
Variable: GLU 
══════════════════════════
  Model selection (AICc):
                                          Model      AIC     AICc     DAICc
Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 46.01566 48.06444  0.000000
Quadratic_TSF only           Quadratic_TSF only 50.78811 51.71835  3.653907
BurnSev only                       BurnSev only 56.29670 57.22693  9.162494
Linear_TSF + BurnSev       Linear_TSF + BurnSev 56.39779 57.82636  9.761926
Null (random only)           Null (random only) 60.25705 60.52372 12.459281
Linear_TSF only                 Linear_TSF only 60.48335 61.02880 12.964365

  Best model: Quadratic_TSF + BurnSev 

  ANOVA (type III, Satterthwaite):
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
  poly1   0.28599 0.28599     1    40  2.2794 0.1389592    
  poly2   1.61143 1.61143     1    40 12.8436 0.0009095 ***
    BurnSev 1.38160 0.69080     2     3  5.5059 0.0990692 .  
  ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  Marginal R² = 0.371 | Conditional R² = 0.388
  
  Post-hoc (Tukey):
    contrast   estimate    SE df t.ratio p.value
  UNB - LOW     0.221 0.138  3   1.600  0.3711
  UNB - HIGH    0.458 0.138  3   3.318  0.0890
  LOW - HIGH    0.237 0.138  3   1.718  0.3335
  
  Degrees-of-freedom method: kenward-roger 
  P value adjustment: tukey method for comparing a family of 3 estimates 
  
  ══════════════════════════
  Variable: UR 
  ══════════════════════════
  Model selection (AICc):
    Model      AIC     AICc     DAICc
  Linear_TSF + BurnSev       Linear_TSF + BurnSev 94.34413 95.77270 0.0000000
  Linear_TSF only                 Linear_TSF only 95.64872 96.19417 0.4214689
  BurnSev only                       BurnSev only 97.35099 98.28122 2.5085206
  Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 96.32998 98.37877 2.6060628
  Quadratic_TSF only           Quadratic_TSF only 97.63544 98.56567 2.7929725
  Null (random only)           Null (random only) 98.32990 98.59656 2.8238593
  
  Best model: Linear_TSF + BurnSev 
  
  ANOVA (type III, Satterthwaite):
    Type III Analysis of Variance Table with Satterthwaite's method
        Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
poly1   1.7178 1.71780     1    44  4.8375 0.03315 *
BurnSev 1.9737 0.98686     2    44  2.7791 0.07302 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  Marginal R² = 0.181 | Conditional R² = NA

  Post-hoc (Tukey):
 contrast   estimate    SE df t.ratio p.value
 UNB - LOW    0.0595 0.211  3   0.282  0.9576
 UNB - HIGH   0.4568 0.211  3   2.168  0.2232
 LOW - HIGH   0.3973 0.211  3   1.886  0.2866

Degrees-of-freedom method: kenward-roger 
P value adjustment: tukey method for comparing a family of 3 estimates 

══════════════════════════
Variable: PHP 
══════════════════════════
  Model selection (AICc):
                                          Model      AIC     AICc     DAICc
Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 73.87260 75.92138  0.000000
BurnSev only                       BurnSev only 79.66604 80.59627  4.674889
Linear_TSF + BurnSev       Linear_TSF + BurnSev 79.88623 81.31480  5.393418
Quadratic_TSF only           Quadratic_TSF only 82.66021 83.59044  7.669058
Null (random only)           Null (random only) 87.32580 87.59246 11.671081
Linear_TSF only                 Linear_TSF only 87.75240 88.29786 12.376473

  Best model: Quadratic_TSF + BurnSev 

  ANOVA (type III, Satterthwaite):
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
  poly1   0.4367 0.43671     1    43  1.9195 0.1730592    
  poly2   1.7776 1.77757     1    43  7.8129 0.0077211 ** 
    BurnSev 4.2918 2.14589     2    43  9.4318 0.0004015 ***
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  Marginal R² = 0.378 | Conditional R² = NA
  
  Post-hoc (Tukey):
    contrast   estimate    SE df t.ratio p.value
  UNB - LOW     0.448 0.169  3   2.658  0.1477
  UNB - HIGH    0.726 0.169  3   4.304  0.0463
  LOW - HIGH    0.278 0.169  3   1.646  0.3559
  
  Degrees-of-freedom method: kenward-roger 
  P value adjustment: tukey method for comparing a family of 3 estimates 
  
  ══════════════════════════
  Variable: DHA 
  ══════════════════════════
  Model selection (AICc):
    Model      AIC     AICc     DAICc
  Linear_TSF + BurnSev       Linear_TSF + BurnSev 293.6008 295.0293  0.000000
  Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 294.8177 296.8665  1.837121
  Linear_TSF only                 Linear_TSF only 301.2923 301.8377  6.808380
  Quadratic_TSF only           Quadratic_TSF only 302.5763 303.5066  8.477251
  BurnSev only                       BurnSev only 346.3382 347.2684 52.239097
  Null (random only)           Null (random only) 348.4382 348.7049 53.675553
  
  Best model: Linear_TSF + BurnSev 
  
  ANOVA (type III, Satterthwaite):
    Type III Analysis of Variance Table with Satterthwaite's method
         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
poly1   2111.48 2111.48     1    44 93.6278 1.826e-12 ***
BurnSev  421.25  210.63     2    44  9.3397 0.0004161 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  Marginal R² = 0.705 | Conditional R² = NA

  Post-hoc (Tukey):
 contrast   estimate   SE df t.ratio p.value
 UNB - LOW      3.92 1.68  3   2.337  0.1930
 UNB - HIGH     7.25 1.68  3   4.317  0.0459
 LOW - HIGH     3.32 1.68  3   1.980  0.2634

Degrees-of-freedom method: kenward-roger 
P value adjustment: tukey method for comparing a family of 3 estimates 

══════════════════════════
Variable: BSR 
══════════════════════════
  Model selection (AICc):
                                          Model      AIC     AICc     DAICc
Quadratic_TSF only           Quadratic_TSF only 146.2848 147.3659  0.000000
Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 148.8316 151.2316  3.865736
Linear_TSF only                 Linear_TSF only 161.9905 162.6221 15.256216
Null (random only)           Null (random only) 165.8117 166.1194 18.753491
Linear_TSF + BurnSev       Linear_TSF + BurnSev 165.0429 166.7096 19.343680
BurnSev only                       BurnSev only 168.9879 170.0690 22.703111

  Best model: Quadratic_TSF only 

  ANOVA (type III, Satterthwaite):
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
  poly1 14.299  14.299     1    39   8.838  0.005036 ** 
    poly2 33.084  33.084     1    39  20.449 5.583e-05 ***
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  Marginal R² = 0.417 | Conditional R² = NA
  
  ══════════════════════════
  Variable: GRSP 
  ══════════════════════════
  Model selection (AICc):
    Model      AIC     AICc      DAICc
  Quadratic_TSF only           Quadratic_TSF only 714.1917 715.1219  0.0000000
  Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 714.0145 716.0632  0.9413376
  Linear_TSF only                 Linear_TSF only 739.5173 740.0627 24.9408071
  Linear_TSF + BurnSev       Linear_TSF + BurnSev 739.6955 741.1240 26.0021399
  Null (random only)           Null (random only) 741.1019 741.3686 26.2466949
  BurnSev only                       BurnSev only 741.4435 742.3737 27.2517866
  
  Best model: Quadratic_TSF only 
  
  ANOVA (type III, Satterthwaite):
    Type III Analysis of Variance Table with Satterthwaite's method
       Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)    
poly1  875858  875858     1    40   6.831  0.01257 *  
poly2 4701572 4701572     1    40  36.668 3.95e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  Marginal R² = 0.426 | Conditional R² = 0.540

══════════════════════════
Variable: Cmic 
══════════════════════════
  Model selection (AICc):
                                          Model      AIC     AICc    DAICc
Null (random only)           Null (random only) 741.6271 741.8938 0.000000
Linear_TSF only                 Linear_TSF only 743.1931 743.7385 1.844712
Quadratic_TSF only           Quadratic_TSF only 743.9988 744.9290 3.035245
BurnSev only                       BurnSev only 744.5766 745.5069 3.613067
Linear_TSF + BurnSev       Linear_TSF + BurnSev 746.1329 747.5615 5.667681
Quadratic_TSF + BurnSev Quadratic_TSF + BurnSev 746.9116 748.9604 7.066633

  Best model: Null (random only) 

  ANOVA (type III, Satterthwaite):
Type III Analysis of Variance Table with Satterthwaite's method
  Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
  Marginal R² = 0.000 | Conditional R² = 0.013
  
  > df_avg <- df %>%
    +   filter(TRAT %in% c("PRE","15 MONTHS AFTER")) %>%
    +   group_by(Monolith, BurnSev, TRAT) %>%
    +   summarise(across(all_of(enzyme_vars), ~mean(.x, na.rm = TRUE)),
                  +             .groups = "drop") %>%
    +   pivot_longer(cols = all_of(enzyme_vars), names_to = "Variable", values_to = "Value")
  > paired_results <- list()
  > for (sev in c("UNB","LOW","HIGH")) {
    +   for (v in enzyme_vars) {
      +     dat_pair <- df_avg %>%
        +       filter(BurnSev == sev, Variable == v) %>%
        +       pivot_wider(names_from = TRAT, values_from = Value) %>%
        +       filter(!is.na(PRE), !is.na(`15 MONTHS AFTER`))
      + 
        +     if (nrow(dat_pair) >= 2) {
          +       t_res  <- t.test(dat_pair$`15 MONTHS AFTER`, dat_pair$PRE,
                                   +                        paired = TRUE, alternative = "two.sided")
          +       d_val  <- (mean(dat_pair$`15 MONTHS AFTER` - dat_pair$PRE, na.rm = TRUE)) /
            +                  sd(dat_pair$`15 MONTHS AFTER` - dat_pair$PRE, na.rm = TRUE)
          + 
            +       paired_results[[paste(sev, v, sep = "_")]] <- tibble(
              +         BurnSev  = sev,
              +         Variable = v,
              +         t_stat   = t_res$statistic,
              +         df       = t_res$parameter,
              +         p_value  = t_res$p.value,
              +         Cohen_d  = d_val,
              +         mean_diff = mean(dat_pair$`15 MONTHS AFTER` - dat_pair$PRE, na.rm = TRUE)
              +       )
            +     }
      +   }
    + }
  > paired_df <- bind_rows(paired_results) %>%
    +   mutate(p_adj    = p.adjust(p_value, method = "BH"),  # Benjamini-Hochberg
               +          sig_label = case_when(
                 +            p_adj < 0.001 ~ "***",
                 +            p_adj < 0.01  ~ "**",
                 +            p_adj < 0.05  ~ "*",
                 +            TRUE          ~ "ns"
                 +          ))
  > print(paired_df)
  # A tibble: 21 × 9
  BurnSev Variable   t_stat    df p_value  Cohen_d mean_diff p_adj sig_label
  <chr>   <chr>       <dbl> <dbl>   <dbl>    <dbl>     <dbl> <dbl> <chr>    
    1 UNB     GLU        1.15       1  0.456    0.813     0.539  0.683 ns       
  2 UNB     UR        -0.312      1  0.807   -0.221    -0.0962 0.892 ns       
  3 UNB     PHP        6.82       1  0.0926   4.83      0.951  0.364 ns       
  4 UNB     DHA      -33.0        1  0.0193 -23.3     -21.9    0.364 ns       
  5 UNB     BSR        0.609      1  0.652    0.430     2.36   0.810 ns       
  6 UNB     GRSP       2.11       1  0.282    1.49    659.     0.591 ns       
  7 UNB     Cmic       0.600      1  0.656    0.424   481.     0.810 ns       
  8 LOW     GLU        0.0989     1  0.937    0.0700    0.0116 0.961 ns       
  9 LOW     UR        -4.77       1  0.132   -3.37     -0.330  0.364 ns       
  10 LOW     PHP        1.31       1  0.416    0.924     0.127  0.683 ns       
  # ℹ 11 more rows
  # ℹ Use `print(n = ...)` to see more rows

    > pca_vars <- c("GLU","UR","PHP","DHA","BSR","GRSP","Cmic","SOM","P","EC","Ca")
> df_pca <- df %>%
+   filter(TRAT != "PRE") %>%
+   select(Monolith, BurnSev, TRAT, TSF_months, all_of(pca_vars)) %>%
+   drop_na()
> pca_vars <- c("GLU","UR","PHP","DHA","BSR","GRSP","Cmic","SOM","P","EC","Ca")
> df_pca <- df %>%
+   filter(TRAT != "PRE") %>%
+   select(Monolith, BurnSev, TRAT, TSF_months, all_of(pca_vars)) %>%
+   drop_na()
> # ── Ejecutar PCA con base R ───────────────────────────────────────────────────
> pca_res <- prcomp(df_pca[, pca_vars], center = TRUE, scale. = TRUE)
> # Varianza explicada
> var_exp  <- summary(pca_res)$importance[2, ] * 100   # % varianza por PC
> cum_var  <- summary(pca_res)$importance[3, ] * 100   # % acumulado
> cat("Variance explained:\n")
Variance explained:
> print(round(var_exp[1:5], 1))
 PC1  PC2  PC3  PC4  PC5 
37.2 23.3 13.4  9.7  7.3 
