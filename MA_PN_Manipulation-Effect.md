Master Thesis Project PN - Manipulation Effects
================
Paula Neubauer
2025-07-31

``` r
df_grip <- df_grip %>%
  filter(id != "19102022")
```

# Subject Level

## distribution of dependent variables

``` r
df_grip <- df_grip %>%
  group_by(id) %>%
  mutate(
    z_rpd = (rpd - mean(rpd, na.rm = TRUE)) / sd(rpd, na.rm = TRUE)
  ) %>%
  ungroup()

df_grip <- df_grip %>%
  group_by(id) %>%
  mutate(
    z_rpd_low = (rpd_low - mean(rpd_low, na.rm = TRUE)) / sd(rpd_low, na.rm = TRUE)
  ) %>%
  ungroup()

df_grip$z_rpd <- as.numeric(df_grip$z_rpd)

par(mfrow = c(3,2), mar = c(4, 4, 2, 1))
hist(df_grip$z_rpd,
     main = "Distribution of z_rpd (500-1500 ms)",
     xlab = "rpd)",
     breaks = 200)
hist(df_grip$z_rpd_low,
     main = "Distribution of z_rpd_low (0-250 ms)",
     xlab = "rpd_low",
     xlim = c(-4, 6),
     breaks = 200)
hist(df_grip$trial_corr_rpd,
     main = "Distribution of Trial Corrected RPD",
     xlab = "Trial Corrected RPD",
     xlim = c(-1.4, 1.7),  
     breaks = 200)          
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Sample

``` r
# New data frame: 1 row per subject containing sample description data (consistent across conditions)
sample_description_df <- df_grip %>%
  group_by(id) %>%
  summarise(across(everything(), ~ {
    vals <- na.omit(.x) # removes NAs so we only check actual values
    # keeps the column if all non-NA values are identical + returns that consistent value (or NA if all were NA)
    if (length(unique(vals)) <= 1) unique(vals)[1] else NA 
  }), .groups = "drop") %>% 
  # Keep columns that have at least some non-NA values + removes columns that are completely NA across all subjects
  select(where(~ any(!is.na(.))))
sample_description_df$block_baseline_mean <- NULL
```

``` r
# Sample size
nrow(sample_description_df)
```

    ## [1] 144

``` r
table(sample_description_df$group.y)
```

    ## 
    ## ASD CON MHC 
    ##  51  53  40

## Sample Description

``` r
# Function returning group mean + sd
fun_return_descriptives <- function(group){
  group_df <- sample_description_df[
    sample_description_df$group.y == group,
    c("gender", "age", "CBCL", "YSR","z_grip_strength", "id",
      "verbal_IQ", "non_verbal_IQ")] 
  # gender
  n <- length(unique(group_df$id))
  male <- (length(which(group_df$gender == "männlich")))
  female <- (length(which(group_df$gender == "weiblich")))
  gender_f_m <- paste(female,"/",male)
  # age
  age_mean <- round(mean(group_df$age, na.rm = TRUE), digits = 1)
  age_sd <- round(sd(group_df$age, na.rm = TRUE), digits = 1)
  age <- paste(age_mean, "(",age_sd,")" )
  # grip strength
  grip_strength_mean <- round(mean(group_df$z_grip_strength, na.rm = TRUE), digits = 1)
  grip_strength_sd <- round(sd(group_df$z_grip_strength, na.rm = TRUE), digits = 1)
  grip_strength <- paste(grip_strength_mean, "(", grip_strength_sd, ")")
  # verbal IQ
  verbal_IQ_mean <- round(mean(group_df$verbal_IQ, na.rm = TRUE), digits = 1)
  verbal_IQ_sd <- round(sd(group_df$verbal_IQ, na.rm = TRUE), digits = 1)
  verbal_IQ <- paste(verbal_IQ_mean, "(", verbal_IQ_sd, ")")
  # non-verbal IQ
  non_verbal_IQ_mean <- round(mean(group_df$non_verbal_IQ, na.rm = TRUE), digits = 1)
  non_verbal_IQ_sd <- round(sd(group_df$non_verbal_IQ, na.rm = TRUE), digits = 1)
  non_verbal_IQ <- paste(non_verbal_IQ_mean, "(", non_verbal_IQ_sd, ")")

  
  group_description <- data.frame(
    n,
    gender_f_m,
    age,
    verbal_IQ,
    non_verbal_IQ,
    grip_strength)
  t(group_description)
}

asd_description <- fun_return_descriptives(group = "ASD")
colnames(asd_description) <- "ASD"
con_description <- fun_return_descriptives(group = "CON")
colnames(con_description) <- "CON"
mhc_description <- fun_return_descriptives(group = "MHC")
colnames(mhc_description) <- "MHC"
```

``` r
## p-values for descriptive statistics
age_anova <- aov(age ~ group.y, data = sample_description_df)
age_anova_p <- summary(age_anova)[[1]][["Pr(>F)"]][[1]]
grip_strength_anova <- aov(z_grip_strength ~ group.y, data = sample_description_df)
grip_strength_anova_p <- summary(grip_strength_anova)[[1]][["Pr(>F)"]][[1]]
verbal_IQ_anova <- aov(verbal_IQ ~ group.y, data = sample_description_df)
verbal_IQ_anova_p <- summary(verbal_IQ_anova)[[1]][["Pr(>F)"]][[1]]
non_verbal_IQ_anova <- aov(non_verbal_IQ ~ group.y, data = sample_description_df)
non_verbal_IQ_anova_p <- summary(non_verbal_IQ_anova)[[1]][["Pr(>F)"]][[1]]
gender_chi2 <- chisq.test(sample_description_df$gender, sample_description_df$group.y)
gender_chi2_p <- gender_chi2$p.value

p_values <- c(
  NA,
  gender_chi2_p,
  age_anova_p,
  grip_strength_anova_p,
  verbal_IQ_anova_p,
  non_verbal_IQ_anova_p)

p_value <- sapply(p_values, function(p) {
  if (is.na(p)) {
    NA
  } else if (p < 0.001) {
    "< 0.001"
  } else {
    format(p, scientific = F)
    round(p, 3)
  }
})

sample_table <- cbind(asd_description, con_description, mhc_description, p_value)
# Redo table row and column names
rownames(sample_table)[rownames(sample_table) == "grip_strength"] <- "grip strength [z]"
rownames(sample_table)[rownames(sample_table) == "gender_f_m"] <- "gender (f/m)"
rownames(sample_table)[rownames(sample_table) == "verbal_IQ"] <- "verbal IQ"
rownames(sample_table)[rownames(sample_table) == "non_verbal_IQ"] <- "non verbal IQ"
colnames(sample_table)[colnames(sample_table) == "p_value"] <- "p value"

sample_description_table <- knitr::kable(sample_table,
                                         caption = "Sample description",
                                         digits = 4,
                                         align = "c")
sample_description_table
```

|                     |      ASD       |      CON       |      MHC       | p value  |
|:--------------------|:--------------:|:--------------:|:--------------:|:--------:|
| n                   |       51       |       53       |       40       |    NA    |
| gender (f/m)        |    10 / 41     |    29 / 24     |    28 / 12     | \< 0.001 |
| age                 |   15.2 ( 2 )   |   15 ( 1.8 )   |  15.9 ( 1.7 )  |  0.051   |
| verbal IQ           | 97.2 ( 15.6 )  | 105.2 ( 11.3 ) |  101 ( 11.1 )  |  0.006   |
| non verbal IQ       | 101.3 ( 15.1 ) | 107.6 ( 11.3 ) | 103.8 ( 11.8 ) |   0.01   |
| grip strength \[z\] |   -2 ( 1.1 )   |  -1.4 ( 1.4 )  |  -1.2 ( 1.2 )  |  0.048   |

Sample description

``` r
# Grip strength
# 1. Descriptive statistics
sample_table <- df_grip %>%
  group_by(group.y) %>%
  summarise(
    N = n(),
    Mean = mean(grip_strength.y, na.rm = TRUE),
    SD = sd(grip_strength.y, na.rm = TRUE)
  )

# 2. Run ANOVA
anova_result <- aov(grip_strength.y ~ group.y, data = df_grip)
p_val <- summary(anova_result)[[1]][["Pr(>F)"]][1]

# 3. Format caption with p-value
caption_text <- paste0("Descriptive statistics for grip strength by group (ANOVA p = ",
                       formatC(p_val, format = "f", digits = 4), ")")

sample_description_table <- knitr::kable(sample_table,
                                         caption = caption_text,
                                         digits = 4,
                                         align = "c")

# Print table
sample_description_table
```

| group.y |  N   |  Mean   |   SD    |
|:-------:|:----:|:-------:|:-------:|
|   ASD   | 1020 | 38.2240 | 15.7467 |
|   CON   | 1060 | 42.8821 | 15.1819 |
|   MHC   | 840  | 46.5476 | 17.9537 |

Descriptive statistics for grip strength by group (ANOVA p = 0.0000)

## Plots

``` r
## Plot: Age
ggplot(sample_description_df, aes(x = group.y, y = age), col = group.y) +
  geom_boxplot(fill = "grey") +
  geom_jitter() +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(19.5, 20, 19)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("Age") +
  theme(plot.title = element_text(face="bold"))
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
## Plot: IQ-verbal
plot_verbal_iq <- ggplot(sample_description_df, aes(x = group.y, y = verbal_IQ), col = group.y) + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(135, 140, 130)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("IQ-verbal") +
  theme(plot.title = element_text(face="bold")) +
  theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))

## Plot: IQ-non-verbal
plot_non_verbal_iq <- ggplot(sample_description_df, aes(x = group.y, y = non_verbal_IQ), col = "group.y") + 
  geom_boxplot(fill = "grey") +
  geom_signif(comparisons = list(c("ASD", "CON"), c("ASD", "MHC"), c("MHC", "CON")),
              map_signif_level=TRUE,
              y_position = c(135, 140, 130)) +
  theme_bw() + 
  theme_classic() +
  ggtitle("IQ-non-verbal") +
  theme(plot.title = element_text(face="bold")) +
  theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  
grid.arrange(plot_verbal_iq, plot_non_verbal_iq, ncol = 2)
```

    ## Warning: Removed 4 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 4 rows containing non-finite outside the scale range
    ## (`stat_signif()`).

    ## Warning in wilcox.test.default(c(90, 110, 95, 87.5, 112.5, 102.5, 107.5, : kann
    ## bei Bindungen keinen exakten p-Wert Berechnen

    ## Warning: Removed 4 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 4 rows containing non-finite outside the scale range
    ## (`stat_signif()`).

    ## Warning in wilcox.test.default(c(107.5, 115, 92.5, 100, 102.5, 90, 132.5, :
    ## kann bei Bindungen keinen exakten p-Wert Berechnen

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
hist(df_grip$trial_corr_rpd, main = "Histogram of Trial Corrected RPD", xlab = "Trial Corrected RPD", col = "lightblue", border = "black")
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# Compute descriptive stats
mean_value <- mean(df_grip$trial_corr_rpd, na.rm = TRUE)
sd_value <- sd(df_grip$trial_corr_rpd, na.rm = TRUE)
median_value <- median(df_grip$trial_corr_rpd, na.rm = TRUE)
min_value <- min(df_grip$trial_corr_rpd, na.rm = TRUE)
max_value <- max(df_grip$trial_corr_rpd, na.rm = TRUE)

# Create a table
descriptives <- data.frame(
  Statistic = c("Mean", "Standard Deviation", "Median", "Min", "Max"),
  Value = round(c(mean_value, sd_value, median_value, min_value, max_value), 3)
)

# Show as table
knitr::kable(descriptives, caption = "Descriptive Statistics for trial_corr_rpd")
```

| Statistic          |  Value |
|:-------------------|-------:|
| Mean               |  0.002 |
| Standard Deviation |  0.396 |
| Median             | -0.041 |
| Min                | -1.440 |
| Max                |  1.603 |

Descriptive Statistics for trial_corr_rpd

## Covariate Check

``` r
# gender as a factor
df_grip$gender <- as.factor(df_grip$gender)

# Fit a mixed model with covariates
lmm_covariates <- lmer(trial_corr_rpd ~ verbal_IQ + non_verbal_IQ + age + gender + grip_strength.y + (1|id), data = df_grip, REML = FALSE)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
anova_table <- anova(lmm_covariates)
knitr::kable(anova_table, caption = "ANOVA Results")
```

|                 |    Sum Sq |   Mean Sq | NumDF | DenDF |   F value |   Pr(\>F) |
|:----------------|----------:|----------:|------:|------:|----------:|----------:|
| verbal_IQ       | 0.0681421 | 0.0681421 |     1 |  2443 | 0.4365761 | 0.5088420 |
| non_verbal_IQ   | 0.0140252 | 0.0140252 |     1 |  2443 | 0.0898571 | 0.7643844 |
| age             | 0.0939354 | 0.0939354 |     1 |  2443 | 0.6018295 | 0.4379559 |
| gender          | 0.0018643 | 0.0018643 |     1 |  2443 | 0.0119444 | 0.9129811 |
| grip_strength.y | 0.3545530 | 0.3545530 |     1 |  2443 | 2.2715668 | 0.1318959 |

ANOVA Results

## Linear Mixed-Effects Model: Effect of Experimental Manipulation

``` r
# Result 1:
df_grip$trial_phase <- factor(
  df_grip$trial_phase,
  levels = c("baseline_pre_squeeze", "squeeze", "baseline_post_squeeze", "relax")
)

lmm_baseline <- lmer(
  trial_corr_rpd ~ trial_phase + (1 | id),
  data = df_grip
)
anova(lmm_baseline)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Sum Sq"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Mean Sq"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["NumDF"],"name":[3],"type":["int"],"align":["right"]},{"label":["DenDF"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["F value"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Pr(>F)"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"120.7409","2":"40.24696","3":"3","4":"2438.447","5":"370.476","6":"3.081842e-198","_rn_":"trial_phase"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
anova_table <- anova(lmm_baseline)
knitr::kable(anova_table, caption = "ANOVA Results")
```

|             |   Sum Sq |  Mean Sq | NumDF |    DenDF | F value | Pr(\>F) |
|:------------|---------:|---------:|------:|---------:|--------:|--------:|
| trial_phase | 120.7409 | 40.24696 |     3 | 2438.447 | 370.476 |       0 |

ANOVA Results

``` r
r2_nakagawa(lmm_baseline)  
```

    ## # R2 for Mixed Models
    ## 
    ##   Conditional R2: 0.309
    ##      Marginal R2: 0.304

``` r
emmeans(lmm_baseline, ~ trial_phase)
```

    ##  trial_phase            emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.0199 0.0132 1290 -0.00603   0.0459
    ##  squeeze                0.3159 0.0124 1144  0.29157   0.3403
    ##  baseline_post_squeeze -0.2182 0.0160 1589 -0.24952  -0.1869
    ##  relax                 -0.1932 0.0124 1145 -0.21761  -0.1689
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

### Post Hoc

``` r
contrast(emmeans(lmm_baseline, ~ trial_phase), "pairwise")
```

    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.2960 0.0179 2409 -16.553
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2381 0.0205 2483  11.615
    ##  baseline_pre_squeeze - relax                   0.2132 0.0179 2409  11.919
    ##  squeeze - baseline_post_squeeze                0.5341 0.0200 2488  26.727
    ##  squeeze - relax                                0.5092 0.0173 2387  29.463
    ##  baseline_post_squeeze - relax                 -0.0249 0.0200 2487  -1.248
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.5962
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
confint(contrast(emmeans(lmm_baseline, ~ trial_phase), "pairwise"))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["chr"],"align":["left"]},{"label":["estimate"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"baseline_pre_squeeze - squeeze","2":"-0.29598100","3":"0.01788069","4":"2408.773","5":"-0.34194914","6":"-0.25001285","_rn_":"1"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"0.23813480","3":"0.02050294","4":"2482.686","5":"0.18542642","6":"0.29084318","_rn_":"2"},{"1":"baseline_pre_squeeze - relax","2":"0.21318829","3":"0.01788663","4":"2409.425","5":"0.16720490","6":"0.25917168","_rn_":"3"},{"1":"squeeze - baseline_post_squeeze","2":"0.53411580","3":"0.01998395","4":"2488.239","5":"0.48274169","6":"0.58548991","_rn_":"4"},{"1":"squeeze - relax","2":"0.50916928","3":"0.01728186","4":"2387.080","5":"0.46474036","6":"0.55359821","_rn_":"5"},{"1":"baseline_post_squeeze - relax","2":"-0.02494651","3":"0.01998851","4":"2487.062","5":"-0.07633237","6":"0.02643934","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
emmeans_result <- emmeans(lmm_baseline, ~ trial_phase)  
plot_data <- as.data.frame(emmeans_result)

# Ensure correct order of trial phases
plot_data$trial_phase <- factor(
  plot_data$trial_phase,
  levels = c("baseline_pre_squeeze", "squeeze", "baseline_post_squeeze", "relax")
)

# plot
emm_plot <- ggplot(plot_data) +
  geom_crossbar(aes(
    x = trial_phase, y = emmean,
    ymin = emmean - SE, ymax = emmean + SE),
    fill = "gray70", color = "black",
    width = 0.3, alpha = 0.8) +
  geom_errorbar(aes(
    x = trial_phase, ymin = lower.CL, ymax = upper.CL),
    width = 0.2, color = "black") +
  scale_x_discrete(labels = c(
    "baseline_pre_squeeze" = "Baseline Pre Grip",
    "squeeze" = "Grip",
    "baseline_post_squeeze" = "Baseline Post Grip",
    "relax" = "Relax"
  )) +
  theme_bw() +
  labs(
    title = "Estimated Marginal Pupil Response by Trial Phase",
    x = "Trial Phase",
    y = "Estimated Marginal Pupil Response"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display plot
print(emm_plot)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Effort Check: Grip strength

``` r
df_grip %>%
  filter(trial_phase == "squeeze") %>%
  summarise(
    mean_grip = mean(grip_strength.y, na.rm = TRUE),
    sd_grip = sd(grip_strength.y, na.rm = TRUE)
  )
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["mean_grip"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["sd_grip"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"42.3361","2":"16.56555"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
lmm_grip <- lmer(trial_corr_rpd ~ grip_strength.y + (1 | id),
              data = df_grip)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
anova_table <- anova(lmm_grip)
knitr::kable(anova_table, caption = "ANOVA Results")
```

|                 |    Sum Sq |   Mean Sq | NumDF | DenDF |  F value |   Pr(\>F) |
|:----------------|----------:|----------:|------:|------:|---------:|----------:|
| grip_strength.y | 0.3580207 | 0.3580207 |     1 |  2508 | 2.277853 | 0.1313594 |

ANOVA Results

``` r
lmm_gripstrength <- lmer(
  trial_corr_rpd ~ trial_phase * grip_strength.y + (1 | id),
  data = df_grip 
)
anova(lmm_gripstrength)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Sum Sq"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Mean Sq"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["NumDF"],"name":[3],"type":["int"],"align":["right"]},{"label":["DenDF"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["F value"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Pr(>F)"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"4.971235","2":"1.657078","3":"3","4":"2413.9406","5":"15.634866","6":"4.518163e-10","_rn_":"trial_phase"},{"1":"0.140108","2":"0.140108","3":"1","4":"268.8517","5":"1.321947","6":"2.512647e-01","_rn_":"grip_strength.y"},{"1":"7.760945","2":"2.586982","3":"3","4":"2418.8881","5":"24.408690","6":"1.475135e-15","_rn_":"trial_phase:grip_strength.y"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
emtrends(lmm_gripstrength, ~ trial_phase, var = 'grip_strength.y', infer=c(TRUE, TRUE))
```

    ##  trial_phase           grip_strength.y.trend       SE   df lower.CL  upper.CL
    ##  baseline_pre_squeeze              -2.10e-03 0.000809 1668 -0.00369 -0.000518
    ##  squeeze                            5.86e-03 0.000743 1516  0.00441  0.007320
    ##  baseline_post_squeeze             -8.88e-05 0.001020 2061 -0.00209  0.001917
    ##  relax                             -1.68e-03 0.000745 1532 -0.00314 -0.000220
    ##  t.ratio p.value
    ##   -2.602  0.0093
    ##    7.890  <.0001
    ##   -0.087  0.9309
    ##   -2.257  0.0242
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

## Group differences

``` r
# Result 2:
df_grip$group.y <- as.factor(df_grip$group.y)

lmm_group <- lmer(
  trial_corr_rpd ~ trial_phase * group.y + (1 | id),
  data = df_grip
)
anova(lmm_group)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Sum Sq"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Mean Sq"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["NumDF"],"name":[3],"type":["int"],"align":["right"]},{"label":["DenDF"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["F value"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Pr(>F)"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"117.9422619","2":"39.3140873","3":"3","4":"2431.1943","5":"364.508410","6":"1.817443e-195","_rn_":"trial_phase"},{"1":"0.2747287","2":"0.1373643","3":"2","4":"142.8635","5":"1.273601","6":"2.829797e-01","_rn_":"group.y"},{"1":"2.5042515","2":"0.4173752","3":"6","4":"2430.4815","5":"3.869778","6":"7.548940e-04","_rn_":"trial_phase:group.y"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
anova_table <- anova(lmm_group)
knitr::kable(anova_table, caption = "ANOVA Results")
```

|  | Sum Sq | Mean Sq | NumDF | DenDF | F value | Pr(\>F) |
|:---|---:|---:|---:|---:|---:|---:|
| trial_phase | 117.9422619 | 39.3140873 | 3 | 2431.1943 | 364.508410 | 0.0000000 |
| group.y | 0.2747287 | 0.1373643 | 2 | 142.8635 | 1.273601 | 0.2829797 |
| trial_phase:group.y | 2.5042515 | 0.4173752 | 6 | 2430.4815 | 3.869778 | 0.0007549 |

ANOVA Results

``` r
r2_nakagawa(lmm_group)
```

    ## # R2 for Mixed Models
    ## 
    ##   Conditional R2: 0.315
    ##      Marginal R2: 0.310

``` r
emmeans(lmm_group, ~ group.y | trial_phase)
```

    ## trial_phase = baseline_pre_squeeze:
    ##  group.y   emmean     SE   df  lower.CL upper.CL
    ##  ASD      0.04326 0.0225 1305 -0.000841   0.0874
    ##  CON      0.01036 0.0217 1291 -0.032226   0.0529
    ##  MHC      0.00421 0.0248 1198 -0.044346   0.0528
    ## 
    ## trial_phase = squeeze:
    ##  group.y   emmean     SE   df  lower.CL upper.CL
    ##  ASD      0.27527 0.0210 1162  0.234093   0.3165
    ##  CON      0.38139 0.0205 1153  0.341145   0.4216
    ##  MHC      0.28218 0.0231 1042  0.236893   0.3275
    ## 
    ## trial_phase = baseline_post_squeeze:
    ##  group.y   emmean     SE   df  lower.CL upper.CL
    ##  ASD     -0.22508 0.0256 1581 -0.275281  -0.1749
    ##  CON     -0.19134 0.0274 1650 -0.245160  -0.1375
    ##  MHC     -0.24124 0.0303 1452 -0.300710  -0.1818
    ## 
    ## trial_phase = relax:
    ##  group.y   emmean     SE   df  lower.CL upper.CL
    ##  ASD     -0.16057 0.0209 1158 -0.201676  -0.1195
    ##  CON     -0.22969 0.0205 1157 -0.270006  -0.1894
    ##  MHC     -0.18702 0.0231 1047 -0.232414  -0.1416
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

### Post Hoc

``` r
# pairwise phase comparison within each group (not about group differences)
contrast(emmeans(lmm_group, ~ trial_phase | group.y), "pairwise")
```

    ## group.y = ASD:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.2320 0.0303 2408  -7.660
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2683 0.0336 2452   7.978
    ##  baseline_pre_squeeze - relax                   0.2038 0.0303 2411   6.736
    ##  squeeze - baseline_post_squeeze                0.5004 0.0327 2453  15.316
    ##  squeeze - relax                                0.4358 0.0292 2382  14.941
    ##  baseline_post_squeeze - relax                 -0.0645 0.0326 2449  -1.976
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.1973
    ## 
    ## group.y = CON:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.3710 0.0294 2397 -12.619
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2017 0.0346 2480   5.834
    ##  baseline_pre_squeeze - relax                   0.2400 0.0294 2395   8.157
    ##  squeeze - baseline_post_squeeze                0.5727 0.0339 2496  16.916
    ##  squeeze - relax                                0.6111 0.0286 2381  21.398
    ##  baseline_post_squeeze - relax                  0.0384 0.0339 2495   1.132
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.6697
    ## 
    ## group.y = MHC:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.2780 0.0333 2403  -8.354
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2454 0.0386 2494   6.356
    ##  baseline_pre_squeeze - relax                   0.1912 0.0333 2404   5.741
    ##  squeeze - baseline_post_squeeze                0.5234 0.0376 2492  13.931
    ##  squeeze - relax                                0.4692 0.0321 2381  14.622
    ##  baseline_post_squeeze - relax                 -0.0542 0.0376 2492  -1.442
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.4734
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
confint(contrast(emmeans(lmm_group, ~ trial_phase | group.y), "pairwise"))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["group.y"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"baseline_pre_squeeze - squeeze","2":"ASD","3":"-0.23201272","4":"0.03028792","5":"2407.590","6":"-0.3098777","7":"-0.15414774","_rn_":"1"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"ASD","3":"0.26833907","4":"0.03363551","5":"2452.212","6":"0.1818691","7":"0.35480905","_rn_":"2"},{"1":"baseline_pre_squeeze - relax","2":"ASD","3":"0.20383489","4":"0.03026206","5":"2411.295","6":"0.1260365","7":"0.28163331","_rn_":"3"},{"1":"squeeze - baseline_post_squeeze","2":"ASD","3":"0.50035179","4":"0.03266886","5":"2452.617","6":"0.4163669","7":"0.58433670","_rn_":"4"},{"1":"squeeze - relax","2":"ASD","3":"0.43584761","4":"0.02917136","5":"2381.495","6":"0.3608525","7":"0.51084270","_rn_":"5"},{"1":"baseline_post_squeeze - relax","2":"ASD","3":"-0.06450418","4":"0.03264088","5":"2449.273","6":"-0.1484172","7":"0.01940887","_rn_":"6"},{"1":"baseline_pre_squeeze - squeeze","2":"CON","3":"-0.37103109","4":"0.02940173","5":"2396.734","6":"-0.4466181","7":"-0.29544409","_rn_":"7"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"CON","3":"0.20169515","4":"0.03457481","5":"2480.475","6":"0.1128111","7":"0.29057917","_rn_":"8"},{"1":"baseline_pre_squeeze - relax","2":"CON","3":"0.24004610","4":"0.02942707","5":"2394.759","6":"0.1643939","7":"0.31569828","_rn_":"9"},{"1":"squeeze - baseline_post_squeeze","2":"CON","3":"0.57272624","4":"0.03385662","5":"2495.584","6":"0.4856889","7":"0.65976359","_rn_":"10"},{"1":"squeeze - relax","2":"CON","3":"0.61107719","4":"0.02855790","5":"2380.519","6":"0.5376592","7":"0.68449518","_rn_":"11"},{"1":"baseline_post_squeeze - relax","2":"CON","3":"0.03835095","4":"0.03387847","5":"2494.693","6":"-0.0487426","7":"0.12544450","_rn_":"12"},{"1":"baseline_pre_squeeze - squeeze","2":"MHC","3":"-0.27796427","4":"0.03327227","5":"2403.482","6":"-0.3635016","7":"-0.19242692","_rn_":"13"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"MHC","3":"0.24544994","4":"0.03861696","5":"2493.690","6":"0.1461748","7":"0.34472505","_rn_":"14"},{"1":"baseline_pre_squeeze - relax","2":"MHC","3":"0.19123684","4":"0.03330923","5":"2403.536","6":"0.1056045","7":"0.27686920","_rn_":"15"},{"1":"squeeze - baseline_post_squeeze","2":"MHC","3":"0.52341421","4":"0.03757257","5":"2491.587","6":"0.4268239","7":"0.62000450","_rn_":"16"},{"1":"squeeze - relax","2":"MHC","3":"0.46920111","4":"0.03208842","5":"2380.607","6":"0.3867067","7":"0.55169554","_rn_":"17"},{"1":"baseline_post_squeeze - relax","2":"MHC","3":"-0.05421309","4":"0.03760630","5":"2492.372","6":"-0.1508901","7":"0.04246390","_rn_":"18"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
# Pairwise group comparisons within each trial_phase
contrast(emmeans(lmm_group, ~ group.y | trial_phase), "pairwise")
```

    ## trial_phase = baseline_pre_squeeze:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  ASD - CON  0.03290 0.0312 1298   1.053  0.5436
    ##  ASD - MHC  0.03905 0.0334 1245   1.168  0.4727
    ##  CON - MHC  0.00614 0.0329 1237   0.187  0.9810
    ## 
    ## trial_phase = squeeze:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  ASD - CON -0.10611 0.0293 1158  -3.616  0.0009
    ##  ASD - MHC -0.00690 0.0312 1095  -0.221  0.9734
    ##  CON - MHC  0.09921 0.0309 1089   3.213  0.0039
    ## 
    ## trial_phase = baseline_post_squeeze:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  ASD - CON -0.03374 0.0375 1618  -0.899  0.6409
    ##  ASD - MHC  0.01616 0.0397 1505   0.407  0.9126
    ##  CON - MHC  0.04990 0.0409 1539   1.220  0.4413
    ## 
    ## trial_phase = relax:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  ASD - CON  0.06912 0.0293 1157   2.355  0.0489
    ##  ASD - MHC  0.02645 0.0312 1095   0.848  0.6735
    ##  CON - MHC -0.04267 0.0309 1094  -1.379  0.3523
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
confint(contrast(emmeans(lmm_group, ~ group.y | trial_phase), "pairwise"))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["trial_phase"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"ASD - CON","2":"baseline_pre_squeeze","3":"0.032904751","4":"0.03124842","5":"1297.936","6":"-0.0404164549","7":"0.10622596","_rn_":"1"},{"1":"ASD - MHC","2":"baseline_pre_squeeze","3":"0.039047767","4":"0.03343532","5":"1244.972","6":"-0.0394086125","7":"0.11750415","_rn_":"2"},{"1":"CON - MHC","2":"baseline_pre_squeeze","3":"0.006143015","4":"0.03291935","5":"1237.197","6":"-0.0711032342","7":"0.08338927","_rn_":"3"},{"1":"ASD - CON","2":"squeeze","3":"-0.106113616","4":"0.02934567","5":"1157.646","6":"-0.1749798201","7":"-0.03724741","_rn_":"4"},{"1":"ASD - MHC","2":"squeeze","3":"-0.006903782","4":"0.03119439","5":"1094.847","6":"-0.0801138540","7":"0.06630629","_rn_":"5"},{"1":"CON - MHC","2":"squeeze","3":"0.099209834","4":"0.03087436","5":"1089.436","6":"0.0267503482","7":"0.17166932","_rn_":"6"},{"1":"ASD - CON","2":"baseline_post_squeeze","3":"-0.033739170","4":"0.03752312","5":"1617.964","6":"-0.1217632834","7":"0.05428494","_rn_":"7"},{"1":"ASD - MHC","2":"baseline_post_squeeze","3":"0.016158635","4":"0.03967673","5":"1504.778","6":"-0.0769240211","7":"0.10924129","_rn_":"8"},{"1":"CON - MHC","2":"baseline_post_squeeze","3":"0.049897805","4":"0.04089135","5":"1538.947","6":"-0.0460322805","7":"0.14582789","_rn_":"9"},{"1":"ASD - CON","2":"relax","3":"0.069115962","4":"0.02934333","5":"1157.404","6":"0.0002552369","7":"0.13797669","_rn_":"10"},{"1":"ASD - MHC","2":"relax","3":"0.026449720","4":"0.03120727","5":"1095.415","6":"-0.0467905058","7":"0.09968995","_rn_":"11"},{"1":"CON - MHC","2":"relax","3":"-0.042666241","4":"0.03093982","5":"1094.283","6":"-0.1152789129","7":"0.02994643","_rn_":"12"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
# Compute estimated marginal means
emm_result <- emmeans(lmm_group, ~ group.y | trial_phase)
emm_data <- as.data.frame(emm_result)

# Ensure correct phase order
emm_data$trial_phase <- factor(
  emm_data$trial_phase,
  levels = c("baseline_pre_squeeze", "squeeze", "baseline_post_squeeze", "relax")
)

# Plot estimated marginal means with confidence intervals
ggplot(emm_data, aes(x = trial_phase, y = emmean, color = group.y)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2,
    size = 1,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ group.y) +
  scale_x_discrete(labels = c(
    "baseline_pre_squeeze" = "Baseline Pre Grip",
    "squeeze" = "Grip",
    "baseline_post_squeeze" = "Baseline Post Grip",
    "relax" = "Relax"
  )) +
  labs(
    title = "Estimated Marginal Pupil Response by Group and Trial Phase",
    x = "Trial Phase",
    y = "Estimated Marginal Pupil Response"
  ) +
  theme_bw() +
  scale_color_manual(values = c("ASD" = "red", "CON" = "#045D5D", "MHC" = "#A05000")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    text = element_text(size = 14)
  )
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
# Plot estimated marginal means by phase (facets), color = group
ggplot(emm_data, aes(x = group.y, y = emmean, color = group.y)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2,
    size = 1,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ trial_phase) +
  scale_color_manual(values = c(
    "ASD" = "#9900CC",   # bright purple
    "CON" = "#00FF00",   # bright green
    "MHC" = "#00FFFF"    # bright cyan
  )) +
  labs(
    title = "Estimated Marginal Pupil Response by Trial Phase and Group",
    x = "Group",
    y = "Estimated Marginal Pupil Response"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    text = element_text(size = 14)
  )
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

## Across Manipulation trials - Habituation Effects

``` r
# Result 3:
df_grip$manipulation_trial <- as.factor(df_grip$manipulation_trial)

lmm_trials <- lmer(
  trial_corr_rpd ~ manipulation_trial * trial_phase + (1 | id),
  data = df_grip
)
anova(lmm_trials)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Sum Sq"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Mean Sq"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["NumDF"],"name":[3],"type":["int"],"align":["right"]},{"label":["DenDF"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["F value"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Pr(>F)"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"8.313045","2":"2.078261","3":"4","4":"2398.186","5":"20.402571","6":"1.515171e-16","_rn_":"manipulation_trial"},{"1":"119.862174","2":"39.954058","3":"3","4":"2418.510","5":"392.234370","6":"1.399885e-207","_rn_":"trial_phase"},{"1":"6.273156","2":"0.522763","3":"12","4":"2392.232","5":"5.132035","6":"1.525822e-08","_rn_":"manipulation_trial:trial_phase"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
anova_table <- anova(lmm_trials)
knitr::kable(anova_table, caption = "ANOVA Results")
```

|  | Sum Sq | Mean Sq | NumDF | DenDF | F value | Pr(\>F) |
|:---|---:|---:|---:|---:|---:|---:|
| manipulation_trial | 8.313045 | 2.078261 | 4 | 2398.186 | 20.402571 | 0 |
| trial_phase | 119.862174 | 39.954058 | 3 | 2418.510 | 392.234370 | 0 |
| manipulation_trial:trial_phase | 6.273156 | 0.522763 | 12 | 2392.232 | 5.132035 | 0 |

ANOVA Results

``` r
r2_nakagawa(lmm_trials) 
```

    ## # R2 for Mixed Models
    ## 
    ##   Conditional R2: 0.354
    ##      Marginal R2: 0.347

### Post Hoc

``` r
# Compare trial phases within each manipulation trial (1 to 5)
contrast(emmeans(lmm_trials, ~ trial_phase | manipulation_trial), "pairwise")
```

    ## manipulation_trial = 1:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.2120 0.0391 2388  -5.427
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.0510 0.0460 2434   1.108
    ##  baseline_pre_squeeze - relax                   0.1527 0.0390 2389   3.917
    ##  squeeze - baseline_post_squeeze                0.2630 0.0447 2430   5.883
    ##  squeeze - relax                                0.3647 0.0374 2370   9.746
    ##  baseline_post_squeeze - relax                  0.1017 0.0446 2431   2.278
    ##  p.value
    ##   <.0001
    ##   0.6843
    ##   0.0005
    ##   <.0001
    ##   <.0001
    ##   0.1034
    ## 
    ## manipulation_trial = 2:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.3569 0.0383 2380  -9.328
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.3114 0.0428 2415   7.277
    ##  baseline_pre_squeeze - relax                   0.2068 0.0383 2380   5.404
    ##  squeeze - baseline_post_squeeze                0.6683 0.0420 2414  15.918
    ##  squeeze - relax                                0.5637 0.0374 2370  15.091
    ##  baseline_post_squeeze - relax                 -0.1046 0.0420 2414  -2.491
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.0615
    ## 
    ## manipulation_trial = 3:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.3590 0.0385 2383  -9.325
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2519 0.0437 2424   5.767
    ##  baseline_pre_squeeze - relax                   0.1963 0.0386 2384   5.089
    ##  squeeze - baseline_post_squeeze                0.6109 0.0427 2419  14.314
    ##  squeeze - relax                                0.5553 0.0374 2371  14.839
    ##  baseline_post_squeeze - relax                 -0.0556 0.0427 2418  -1.302
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.5618
    ## 
    ## manipulation_trial = 4:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.3232 0.0386 2384  -8.376
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2828 0.0437 2421   6.465
    ##  baseline_pre_squeeze - relax                   0.2183 0.0387 2384   5.639
    ##  squeeze - baseline_post_squeeze                0.6060 0.0427 2419  14.199
    ##  squeeze - relax                                0.5414 0.0375 2371  14.443
    ##  baseline_post_squeeze - relax                 -0.0646 0.0428 2419  -1.509
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.4321
    ## 
    ## manipulation_trial = 5:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze                -0.2286 0.0392 2389  -5.827
    ##  baseline_pre_squeeze - baseline_post_squeeze   0.2792 0.0460 2434   6.070
    ##  baseline_pre_squeeze - relax                   0.2899 0.0392 2390   7.401
    ##  squeeze - baseline_post_squeeze                0.5078 0.0445 2431  11.405
    ##  squeeze - relax                                0.5185 0.0374 2370  13.855
    ##  baseline_post_squeeze - relax                  0.0107 0.0445 2430   0.240
    ##  p.value
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.9951
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
confint(contrast(emmeans(lmm_trials, ~ trial_phase | manipulation_trial), "pairwise"))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["manipulation_trial"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"baseline_pre_squeeze - squeeze","2":"1","3":"-0.21195431","4":"0.03905529","5":"2387.772","6":"-0.31235927","7":"-0.111549351","_rn_":"1"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"1","3":"0.05101541","4":"0.04602277","5":"2433.853","6":"-0.06730027","7":"0.169331090","_rn_":"2"},{"1":"baseline_pre_squeeze - relax","2":"1","3":"0.15274097","4":"0.03899428","5":"2388.517","6":"0.05249287","7":"0.252989062","_rn_":"3"},{"1":"squeeze - baseline_post_squeeze","2":"1","3":"0.26296972","4":"0.04470126","5":"2430.215","6":"0.14805128","7":"0.377888170","_rn_":"4"},{"1":"squeeze - relax","2":"1","3":"0.36469528","4":"0.03741981","5":"2370.503","6":"0.26849438","7":"0.460896178","_rn_":"5"},{"1":"baseline_post_squeeze - relax","2":"1","3":"0.10172556","4":"0.04464800","5":"2430.784","6":"-0.01305595","7":"0.216507057","_rn_":"6"},{"1":"baseline_pre_squeeze - squeeze","2":"2","3":"-0.35693693","4":"0.03826615","5":"2380.407","6":"-0.45531336","7":"-0.258560507","_rn_":"7"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"2","3":"0.31137507","4":"0.04279152","5":"2415.428","6":"0.20136575","7":"0.421384385","_rn_":"8"},{"1":"baseline_pre_squeeze - relax","2":"2","3":"0.20678755","4":"0.03826615","5":"2380.407","6":"0.10841112","7":"0.305163970","_rn_":"9"},{"1":"squeeze - baseline_post_squeeze","2":"2","3":"0.66831200","4":"0.04198380","5":"2413.905","6":"0.56037912","7":"0.776244875","_rn_":"10"},{"1":"squeeze - relax","2":"2","3":"0.56372448","4":"0.03735478","5":"2369.698","6":"0.46769074","7":"0.659758215","_rn_":"11"},{"1":"baseline_post_squeeze - relax","2":"2","3":"-0.10458752","4":"0.04198380","5":"2413.905","6":"-0.21252040","7":"0.003345354","_rn_":"12"},{"1":"baseline_pre_squeeze - squeeze","2":"3","3":"-0.35901252","4":"0.03849882","5":"2382.876","6":"-0.45798703","7":"-0.260038007","_rn_":"13"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"3","3":"0.25189501","4":"0.04368237","5":"2423.715","6":"0.13959572","7":"0.364194291","_rn_":"14"},{"1":"baseline_pre_squeeze - relax","2":"3","3":"0.19625697","4":"0.03856208","5":"2383.760","6":"0.09711984","7":"0.295394095","_rn_":"15"},{"1":"squeeze - baseline_post_squeeze","2":"3","3":"0.61090752","4":"0.04267986","5":"2418.865","6":"0.50118536","7":"0.720629693","_rn_":"16"},{"1":"squeeze - relax","2":"3","3":"0.55526949","4":"0.03741981","5":"2370.518","6":"0.45906858","7":"0.651470391","_rn_":"17"},{"1":"baseline_post_squeeze - relax","2":"3","3":"-0.05563804","4":"0.04273558","5":"2418.221","6":"-0.16550347","7":"0.054227390","_rn_":"18"},{"1":"baseline_pre_squeeze - squeeze","2":"4","3":"-0.32315189","4":"0.03857849","5":"2383.889","6":"-0.42233119","7":"-0.223972587","_rn_":"19"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"4","3":"0.28284410","4":"0.04374910","5":"2421.171","6":"0.17037320","7":"0.395314994","_rn_":"20"},{"1":"baseline_pre_squeeze - relax","2":"4","3":"0.21826474","4":"0.03870403","5":"2383.923","6":"0.11876269","7":"0.317766795","_rn_":"21"},{"1":"squeeze - baseline_post_squeeze","2":"4","3":"0.60599598","4":"0.04267984","5":"2418.788","6":"0.49627387","7":"0.715718098","_rn_":"22"},{"1":"squeeze - relax","2":"4","3":"0.54141663","4":"0.03748563","5":"2371.329","6":"0.44504654","7":"0.637786723","_rn_":"23"},{"1":"baseline_post_squeeze - relax","2":"4","3":"-0.06457935","4":"0.04279389","5":"2419.112","6":"-0.17459467","7":"0.045435961","_rn_":"24"},{"1":"baseline_pre_squeeze - squeeze","2":"5","3":"-0.22857099","4":"0.03922815","5":"2389.119","6":"-0.32942031","7":"-0.127721672","_rn_":"25"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"5","3":"0.27920110","4":"0.04599361","5":"2434.289","6":"0.16096041","7":"0.397441798","_rn_":"26"},{"1":"baseline_pre_squeeze - relax","2":"5","3":"0.28988380","4":"0.03916743","5":"2389.864","6":"0.18919061","7":"0.390576985","_rn_":"27"},{"1":"squeeze - baseline_post_squeeze","2":"5","3":"0.50777209","4":"0.04452104","5":"2430.963","6":"0.39331699","7":"0.622227196","_rn_":"28"},{"1":"squeeze - relax","2":"5","3":"0.51845478","4":"0.03741981","5":"2370.503","6":"0.42225388","7":"0.614655684","_rn_":"29"},{"1":"baseline_post_squeeze - relax","2":"5","3":"0.01068269","4":"0.04446544","5":"2429.886","6":"-0.10362951","7":"0.124994896","_rn_":"30"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
# Get estimated marginal means
emms <- emmeans(lmm_trials, ~ trial_phase | manipulation_trial)

# Convert to a dataframe
emms_df <- as.data.frame(emms)

# Relabel trial_phase and manipulation_trial
emms_df <- emms_df %>%
  mutate(
    trial_phase = factor(trial_phase,
                         levels = c("baseline_pre_squeeze", "squeeze", "baseline_post_squeeze", "relax"),
                         labels = c("Baseline Pre-Grip", "Grip", "Baseline Post-Grip", "Relax")),
    trial = as.numeric(as.character(manipulation_trial))
  )

# Plot
ggplot(emms_df, aes(x = trial, 
                    y = emmean, 
                    color = trial_phase,
                    group = trial_phase)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  labs(
    title = "Estimated Marginal Means of RPD by Trial Phase and Trial",
    x = "Trial",
    y = "Estimated Relative Pupil Dilation (RPD)",
    color = "Trial Phase"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
# Compare manipulation trial (1 to 5) within each trial phase
contrast(emmeans(lmm_trials, ~ manipulation_trial | trial_phase), "pairwise")
```

    ## trial_phase = baseline_pre_squeeze:
    ##  contrast                                  estimate     SE   df t.ratio p.value
    ##  manipulation_trial1 - manipulation_trial2 -0.18845 0.0399 2394  -4.727  <.0001
    ##  manipulation_trial1 - manipulation_trial3 -0.17041 0.0401 2396  -4.251  0.0002
    ##  manipulation_trial1 - manipulation_trial4 -0.19730 0.0402 2398  -4.912  <.0001
    ##  manipulation_trial1 - manipulation_trial5 -0.23991 0.0407 2406  -5.889  <.0001
    ##  manipulation_trial2 - manipulation_trial3  0.01804 0.0394 2390   0.458  0.9909
    ##  manipulation_trial2 - manipulation_trial4 -0.00885 0.0395 2392  -0.224  0.9994
    ##  manipulation_trial2 - manipulation_trial5 -0.05146 0.0400 2393  -1.285  0.7003
    ##  manipulation_trial3 - manipulation_trial4 -0.02690 0.0397 2390  -0.678  0.9613
    ##  manipulation_trial3 - manipulation_trial5 -0.06950 0.0403 2396  -1.726  0.4178
    ##  manipulation_trial4 - manipulation_trial5 -0.04260 0.0403 2392  -1.056  0.8287
    ## 
    ## trial_phase = squeeze:
    ##  contrast                                  estimate     SE   df t.ratio p.value
    ##  manipulation_trial1 - manipulation_trial2 -0.33343 0.0374 2370  -8.911  <.0001
    ##  manipulation_trial1 - manipulation_trial3 -0.31747 0.0374 2370  -8.484  <.0001
    ##  manipulation_trial1 - manipulation_trial4 -0.30850 0.0374 2370  -8.244  <.0001
    ##  manipulation_trial1 - manipulation_trial5 -0.25652 0.0375 2371  -6.843  <.0001
    ##  manipulation_trial2 - manipulation_trial3  0.01597 0.0374 2370   0.427  0.9930
    ##  manipulation_trial2 - manipulation_trial4  0.02493 0.0374 2370   0.667  0.9633
    ##  manipulation_trial2 - manipulation_trial5  0.07691 0.0374 2370   2.055  0.2402
    ##  manipulation_trial3 - manipulation_trial4  0.00897 0.0374 2370   0.240  0.9993
    ##  manipulation_trial3 - manipulation_trial5  0.06094 0.0374 2370   1.629  0.4791
    ##  manipulation_trial4 - manipulation_trial5  0.05198 0.0374 2370   1.389  0.6348
    ## 
    ## trial_phase = baseline_post_squeeze:
    ##  contrast                                  estimate     SE   df t.ratio p.value
    ##  manipulation_trial1 - manipulation_trial2  0.07191 0.0486 2431   1.481  0.5751
    ##  manipulation_trial1 - manipulation_trial3  0.03047 0.0492 2438   0.620  0.9720
    ##  manipulation_trial1 - manipulation_trial4  0.03453 0.0492 2436   0.702  0.9560
    ##  manipulation_trial1 - manipulation_trial5 -0.01172 0.0507 2437  -0.231  0.9994
    ##  manipulation_trial2 - manipulation_trial3 -0.04144 0.0468 2420  -0.886  0.9020
    ##  manipulation_trial2 - manipulation_trial4 -0.03738 0.0468 2432  -0.799  0.9309
    ##  manipulation_trial2 - manipulation_trial5 -0.08363 0.0484 2433  -1.728  0.4168
    ##  manipulation_trial3 - manipulation_trial4  0.00405 0.0474 2414   0.086  1.0000
    ##  manipulation_trial3 - manipulation_trial5 -0.04219 0.0490 2429  -0.861  0.9110
    ##  manipulation_trial4 - manipulation_trial5 -0.04625 0.0490 2431  -0.944  0.8797
    ## 
    ## trial_phase = relax:
    ##  contrast                                  estimate     SE   df t.ratio p.value
    ##  manipulation_trial1 - manipulation_trial2 -0.13440 0.0374 2370  -3.598  0.0030
    ##  manipulation_trial1 - manipulation_trial3 -0.12689 0.0374 2371  -3.391  0.0064
    ##  manipulation_trial1 - manipulation_trial4 -0.13178 0.0375 2371  -3.515  0.0041
    ##  manipulation_trial1 - manipulation_trial5 -0.10276 0.0374 2370  -2.751  0.0472
    ##  manipulation_trial2 - manipulation_trial3  0.00751 0.0374 2371   0.201  0.9996
    ##  manipulation_trial2 - manipulation_trial4  0.00262 0.0375 2371   0.070  1.0000
    ##  manipulation_trial2 - manipulation_trial5  0.03164 0.0374 2370   0.847  0.9158
    ##  manipulation_trial3 - manipulation_trial4 -0.00489 0.0376 2372  -0.130  0.9999
    ##  manipulation_trial3 - manipulation_trial5  0.02413 0.0374 2371   0.645  0.9676
    ##  manipulation_trial4 - manipulation_trial5  0.02901 0.0375 2371   0.774  0.9381
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 5 estimates

``` r
confint(contrast(emmeans(lmm_trials, ~ manipulation_trial| trial_phase), "pairwise"))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["trial_phase"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"manipulation_trial1 - manipulation_trial2","2":"baseline_pre_squeeze","3":"-0.188449428","4":"0.03986503","5":"2393.836","6":"-0.29727521","7":"-0.0796236437","_rn_":"1"},{"1":"manipulation_trial1 - manipulation_trial3","2":"baseline_pre_squeeze","3":"-0.170407353","4":"0.04008912","5":"2396.472","6":"-0.27984479","7":"-0.0609699198","_rn_":"2"},{"1":"manipulation_trial1 - manipulation_trial4","2":"baseline_pre_squeeze","3":"-0.197302673","4":"0.04016589","5":"2397.583","6":"-0.30694964","7":"-0.0876557095","_rn_":"3"},{"1":"manipulation_trial1 - manipulation_trial5","2":"baseline_pre_squeeze","3":"-0.239907368","4":"0.04073546","5":"2405.906","6":"-0.35110888","7":"-0.1287058599","_rn_":"4"},{"1":"manipulation_trial2 - manipulation_trial3","2":"baseline_pre_squeeze","3":"0.018042075","4":"0.03938120","5":"2389.521","6":"-0.08946308","7":"0.1255472306","_rn_":"5"},{"1":"manipulation_trial2 - manipulation_trial4","2":"baseline_pre_squeeze","3":"-0.008853245","4":"0.03946077","5":"2392.050","6":"-0.11657553","7":"0.0988690400","_rn_":"6"},{"1":"manipulation_trial2 - manipulation_trial5","2":"baseline_pre_squeeze","3":"-0.051457941","4":"0.04003272","5":"2393.340","6":"-0.16074151","7":"0.0578256316","_rn_":"7"},{"1":"manipulation_trial3 - manipulation_trial4","2":"baseline_pre_squeeze","3":"-0.026895320","4":"0.03968210","5":"2389.628","6":"-0.13522189","7":"0.0814312448","_rn_":"8"},{"1":"manipulation_trial3 - manipulation_trial5","2":"baseline_pre_squeeze","3":"-0.069500016","4":"0.04025587","5":"2396.013","6":"-0.17939267","7":"0.0403926350","_rn_":"9"},{"1":"manipulation_trial4 - manipulation_trial5","2":"baseline_pre_squeeze","3":"-0.042604695","4":"0.04032719","5":"2391.811","6":"-0.15269218","7":"0.0674827914","_rn_":"10"},{"1":"manipulation_trial1 - manipulation_trial2","2":"squeeze","3":"-0.333432047","4":"0.03741981","5":"2370.503","6":"-0.43558350","7":"-0.2312805892","_rn_":"11"},{"1":"manipulation_trial1 - manipulation_trial3","2":"squeeze","3":"-0.317465557","4":"0.03741981","5":"2370.503","6":"-0.41961702","7":"-0.2153140996","_rn_":"12"},{"1":"manipulation_trial1 - manipulation_trial4","2":"squeeze","3":"-0.308500250","4":"0.03741981","5":"2370.503","6":"-0.41065171","7":"-0.2063487924","_rn_":"13"},{"1":"manipulation_trial1 - manipulation_trial5","2":"squeeze","3":"-0.256524045","4":"0.03748474","5":"2371.311","6":"-0.35885272","7":"-0.1541953708","_rn_":"14"},{"1":"manipulation_trial2 - manipulation_trial3","2":"squeeze","3":"0.015966490","4":"0.03735478","5":"2369.698","6":"-0.08600747","7":"0.1179404482","_rn_":"15"},{"1":"manipulation_trial2 - manipulation_trial4","2":"squeeze","3":"0.024931797","4":"0.03735478","5":"2369.698","6":"-0.07704216","7":"0.1269057554","_rn_":"16"},{"1":"manipulation_trial2 - manipulation_trial5","2":"squeeze","3":"0.076908002","4":"0.03741981","5":"2370.503","6":"-0.02524346","7":"0.1790594600","_rn_":"17"},{"1":"manipulation_trial3 - manipulation_trial4","2":"squeeze","3":"0.008965307","4":"0.03735478","5":"2369.698","6":"-0.09300865","7":"0.1109392658","_rn_":"18"},{"1":"manipulation_trial3 - manipulation_trial5","2":"squeeze","3":"0.060941513","4":"0.03741981","5":"2370.503","6":"-0.04120995","7":"0.1630929704","_rn_":"19"},{"1":"manipulation_trial4 - manipulation_trial5","2":"squeeze","3":"0.051976205","4":"0.03741981","5":"2370.503","6":"-0.05017525","7":"0.1541276633","_rn_":"20"},{"1":"manipulation_trial1 - manipulation_trial2","2":"baseline_post_squeeze","3":"0.071910229","4":"0.04855914","5":"2431.279","6":"-0.06064766","7":"0.2044681189","_rn_":"21"},{"1":"manipulation_trial1 - manipulation_trial3","2":"baseline_post_squeeze","3":"0.030472245","4":"0.04916743","5":"2437.697","6":"-0.10374590","7":"0.1646903906","_rn_":"22"},{"1":"manipulation_trial1 - manipulation_trial4","2":"baseline_post_squeeze","3":"0.034526012","4":"0.04916436","5":"2435.577","6":"-0.09968384","7":"0.1687358689","_rn_":"23"},{"1":"manipulation_trial1 - manipulation_trial5","2":"baseline_post_squeeze","3":"-0.011721675","4":"0.05071440","5":"2436.941","6":"-0.15016280","7":"0.1267194459","_rn_":"24"},{"1":"manipulation_trial2 - manipulation_trial3","2":"baseline_post_squeeze","3":"-0.041437984","4":"0.04675441","5":"2420.021","6":"-0.16906973","7":"0.0861937607","_rn_":"25"},{"1":"manipulation_trial2 - manipulation_trial4","2":"baseline_post_squeeze","3":"-0.037384217","4":"0.04677020","5":"2431.783","6":"-0.16505859","7":"0.0902901557","_rn_":"26"},{"1":"manipulation_trial2 - manipulation_trial5","2":"baseline_post_squeeze","3":"-0.083631904","4":"0.04839523","5":"2432.662","6":"-0.21574229","7":"0.0484784817","_rn_":"27"},{"1":"manipulation_trial3 - manipulation_trial4","2":"baseline_post_squeeze","3":"0.004053767","4":"0.04736732","5":"2413.654","6":"-0.12525138","7":"0.1333589172","_rn_":"28"},{"1":"manipulation_trial3 - manipulation_trial5","2":"baseline_post_squeeze","3":"-0.042193920","4":"0.04899068","5":"2428.846","6":"-0.17592993","7":"0.0915420897","_rn_":"29"},{"1":"manipulation_trial4 - manipulation_trial5","2":"baseline_post_squeeze","3":"-0.046247687","4":"0.04899366","5":"2430.852","6":"-0.17999177","7":"0.0874963912","_rn_":"30"},{"1":"manipulation_trial1 - manipulation_trial2","2":"relax","3":"-0.134402848","4":"0.03735478","5":"2369.698","6":"-0.23637681","7":"-0.0324288891","_rn_":"31"},{"1":"manipulation_trial1 - manipulation_trial3","2":"relax","3":"-0.126891349","4":"0.03741981","5":"2370.518","6":"-0.22904281","7":"-0.0247398858","_rn_":"32"},{"1":"manipulation_trial1 - manipulation_trial4","2":"relax","3":"-0.131778896","4":"0.03748563","5":"2371.329","6":"-0.23411001","7":"-0.0294477839","_rn_":"33"},{"1":"manipulation_trial1 - manipulation_trial5","2":"relax","3":"-0.102764538","4":"0.03735478","5":"2369.698","6":"-0.20473850","7":"-0.0007905796","_rn_":"34"},{"1":"manipulation_trial2 - manipulation_trial3","2":"relax","3":"0.007511498","4":"0.03741981","5":"2370.518","6":"-0.09463997","7":"0.1096629619","_rn_":"35"},{"1":"manipulation_trial2 - manipulation_trial4","2":"relax","3":"0.002623952","4":"0.03748563","5":"2371.329","6":"-0.09970716","7":"0.1049550638","_rn_":"36"},{"1":"manipulation_trial2 - manipulation_trial5","2":"relax","3":"0.031638310","4":"0.03735478","5":"2369.698","6":"-0.07033565","7":"0.1336122682","_rn_":"37"},{"1":"manipulation_trial3 - manipulation_trial4","2":"relax","3":"-0.004887547","4":"0.03755046","5":"2372.156","6":"-0.10739560","7":"0.0976205020","_rn_":"38"},{"1":"manipulation_trial3 - manipulation_trial5","2":"relax","3":"0.024126811","4":"0.03741981","5":"2370.518","6":"-0.07802465","7":"0.1262782746","_rn_":"39"},{"1":"manipulation_trial4 - manipulation_trial5","2":"relax","3":"0.029014358","4":"0.03748563","5":"2371.329","6":"-0.07331675","7":"0.1313454699","_rn_":"40"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
# Get EMMs by trial_phase within each manipulation_trial
emm_data <- emmip(
  lmm_trials,
  ~ trial_phase | manipulation_trial,
  CIs = TRUE,
  plotit = FALSE
)

# Ensure trial phase order
emm_data$trial_phase <- factor(
  emm_data$trial_phase,
  levels = c("baseline_pre_squeeze", "squeeze", "baseline_post_squeeze", "relax")
)

# Plot
ggplot(emm_data, aes(x = trial_phase, y = yvar, color = as.factor(manipulation_trial), group = manipulation_trial)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  #geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2) +
  labs(
    title = "Effect of Trial Phase on Pupil Response by Manipulation Trial",
    x = "Trial Phase",
    y = "Estimated Marginal Pupil Response",
    color = "Manipulation Trial"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## Across Manipulation trials with groups - Habituation Effects

``` r
# Result 4
lmm_trials_group <- lmer(
  trial_corr_rpd ~ group.y * manipulation_trial * trial_phase + (1 | id),
  data = df_grip
)
anova(lmm_trials_group)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Sum Sq"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Mean Sq"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["NumDF"],"name":[3],"type":["int"],"align":["right"]},{"label":["DenDF"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["F value"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Pr(>F)"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"0.2547551","2":"0.12737757","3":"2","4":"142.8629","5":"1.2546423","6":"2.883007e-01","_rn_":"group.y"},{"1":"8.0322001","2":"2.00805002","3":"4","4":"2358.0847","5":"19.7788707","6":"4.969356e-16","_rn_":"manipulation_trial"},{"1":"116.8698260","2":"38.95660867","3":"3","4":"2379.9064","5":"383.7144074","6":"2.858062e-203","_rn_":"trial_phase"},{"1":"0.1134831","2":"0.01418539","3":"8","4":"2358.0868","5":"0.1397231","6":"9.973795e-01","_rn_":"group.y:manipulation_trial"},{"1":"2.4886897","2":"0.41478162","3":"6","4":"2379.2287","5":"4.0855117","6":"4.392493e-04","_rn_":"group.y:trial_phase"},{"1":"6.4601257","2":"0.53834381","3":"12","4":"2352.7388","5":"5.3025734","6":"6.595825e-09","_rn_":"manipulation_trial:trial_phase"},{"1":"1.9601540","2":"0.08167308","3":"24","4":"2352.5445","5":"0.8044627","6":"7.348206e-01","_rn_":"group.y:manipulation_trial:trial_phase"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
anova_table <- anova(lmm_trials_group)
knitr::kable(anova_table, caption = "ANOVA Results")
```

|  | Sum Sq | Mean Sq | NumDF | DenDF | F value | Pr(\>F) |
|:---|---:|---:|---:|---:|---:|---:|
| group.y | 0.2547551 | 0.1273776 | 2 | 142.8629 | 1.2546423 | 0.2883007 |
| manipulation_trial | 8.0322001 | 2.0080500 | 4 | 2358.0847 | 19.7788707 | 0.0000000 |
| trial_phase | 116.8698260 | 38.9566087 | 3 | 2379.9064 | 383.7144074 | 0.0000000 |
| group.y:manipulation_trial | 0.1134831 | 0.0141854 | 8 | 2358.0868 | 0.1397231 | 0.9973795 |
| group.y:trial_phase | 2.4886897 | 0.4147816 | 6 | 2379.2287 | 4.0855117 | 0.0004392 |
| manipulation_trial:trial_phase | 6.4601257 | 0.5383438 | 12 | 2352.7388 | 5.3025734 | 0.0000000 |
| group.y:manipulation_trial:trial_phase | 1.9601540 | 0.0816731 | 24 | 2352.5445 | 0.8044627 | 0.7348206 |

ANOVA Results

``` r
r2_nakagawa(lmm_trials_group) 
```

    ## # R2 for Mixed Models
    ## 
    ##   Conditional R2: 0.363
    ##      Marginal R2: 0.356

``` r
emmeans(lmm_trials_group, ~ trial_phase | group.y * manipulation_trial)
```

    ## group.y = ASD, manipulation_trial = 1:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze  -0.113753 0.0494 2468 -0.21070  -0.0168
    ##  squeeze                0.064348 0.0453 2465 -0.02450   0.1532
    ##  baseline_post_squeeze -0.219991 0.0575 2470 -0.33283  -0.1072
    ##  relax                 -0.232878 0.0449 2465 -0.32086  -0.1449
    ## 
    ## group.y = CON, manipulation_trial = 1:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze  -0.138383 0.0462 2466 -0.22907  -0.0477
    ##  squeeze                0.126480 0.0440 2465  0.04018   0.2128
    ##  baseline_post_squeeze -0.213811 0.0605 2470 -0.33254  -0.0951
    ##  relax                 -0.320260 0.0440 2465 -0.40656  -0.2340
    ## 
    ## group.y = MHC, manipulation_trial = 1:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze  -0.172356 0.0550 2466 -0.28013  -0.0646
    ##  squeeze                0.013964 0.0495 2459 -0.08303   0.1110
    ##  baseline_post_squeeze -0.111864 0.0717 2467 -0.25247   0.0287
    ##  relax                 -0.328923 0.0495 2459 -0.42592  -0.2319
    ## 
    ## group.y = ASD, manipulation_trial = 2:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.023009 0.0462 2466 -0.06767   0.1137
    ##  squeeze                0.368390 0.0449 2465  0.28041   0.4564
    ##  baseline_post_squeeze -0.236849 0.0534 2469 -0.34156  -0.1321
    ##  relax                 -0.128991 0.0449 2465 -0.21697  -0.0410
    ## 
    ## group.y = CON, manipulation_trial = 2:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.034462 0.0458 2466 -0.05529   0.1242
    ##  squeeze                0.454682 0.0440 2465  0.36838   0.5410
    ##  baseline_post_squeeze -0.227493 0.0549 2469 -0.33524  -0.1197
    ##  relax                 -0.183744 0.0440 2465 -0.27005  -0.0974
    ## 
    ## group.y = MHC, manipulation_trial = 2:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.103640 0.0534 2462 -0.00113   0.2084
    ##  squeeze                0.389887 0.0495 2459  0.29289   0.4869
    ##  baseline_post_squeeze -0.343087 0.0629 2466 -0.46639  -0.2198
    ##  relax                 -0.160131 0.0495 2459 -0.25713  -0.0631
    ## 
    ## group.y = ASD, manipulation_trial = 3:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.041101 0.0494 2468 -0.05584   0.1380
    ##  squeeze                0.335861 0.0449 2465  0.24788   0.4238
    ##  baseline_post_squeeze -0.205267 0.0542 2469 -0.31146  -0.0991
    ##  relax                 -0.113939 0.0453 2465 -0.20279  -0.0251
    ## 
    ## group.y = CON, manipulation_trial = 3:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.045869 0.0453 2466 -0.04298   0.1347
    ##  squeeze                0.452420 0.0440 2465  0.36612   0.5387
    ##  baseline_post_squeeze -0.201975 0.0595 2470 -0.31864  -0.0853
    ##  relax                 -0.233490 0.0440 2465 -0.31979  -0.1472
    ## 
    ## group.y = MHC, manipulation_trial = 3:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.000125 0.0520 2461 -0.10185   0.1021
    ##  squeeze                0.376738 0.0495 2459  0.27974   0.4737
    ##  baseline_post_squeeze -0.261709 0.0617 2465 -0.38270  -0.1407
    ##  relax                 -0.140519 0.0495 2459 -0.23751  -0.0435
    ## 
    ## group.y = ASD, manipulation_trial = 4:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.115331 0.0483 2467  0.02062   0.2100
    ##  squeeze                0.349122 0.0449 2465  0.26114   0.4371
    ##  baseline_post_squeeze -0.286554 0.0549 2469 -0.39430  -0.1788
    ##  relax                 -0.137049 0.0449 2465 -0.22503  -0.0491
    ## 
    ## group.y = CON, manipulation_trial = 4:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.007307 0.0462 2466 -0.08338   0.0980
    ##  squeeze                0.444809 0.0440 2465  0.35851   0.5311
    ##  baseline_post_squeeze -0.158991 0.0595 2470 -0.27565  -0.0423
    ##  relax                 -0.167423 0.0444 2465 -0.25455  -0.0803
    ## 
    ## group.y = MHC, manipulation_trial = 4:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.055146 0.0527 2461 -0.04820   0.1585
    ##  squeeze                0.339075 0.0495 2459  0.24208   0.4361
    ##  baseline_post_squeeze -0.218981 0.0606 2465 -0.33779  -0.1002
    ##  relax                 -0.180679 0.0501 2460 -0.27885  -0.0825
    ## 
    ## group.y = ASD, manipulation_trial = 5:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.146594 0.0483 2467  0.05188   0.2413
    ##  squeeze                0.253992 0.0453 2465  0.16514   0.3428
    ##  baseline_post_squeeze -0.172612 0.0558 2469 -0.28198  -0.0632
    ##  relax                 -0.189077 0.0449 2465 -0.27705  -0.1011
    ## 
    ## group.y = CON, manipulation_trial = 5:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.117238 0.0500 2468  0.01912   0.2154
    ##  squeeze                0.428539 0.0440 2465  0.34224   0.5148
    ##  baseline_post_squeeze -0.141158 0.0617 2470 -0.26206  -0.0203
    ##  relax                 -0.242301 0.0440 2465 -0.32860  -0.1560
    ## 
    ## group.y = MHC, manipulation_trial = 5:
    ##  trial_phase              emmean     SE   df lower.CL upper.CL
    ##  baseline_pre_squeeze   0.026790 0.0527 2462 -0.07655   0.1301
    ##  squeeze                0.292658 0.0495 2459  0.19566   0.3897
    ##  baseline_post_squeeze -0.239418 0.0717 2469 -0.37996  -0.0989
    ##  relax                 -0.123257 0.0495 2459 -0.22025  -0.0263
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

### Post Hoc

``` r
contrast(emmeans(lmm_trials_group, ~ trial_phase | group.y * manipulation_trial), "pairwise")
```

    ## group.y = ASD, manipulation_trial = 1:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.17810 0.0667 2351  -2.669
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.10624 0.0756 2384   1.406
    ##  baseline_pre_squeeze - relax                  0.11912 0.0664 2353   1.793
    ##  squeeze - baseline_post_squeeze               0.28434 0.0729 2380   3.899
    ##  squeeze - relax                               0.29723 0.0634 2334   4.687
    ##  baseline_post_squeeze - relax                 0.01289 0.0727 2381   0.177
    ##  p.value
    ##   0.0383
    ##   0.4956
    ##   0.2768
    ##   0.0006
    ##   <.0001
    ##   0.9980
    ## 
    ## group.y = CON, manipulation_trial = 1:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.26486 0.0635 2343  -4.171
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.07543 0.0759 2390   0.994
    ##  baseline_pre_squeeze - relax                  0.18188 0.0635 2343   2.864
    ##  squeeze - baseline_post_squeeze               0.34029 0.0746 2393   4.564
    ##  squeeze - relax                               0.44674 0.0619 2332   7.218
    ##  baseline_post_squeeze - relax                 0.10645 0.0746 2393   1.428
    ##  p.value
    ##   0.0002
    ##   0.7529
    ##   0.0219
    ##   <.0001
    ##   <.0001
    ##   0.4821
    ## 
    ## group.y = MHC, manipulation_trial = 1:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.18632 0.0736 2356  -2.533
    ##  baseline_pre_squeeze - baseline_post_squeeze -0.06049 0.0900 2414  -0.672
    ##  baseline_pre_squeeze - relax                  0.15657 0.0736 2356   2.129
    ##  squeeze - baseline_post_squeeze               0.12583 0.0867 2404   1.451
    ##  squeeze - relax                               0.34289 0.0695 2332   4.931
    ##  baseline_post_squeeze - relax                 0.21706 0.0867 2404   2.502
    ##  p.value
    ##   0.0552
    ##   0.9077
    ##   0.1443
    ##   0.4679
    ##   <.0001
    ##   0.0598
    ## 
    ## group.y = ASD, manipulation_trial = 2:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.34538 0.0641 2338  -5.389
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.25986 0.0703 2367   3.696
    ##  baseline_pre_squeeze - relax                  0.15200 0.0641 2338   2.372
    ##  squeeze - baseline_post_squeeze               0.60524 0.0694 2368   8.718
    ##  squeeze - relax                               0.49738 0.0631 2332   7.883
    ##  baseline_post_squeeze - relax                -0.10786 0.0694 2368  -1.554
    ##  p.value
    ##   <.0001
    ##   0.0013
    ##   0.0828
    ##   <.0001
    ##   <.0001
    ##   0.4056
    ## 
    ## group.y = CON, manipulation_trial = 2:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.42022 0.0632 2341  -6.653
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.26195 0.0712 2379   3.679
    ##  baseline_pre_squeeze - relax                  0.21821 0.0632 2341   3.455
    ##  squeeze - baseline_post_squeeze               0.68218 0.0701 2377   9.732
    ##  squeeze - relax                               0.63843 0.0619 2332  10.314
    ##  baseline_post_squeeze - relax                -0.04375 0.0701 2377  -0.624
    ##  p.value
    ##   <.0001
    ##   0.0014
    ##   0.0031
    ##   <.0001
    ##   <.0001
    ##   0.9244
    ## 
    ## group.y = MHC, manipulation_trial = 2:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.28625 0.0724 2349  -3.954
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.44673 0.0821 2385   5.440
    ##  baseline_pre_squeeze - relax                  0.26377 0.0724 2349   3.643
    ##  squeeze - baseline_post_squeeze               0.73297 0.0796 2381   9.206
    ##  squeeze - relax                               0.55002 0.0695 2332   7.910
    ##  baseline_post_squeeze - relax                -0.18296 0.0796 2381  -2.298
    ##  p.value
    ##   0.0005
    ##   <.0001
    ##   0.0016
    ##   <.0001
    ##   <.0001
    ##   0.0988
    ## 
    ## group.y = ASD, manipulation_trial = 3:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.29476 0.0664 2353  -4.437
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.24637 0.0730 2382   3.374
    ##  baseline_pre_squeeze - relax                  0.15504 0.0667 2355   2.323
    ##  squeeze - baseline_post_squeeze               0.54113 0.0700 2370   7.729
    ##  squeeze - relax                               0.44980 0.0634 2334   7.093
    ##  baseline_post_squeeze - relax                -0.09133 0.0703 2369  -1.299
    ##  p.value
    ##   0.0001
    ##   0.0042
    ##   0.0930
    ##   <.0001
    ##   <.0001
    ##   0.5634
    ## 
    ## group.y = CON, manipulation_trial = 3:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.40655 0.0628 2338  -6.471
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.24784 0.0745 2390   3.327
    ##  baseline_pre_squeeze - relax                  0.27936 0.0628 2338   4.446
    ##  squeeze - baseline_post_squeeze               0.65440 0.0737 2391   8.878
    ##  squeeze - relax                               0.68591 0.0619 2332  11.082
    ##  baseline_post_squeeze - relax                 0.03152 0.0737 2391   0.428
    ##  p.value
    ##   <.0001
    ##   0.0049
    ##   0.0001
    ##   <.0001
    ##   <.0001
    ##   0.9738
    ## 
    ## group.y = MHC, manipulation_trial = 3:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.37661 0.0714 2343  -5.278
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.26183 0.0803 2381   3.261
    ##  baseline_pre_squeeze - relax                  0.14064 0.0714 2343   1.971
    ##  squeeze - baseline_post_squeeze               0.63845 0.0787 2378   8.113
    ##  squeeze - relax                               0.51726 0.0695 2332   7.439
    ##  baseline_post_squeeze - relax                -0.12119 0.0787 2378  -1.540
    ##  p.value
    ##   <.0001
    ##   0.0062
    ##   0.1993
    ##   <.0001
    ##   <.0001
    ##   0.4136
    ## 
    ## group.y = ASD, manipulation_trial = 4:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.23379 0.0656 2348  -3.565
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.40188 0.0728 2374   5.518
    ##  baseline_pre_squeeze - relax                  0.25238 0.0656 2348   3.848
    ##  squeeze - baseline_post_squeeze               0.63568 0.0706 2373   9.001
    ##  squeeze - relax                               0.48617 0.0631 2332   7.705
    ##  baseline_post_squeeze - relax                -0.14950 0.0706 2373  -2.117
    ##  p.value
    ##   0.0021
    ##   <.0001
    ##   0.0007
    ##   <.0001
    ##   <.0001
    ##   0.1479
    ## 
    ## group.y = CON, manipulation_trial = 4:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.43750 0.0635 2343  -6.889
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.16630 0.0751 2387   2.216
    ##  baseline_pre_squeeze - relax                  0.17473 0.0638 2341   2.739
    ##  squeeze - baseline_post_squeeze               0.60380 0.0737 2391   8.191
    ##  squeeze - relax                               0.61223 0.0622 2334   9.844
    ##  baseline_post_squeeze - relax                 0.00843 0.0740 2389   0.114
    ##  p.value
    ##   <.0001
    ##   0.1192
    ##   0.0315
    ##   <.0001
    ##   <.0001
    ##   0.9995
    ## 
    ## group.y = MHC, manipulation_trial = 4:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.28393 0.0719 2346  -3.951
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.27413 0.0799 2386   3.430
    ##  baseline_pre_squeeze - relax                  0.23582 0.0723 2349   3.263
    ##  squeeze - baseline_post_squeeze               0.55806 0.0778 2374   7.171
    ##  squeeze - relax                               0.51975 0.0700 2334   7.430
    ##  baseline_post_squeeze - relax                -0.03830 0.0782 2378  -0.490
    ##  p.value
    ##   0.0005
    ##   0.0034
    ##   0.0062
    ##   <.0001
    ##   <.0001
    ##   0.9614
    ## 
    ## group.y = ASD, manipulation_trial = 5:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.10740 0.0659 2346  -1.630
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.31921 0.0735 2377   4.345
    ##  baseline_pre_squeeze - relax                  0.33567 0.0656 2348   5.118
    ##  squeeze - baseline_post_squeeze               0.42660 0.0716 2379   5.962
    ##  squeeze - relax                               0.44307 0.0634 2334   6.987
    ##  baseline_post_squeeze - relax                 0.01647 0.0713 2376   0.231
    ##  p.value
    ##   0.3618
    ##   0.0001
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.9957
    ## 
    ## group.y = CON, manipulation_trial = 5:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.31130 0.0663 2359  -4.694
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.25840 0.0791 2396   3.266
    ##  baseline_pre_squeeze - relax                  0.35954 0.0663 2359   5.422
    ##  squeeze - baseline_post_squeeze               0.56970 0.0755 2396   7.549
    ##  squeeze - relax                               0.67084 0.0619 2332  10.838
    ##  baseline_post_squeeze - relax                 0.10114 0.0755 2396   1.340
    ##  p.value
    ##   <.0001
    ##   0.0061
    ##   <.0001
    ##   <.0001
    ##   <.0001
    ##   0.5373
    ## 
    ## group.y = MHC, manipulation_trial = 5:
    ##  contrast                                     estimate     SE   df t.ratio
    ##  baseline_pre_squeeze - squeeze               -0.26587 0.0719 2346  -3.699
    ##  baseline_pre_squeeze - baseline_post_squeeze  0.26621 0.0887 2416   3.003
    ##  baseline_pre_squeeze - relax                  0.15005 0.0719 2346   2.088
    ##  squeeze - baseline_post_squeeze               0.53208 0.0867 2403   6.134
    ##  squeeze - relax                               0.41592 0.0695 2332   5.982
    ##  baseline_post_squeeze - relax                -0.11616 0.0867 2403  -1.339
    ##  p.value
    ##   0.0013
    ##   0.0144
    ##   0.1573
    ##   <.0001
    ##   <.0001
    ##   0.5380
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

``` r
confint(contrast(emmeans(lmm_trials_group, ~ trial_phase | group.y * manipulation_trial), "pairwise"))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["group.y"],"name":[2],"type":["fct"],"align":["left"]},{"label":["manipulation_trial"],"name":[3],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[8],"type":["dbl"],"align":["right"]}],"data":[{"1":"baseline_pre_squeeze - squeeze","2":"ASD","3":"1","4":"-0.178101184","5":"0.06672291","6":"2350.643","7":"-0.349637121","8":"-0.006565246","_rn_":"1"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"ASD","3":"1","4":"0.106237164","5":"0.07555297","6":"2384.363","7":"-0.087997743","8":"0.300472071","_rn_":"2"},{"1":"baseline_pre_squeeze - relax","2":"ASD","3":"1","4":"0.119124072","5":"0.06642714","6":"2352.740","7":"-0.051651369","8":"0.289899513","_rn_":"3"},{"1":"squeeze - baseline_post_squeeze","2":"ASD","3":"1","4":"0.284338348","5":"0.07293097","6":"2379.582","7":"0.096843928","8":"0.471832767","_rn_":"4"},{"1":"squeeze - relax","2":"ASD","3":"1","4":"0.297225256","5":"0.06341624","6":"2333.844","7":"0.134189504","8":"0.460261007","_rn_":"5"},{"1":"baseline_post_squeeze - relax","2":"ASD","3":"1","4":"0.012886908","5":"0.07266047","6":"2381.328","7":"-0.173912002","8":"0.199685818","_rn_":"6"},{"1":"baseline_pre_squeeze - squeeze","2":"CON","3":"1","4":"-0.264863005","5":"0.06350463","6":"2342.764","7":"-0.428125544","8":"-0.101600465","_rn_":"7"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"CON","3":"1","4":"0.075427893","5":"0.07588767","6":"2389.669","7":"-0.119667180","8":"0.270522967","_rn_":"8"},{"1":"baseline_pre_squeeze - relax","2":"CON","3":"1","4":"0.181877381","5":"0.06350463","6":"2342.764","7":"0.018614841","8":"0.345139921","_rn_":"9"},{"1":"squeeze - baseline_post_squeeze","2":"CON","3":"1","4":"0.340290898","5":"0.07456492","6":"2393.052","7":"0.148596614","8":"0.531985182","_rn_":"10"},{"1":"squeeze - relax","2":"CON","3":"1","4":"0.446740386","5":"0.06189614","6":"2331.570","7":"0.287612510","8":"0.605868261","_rn_":"11"},{"1":"baseline_post_squeeze - relax","2":"CON","3":"1","4":"0.106449488","5":"0.07456492","6":"2393.052","7":"-0.085244796","8":"0.298143772","_rn_":"12"},{"1":"baseline_pre_squeeze - squeeze","2":"MHC","3":"1","4":"-0.186320192","5":"0.07355600","6":"2356.271","7":"-0.375422801","8":"0.002782417","_rn_":"13"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"MHC","3":"1","4":"-0.060492801","5":"0.09001677","6":"2414.383","7":"-0.291909861","8":"0.170924259","_rn_":"14"},{"1":"baseline_pre_squeeze - relax","2":"MHC","3":"1","4":"0.156566481","5":"0.07355600","6":"2356.271","7":"-0.032536128","8":"0.345669090","_rn_":"15"},{"1":"squeeze - baseline_post_squeeze","2":"MHC","3":"1","4":"0.125827391","5":"0.08674509","6":"2403.592","7":"-0.097179458","8":"0.348834240","_rn_":"16"},{"1":"squeeze - relax","2":"MHC","3":"1","4":"0.342886673","5":"0.06953074","6":"2331.570","7":"0.164131115","8":"0.521642231","_rn_":"17"},{"1":"baseline_post_squeeze - relax","2":"MHC","3":"1","4":"0.217059282","5":"0.08674509","6":"2403.592","7":"-0.005947567","8":"0.440066131","_rn_":"18"},{"1":"baseline_pre_squeeze - squeeze","2":"ASD","3":"2","4":"-0.345380241","5":"0.06408694","6":"2338.310","7":"-0.510140054","8":"-0.180620427","_rn_":"19"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"ASD","3":"2","4":"0.259858245","5":"0.07031454","6":"2366.880","7":"0.079089590","8":"0.440626901","_rn_":"20"},{"1":"baseline_pre_squeeze - relax","2":"ASD","3":"2","4":"0.152000272","5":"0.06408694","6":"2338.310","7":"-0.012759542","8":"0.316760085","_rn_":"21"},{"1":"squeeze - baseline_post_squeeze","2":"ASD","3":"2","4":"0.605238486","5":"0.06942467","6":"2367.908","7":"0.426757609","8":"0.783719363","_rn_":"22"},{"1":"squeeze - relax","2":"ASD","3":"2","4":"0.497380512","5":"0.06309812","6":"2331.570","7":"0.335162487","8":"0.659598538","_rn_":"23"},{"1":"baseline_post_squeeze - relax","2":"ASD","3":"2","4":"-0.107857974","5":"0.06942467","6":"2367.908","7":"-0.286338851","8":"0.070622904","_rn_":"24"},{"1":"baseline_pre_squeeze - squeeze","2":"CON","3":"2","4":"-0.420220919","5":"0.06316021","6":"2340.724","7":"-0.582598119","8":"-0.257843719","_rn_":"25"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"CON","3":"2","4":"0.261954836","5":"0.07120753","6":"2379.339","7":"0.078891111","8":"0.445018561","_rn_":"26"},{"1":"baseline_pre_squeeze - relax","2":"CON","3":"2","4":"0.218205664","5":"0.06316021","6":"2340.724","7":"0.055828464","8":"0.380582864","_rn_":"27"},{"1":"squeeze - baseline_post_squeeze","2":"CON","3":"2","4":"0.682175755","5":"0.07009414","6":"2377.274","7":"0.501974278","8":"0.862377232","_rn_":"28"},{"1":"squeeze - relax","2":"CON","3":"2","4":"0.638426583","5":"0.06189614","6":"2331.570","7":"0.479298708","8":"0.797554459","_rn_":"29"},{"1":"baseline_post_squeeze - relax","2":"CON","3":"2","4":"-0.043749172","5":"0.07009414","6":"2377.274","7":"-0.223950649","8":"0.136452305","_rn_":"30"},{"1":"baseline_pre_squeeze - squeeze","2":"MHC","3":"2","4":"-0.286246817","5":"0.07240043","6":"2348.696","7":"-0.472379024","8":"-0.100114611","_rn_":"31"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"MHC","3":"2","4":"0.446726656","5":"0.08212288","6":"2385.393","7":"0.235601602","8":"0.657851711","_rn_":"32"},{"1":"baseline_pre_squeeze - relax","2":"MHC","3":"2","4":"0.263771244","5":"0.07240043","6":"2348.696","7":"0.077639038","8":"0.449903451","_rn_":"33"},{"1":"squeeze - baseline_post_squeeze","2":"MHC","3":"2","4":"0.732973474","5":"0.07961745","6":"2381.373","7":"0.528289227","8":"0.937657720","_rn_":"34"},{"1":"squeeze - relax","2":"MHC","3":"2","4":"0.550018061","5":"0.06953074","6":"2331.570","7":"0.371262503","8":"0.728773619","_rn_":"35"},{"1":"baseline_post_squeeze - relax","2":"MHC","3":"2","4":"-0.182955412","5":"0.07961745","6":"2381.373","7":"-0.387639659","8":"0.021728835","_rn_":"36"},{"1":"baseline_pre_squeeze - squeeze","2":"ASD","3":"3","4":"-0.294759845","5":"0.06642717","6":"2352.850","7":"-0.465535365","8":"-0.123984325","_rn_":"37"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"ASD","3":"3","4":"0.246368146","5":"0.07302160","6":"2382.125","7":"0.058640870","8":"0.434095422","_rn_":"38"},{"1":"baseline_pre_squeeze - relax","2":"ASD","3":"3","4":"0.155040096","5":"0.06673084","6":"2355.479","7":"-0.016515966","8":"0.326596158","_rn_":"39"},{"1":"squeeze - baseline_post_squeeze","2":"ASD","3":"3","4":"0.541127991","5":"0.07000874","6":"2370.475","7":"0.361145693","8":"0.721110288","_rn_":"40"},{"1":"squeeze - relax","2":"ASD","3":"3","4":"0.449799941","5":"0.06341625","6":"2333.887","7":"0.286764159","8":"0.612835723","_rn_":"41"},{"1":"baseline_post_squeeze - relax","2":"ASD","3":"3","4":"-0.091328050","5":"0.07028942","6":"2368.537","7":"-0.272032057","8":"0.089375957","_rn_":"42"},{"1":"baseline_pre_squeeze - squeeze","2":"CON","3":"3","4":"-0.406550920","5":"0.06282761","6":"2338.056","7":"-0.568073159","8":"-0.245028680","_rn_":"43"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"CON","3":"3","4":"0.247844235","5":"0.07449007","6":"2390.241","7":"0.056342210","8":"0.439346260","_rn_":"44"},{"1":"baseline_pre_squeeze - relax","2":"CON","3":"3","4":"0.279359563","5":"0.06282761","6":"2338.056","7":"0.117837323","8":"0.440881803","_rn_":"45"},{"1":"squeeze - baseline_post_squeeze","2":"CON","3":"3","4":"0.654395154","5":"0.07371275","6":"2390.577","7":"0.464891502","8":"0.843898807","_rn_":"46"},{"1":"squeeze - relax","2":"CON","3":"3","4":"0.685910483","5":"0.06189614","6":"2331.570","7":"0.526782607","8":"0.845038358","_rn_":"47"},{"1":"baseline_post_squeeze - relax","2":"CON","3":"3","4":"0.031515328","5":"0.07371275","6":"2390.577","7":"-0.157988324","8":"0.221018981","_rn_":"48"},{"1":"baseline_pre_squeeze - squeeze","2":"MHC","3":"3","4":"-0.376612910","5":"0.07135655","6":"2342.879","7":"-0.560061758","8":"-0.193164062","_rn_":"49"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"MHC","3":"3","4":"0.261834118","5":"0.08030202","6":"2380.535","7":"0.055389901","8":"0.468278334","_rn_":"50"},{"1":"baseline_pre_squeeze - relax","2":"MHC","3":"3","4":"0.140643798","5":"0.07135655","6":"2342.879","7":"-0.042805051","8":"0.324092646","_rn_":"51"},{"1":"squeeze - baseline_post_squeeze","2":"MHC","3":"3","4":"0.638447027","5":"0.07869332","6":"2378.180","7":"0.436138392","8":"0.840755663","_rn_":"52"},{"1":"squeeze - relax","2":"MHC","3":"3","4":"0.517256707","5":"0.06953074","6":"2331.570","7":"0.338501149","8":"0.696012265","_rn_":"53"},{"1":"baseline_post_squeeze - relax","2":"MHC","3":"3","4":"-0.121190320","5":"0.07869332","6":"2378.180","7":"-0.323498955","8":"0.081118315","_rn_":"54"},{"1":"baseline_pre_squeeze - squeeze","2":"ASD","3":"4","4":"-0.233791612","5":"0.06558591","6":"2348.092","7":"-0.402404596","8":"-0.065178627","_rn_":"55"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"ASD","3":"4","4":"0.401884376","5":"0.07283625","6":"2373.850","7":"0.214633153","8":"0.589135600","_rn_":"56"},{"1":"baseline_pre_squeeze - relax","2":"ASD","3":"4","4":"0.252379573","5":"0.06558591","6":"2348.092","7":"0.083766589","8":"0.420992558","_rn_":"57"},{"1":"squeeze - baseline_post_squeeze","2":"ASD","3":"4","4":"0.635675988","5":"0.07062194","6":"2373.172","7":"0.454117372","8":"0.817234604","_rn_":"58"},{"1":"squeeze - relax","2":"ASD","3":"4","4":"0.486171185","5":"0.06309812","6":"2331.570","7":"0.323953160","8":"0.648389211","_rn_":"59"},{"1":"baseline_post_squeeze - relax","2":"ASD","3":"4","4":"-0.149504803","5":"0.07062194","6":"2373.172","7":"-0.331063419","8":"0.032053813","_rn_":"60"},{"1":"baseline_pre_squeeze - squeeze","2":"CON","3":"4","4":"-0.437502116","5":"0.06350462","6":"2342.769","7":"-0.600764636","8":"-0.274239595","_rn_":"61"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"CON","3":"4","4":"0.166297600","5":"0.07505010","6":"2386.865","7":"-0.026644370","8":"0.359239570","_rn_":"62"},{"1":"baseline_pre_squeeze - relax","2":"CON","3":"4","4":"0.174730059","5":"0.06379083","6":"2340.528","7":"0.010731612","8":"0.338728506","_rn_":"63"},{"1":"squeeze - baseline_post_squeeze","2":"CON","3":"4","4":"0.603799716","5":"0.07371273","6":"2390.532","7":"0.414296116","8":"0.793303317","_rn_":"64"},{"1":"squeeze - relax","2":"CON","3":"4","4":"0.612232175","5":"0.06219626","6":"2333.888","7":"0.452332847","8":"0.772131503","_rn_":"65"},{"1":"baseline_post_squeeze - relax","2":"CON","3":"4","4":"0.008432459","5":"0.07395945","6":"2388.897","7":"-0.181705517","8":"0.198570435","_rn_":"66"},{"1":"baseline_pre_squeeze - squeeze","2":"MHC","3":"4","4":"-0.283928778","5":"0.07186633","6":"2345.946","7":"-0.468688052","8":"-0.099169504","_rn_":"67"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"MHC","3":"4","4":"0.274127002","5":"0.07992583","6":"2385.895","7":"0.068650234","8":"0.479603771","_rn_":"68"},{"1":"baseline_pre_squeeze - relax","2":"MHC","3":"4","4":"0.235824757","5":"0.07228087","6":"2348.826","7":"0.049999916","8":"0.421649597","_rn_":"69"},{"1":"squeeze - baseline_post_squeeze","2":"MHC","3":"4","4":"0.558055780","5":"0.07782510","6":"2374.416","7":"0.357978988","8":"0.758132573","_rn_":"70"},{"1":"squeeze - relax","2":"MHC","3":"4","4":"0.519753535","5":"0.06995798","6":"2334.238","7":"0.339899755","8":"0.699607315","_rn_":"71"},{"1":"baseline_post_squeeze - relax","2":"MHC","3":"4","4":"-0.038302246","5":"0.07821081","6":"2377.635","7":"-0.239370441","8":"0.162765949","_rn_":"72"},{"1":"baseline_pre_squeeze - squeeze","2":"ASD","3":"5","4":"-0.107398121","5":"0.06588544","6":"2345.912","7":"-0.276781288","8":"0.061985045","_rn_":"73"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"ASD","3":"5","4":"0.319205782","5":"0.07346244","6":"2377.048","7":"0.130344889","8":"0.508066674","_rn_":"74"},{"1":"baseline_pre_squeeze - relax","2":"ASD","3":"5","4":"0.335671246","5":"0.06558590","6":"2348.054","7":"0.167058288","8":"0.504284203","_rn_":"75"},{"1":"squeeze - baseline_post_squeeze","2":"ASD","3":"5","4":"0.426603903","5":"0.07155170","6":"2378.807","7":"0.242655336","8":"0.610552470","_rn_":"76"},{"1":"squeeze - relax","2":"ASD","3":"5","4":"0.443069367","5":"0.06341624","6":"2333.845","7":"0.280033615","8":"0.606105119","_rn_":"77"},{"1":"baseline_post_squeeze - relax","2":"ASD","3":"5","4":"0.016465464","5":"0.07126657","6":"2375.916","7":"-0.166750242","8":"0.199681170","_rn_":"78"},{"1":"baseline_pre_squeeze - squeeze","2":"CON","3":"5","4":"-0.311301043","5":"0.06631637","6":"2359.370","7":"-0.481791378","8":"-0.140810709","_rn_":"79"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"CON","3":"5","4":"0.258395986","5":"0.07910875","6":"2395.590","7":"0.055020394","8":"0.461771579","_rn_":"80"},{"1":"baseline_pre_squeeze - relax","2":"CON","3":"5","4":"0.359538550","5":"0.06631637","6":"2359.370","7":"0.189048216","8":"0.530028885","_rn_":"81"},{"1":"squeeze - baseline_post_squeeze","2":"CON","3":"5","4":"0.569697030","5":"0.07546966","6":"2395.839","7":"0.375676963","8":"0.763717096","_rn_":"82"},{"1":"squeeze - relax","2":"CON","3":"5","4":"0.670839594","5":"0.06189614","6":"2331.570","7":"0.511711718","8":"0.829967469","_rn_":"83"},{"1":"baseline_post_squeeze - relax","2":"CON","3":"5","4":"0.101142564","5":"0.07546966","6":"2395.839","7":"-0.092877502","8":"0.295162630","_rn_":"84"},{"1":"baseline_pre_squeeze - squeeze","2":"MHC","3":"5","4":"-0.265868090","5":"0.07186621","6":"2345.642","7":"-0.450627068","8":"-0.081109112","_rn_":"85"},{"1":"baseline_pre_squeeze - baseline_post_squeeze","2":"MHC","3":"5","4":"0.266208025","5":"0.08865353","6":"2416.227","7":"0.038295734","8":"0.494120317","_rn_":"86"},{"1":"baseline_pre_squeeze - relax","2":"MHC","3":"5","4":"0.150047270","5":"0.07186621","6":"2345.642","7":"-0.034711707","8":"0.334806248","_rn_":"87"},{"1":"squeeze - baseline_post_squeeze","2":"MHC","3":"5","4":"0.532076115","5":"0.08674346","6":"2403.140","7":"0.309073429","8":"0.755078801","_rn_":"88"},{"1":"squeeze - relax","2":"MHC","3":"5","4":"0.415915360","5":"0.06953074","6":"2331.570","7":"0.237159802","8":"0.594670918","_rn_":"89"},{"1":"baseline_post_squeeze - relax","2":"MHC","3":"5","4":"-0.116160755","5":"0.08674346","6":"2403.140","7":"-0.339163441","8":"0.106841931","_rn_":"90"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
# Get estimated marginal means for the interaction
emm_data <- emmip(
  lmm_trials_group,
  ~ trial_phase | group.y * manipulation_trial,
  CIs = TRUE,
  plotit = FALSE
)

# Ensure correct order of trial phases
emm_data$trial_phase <- factor(
  emm_data$trial_phase,
  levels = c("baseline_pre_squeeze", "squeeze", "baseline_post_squeeze", "relax")
)

# Plot
ggplot(emm_data, aes(x = trial_phase, y = yvar, group = manipulation_trial, color = as.factor(manipulation_trial))) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  #geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2) +
  facet_wrap(~ group.y) +
  scale_x_discrete(labels = c(
    "baseline_pre_squeeze" = "Baseline Pre Grip",
    "squeeze" = "Grip",
    "baseline_post_squeeze" = "Baseline Post Grip",
    "relax" = "Relax"
  )) +
  labs(
    title = "Group × Manipulation Trial × Phase Interaction on Pupil Response",
    x = "Trial Phase",
    y = "Estimated Marginal Pupil Response",
    color = "Manipulation Trial"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

## Assumption checks

### Baseline model - Effect of Experimental Manipulation

``` r
check_model(lmm_baseline)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
performance::check_model(lmm_baseline)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-35-2.png)<!-- -->

### Group

``` r
check_model(lmm_group)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
performance::check_model(lmm_group)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-36-2.png)<!-- -->

### Across trials - Habituation effects

``` r
check_model(lmm_trials)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
performance::check_model(lmm_trials)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-37-2.png)<!-- -->

### Habituation effects with groups

``` r
check_model(lmm_trials_group)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

``` r
performance::check_model(lmm_trials_group)
```

![](MA_PN_Manipulation-Effect_files/figure-gfm/unnamed-chunk-38-2.png)<!-- -->
