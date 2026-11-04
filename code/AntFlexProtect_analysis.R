### Beyond a Reflex: Intraspecific Differences in the Flexible Self-Protection of Ants ###
### Last edited 29/10/2026 (DMY) ###

# Clear workspace and load libraries
rm(list = ls())
pacman::p_load(tidyverse, readr, emmeans, betareg, mgcv, cowplot, patchwork, pwr, lmtest)

# =========================================================================
# MAIN STUDY: Dataframe creation + treatment randomisation
# =========================================================================
# Generate treatment assignments 
set.seed(1996)
temp <- tibble(
  AntID = 1:192,
  Treatment = sample(rep(
    c("Worker_Cont", "Worker_Proc", "Worker_Treat",
      "Soldier_Cont", "Soldier_Proc", "Soldier_Treat"),
    each = 32
  ))
)

# Counter-balanced antenna side assignment per treatment group
AntID <- temp %>%
  filter(Treatment == "Worker_Cont") %>%
  mutate(Treat_which_side = sample(rep(c("Left", "Right"), each = 16)))

AntID <- temp %>%
  filter(Treatment == "Worker_Proc") %>%
  mutate(Treat_which_side = sample(rep(c("Left", "Right"), each = 16))) %>%
  add_row(AntID, .)

AntID <- temp %>%
  filter(Treatment == "Worker_Treat") %>%
  mutate(Treat_which_side = sample(rep(c("Left", "Right"), each = 16))) %>%
  add_row(AntID, .)

AntID <- temp %>%
  filter(Treatment == "Soldier_Cont") %>%
  mutate(Treat_which_side = sample(rep(c("Left", "Right"), each = 16))) %>%
  add_row(AntID, .)

AntID <- temp %>%
  filter(Treatment == "Soldier_Proc") %>%
  mutate(Treat_which_side = sample(rep(c("Left", "Right"), each = 16))) %>%
  add_row(AntID, .)

AntID <- temp %>%
  filter(Treatment == "Soldier_Treat") %>%
  mutate(Treat_which_side = sample(rep(c("Left", "Right"), each = 16))) %>%
  add_row(AntID, .) %>%
  .[sample(1:nrow(.)), ]

AntID <- AntID %>%
  arrange(AntID) %>%
  mutate(AntID = paste0("ant_", AntID))

# =========================================================================
# MAIN STUDY: Analysis
# =========================================================================
# Import grooming data
groom_raw <- read_csv("data/grooming_data.csv")

# Join grooming data with treatment assignments
groom_trim <- groom_raw %>%
  full_join(AntID, by = "AntID") %>%
  select(AntID, Treatment, Treat_which_side, Behaviour, Start_sec, Stop_sec, Duration_sec) %>%
  na.omit() %>%
  mutate(Frequency = 1,
         Minute = floor(Stop_sec / 60))

# Assign grooming to treated vs untreated antenna
groom_trim <- groom_trim %>%
  mutate(
    Groom_treated_freq = case_when(
      Treat_which_side == "Left"  & Behaviour == "Left_antenna"  ~ 1,
      Treat_which_side == "Right" & Behaviour == "Right_antenna" ~ 1,
      TRUE ~ NA_real_
    ),
    Groom_untreated_freq = case_when(
      Treat_which_side == "Left"  & Behaviour == "Right_antenna" ~ 1,
      Treat_which_side == "Right" & Behaviour == "Left_antenna"  ~ 1,
      TRUE ~ NA_real_
    ),
    Groom_treated_time = case_when(
      Treat_which_side == "Left"  & Behaviour == "Left_antenna"  ~ Duration_sec,
      Treat_which_side == "Right" & Behaviour == "Right_antenna" ~ Duration_sec,
      TRUE ~ NA_real_
    ),
    Groom_untreated_time = case_when(
      Treat_which_side == "Left"  & Behaviour == "Right_antenna" ~ Duration_sec,
      Treat_which_side == "Right" & Behaviour == "Left_antenna"  ~ Duration_sec,
      TRUE ~ NA_real_
    )
  ) %>%
  select(AntID, Treatment, Minute, Groom_treated_freq, Groom_untreated_freq,
         Groom_treated_time, Groom_untreated_time)

# Aggregate per ant
groom_total <- groom_trim %>%
  group_by(AntID, Treatment) %>%
  summarise(
    Groom_treated_freq = sum(Groom_treated_freq, na.rm = TRUE),
    Groom_untreated_freq = sum(Groom_untreated_freq, na.rm = TRUE),
    Groom_treated_time = sum(Groom_treated_time, na.rm = TRUE),
    Groom_untreated_time = sum(Groom_untreated_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total_groom_freq = Groom_treated_freq + Groom_untreated_freq,
    Total_groom_time = Groom_treated_time + Groom_untreated_time,
    Prop_treated_freq = Groom_treated_freq / Total_groom_freq,
    Prop_treated_time = Groom_treated_time / Total_groom_time
  )

# Split Treatment into Caste and Condition factors
groom_total <- groom_total %>%
  mutate(
    Caste = if_else(str_detect(Treatment, "Soldier"), "Soldier", "Worker"),
    Condition = str_remove(Treatment, "^(Soldier|Worker)_"),
    Caste = factor(Caste),
    Condition = factor(Condition, levels = c("Cont", "Proc", "Treat"))
  )

# Treatment labels for plots
treat_labels <- c("Soldier\nNo-Touch Control", "Soldier\nInnocuous Control", "Soldier\nNoxious Treatment",
                  "Forager\nNo-Touch Control", "Forager\nInnocuous Control", "Forager\nNoxious Treatment")

# =========================================================================
# MAIN: Beta regression on proportional grooming
# =========================================================================

# Squeeze proportions away from 0 and 1
n <- nrow(groom_total)
groom_total <- groom_total %>%
  mutate(
    Prop_treated_freq_beta = (Prop_treated_freq * (n - 1) + 0.5) / n,
    Prop_treated_time_beta = (Prop_treated_time * (n - 1) + 0.5) / n
  )

# Proportion frequency
groom_total$Caste = relevel(groom_total$Caste, ref = "Soldier")
freq_model <- betareg(Prop_treated_freq_beta ~ Caste * Condition, data = groom_total)

summary(freq_model)  # MAIN RESULT: REPORT IN A TABLE (combined with below in a stacked table)

# Frequency diagnostics. All good.
par(mfrow = c(2, 2))
plot(freq_model)
par(mfrow = c(1, 1))

# Frequency post-hoc
freq_emm <- emmeans(freq_model, pairwise ~ Condition | Caste, adjust = "tukey")
freq_emm

# Proportion duration
groom_total$Condition = as.factor(groom_total$Condition)
groom_total$Condition = relevel(groom_total$Condition, ref = "Cont")
groom_total$Caste = relevel(groom_total$Caste, ref = "Soldier")
time_model <- betareg(Prop_treated_time_beta ~ Caste * Condition, data = groom_total)
summary(time_model)  # MAIN RESULT: REPORT IN A TABLE

# Duration diagnostics
par(mfrow = c(2, 2))
plot(time_model)
par(mfrow = c(1, 1))

# Duration post-hoc
time_emm <- emmeans(time_model, pairwise ~ Condition | Caste, adjust = "tukey")
time_emm


# =========================================================================
# FIGURE 1: Proportion boxplots
# =========================================================================

graph_prop_freq <- ggplot(groom_total, aes(x = Treatment, y = Prop_treated_freq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitter(width = 0.07), alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  scale_x_discrete(labels = treat_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = NULL, y = "Proportion of Grooming Frequency\nDirected at Treated Antenna") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        text = element_text(family = "Times New Roman", size = 10), 
        axis.text = element_text(color = "black"))

graph_prop_time <- ggplot(groom_total, aes(x = Treatment, y = Prop_treated_time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitter(width = 0.07), alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  scale_x_discrete(labels = treat_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Treatment Group", y = "Proportion of Grooming Duration\nDirected at Treated Antenna") +
  theme_classic() +
  theme(text = element_text(family = "Times New Roman", size = 10), 
        axis.text = element_text(color = "black"))

fig_proportion_boxplots <- graph_prop_freq / graph_prop_time
fig_proportion_boxplots
ggsave("output/Fig1_proportion_boxplots.png", fig_proportion_boxplots, width = 16, height = 16, units = "cm")

# =========================================================================
# MAIN: In text descriptive statistics (raw mean and SD of proportions)
# =========================================================================
desc <- groom_total %>%
  group_by(Caste, Condition) %>%
  summarise(n = n(),
            freq_mean = mean(Prop_treated_freq), freq_sd = sd(Prop_treated_freq),
            time_mean = mean(Prop_treated_time), time_sd = sd(Prop_treated_time),
            .groups = "drop")
print(desc)
write_csv(desc, "output/intext_descriptives.csv")

# =========================================================================
# MAIN: Temporal GAM
# =========================================================================

ts_counts <- groom_trim %>%
  filter(Minute <= 9) %>%
  group_by(AntID, Minute) %>%
  summarise(Treated = sum(Groom_treated_freq, na.rm = TRUE),
            Untreated = sum(Groom_untreated_freq, na.rm = TRUE),
            .groups = "drop")

ants <- groom_total %>%
  transmute(AntID,
            Caste = factor(if_else(Caste == "Soldier", "Soldier", "Worker"),
                           levels = c("Worker", "Soldier")),
            Condition)

ts_full <- ants %>%
  crossing(Minute = 0:9) %>%
  left_join(ts_counts, by = c("AntID", "Minute")) %>%
  mutate(across(c(Treated, Untreated), ~ replace_na(.x, 0))) %>%
  pivot_longer(c(Treated, Untreated), names_to = "Antenna", values_to = "Bouts") %>%
  mutate(Antenna = factor(Antenna, levels = c("Untreated", "Treated")),
         AntID = factor(AntID),
         grp = interaction(Caste, Condition, Antenna))

cat("Rows:", nrow(ts_full), " proportion zero:", round(mean(ts_full$Bouts == 0), 3), "\n")

gam_rev <- bam(Bouts ~ Caste * Condition * Antenna +
                 s(Minute, by = grp, k = 6) + s(AntID, bs = "re"),
               family = nb(), data = ts_full, method = "fREML", discrete = TRUE)
gam_rev_null <- bam(Bouts ~ Antenna + s(Minute, by = Antenna, k = 6) + s(AntID, bs = "re"),
                    family = nb(), data = ts_full, method = "fREML", discrete = TRUE)
gam_rev_s <- summary(gam_rev)
print(gam_rev_s)
cat("theta:", gam_rev$family$getTheta(TRUE), "\n")
cat("AIC full:", AIC(gam_rev), " null:", AIC(gam_rev_null),
    " dAIC:", AIC(gam_rev_null) - AIC(gam_rev), "\n")

# Robustness: Poisson dispersion; k = 10; antenna-within-ant random effect
gam_pois <- bam(Bouts ~ Caste * Condition * Antenna + s(Minute, by = grp, k = 6) + s(AntID, bs = "re"),
                family = poisson(), data = ts_full, method = "fREML", discrete = TRUE)
cat("Poisson dispersion:", sum(residuals(gam_pois, type = "pearson")^2) / df.residual(gam_pois), "\n")
gam_k10 <- update(gam_rev, . ~ . - s(Minute, by = grp, k = 6) + s(Minute, by = grp, k = 10))
gam_re2 <- update(gam_rev, . ~ . + s(AntID, Antenna, bs = "re"))
print(lapply(list(k10 = gam_k10, re2 = gam_re2), \(m) round(summary(m)$p.table, 3)))

write_csv(as_tibble(gam_rev_s$p.table, rownames = "term"), "output/Table2_gam_parametric.csv")
write_csv(as_tibble(gam_rev_s$s.table, rownames = "term"), "output/Table3_gam_smooths.csv")

# Contrasts (response scale = ratios of bout rates)
emm_ant <- emmeans(gam_rev, ~ Antenna | Caste * Condition, type = "response")
print(confint(contrast(emm_ant, "revpairwise")))          # treated / untreated
print(contrast(emm_ant, "revpairwise"))
emm_int <- emmeans(gam_rev, ~ Antenna * Condition | Caste)
print(contrast(emm_int, interaction = c("revpairwise", "trt.vs.ctrl"), type = "response"))
emm_cond <- emmeans(gam_rev, ~ Condition | Caste * Antenna, type = "response")
print(contrast(emm_cond, "trt.vs.ctrl"))                  # vs no-touch, per antenna
emm_caste <- emmeans(gam_rev, ~ Caste | Condition * Antenna, type = "response")
print(contrast(emm_caste, "revpairwise"))                 # soldier / worker
print(contrast(emmeans(gam_rev, ~ Caste, type = "response"), "revpairwise"))

# Treated / untreated ratio by minute, noxious groups
emm_min <- emmeans(gam_rev, ~ Antenna | Caste * Condition * Minute,
                   at = list(Minute = 0:9), type = "response")
by_min <- as_tibble(contrast(emm_min, "revpairwise")) %>% filter(Condition == "Treat")
print(by_min, n = 20)

# =========================================================================
# MAIN: Figure 2
# =========================================================================
nd <- expand.grid(Minute = seq(0, 9, length.out = 100),
                  Caste = levels(ts_full$Caste),
                  Condition = levels(ts_full$Condition),
                  Antenna = levels(ts_full$Antenna)) %>%
  mutate(grp = interaction(Caste, Condition, Antenna),
         AntID = ts_full$AntID[1])
pr <- predict(gam_rev, nd, exclude = "s(AntID)", se.fit = TRUE)
nd <- nd %>% mutate(fit = exp(pr$fit),
                    lower = exp(pr$fit - 1.96 * pr$se.fit),
                    upper = exp(pr$fit + 1.96 * pr$se.fit),
                    Caste = recode(Caste, Worker = "Forager"),
                    Antenna = factor(Antenna, levels = c("Treated", "Untreated")))

fig_rev <- ggplot(nd, aes(Minute, colour = Antenna, fill = Antenna)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit)) +
  facet_grid(Caste ~ Condition, labeller = labeller(
    Condition = c(Cont = "No-touch control", Proc = "Innocuous control", Treat = "Noxious treatment"))) +
  scale_x_continuous(breaks = 0:9) +
  scale_colour_manual(values = c(Treated = "#E69F00", Untreated = "#56B4E9")) +
  scale_fill_manual(values = c(Treated = "#E69F00", Untreated = "#56B4E9")) +
  theme_classic() +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  labs(x = "Minute", y = "Grooming bouts per minute", colour = "Antenna", fill = "Antenna")
ggsave("output/Fig2_Timeseries.png", fig_rev, width = 3200, height = 2400, units = "px")

# =========================================================================
# MAIN: Variable-dispersion beta regression
# =========================================================================
# Test whether response variability (precision) differs across groups.
# Higher precision = tighter spread; a group effect here means the
# treatment doesn't just shift the mean but changes how variable
# individuals are.

# Frequency: variable-dispersion model
freq_model_vd <- betareg(Prop_treated_freq_beta ~ Caste * Condition | Caste * Condition,
                         data = groom_total)
summary(freq_model_vd)  # Precision sub-model coefficients

# LR test: does allowing group-level dispersion improve fit?
lrtest(freq_model, freq_model_vd)

# Duration: variable-dispersion model
time_model_vd <- betareg(Prop_treated_time_beta ~ Caste * Condition | Caste * Condition,
                         data = groom_total)
summary(time_model_vd)  # Precision sub-model coefficients

# LR test
lrtest(time_model, time_model_vd)

# =========================================================================
# Supplementary Figure 1: Temperature calibration
# =========================================================================
pilot_data = read_csv("data/pilot_data.csv")
print(count(pilot_data, caste, temperature_deg))

pilot_means = pilot_data %>%
  select(caste,temperature_deg,antenna_use) %>%
  na.omit() %>%
  group_by(caste,temperature_deg) %>%
  summarise(mean = mean(antenna_use),
            sd = sd(antenna_use))

pilot_graph = ggplot(data = pilot_means, aes(x = temperature_deg, y = mean, colour = caste)) +
  geom_point(data = pilot_means, aes(x = temperature_deg, y = mean), size = 3, position = position_jitter(width = 0.5, height = 0.02, seed = 1)) +
  geom_errorbar(data = pilot_means, aes(ymin = mean - sd, ymax = mean + sd, width = 2), position = position_jitter(width = 0.5, height = 0.02, seed = 1)) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_point(data = pilot_data, aes(x = temperature_deg, y = antenna_use), size = 1, position = position_jitter(width = 1, height = 0.02, seed = 1), alpha = 0.3) +
  coord_cartesian(ylim = c(-0.1, 1.1), xlim = c(15, 75), expand = FALSE) + # Set limits and prevent expansion
  theme_classic() +
  scale_colour_brewer(palette = "Set2", labels = c("Forager", "Soldier")) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(20, 70, by = 10)) +
  theme(text = element_text(family = "Times New Roman", size = 10), axis.text = element_text(color = "black")) +
  geom_vline(xintercept = 55, linetype = "dashed", colour = "red") +
  theme(legend.position = c(0.11, 0.7)) +
  labs(x = "Probe Temperature (°C)", y = "Antennal Usage", color = "Ant Caste")

pilot_graph
ggsave(filename = "output/Fig_S1_probe_calibration.png", plot = pilot_graph, width = 16, height = 12, units = "cm")

# =========================================================================
# Supplementary Figure 2: INSTAR reporting
# =========================================================================
# Render the INSTAR summary figure (Supplementary Figure 2) from the
# completed reporting sheet.
# install.packages("remotes"); remotes::install_github("thomased/instarreport")
library(instarreport)

items  <- read_items("data/INSTAR_Allman_et_al.csv")
report <- instar_report(items)   # title/authors/journal come from the sheet's reserved rows
print(report)
subset(summary(report), status != "reported")   # anything not yet filled

save_figure(report, "output/INSTAR_Allman_et_al.pdf")
save_figure(report, "output/INSTAR_Allman_et_al.png", dpi = 300)

# =========================================================================
# Supplementary Figure 3: Proportion time series
# =========================================================================
graph_ts_prop_freq <- ggplot(ts_full, aes(x = Minute, y = Prop_treated_freq, colour = Treatment)) +
  geom_smooth(method = "gam", alpha = 0.25, formula = y ~ s(x, bs = "tp"),
              method.args = list(family = "quasibinomial"), level = 0.95) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(xlim = c(1, 8), ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman", size = 10), axis.text = element_text(color = "black")) +
  labs(x = NULL, y = "Proportion of Grooming Frequency\nDirected at Treated Antenna")

graph_ts_prop_time <- ggplot(time_series, aes(x = Minute, y = Prop_treated_time, colour = Treatment)) +
  geom_smooth(method = "gam", alpha = 0.25, formula = y ~ s(x, bs = "tp"),
              method.args = list(family = "quasibinomial"), level = 0.95) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(xlim = c(1, 8), ylim = c(0, 1)) +
  guides(colour = guide_legend(ncol = 2)) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(family = "Times New Roman", size = 10), axis.text = element_text(color = "black")) +
  scale_color_discrete(labels = c("Soldier No-Touch", "Soldier Innocuous Control", "Soldier Noxious Treatment", "Forager No-Touch", "Forager Innocuous Control", "Forager Noxious Treatment")) +
  labs(x = "Minute", y = "Proportion of Grooming Duration\nDirected at Treated Antenna")

fig_supp_proportion_ts <- graph_ts_prop_freq / graph_ts_prop_time
ggsave("output/Supplementary_2_proportional_timeseries.png", fig_supp_proportion_ts, width = 16, height = 20, units = "cm")


