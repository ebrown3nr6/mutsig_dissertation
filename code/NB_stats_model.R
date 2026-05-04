required_packages <- c("MASS", "tidyverse", "broom", "emmeans", "car")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(MASS)
library(tidyverse)
library(broom)
library(emmeans)
library(car)


#find paths and load information from csv
output_folder <- "~/shared-team/2025-masters-project/people/eleanor/original_and_processed_files/processed_c/charts/mut_signatures"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

summary_csv <- "~/shared-team/2025-masters-project/people/eleanor/original_and_processed_files/processed_c/mutation_summary_final_c.csv"

REFERENCE_SPECIES <- "styphimurium"

df <- read_csv(summary_csv, show_col_types = FALSE) %>%
  mutate(
    mutation_class = factor(mutation_class),
    species        = factor(species),
    species        = relevel(species, ref = REFERENCE_SPECIES),
    is_intergenic  = as.integer(is_intergenic)
  )


# check for 0 opportunities
n_zero_opp <- sum(df$opportunities == 0, na.rm = TRUE)
if (n_zero_opp > 0) {
  cat("\nWARNING:", n_zero_opp, "rows have opportunities = 0 and will be removed",
      "(log(0) is undefined)\n")
  df <- df %>% filter(opportunities > 0)
  cat("Rows remaining:", nrow(df), "\n")
} else {
  cat("\nNo zero opportunities found \n")
}

# fitting binomial model

m1 <- glm.nb(
  no_mutations ~ mutation_class * species + is_intergenic + offset(log(opportunities)),
  data = df
)

# find dispersion with alpha
m1_theta <- m1$theta
m1_alpha <- 1 / m1_theta

print(summary(m1))


# make coefficients table easier to read

m1_coefs <- tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(IRR = estimate, lower_CI = conf.low, upper_CI = conf.high) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4)),
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(m1_coefs, n = Inf)
write_csv(m1_coefs, file.path(output_folder, "model1_coefficients.csv"))
cat("\nSaved: model1_coefficients.csv\n")

# estimated marginal means model

m1_emm <- emmeans(m1, ~ species | mutation_class, type = "response")
print(m1_emm)

# pairwise contrast using predicted rates from emmeans, Bonf corrected
m1_contrasts <- contrast(m1_emm, method = "pairwise", adjust = "bonferroni")

m1_contrasts_df <- as.data.frame(summary(m1_contrasts)) %>%
  rename(p_value_bonferroni = p.value) %>%
  mutate(
    significant = p_value_bonferroni < 0.05,
    sig_label   = case_when(
      p_value_bonferroni < 0.001 ~ "***",
      p_value_bonferroni < 0.01  ~ "**",
      p_value_bonferroni < 0.05  ~ "*",
      TRUE                       ~ "ns"
    )
  )

print(m1_contrasts_df)
write_csv(m1_contrasts_df, file.path(output_folder, "model1_species_contrasts.csv"))
cat("\nSaved: model1_species_contrasts.csv\n")


#fligner test of variation for non parametric data
mut_classes <- levels(m1_pred$mutation_class)
pairs <- combn(mut_classes, 2, simplify = FALSE)

pair_results <- lapply(pairs, function(pair) {
  sub <- m1_pred %>% filter(mutation_class %in% pair)
  test <- fligner.test(predicted ~ mutation_class, data = sub)
  data.frame(
    class_1  = pair[1],
    class_2  = pair[2],
    statistic = round(test$statistic, 4),
    p_raw    = test$p.value
  )
})

pair_df <- bind_rows(pair_results)

# Bonferroni correction
pair_df$p_corrected <- p.adjust(pair_df$p_raw, method = "bonferroni")
pair_df$significant <- pair_df$p_corrected < 0.05

print(pair_df %>% arrange(p_corrected))

write_csv(as.data.frame(m1_emm), file.path(output_folder, "model1_predicted_rates.csv"))

