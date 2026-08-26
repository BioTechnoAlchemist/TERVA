setwd("/Users/limikk/Desktop/TERVA/TERVA2DATA/VALIDATIONS_2023/Markus_morphometric")

library(tidyverse)

# Load data ####

#Slide 13

Section13S_1 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_13S_1/13S_section_1/DigitalSlide_C1M_13S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section13S_2 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_13S_1/13S_section_2/DigitalSlide_C1M_13S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)
Section13S_3 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_13S_1/13S_section_3/DigitalSlide_C1M_13S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section13S_4 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_13S_1/13S_section_4/DigitalSlide_C1M_13S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section13S_5 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_13S_1/13S_section_5/DigitalSlide_C1M_13S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

#Slide 18 (different name in new data)

Section18S_1 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_18S_1/18S_section_1/DigitalSlide_C1M_18S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section18S_2 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_18S_1/18S_section_2/DigitalSlide_C1M_18S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section18S_3 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_18S_1/18S_section_3/DigitalSlide_C1M_18S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section18S_4 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_18S_1/18S_section_4/DigitalSlide_C1M_18S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

Section18S_5 <- read_delim("Exports/2 FINAL analysis working folder - lipid analysis/DigitalSlide_C1M_18S_1/18S_section_5/DigitalSlide_C1M_18S_1_objects.tsv", 
                           delim = "\t", escape_double = FALSE, col_types = cols(`Study level 1` = col_skip(), 
                                                                                 `Study level 2` = col_skip(), 
                                                                                 `Study level 3` = col_skip(), 
                                                                                 Image = col_skip(), 
                                                                                 LayerData = col_skip(), 
                                                                                 `Area of each single lipid` = col_number()), 
                           trim_ws = TRUE)

# Combine data ####


list <- list(Section13S_1=Section13S_1,Section13S_2=Section13S_2,Section13S_3=Section13S_3,Section13S_4=Section13S_4,Section13S_5=Section13S_5,
             Section18S_1=Section18S_1,Section18S_2=Section18S_2,Section18S_3=Section18S_3,Section18S_4=Section18S_4,Section18S_5=Section18S_5)

list_sums <- lapply(list, function(x) {
  x <- x %>% mutate(sum_area = sum(`Area of each single lipid`), avg_area = mean(`Area of each single lipid`), sd_area = sd(`Area of each single lipid`))
})

list2env(list_sums, globalenv())

#Add status

Section13S_1$status <- as.factor("Obese")
Section13S_2$status <- as.factor("Obese")
Section13S_3$status <- as.factor("Obese")
Section13S_4$status <- as.factor("Non-obese")
Section13S_5$status <- as.factor("Non-obese")
Section18S_1$status <- as.factor("Obese")
Section18S_2$status <- as.factor("Obese")
Section18S_3$status <- as.factor("Non-obese")
Section18S_4$status <- as.factor("Non-obese")
Section18S_5$status <- as.factor("Non-obese")

# Load necessary libraries
library(ggplot2)
library(dplyr)

prepare_list <- function() {
  all_objects <- ls(envir = .GlobalEnv)
  X_names <- grep("^Section[0-9]+S_[0-9]+$", all_objects, value = TRUE)
  
  X_list <- mget(X_names, envir = .GlobalEnv)
  
  return(list(objects = X_list, names = X_names))
}

object_list <- prepare_list()

df <- bind_rows(object_list)
df_combined <- bind_rows(df$objects, .id = "sample_id") #ensures the data is combined into one big tibble, not with subclasses by samples.

df_summary <- df_combined %>%
  group_by(sample_id, status) %>%
  summarise(avg_area = mean(avg_area, na.rm = TRUE), .groups = "drop") #Make a summary to test for normality within sample.

# Normality test: Apply Shapiro-Wilk test for each sample
shapiro_results <- df_summary %>%
  group_by(status) %>%
  summarise(p_value = shapiro.test(avg_area)$p.value, .groups = "drop")

print(shapiro_results)  # Print normality test results

# Determine which test to use
if (all(shapiro_results$p_value > 0.05)) {
  # If both groups are normal, use t-test
  test_result <- t.test(avg_area ~ status, data = df_summary, var.equal = FALSE, alternative = "greater")
  test_type <- "T-Test (parametric)"
} else {
  # If not normal, use Wilcoxon test
  test_result <- wilcox.test(avg_area ~ status, data = df_summary)
  test_type <- "Wilcoxon Rank-Sum Test (non-parametric)"
}


# Print the test result and type
print(test_type)
print(test_result)

# Visualization: QQ plot for normality check
ggplot(df_summary, aes(sample = avg_area)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~status) +
  theme_minimal() +
  ggtitle("QQ Plot for Normality Check")

# Boxplot
p <- ggplot(df_summary, aes(x = status, y = avg_area, fill = status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.05, alpha = 0.5) +
  scale_fill_manual(values = c("mediumblue","goldenrod")) +
  theme_minimal() +
  labs(title = "Comparison of average lipid sizes between groups, N = 10",
       subtitle = "One-way Welch t-test: avg_area ~ status",
       x = "Group",
       y = "Average lipid area") +
  theme(legend.position = "none") +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 4, color = "black")

p + annotate("text", x = 1.5, y = max(df_summary$avg_area, na.rm = TRUE), 
             label = paste("P-value =", format(test_result$p.value, digits = 3)), 
             size = 4, color = "black")

# Bin level analysis 

library(rstatix)
library(ggplot2)
library(dplyr)
library(tidyr)

# Clean object-level data
df_lipid <- df_combined %>%
  select(sample_id, status, area = `Area of each single lipid`) %>%
  filter(!is.na(area), area > 0) %>%
  mutate(status = factor(status, levels = c("Non-obese", "Obese")))

# choose sensible bins
print(summary(df_lipid$area))
cat("\nQuantiles:\n")
print(quantile(df_lipid$area, probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)))

#    Median ≈ 1.3, Q75 ≈ 5.8, Q90 ≈ 30, Q99 ≈ 564
bin_breaks <- c(0, 0.5, 1, 3, 10, 50, 200, 1000, Inf)
bin_labels <- c("<0.5", "0.5–1", "1–3", "3–10",
                "10–50", "50–200", "200–1000", ">1000")

df_lipid_binned <- df_lipid %>%
  mutate(
    size_bin = cut(
      area,
      breaks         = bin_breaks,
      labels         = bin_labels,
      right          = FALSE,
      include.lowest = TRUE
    )
  )

# Frequency per biological replicate
# use nesting() so sample_id stays with its correct status
df_bin_freq <- df_lipid_binned %>%
  count(sample_id, status, size_bin, name = "n_cells") %>%
  complete(
    nesting(sample_id, status),
    size_bin,
    fill = list(n_cells = 0)
  ) %>%
  group_by(sample_id, status) %>%
  mutate(
    total_cells = sum(n_cells),
    percent     = ifelse(total_cells > 0, (n_cells / total_cells) * 100, 0)
  ) %>%
  ungroup()

# Verify correct sample sizes
cat("── Samples per group ──\n")
print(df_bin_freq %>% distinct(sample_id, status) %>% count(status))

# Bar plot
p_bins <- ggplot(df_bin_freq, aes(x = size_bin, y = percent, fill = status)) +
  stat_summary(
    fun      = mean,
    geom     = "bar",
    position = position_dodge(width = 0.8),
    alpha    = 0.85, width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom     = "errorbar",
    position = position_dodge(width = 0.8),
    width    = 0.25
  ) +
  scale_fill_manual(values = c("Non-obese" = "mediumblue", "Obese" = "goldenrod")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Lipid droplet size frequency distribution",
    x     = expression("Lipid area bin ("*mu*"m"^2*")"),
    y     = "Lipid droplets per bin (%, mean ± SE)",
    fill  = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_bins)

# Bin-wise Wilcoxon tests
bin_tests <- df_bin_freq %>%
  group_by(size_bin) %>%
  summarise(
    n_non_obese    = sum(status == "Non-obese"),
    n_obese        = sum(status == "Obese"),
    mean_non_obese = mean(percent[status == "Non-obese"]),
    mean_obese     = mean(percent[status == "Obese"]),
    p_value = {
      vals_no <- percent[status == "Non-obese"]
      vals_ob <- percent[status == "Obese"]
      if (sd(vals_no) == 0 & sd(vals_ob) == 0) {
        NA_real_
      } else {
        wilcox.test(vals_no, vals_ob, exact = FALSE)$p.value
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      is.na(p_adj)  ~ "NT",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )
print(as.data.frame(bin_tests))

# Sample-level distribution features
# Compute per-sample distributional summaries beyond the mean
library(moments)  # for skewness/kurtosis

df_dist_features <- df_lipid %>%
  group_by(sample_id, status) %>%
  summarise(
    n_droplets = n(),
    mean_area  = mean(area),
    median_area = median(area),
    iqr_area    = IQR(area),
    sd_area     = sd(area),
    cv_area     = sd(area) / mean(area),          # coefficient of variation
    q90         = quantile(area, 0.90),
    q95         = quantile(area, 0.95),
    q99         = quantile(area, 0.99),
    skewness    = skewness(area),                  # tail heaviness
    pct_above_50  = mean(area > 50) * 100,         # % large droplets
    pct_above_200 = mean(area > 200) * 100,        # % very large droplets
    .groups = "drop"
  )

cat("── Per-sample distribution features ──\n")
print(as.data.frame(df_dist_features))

# Test each distributional feature
features_to_test <- c("mean_area", "median_area", "iqr_area", "cv_area",
                      "q90", "q95", "q99", "skewness",
                      "pct_above_50", "pct_above_200")

feature_tests <- lapply(features_to_test, function(feat) {
  vals_no <- df_dist_features[[feat]][df_dist_features$status == "Non-obese"]
  vals_ob <- df_dist_features[[feat]][df_dist_features$status == "Obese"]
  
  wt <- wilcox.test(vals_no, vals_ob, exact = FALSE)
  
  tibble(
    feature        = feat,
    mean_non_obese = mean(vals_no),
    mean_obese     = mean(vals_ob),
    fold_change    = mean(vals_ob) / mean(vals_no),
    p_value        = wt$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

print(as.data.frame(feature_tests))

# Per-sample ECDF overlay
p_ecdf <- ggplot(df_lipid, aes(x = area, colour = status, group = sample_id)) +
  stat_ecdf(linewidth = 0.5, alpha = 0.6) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(values = c("Non-obese" = "mediumblue", "Obese" = "goldenrod")) +
  annotation_logticks(sides = "b") +
  theme_minimal(base_size = 13) +
  labs(
    title  = "Per-sample cumulative distributions of lipid droplet areas",
    x      = expression("Lipid area ("*mu*"m"^2*", log scale)"),
    y      = "Cumulative fraction",
    colour = "Group"
  )

print(p_ecdf)

# Dot plots of key distributional features
library(patchwork)

make_dot_plot <- function(data, yvar, ylabel, pval) {
  ggplot(data, aes(x = status, y = .data[[yvar]], colour = status)) +
    geom_jitter(width = 0.05, size = 3, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.5,
                 colour = "black") +
    scale_colour_manual(values = c("Non-obese" = "mediumblue", "Obese" = "goldenrod")) +
    theme_minimal(base_size = 12) +
    labs(y = ylabel, x = NULL,
         subtitle = paste0("p(adj) = ", format(pval, digits = 2))) +
    theme(legend.position = "none")
}

p1 <- make_dot_plot(df_dist_features, "median_area", 
                    expression("Median area ("*mu*"m"^2*")"),
                    feature_tests$p_adj[feature_tests$feature == "median_area"])

p2 <- make_dot_plot(df_dist_features, "cv_area", 
                    "Coefficient of variation",
                    feature_tests$p_adj[feature_tests$feature == "cv_area"])

p3 <- make_dot_plot(df_dist_features, "q95", 
                    expression("95th percentile ("*mu*"m"^2*")"),
                    feature_tests$p_adj[feature_tests$feature == "q95"])

p4 <- make_dot_plot(df_dist_features, "pct_above_200", 
                    "% droplets > 200",
                    feature_tests$p_adj[feature_tests$feature == "pct_above_200"])

p_features <- (p1 | p2 | p3 | p4) +
  plot_annotation(
    title = "Sample-level distributional features of lipid droplet size",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(p_features)

# Final annotated bar plot
p_final <- p_bins +
  labs(caption = paste0(
    "Bin-wise: Wilcoxon rank-sum (N = 5 vs 5), BH-adjusted\n",
    "No significant bin-wise differences after multiple testing correction"
  ))

print(p_final)

# Sensitivity analysis: exclude Section18S_2
df_dist_sensitivity <- df_dist_features %>%
  filter(sample_id != "Section18S_2")

cat("── Without Section18S_2 (N=5 vs 4) ──\n")
sensitivity_tests <- lapply(features_to_test, function(feat) {
  vals_no <- df_dist_sensitivity[[feat]][df_dist_sensitivity$status == "Non-obese"]
  vals_ob <- df_dist_sensitivity[[feat]][df_dist_sensitivity$status == "Obese"]
  
  wt <- wilcox.test(vals_no, vals_ob, exact = FALSE)
  
  tibble(
    feature        = feat,
    mean_non_obese = mean(vals_no),
    mean_obese     = mean(vals_ob),
    fold_change    = mean(vals_ob) / mean(vals_no),
    p_value        = wt$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

print(as.data.frame(sensitivity_tests))


library(patchwork)

supp_figure <- ((p_bins + labs(tag = "A")) + 
                  (p_ecdf + labs(tag = "B"))) /
  (wrap_elements(p_features) + labs(tag = "C")) +
  plot_layout(heights = c(1, 0.8)) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("Supplementary_Figure_LipidDistribution.tif",
       supp_figure,
       width = 14, height = 10,
       dpi = 300,
       device = "tiff",
       compression = "lzw")
