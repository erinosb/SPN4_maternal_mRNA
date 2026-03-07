Fig8_chs-1_smFOX
================
Naly Torres
2025-07-22

- [1. Read input counts data](#1-read-input-counts-data)
- [2. Plot total mRNA abundance per
  strain](#2-plot-total-mrna-abundance-per-strain)
- [3. Plot mean abundace of groups across cell
  stages](#3-plot-mean-abundace-of-groups-across-cell-stages)
  - [4. Statistical analysis](#4-statistical-analysis)
- [4A. Normally distributed data](#4a-normally-distributed-data)
  - [more statistics - Check for normality. try other
    tests](#more-statistics---check-for-normality-try-other-tests)

``` r
library(knitr)
library(ggplot2)
library(dplyr)
library(rstatix)
```

### 1. Read input counts data

``` r
# Define relative input folder path
input_path <- "../01_input/Fig6_chs-1-abundance_combined_quantification.csv"

# Read the CSV file
df <- read.csv(input_path, stringsAsFactors = FALSE)

# Check data
#head(df)
#str(df)
# knitr::kable(df)
```

### 2. Plot total mRNA abundance per strain

``` r
# Define the custom order for subdirectories
custom_order <- c("2-cell_N2","2-cell_DG5972", "4-cell_N2", "4-cell_DG5972")

# Convert subdirectory to a factor with custom order
df$subdirectory <- factor(df$subdirectory, levels = custom_order)

# Define custom fill colors for the subdirectories
fill_colors <- c(
  "2-cell_N2" = "gray", 
  "2-cell_DG5972" = "purple", 
  "4-cell_N2" = "gray", 
  "4-cell_DG5972" = "purple"
)

# Create a combined violin and boxplot (flipped vertically)
p <- ggplot(df, aes(x = factor(subdirectory, levels = custom_order), y = `chs_1_mRNA_total_molecules`, fill = subdirectory)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.4, width = 0.8, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black") +
  geom_jitter(height = 0, width = 0.2, color = "black", size = 0.5) +
  labs(title = "Abundance of chs-1 mRNA in FOX mutants",
       x = "3'UTR deletion strain",
       y = "chs-1 mRNA abundance \n (total molecules)") +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),  # Rotate labels on x-axis
        legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 12),
        axis.ticks = element_blank()) +
  coord_cartesian(ylim = c(0, 15000))

# Display the plot
print(p)
```

![](Fig6_chs-1_abundance_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Save the plot as PNG with the date in the filename
ggsave(("../03_output/Fig5_chs-1-smFOX-mRNA-abundance-violin-plot.png"),
       p, width = 10, height = 6)
```

### 3. Plot mean abundace of groups across cell stages

``` r
# Calculate the mean abundance for each combination of cell stage and strain
mean_data <- df %>%
  group_by(Cell_Stage, Strain) %>%
  summarise(mean_abundance = mean(`chs_1_mRNA_total_molecules`)) %>%
  arrange(Strain)  # Arrange the data by strain for easier segmentation

# Create a line plot
p <- ggplot(mean_data, aes(x = Cell_Stage, y = mean_abundance, color = Strain)) +
  geom_line() +
  geom_point() +
  # Add horizontal segments to connect points of the same strain
  geom_segment(data = mean_data %>% group_by(Strain) %>% slice(1) %>% filter(!is.na(lead(Cell_Stage))),
               aes(x = Cell_Stage, xend = lead(Cell_Stage), y = mean_abundance, yend = mean_abundance, color = Strain)) +
  labs(title = "Mean Abundance of chs-1 mRNA per Cell Stage",
       x = "Cell Stage",
       y = "Mean chs-1 mRNA abundance",
       color = "Strain") +
  theme_minimal() +
  theme(legend.position = "right")

# Display the plot
print(p)
```

![](Fig6_chs-1_abundance_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Save the plot as PNG with the date in the filename
ggsave(("../03_output/Fig5_chs-1-smFOX-mRNA-abundance-line-plot.png"),
       p, width = 10, height = 6)
```

#### 4. Statistical analysis

### 4A. Normally distributed data

ANOVA: is there a difference amongst the groups? Tukey test: Which
groups are different?

``` r
# Test significantly differences in chs-1 mRNA abundance between each "subdirectory" ( dev. stage + strain)
# Perform ANOVA
anova_result <- aov(`chs_1_mRNA_total_molecules` ~ subdirectory, data = df)

# Get ANOVA summary as a table
anova_summary <- summary(anova_result)[[1]]

# Convert to data frame
anova_df <- as.data.frame(anova_summary)

# Add row names as a new column (e.g., Factor names)
anova_df$Term <- rownames(anova_df)

# Reorder columns to put Term first
anova_df <- anova_df[, c("Term", setdiff(names(anova_df), "Term"))]

# Save to CSV
output_path_anova <- "../03_output/Fig5_chs-1-smFOX-mRNA-abundance_anova_summary.csv"
write.csv(anova_df, file = output_path_anova, row.names = FALSE)
# knitr::kable(anova_df)


# Perform Tukey's post hoc test
tukey_result <- TukeyHSD(anova_result)

# Extract Tukey results into a data frame
tukey_df <- as.data.frame(tukey_result$subdirectory)

# Add row names as a new column to identify the comparisons
tukey_df$Comparison <- rownames(tukey_df)

# Reorder columns to have 'Comparison' first
tukey_df <- tukey_df[, c("Comparison", setdiff(names(tukey_df), "Comparison"))]

# Save to CSV
output_path <- "../03_output/Fig5_chs-1-smFOX-mRNA-abundance_tukey_posthoc_results.csv"
write.csv(tukey_df, file = output_path, row.names = FALSE)
# knitr::kable(tukey_df)
```

#### more statistics - Check for normality. try other tests

Shapiro-Wilk: To check if your data in each group is normally
distributed. If norma you might use ANOVA. Kruskal-Wallis
(non-parametric version of ANOVA): Check if any of the groups are
significantly different from each other. Dunn’s test (non-parametric
post-hoc test): Compare all pairs of groups to find out what groups are
different.

``` r
# Run the Shapiro-Wilk test for normality on each group
shapiro_df <- df %>%
  group_by(subdirectory) %>%
  summarise(shapiro_p = shapiro_test(chs_1_mRNA_total_molecules)$p.value)

# Save Shapiro-Wilk results
shapiro_path <- "../03_output/Fig5_chs-1-smFOX-mRNA-abundance_shapiro_results.csv"
write.csv(shapiro_df, shapiro_path, row.names = FALSE)
# knitr::kable(shapiro_df)


# If normality assumptions are violated, use non-parametric tests
# Kruskal-Wallis test
kruskal_test_results <- df %>%
  kruskal_test(chs_1_mRNA_total_molecules ~ subdirectory)

# Convert kruskal test to data frame for saving
kruskal_df <- as.data.frame(kruskal_test_results)

kruskal_path <- "../03_output/Fig5_chs-1-smFOX-mRNA-abundance_kruskal_test_results.csv"
write.csv(kruskal_df, kruskal_path, row.names = FALSE)
# knitr::kable(kruskal_df)

# Perform Dunn's post-hoc test if Kruskal-Wallis is significant
if (kruskal_test_results$p < 0.05) {
  dunn_df <- df %>%
    dunn_test(chs_1_mRNA_total_molecules ~ subdirectory, p.adjust.method = "bonferroni")
  
  dunn_path <- "../03_output/Fig5_chs-1-smFOX-mRNA-abundance_dunn_posthoc_results.csv"
  write.csv(dunn_df, dunn_path, row.names = FALSE)
} else {
  cat("No significant differences among groups based on the Kruskal-Wallis test.\n")
}
# knitr::kable(dunn_df)
```

Based Shapiro-Wilk p-value, groups 2-cell_N2, 2-cell_DG3913, and
2-cell_DG5410 are not normal.

kruskal_test_results p \< 0.0001 show statistically significant
difference between the groups in terms of mRNA molecule counts.
