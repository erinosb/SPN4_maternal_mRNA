Assessing the single-cell specific abundance of lin-41 and set-3 from
early embryonic cells
================
Erin Nishimura
February 9, 2026

- [Load required packages](#load-required-packages)
- [Import smFISH abundance data](#import-smfish-abundance-data)
- [Reformat the smFISH abundance
  data](#reformat-the-smfish-abundance-data)
- [Import the cell size metrics](#import-the-cell-size-metrics)
- [Filter by lin-41 counts and plot](#filter-by-lin-41-counts-and-plot)
- [Filter by set-3 counts and plot](#filter-by-set-3-counts-and-plot)
- [Print and save relevant plots](#print-and-save-relevant-plots)
- [Session info](#session-info)

The purpose of this script is to assess the abundance of lin-41 mRNA and
set-3 mRNA, not on a whole embryo level, but on a cell-specific level.
The wormLib package currently tabulates smFISH microscopy images for
abundance based on individual blastomeres for 2-cell and 4-cell stage
embryos. Beyond this, the cells are too difficult to differentiate
without additional markers.

This project is being conducted in the revisions of the Spike et al
manuscript.

This analysis will appear in Figure 2 of Spike et al.

# Load required packages

- tidyverse
- gridExtra

# Import smFISH abundance data

- This is output from wormLib (bigFISH)

``` r
# read in the mRNA abundance data from smFISH quantifications:
FISH <- read.csv(file = "../01_input/combined_cell_data.csv")

# peek at the data:
head(FISH)
```

    ##        ImageID set3_mRNA lin41_mRNA label confidence
    ## 1 230521_N2_01       330        915    P1       0.59
    ## 2 230521_N2_01       328        585    AB       0.89
    ## 3 230521_N2_03       543       1377    P1         NA
    ## 4 230521_N2_03       540        936    AB         NA
    ## 5 230521_N2_04       759       2173    AB       0.91
    ## 6 230521_N2_04       502       1729    P1       0.86
    ##                                  csv_name                subdirectory
    ## 1 quantification_cell_counts_230521_N2_01 2-cell_230521_N2_lin41_set3
    ## 2 quantification_cell_counts_230521_N2_01 2-cell_230521_N2_lin41_set3
    ## 3 quantification_cell_counts_230521_N2_03 2-cell_230521_N2_lin41_set3
    ## 4 quantification_cell_counts_230521_N2_03 2-cell_230521_N2_lin41_set3
    ## 5 quantification_cell_counts_230521_N2_04 2-cell_230521_N2_lin41_set3
    ## 6 quantification_cell_counts_230521_N2_04 2-cell_230521_N2_lin41_set3
    ##   confidenceInWrongLabel confidenceInWrongPrediction erm1_mRNA
    ## 1                     NA                          NA        NA
    ## 2                     NA                          NA        NA
    ## 3                   0.87                          NA        NA
    ## 4                   0.85                          NA        NA
    ## 5                     NA                          NA        NA
    ## 6                     NA                          NA        NA
    ##                                       rep cell_stage
    ## 1 quantification_cell_counts_230521_N2_01     2-cell
    ## 2 quantification_cell_counts_230521_N2_01     2-cell
    ## 3 quantification_cell_counts_230521_N2_03     2-cell
    ## 4 quantification_cell_counts_230521_N2_03     2-cell
    ## 5 quantification_cell_counts_230521_N2_04     2-cell
    ## 6 quantification_cell_counts_230521_N2_04     2-cell

``` r
# data structure:
str(FISH)
```

    ## 'data.frame':    108 obs. of  12 variables:
    ##  $ ImageID                    : chr  "230521_N2_01" "230521_N2_01" "230521_N2_03" "230521_N2_03" ...
    ##  $ set3_mRNA                  : int  330 328 543 540 759 502 539 461 578 803 ...
    ##  $ lin41_mRNA                 : int  915 585 1377 936 2173 1729 1082 1407 1754 1210 ...
    ##  $ label                      : chr  "P1" "AB" "P1" "AB" ...
    ##  $ confidence                 : num  0.59 0.89 NA NA 0.91 0.86 NA NA 0.77 0.88 ...
    ##  $ csv_name                   : chr  "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_03" "quantification_cell_counts_230521_N2_03" ...
    ##  $ subdirectory               : chr  "2-cell_230521_N2_lin41_set3" "2-cell_230521_N2_lin41_set3" "2-cell_230521_N2_lin41_set3" "2-cell_230521_N2_lin41_set3" ...
    ##  $ confidenceInWrongLabel     : num  NA NA 0.87 0.85 NA NA 0.86 0.88 NA NA ...
    ##  $ confidenceInWrongPrediction: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ erm1_mRNA                  : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ rep                        : chr  "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_03" "quantification_cell_counts_230521_N2_03" ...
    ##  $ cell_stage                 : chr  "2-cell" "2-cell" "2-cell" "2-cell" ...

# Reformat the smFISH abundance data

1.  There is an issue where ABp cells are listed as ABb cells. Fix it.
2.  Pivot longer to unify mRNA columns
3.  Add lineage info
4.  Select only desired columns
5.  Set factors

``` r
# === 1. There is an issue where ABp cells are listed as ABb cells. Fix it. ===
FISH <- FISH %>%
  mutate(label = case_when(label == "ABb" ~ "ABp" , .default = label))
unique(FISH$label)
```

    ## [1] "P1"  "AB"  "ABa" "EMS" "ABp" "P2"

``` r
# === 2. Pivot longer to unify mRNA columns ===
FISH_long <- FISH %>%
  pivot_longer(
    cols = ends_with("_mRNA"),
    names_to = "mRNA_name",
    values_to = "abundance"
  )

# === 3. Add lineage info ===
FISH_long <- FISH_long %>%
  mutate(
    lineage = case_when(
      label %in% c("AB", "ABa", "ABp", "EMS") ~ "Somatic",
      label %in% c("P1", "P2", "P3", "P4") ~ "Germline",
      TRUE ~ NA_character_
    ),
    lineage = factor(lineage, levels = c("Somatic", "Germline"))  # Force order
  )

# === 4. Select only desired columns ===
FISH_filtered <- FISH_long %>%
  select(ImageID, label, rep, cell_stage, mRNA_name, abundance, lineage) %>%
  filter(!is.na(abundance))

# === 5. Set factors ===
FISH_filtered$ImageID <- factor(FISH_filtered$ImageID)
FISH_filtered$label <- factor(FISH_filtered$label, level = c("AB", "P1", "ABa", "ABp", "EMS", "P2"))
FISH_filtered$mRNA_name <- factor(FISH_filtered$mRNA_name)
FISH_filtered$cell_stage <- factor(FISH_filtered$cell_stage)

# Final dataset
head(FISH_filtered)
```

    ## # A tibble: 6 × 7
    ##   ImageID      label rep                  cell_stage mRNA_name abundance lineage
    ##   <fct>        <fct> <chr>                <fct>      <fct>         <int> <fct>  
    ## 1 230521_N2_01 P1    quantification_cell… 2-cell     set3_mRNA       330 Germli…
    ## 2 230521_N2_01 P1    quantification_cell… 2-cell     lin41_mR…       915 Germli…
    ## 3 230521_N2_01 AB    quantification_cell… 2-cell     set3_mRNA       328 Somatic
    ## 4 230521_N2_01 AB    quantification_cell… 2-cell     lin41_mR…       585 Somatic
    ## 5 230521_N2_03 P1    quantification_cell… 2-cell     set3_mRNA       543 Germli…
    ## 6 230521_N2_03 P1    quantification_cell… 2-cell     lin41_mR…      1377 Germli…

``` r
str(FISH_filtered)
```

    ## tibble [214 × 7] (S3: tbl_df/tbl/data.frame)
    ##  $ ImageID   : Factor w/ 42 levels "230521_N2_01",..: 1 1 1 1 2 2 2 2 3 3 ...
    ##  $ label     : Factor w/ 6 levels "AB","P1","ABa",..: 2 2 1 1 2 2 1 1 1 1 ...
    ##  $ rep       : chr [1:214] "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_01" "quantification_cell_counts_230521_N2_01" ...
    ##  $ cell_stage: Factor w/ 2 levels "2-cell","4-cell": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ mRNA_name : Factor w/ 3 levels "erm1_mRNA","lin41_mRNA",..: 3 2 3 2 3 2 3 2 3 2 ...
    ##  $ abundance : int [1:214] 330 915 328 585 543 1377 540 936 759 2173 ...
    ##  $ lineage   : Factor w/ 2 levels "Somatic","Germline": 2 2 1 1 2 2 1 1 1 1 ...

# Import the cell size metrics

- These were part of the WormLib (bigFISH) output
- Import cell sizes
- Calculate mean cell volume per cell type
- Merge cell volume data with the FISH Abundance data

``` r
# Read in cell size estimates
cellsizes <- read.table(file = "../01_input/260209_merged_output_cellsizes_shorter.txt", header = TRUE)
head(cellsizes)
```

    ##        ImageID label perimeter major_axis_length minor_axis_length
    ## 1 230521_N2_01    P1  875.5950          306.4911          185.9879
    ## 2 230521_N2_01    AB 1009.4356          322.1181          262.1106
    ## 3 230521_N2_03    P1  813.9209          258.1622          216.2772
    ## 4 230521_N2_03    AB  979.5605          291.6551          263.8667
    ## 5 230521_N2_04    AB  866.5067          262.2477          242.3268
    ## 6 230521_N2_04    P1  662.8843          197.5687          169.9625

``` r
# Merge the FISH Abundance with the cell sizes data:
FISHcounts_with_sizes <- left_join(FISH_filtered, cellsizes)

#print(FISHcounts_with_sizes, n = 214, keep = TRUE)

# Calculate the mean size of each cell type
size_per_cell <- FISHcounts_with_sizes %>%
  group_by(label) %>%
  summarise(circumference = mean(perimeter, na.rm = TRUE)) %>%
  mutate(volume = circumference^3 / 6*(pi^2) )

# append the mean cell circumferences and volumes associated with each cell type back into the FISH_counts_with_sizes data:
FISHcounts_total_sizes <- left_join(FISHcounts_with_sizes, size_per_cell)
FISHcounts_total_sizes$ImageID <- factor(FISHcounts_total_sizes$ImageID)
FISHcounts_total_sizes$label <- factor(FISHcounts_total_sizes$label, level = c("AB", "P1", "ABa", "ABp", "EMS", "P2"))

head(FISHcounts_total_sizes)
```

    ## # A tibble: 6 × 12
    ##   ImageID      label rep        cell_stage mRNA_name abundance lineage perimeter
    ##   <fct>        <fct> <chr>      <fct>      <fct>         <int> <fct>       <dbl>
    ## 1 230521_N2_01 P1    quantific… 2-cell     set3_mRNA       330 Germli…      876.
    ## 2 230521_N2_01 P1    quantific… 2-cell     lin41_mR…       915 Germli…      876.
    ## 3 230521_N2_01 AB    quantific… 2-cell     set3_mRNA       328 Somatic     1009.
    ## 4 230521_N2_01 AB    quantific… 2-cell     lin41_mR…       585 Somatic     1009.
    ## 5 230521_N2_03 P1    quantific… 2-cell     set3_mRNA       543 Germli…      814.
    ## 6 230521_N2_03 P1    quantific… 2-cell     lin41_mR…      1377 Germli…      814.
    ## # ℹ 4 more variables: major_axis_length <dbl>, minor_axis_length <dbl>,
    ## #   circumference <dbl>, volume <dbl>

``` r
dim(FISHcounts_total_sizes)
```

    ## [1] 214  12

# Filter by lin-41 counts and plot

``` r
# Filter by lin-41
lin41_selection <- FISHcounts_total_sizes %>%
  filter(mRNA_name == "lin41_mRNA")

# Are there any unmatched entries?
unmatched <- lin41_selection %>%
  filter(is.na(perimeter))
length(unique(unmatched$ImageID))
```

    ## [1] 0

``` r
# 0

# Divide the abundance by the volume:
scaled_lin41_selection <- lin41_selection %>%
  mutate(abundance_by_vol = abundance / volume)

print(scaled_lin41_selection, width = Inf)
```

    ## # A tibble: 66 × 13
    ##    ImageID      label rep                                     cell_stage
    ##    <fct>        <fct> <chr>                                   <fct>     
    ##  1 230521_N2_01 P1    quantification_cell_counts_230521_N2_01 2-cell    
    ##  2 230521_N2_01 AB    quantification_cell_counts_230521_N2_01 2-cell    
    ##  3 230521_N2_03 P1    quantification_cell_counts_230521_N2_03 2-cell    
    ##  4 230521_N2_03 AB    quantification_cell_counts_230521_N2_03 2-cell    
    ##  5 230521_N2_04 AB    quantification_cell_counts_230521_N2_04 2-cell    
    ##  6 230521_N2_04 P1    quantification_cell_counts_230521_N2_04 2-cell    
    ##  7 230521_N2_07 AB    quantification_cell_counts_230521_N2_07 2-cell    
    ##  8 230521_N2_07 P1    quantification_cell_counts_230521_N2_07 2-cell    
    ##  9 230521_N2_12 P1    quantification_cell_counts_230521_N2_12 2-cell    
    ## 10 230521_N2_12 AB    quantification_cell_counts_230521_N2_12 2-cell    
    ##    mRNA_name  abundance lineage  perimeter major_axis_length minor_axis_length
    ##    <fct>          <int> <fct>        <dbl>             <dbl>             <dbl>
    ##  1 lin41_mRNA       915 Germline      876.              306.              186.
    ##  2 lin41_mRNA       585 Somatic      1009.              322.              262.
    ##  3 lin41_mRNA      1377 Germline      814.              258.              216.
    ##  4 lin41_mRNA       936 Somatic       980.              292.              264.
    ##  5 lin41_mRNA      2173 Somatic       867.              262.              242.
    ##  6 lin41_mRNA      1729 Germline      663.              198.              170.
    ##  7 lin41_mRNA      1082 Somatic       953.              295.              230.
    ##  8 lin41_mRNA      1407 Germline      804.              264.              202.
    ##  9 lin41_mRNA      1754 Germline      790.              256.              204.
    ## 10 lin41_mRNA      1210 Somatic       921.              262.              251.
    ##    circumference      volume abundance_by_vol
    ##            <dbl>       <dbl>            <dbl>
    ##  1          790.  812023187.      0.00000113 
    ##  2          960. 1453316899.      0.000000403
    ##  3          790.  812023187.      0.00000170 
    ##  4          960. 1453316899.      0.000000644
    ##  5          960. 1453316899.      0.00000150 
    ##  6          790.  812023187.      0.00000213 
    ##  7          960. 1453316899.      0.000000745
    ##  8          790.  812023187.      0.00000173 
    ##  9          790.  812023187.      0.00000216 
    ## 10          960. 1453316899.      0.000000833
    ## # ℹ 56 more rows

``` r
# Plot abundance scaled by estimated cell type volume 

u <- ggplot(scaled_lin41_selection, aes(x=label, y=abundance_by_vol, fill=lineage)) + 
    geom_boxplot() +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 1) +
    labs(title = "lin-41 mRNA abundance by volume",
       x = "Cell Identity",
       y = "lin-41 mRNA, relative concentration per volume")
u
```

![](260209_blastomere_specific__files/figure-gfm/lin41-1.png)<!-- -->

# Filter by set-3 counts and plot

``` r
# Filter by lin-41
set3_selection <- FISHcounts_total_sizes %>%
  filter(mRNA_name == "set3_mRNA")

#dim(set3_selection)
# 106 rows

# Divide the abundance by the volume:
scaled_set3_selection <- set3_selection %>%
  mutate(abundance_by_vol = abundance / volume)

#print(scaled_set3_selection, width = Inf)

# Plot abundance scaled by estimated cell type volume 

v <- ggplot(scaled_set3_selection, aes(x=label, y=abundance_by_vol, fill=lineage)) + 
    geom_boxplot() +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 1) +
    labs(title = "set-3 mRNA abundance by volume",
       x = "Cell Identity",
       y = "set-3 mRNA, relative concentration per volume")
v
```

![](260209_blastomere_specific__files/figure-gfm/set3-1.png)<!-- -->

# Print and save relevant plots

``` r
# Set today's date
today <- format(Sys.Date(), "%y%m%d")

# set plotname
filename <- paste("../03_output/", today, "_scaled_boxplots.pdf", sep = "")

# save pdf
pdf(filename, width = 10, height = 6)

grid.arrange(u, v, ncol = 2)

dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Session info

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS 26.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Denver
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3   lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
    ##  [5] dplyr_1.1.4     purrr_1.1.0     readr_2.1.5     tidyr_1.3.1    
    ##  [9] tibble_3.3.0    ggplot2_3.5.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.3.1     tidyselect_1.2.1   scales_1.4.0      
    ##  [5] yaml_2.3.10        fastmap_1.2.0      R6_2.6.1           labeling_0.4.3    
    ##  [9] generics_0.1.4     knitr_1.50         pillar_1.11.0      RColorBrewer_1.1-3
    ## [13] tzdb_0.5.0         rlang_1.1.6        utf8_1.2.6         stringi_1.8.7     
    ## [17] xfun_0.52          timechange_0.3.0   cli_3.6.5          withr_3.0.2       
    ## [21] magrittr_2.0.3     digest_0.6.37      grid_4.3.1         rstudioapi_0.17.1 
    ## [25] hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5        evaluate_1.0.4    
    ## [29] glue_1.8.0         farver_2.1.2       rmarkdown_2.29     tools_4.3.1       
    ## [33] pkgconfig_2.0.3    htmltools_0.5.8.1
