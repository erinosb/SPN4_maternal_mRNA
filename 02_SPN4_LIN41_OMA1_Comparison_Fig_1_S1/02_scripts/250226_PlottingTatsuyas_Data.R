# Plotting Tatusya's Data
#---
#  title: "250709_PlottingTatsuyas_Data.R"
#author: "Erin Nishimura"
#date: "07/09/2025"
#
#---

########### INSTALL PACKAGES ####################

#install.package(sm)
#install.packages("hrbrthemes")


###########  LOAD PACKAGES  #####################

library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(apeglm)
library(DESeq2)
library(extrafont)
library(UpSetR)

###########  READ IN THE DATA  #####################

# import the counts data
getwd()
setwd("~/Library/CloudStorage/Dropbox/github/SPN4_maternal_mRNA/02_SPN4_LIN41_OMA1_Comparison_Fig_1_S1")

#Import the combined "Enrichment" and "Probe" information. There are some duplicates
getwd()
unfiltered_SPN4data <- read.table(file = "01_input/Tsukamoto_SuppData_SPN4_LIN41_OMA1_IP_RNAseq2.txt", header = TRUE, fill = TRUE) 

# EDA
dim(unfiltered_SPN4data)
str(unfiltered_SPN4data)

# There are 46122 genes
head(unfiltered_SPN4data)

# Remove miRNAs, cTel, poorly annotated genes:
colnames(unfiltered_SPN4data)
MostGenesSPN4 <- unfiltered_SPN4data %>%
  filter(!str_detect(name, pattern = "cTel")) %>%
  filter(!str_detect(name, pattern = "21u")) %>%
  filter(!is.na(spn4_enrichment))

# Remove low expression genes:
AllGenes_SPN4 <- MostGenesSPN4 %>%
  mutate(mean = (Log2_SPN.4_LYSATE_FPKM + Log2_OMA.1_LYSATE_FPKM +  Log2_LIN.41_LYSATE_FPKM) / 3) %>%
  filter( mean > -2.5)

str(AllGenes_SPN4)


# Now there are 16,917 genes in the AllGenes_SPN4 object
str(AllGenes_SPN4)


###########  EXPLORATORY DATA ANALYSIS  #####################

# Data Structure and organization
head(AllGenes_SPN4)
colnames(AllGenes_SPN4)

# Range of enrichment levels
summary(AllGenes_SPN4$spn4_enrichment)
summary(AllGenes_SPN4$oma1_enrichment)
summary(AllGenes_SPN4$lin41_enrichment)

# Create a wide version of the data with only the enrichment values
wide_SPN4 <- AllGenes_SPN4 %>% 
  select(gene_ID, name, spn4_enrichment, oma1_enrichment, lin41_enrichment)

# Create a long version of the data with only the enrichment values
long_SPN4 <- pivot_longer(wide_SPN4,
             cols = c(spn4_enrichment, oma1_enrichment, lin41_enrichment), 
             names_to = "enrichment")


long_SPN4

# EDA - XY-Scatter Plot

## Calculate the mean 
AllGenes_SPN4_PlusMean <- cbind(AllGenes_SPN4, MEAN_Log2_SPN4_FPKM =  apply(AllGenes_SPN4[, 9:10], 1, mean))

## MA-Plot for Figure 1:
plot(AllGenes_SPN4_PlusMean$MEAN_Log2_SPN4_FPKM, AllGenes_SPN4_PlusMean$spn4_enrichment,
     ylim = c(-8, 8),
     col = ifelse(AllGenes_SPN4_PlusMean$SPN4 == TRUE,'#FF000040','#00000020'), 
     pch=20,  
     xlab = "Mean Intensity (Log2 FPKM)",
     ylab = "Fold Enrichment (SPN-4 IP over input)", 
     las = 1)

abline(h=0, col = rgb(255, 0, 0, max = 255, alpha = 100), lwd = 2)
abline(h=-2, col = rgb(0, 0, 255, max = 255, alpha = 100), lty = "dotted", lwd = 2)
abline(h=2, col = rgb(0, 0, 255, max = 255, alpha = 100), lty = "dotted", lwd = 2)

# Annotate the assayed genes on the plot
assayedGenes <- c("chs-1", "lin-41", "cpg-2", "oma-2", "car-1", "mex-5", "C04B4.2", "nos-1", "egg-3", "Y37H2A.12", "ZK666.4", "pigv-1", "npr-35", "nasp-2", "R05H11.1", "cgh-1", "mcm-2", "vab-2", "Y19D2B.2")

length(assayedGenes)
length(intersect(AllGenes_SPN4_PlusMean$name, assayedGenes) )

filtered_SPN4_PlusMean <- AllGenes_SPN4_PlusMean %>%
  filter(AllGenes_SPN4_PlusMean$name %in% assayedGenes)

points(filtered_SPN4_PlusMean$MEAN_Log2_SPN4_FPKM, 
       filtered_SPN4_PlusMean$spn4_enrichment,
       col = '#000000')

text(filtered_SPN4_PlusMean$MEAN_Log2_SPN4_FPKM, 
     filtered_SPN4_PlusMean$spn4_enrichment, 
     filtered_SPN4_PlusMean$name,
     pos = 4, 
     cex = 0.75)

# Save a plot 
today <- format(Sys.Date(),"%Y%m%d")
file1 <- paste("03_figures/", today, "_MA_plot.pdf", sep = "")
pdf(file1, height = 6, width = 6)
par(mfrow=c(1,1), pty="s")

plot(AllGenes_SPN4_PlusMean$MEAN_Log2_SPN4_FPKM, AllGenes_SPN4_PlusMean$spn4_enrichment,
     ylim = c(-8, 8),
     col = ifelse(AllGenes_SPN4_PlusMean$SPN4 == TRUE,'#FF000040','#00000020'), 
     pch=20,  
     xlab = "Mean Intensity (Log2 FPKM)",
     ylab = "Fold Enrichment (SPN-4 IP over input)", 
     las = 1)

abline(h=0, col = rgb(255, 0, 0, max = 255, alpha = 100), lwd = 2)
abline(h=-2, col = rgb(0, 0, 255, max = 255, alpha = 100), lty = "dotted", lwd = 2)
abline(h=2, col = rgb(0, 0, 255, max = 255, alpha = 100), lty = "dotted", lwd = 2)

# Annotate the assayed genes on the plot
filtered_SPN4_PlusMean <- AllGenes_SPN4_PlusMean %>%
  filter(AllGenes_SPN4_PlusMean$name %in% assayedGenes)

points(filtered_SPN4_PlusMean$MEAN_Log2_SPN4_FPKM, 
       filtered_SPN4_PlusMean$spn4_enrichment,
       col = '#000000')

text(filtered_SPN4_PlusMean$MEAN_Log2_SPN4_FPKM, 
     filtered_SPN4_PlusMean$spn4_enrichment, 
     filtered_SPN4_PlusMean$name,
     pos = 4, 
     cex = 0.75)

dev.off()

# EDA - Create histograms

# These are useful to ensure that our data is "well behaved". There are no associated figures or supplemental figures for this analysis.

par(mfrow=c(1,1))
histo1 <- ggplot(data=long_SPN4, aes(x=value, group=enrichment, fill=enrichment)) +
  geom_density(alpha=.3, adjust = 2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  ylab("density") +
  xlab("enrichment")

histo1


histo2 <- ggplot(data=long_SPN4, aes(x=value, group=enrichment, fill=enrichment)) +
  geom_density(alpha=.2, adjust = 2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  facet_wrap(~enrichment) +
  ylab("density") +
  xlab("enrichment")

histo2

############### CORR Matrix ############################

# SUPPLEMENTAL FIGURE 1
# This is a nice correlation matrix 
# Create a matrix

head(AllGenes_SPN4)
str(AllGenes_SPN4)
dim(AllGenes_SPN4)

# Select for enrichments over -2
colnames(AllGenes_SPN4)
SPN4_wide_2cutoff_2FPKM <- AllGenes_SPN4 %>%
  filter(spn4_enrichment > -2) %>%
  filter(Log2_SPN.4_LYSATE_FPKM > 2)

# Create rownames
rownames(SPN4_wide_2cutoff_2FPKM) <- SPN4_wide_2cutoff_2FPKM$gene_ID



# Trim and create a matrix:
SPN4_widematrix2 <- SPN4_wide_2cutoff_2FPKM[,6:8]
SPN4_widematrix2
SPN4_widematrix3 <- as.matrix(SPN4_widematrix2)
head(SPN4_widematrix3)


# Calculate the distances between each sample
sampleDists <- dist(t(SPN4_widematrix3)) 
sampleDists

sampleDistMatrix <- as.matrix(sampleDists) #convert from data.frame -> matrix
rownames(sampleDistMatrix) <- colnames(SPN4_widematrix3) # Add some labels
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #Pick some pretty colors

# Draw the heatmap
par(mfrow=c(1,1))
p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors) # Plot the heatmap

p

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file2 <- paste("03_figures/", today, "_corr_matrix_plots.pdf", sep = "")

pdf(file2, height = 6, width = 7)
par(mfrow=c(1,1))
p

dev.off()

################# HEATMAP #################

## SUPPLEMENTAL FIGURE 1

## This is a pretty restricted heatmap. It shows only mRNAs with dramatically changing dynamics between SPN-4, LIN-41, and OMA-1. We didn't include this in the paper. Instead, we created a more expansive heatmap (see below) where there wasn't a strict changing dynamic used for filtering. 
## A nice heatmap showing enriched genes and how their enrichment differs in SPN-4, LIN-41, and OMA-1 IPs.

## REVIEW:
# head(AllGenes_SPN4)
# str(AllGenes_SPN4)
#dim(AllGenes_SPN4)

# Calculate standard deviations
SPN4_wide_stringentCutoff1 <- cbind(AllGenes_SPN4, stdev =  apply(AllGenes_SPN4[, 6:8], 1, sd))
colnames(SPN4_wide_stringentCutoff1)
head(SPN4_wide_stringentCutoff1)
summary(SPN4_wide_stringentCutoff1$stdev)
SPN4_wide_stringentCutoff1

# Make the cutoffs
## Select for enrichments over -2 and more than Log2(2)
## Select for expression levels greater than e^2(4)
## Select for stdev > 0.6
SPN4_wide_stringentCutoff2 <- SPN4_wide_stringentCutoff1 %>%
  filter(spn4_enrichment > -2) %>%
  filter(Log2_SPN.4_LYSATE_FPKM > 4) %>%
  filter(stdev > 0.6)


# Full heatmap of high-ish expression genes:
dim(SPN4_wide_stringentCutoff2)
head(SPN4_wide_stringentCutoff2)
summary(SPN4_wide_stringentCutoff2$spn4_enrichment)

SPN4_wide_stringentCutoff2

stringentMatrix <- as.matrix(SPN4_wide_stringentCutoff2[,c(6,7,8)])
rownames(stringentMatrix) <- SPN4_wide_stringentCutoff2$gene_ID

stringentMatrix
q <- pheatmap(stringentMatrix, 
              scale="row", 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE)

q


# What about clustering? How many clusters are there?
# Note - this isn't a figure in the paper

r <- pheatmap(stringentMatrix, 
              scale="row", 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              cutree_rows = 7,
              show_rownames = FALSE)

r


# append clusters into the same data frame
cl = cutree(r$tree_row,k = 7)
cl
ann = data.frame(cl)
stringentMatrix

rownames(ann) <- rownames(stringentMatrix)
SPN4_wide_2cutoff_2FPKM <- cbind(SPN4_wide_stringentCutoff2, 
                      cluster = cl)

SPN4_wide_2cutoff_2FPKM

s <- pheatmap(stringentMatrix, 
              scale="row", 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              cutree_rows = 7,
              show_rownames = FALSE,
              annotation_row=ann)

s

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file3 <- paste("03_figures/", today, "_heatmap_strictly_changing.pdf", sep = "")
pdf(file3, height = 9, width = 6)

s

dev.off()


################# HEATMAP OF GENES WITH AT LEAST ONE BINDING SITE #################

# The above heatmap was a little too stringent because it required that there be changing dynamics between the enrichment IPs. What we really want to see instead is all mRNAs represented. So, the filter will be - is a transcript IP'd by at least one of these RBPs. Also, we're going to plot this so that there is no z-scaling. This will help us better visualize when a gene's enrichment is equal across all columns.

## Supplemental Figure 1

# Calculate enrichment sum (Sum of SPN4 == TRUE +  OMA1 == TRUE LIN41 == TRUE)
# Start with the raw data unfiltered_SPN4data instead of AllGenes_SPN4 because the unfiltered_SPN4data isn't filtered by gene expressio level

unfiltered_SPN4data

unfiltered_SPN4data$LIN41 <- as.logical(unfiltered_SPN4data$LIN41)
SPN4_SUM <- cbind(unfiltered_SPN4data, sum =  apply(unfiltered_SPN4data[, 3:5], 1, sum))
colnames(SPN4_SUM)
head(SPN4_SUM)
table(SPN4_SUM$sum)

ANY_bound <- SPN4_SUM %>%
  filter(sum >= 1)
  
dim(ANY_bound)
table(ANY_bound$sum)
table(SPN4_SUM$sum)

anyBoundMatrix <- as.matrix(ANY_bound[,c(6,7,8)])
rownames(anyBoundMatrix) <- ANY_bound$name

dim(anyBoundMatrix)
head(anyBoundMatrix)

u <- pheatmap(anyBoundMatrix, 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE)

u

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file4 <- paste("03_figures/", today, "_Heatmap_BoundTargets_noScale.pdf", sep = "")

pdf(file4, height = 9, width = 6)
u
dev.off()



################# UPSET PLOT ##################

# I would like to show an upset plot to illustrate the overlap between different categories (SPN-4 enrichment, OMA-1 enrichment, LIN-41 enrichment)

## Figure 1

## Note - uses the UpSetR package 
# Start with the raw data unfiltered_SPN4data instead of AllGenes_SPN4 because the unfiltered_SPN4data isn't filtered by gene expressio level

unfiltered_SPN4data


allMatrix <- unfiltered_SPN4data[,c(3,4,5)]
rownames(allMatrix) <- unfiltered_SPN4data$name

allNumericMatrix <- 1*allMatrix

AllGenes_SPN4[(duplicated(AllGenes_SPN4$name)),]

upset(allNumericMatrix, text.scale=2)

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file5 <- paste("03_figures/", today, "_upsetPlot.pdf", sep = "")

pdf(file5, height = 4, width = 8)
upset(allNumericMatrix, text.scale=2)
dev.off()


### SAVE a bunch of lists for GO Ontology Analysis #############

# These are lists of all the sets of genes (WBGENE IDs) from all the different categories

today <- format(Sys.Date(),"%Y%m%d")
length(unfiltered_SPN4data$gene_ID)
write(unfiltered_SPN4data$gene_ID, file = paste("04_output_data/", today, "_ALLlist.txt"))

# Gene IDs - bound by a single RBP independent of any other binding
OMA1_list <- unfiltered_SPN4data %>%
  filter(OMA1 == TRUE)
OMA1_list$gene_ID
length(OMA1_list$gene_ID)
write(OMA1_list$gene_ID, file = paste("04_output_data/", today, "_OMA1_list.txt"))

SPN4_list <- unfiltered_SPN4data %>%
  filter(SPN4 == TRUE)
SPN4_list$gene_ID
length(SPN4_list$gene_ID)
write(SPN4_list$gene_ID, file = paste("04_output_data/", today, "_SPN4_list.txt"))

LIN41_list <- unfiltered_SPN4data %>%
  filter(LIN41 == TRUE)
LIN41_list$gene_ID
length(LIN41_list$gene_ID)
write(LIN41_list$gene_ID, file = paste("04_output_data/", today, "_LIN41_list.txt"))


#Gene IDs - Single RBP binding only (no overlap with another RBP)
OMA1_only <- unfiltered_SPN4data %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == FALSE)
OMA1_only$gene_ID
length(OMA1_only$gene_ID)
write(OMA1_only$gene_ID, file = paste("04_output_data/", today, "_OMA1_ONLY_list.txt"))

SPN4_only <- unfiltered_SPN4data %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == FALSE)
SPN4_only$gene_ID
length(SPN4_only$gene_ID)
write(SPN4_only$gene_ID, file = paste("04_output_data/", today, "_SPN4_ONLY_list.txt"))

LIN41_only <- unfiltered_SPN4data %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == TRUE)
LIN41_only$gene_ID
length(LIN41_only$gene_ID)
write(LIN41_only$gene_ID, file = paste("04_output_data/", today, "_LIN41_ONLY_list.txt"))

# Overlapping sets between two RBPs
LIN41_SPN4 <- unfiltered_SPN4data %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == TRUE)
LIN41_SPN4$gene_ID
length(LIN41_SPN4$gene_ID)
write(LIN41_SPN4$gene_ID, file = paste("04_output_data/", today, "_LIN41_and_SPN4_list.txt"))


OMA1_SPN4 <- unfiltered_SPN4data %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == FALSE)
OMA1_SPN4$gene_ID
length(OMA1_SPN4$gene_ID)
write(OMA1_SPN4$gene_ID, file = paste("04_output_data/", today, "_OMA1_SPN4_list.txt"))

OMA1_LIN41 <- unfiltered_SPN4data %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == TRUE)
OMA1_LIN41$gene_ID
length(OMA1_LIN41$gene_ID)
write(OMA1_LIN41$gene_ID, file = paste("04_output_data/", today, "_OMA1_LIN41_list.txt"))

OMA1_SPN4_LIN41 <- unfiltered_SPN4data %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == TRUE)
OMA1_SPN4_LIN41$gene_ID
length(OMA1_SPN4_LIN41$gene_ID)
write(OMA1_SPN4_LIN41$gene_ID, file = paste("04_output_data/", today, "_OMA1_SPN4_LIN41_list.txt"))

################# SELECTION OF STATISTICS  #################

# From Paper draft:

# Overall, the overlap between these mRNA cohorts was greatest between LIN-41 and OMA-1 (Figure 1F, Supplementary Figure S1A, Supplementary Figure S1B). 30 % of the OMA-1 targets were shared with LIN-41 and, inversely, 40 % of LIN-41 targets were shared with OMA-1. In contrast, 14 – 25 % of mRNA sets were overlapping in combinations comparing SPN-4 associated mRNAs to either LIN-41 or OMA-1 sets. 

# How many OMA-1 targets were shared with LIN-41?
round((length(OMA1_LIN41$gene_ID) / length(OMA1_list$gene_ID ))*100, 2)

# How many LIN-411 targets were shared withOMA-1?
round((length(OMA1_LIN41$gene_ID) / length(LIN41_list$gene_ID ))*100, 2)

# How many SPN-4 targets are shared with LIN-41?
round((length(LIN41_SPN4$gene_ID) / length(SPN4_list$gene_ID ))*100, 2)

# How many SPN-4 targets are shared with OMA-1?
round((length(OMA1_SPN4$gene_ID) / length(SPN4_list$gene_ID ))*100, 2)

# How many LIN-41 targets are shared with SPN-4?
round((length(LIN41_SPN4$gene_ID) / length(LIN41_list$gene_ID ))*100, 2)

# How many OMA-1 targets are shared with SPN-4?
round((length(LIN41_SPN4$gene_ID) / length(OMA1_list$gene_ID ))*100, 2)

################# HEATMAP OF SELECT GENES #################

## I want to make a heatmap of select genes that are assayed later in Figure 3

## No figure. Built into a Figure for Figure 3 in section below.

#Heatmap of key assayed genes:

assayedGenes <- c("chs-1", "lin-41", "cpg-2", "oma-2", "car-1", "mex-5", "C04B4.2", "nos-1", "egg-3", "Y37H2A.12", "ZK666.4", "pigv-1", "npr-35", "nasp-2", "R05H11.1", "cgh-1", "mcm-2", "vab-2", "Y19D2B.2")

length(assayedGenes)
length(intersect(unfiltered_SPN4data$name, assayedGenes) )

filtered_SPN4 <- unfiltered_SPN4data %>%
  filter(unfiltered_SPN4data$name %in% assayedGenes)

filtered_SPN4

# Input smFISH Fold Change data. It should be in the input folder
getwd()
smFISHFoldChange <- read.table(file = "01_input/240618_smFISH_foldChangeData.txt", header = TRUE) 
smFISHFoldChange


# Merge assayed Genes data frame (RIP-seq) and fold change data frame (smFISH)
joined_filtered_matrix1 <- left_join(smFISHFoldChange, filtered_SPN4, by = c("transcript" = "name"))
joined_filtered_matrix1 <- unique(joined_filtered_matrix1)

# Double Check the ordering is correct
joined_filtered_matrix1 <- joined_filtered_matrix1 %>%
  arrange(desc(foldChange))

# Save transcript as a rowname
rownames(joined_filtered_matrix1) <- joined_filtered_matrix1$transcript
joined_filtered_matrix1
# Subset
joined_filtered_matrix2 <- joined_filtered_matrix1[,7:9]
joined_filtered_matrix2

v <- pheatmap(joined_filtered_matrix2[,c(3,2,1)], 
              scale="none", 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              border_color = "white", 
              show_rownames = TRUE)

v

# Save the plot:
today <- format(Sys.Date(),"%Y%m%d")
file6 <- paste("03_figures/", today, "_heatmap_smFISH_selecteGenes.pdf", sep = "")

pdf(file6, height = 9, width = 6)
v
dev.off()

########## Fold CHange V. SPN Enrich V. Express. Level #######

## Figure 3

#Heatmap of key assayed genes plotted to include the expression level as well as their expression levels as this seems to make a difference

# Subset
joined_filtered_matrix2 <- joined_filtered_matrix1 %>%
  arrange(desc(foldChange))

str(joined_filtered_matrix2)

w <- pheatmap(joined_filtered_matrix2[,c(10)], 
              scale="none", 
              color = colorRampPalette(c("darkgreen", "white", "orange"), space = "Lab")(100),
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              border_color = "white", 
              show_rownames=TRUE)

w
help(pheatmap)

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file7 <- paste("03_figures/", today, "_heatmap_smFISH_over2lysate_PlusExp.pdf", sep = "")

pdf(file7, height = 9, width = 6)
w
dev.off()

######### Compare this dataset to the Tsukamoto et al., Genetics 2017 LIN-41 and OMA-1 IP + RNAseq dataset ###########

## Supplemental Figure S1 

# Merge the 2017 dataset with data from this manuscript
# import the supplemental data from Tsukamoto et al., Genetics, 2017. file3.txt: 
tsukamoto2017_data <- read.table(file = "01_input/files3_Tsukamoto2017_RNAseq.txt", header = TRUE, fill = TRUE) 

# EDA
head(tsukamoto2017_data$Oma.1_042310)
colnames(unfiltered_SPN4data)
colnames(tsukamoto2017_data)

# Merge with the data from this manuscript
joined_2017_and_2025 <- left_join(unfiltered_SPN4data, tsukamoto2017_data, by = c("gene_ID" = "Ensembl_ID"))

# EDA
dim(unfiltered_SPN4data)
dim(tsukamoto2017_data)
dim(joined_2017_and_2025)

################ OMA1 - calculate correlation:
# Calculate Correlation
cor_oma1_2017_v_2025 <- cor(joined_2017_and_2025$oma1_enrichment, joined_2017_and_2025$Oma.1_IP, , use = "pairwise.complete.obs")
cor_oma1_2017_v_2025
text_oma1 <- paste("r2 = ", round(cor_oma1_2017_v_2025, 2), sep = "")
text_oma1

# Subset OMA-1 significantly enriched transcripts from 2025 IP. Do they have similarly high levels in OMA-1 IP from 2017?
subset_oma1_true <- joined_2017_and_2025[which(joined_2017_and_2025$OMA1 == TRUE), ]
dim(subset_oma1_true)

# Plot OMA-1 correlation
par(pty="s")
plot(joined_2017_and_2025$oma1_enrichment, joined_2017_and_2025$Oma.1_IP, pch = 20, col = "#838B8B22",
     xlim = c(-2.5, 10),
     ylim = c(-2.5,10))
points(subset_oma1_true$oma1_enrichment, subset_oma1_true$Oma.1_IP, pch = 20, col = "#1E90FF44")
text(-1, 9, labels = text_oma1, col = "blue")
abline(lm(Oma.1_IP ~ oma1_enrichment, data = joined_2017_and_2025), col = "blue")


################ LIN-41 - calculate correlation:
# Calculate Correlation
cor_lin41_2017_v_2025 <- cor(joined_2017_and_2025$lin41_enrichment, joined_2017_and_2025$Lin.41_IP, , use = "pairwise.complete.obs")
cor_lin41_2017_v_2025
text_lin41 <- paste("r2 = ", round(cor_lin41_2017_v_2025, 2), sep = "")
text_lin41

# Subset LIN-41 significantly enriched transcripts from 2025 IP. Do they have similarly high levels in LIN-41 IP from 2017?
subset_lin41_true <- joined_2017_and_2025[which(joined_2017_and_2025$LIN41 == TRUE), ]
dim(subset_lin41_true)

# Plot LIN-41 correlation
par(pty="s")
plot(joined_2017_and_2025$lin41_enrichment, joined_2017_and_2025$Lin.41_IP, pch = 20, col = "#838B8B22",
     xlim = c(-2.5, 10),
     ylim = c(-2.5,10))
points(subset_lin41_true$lin41_enrichment, subset_lin41_true$Lin.41_IP, pch = 20, col = "#DA70D644")
text(-1, 9, labels = text_lin41, col = "orchid4")
abline(lm(Lin.41_IP ~ lin41_enrichment, data = joined_2017_and_2025), col = "orchid4")


# Create plot 
graphics.off()
today <- format(Sys.Date(),"%Y%m%d")
file8 <- paste("03_figures/", today, "_Compare_to_2017.pdf", sep = "")

pdf(file8, height = 4, width = 8)

par(mfrow = c(1,2), pty = "s")


plot(joined_2017_and_2025$oma1_enrichment, joined_2017_and_2025$Oma.1_IP, pch = 20, col = "#838B8B22",
     xlim = c(-2.5, 10),
     ylim = c(-2.5,10))
points(subset_oma1_true$oma1_enrichment, subset_oma1_true$Oma.1_IP, pch = 20, col = "#1E90FF44")
text(0, 9, labels = text_oma1, col = "blue")
abline(lm(Oma.1_IP ~ oma1_enrichment, data = joined_2017_and_2025), col = "blue")

plot(joined_2017_and_2025$lin41_enrichment, joined_2017_and_2025$Lin.41_IP, pch = 20, col = "#838B8B22",
     xlim = c(-2.5, 10),
     ylim = c(-2.5,10))
points(subset_lin41_true$lin41_enrichment, subset_lin41_true$Lin.41_IP, pch = 20, col = "#DA70D644")
text(0, 9, labels = text_lin41, col = "orchid4")
abline(lm(Lin.41_IP ~ lin41_enrichment, data = joined_2017_and_2025), col = "orchid4")

dev.off()


######### SAVE SESSION INFO ###############

today <- format(Sys.Date(),"%Y%m%d")
file8 <- paste("04_output_data/", today, "_sessionInfo.txt", sep = "")




writeLines(capture.output(sessionInfo()), file8)

################################
##         END SCRIPT         ##
################################
