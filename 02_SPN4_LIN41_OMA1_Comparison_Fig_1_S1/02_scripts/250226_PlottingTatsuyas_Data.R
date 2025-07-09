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
enriched_SPN4 <- read.table(file = "01_input/RenamedCOLS_combined_results_w_kimbleData_and_stoeckiusData_Tue_Jun_18_2024_1953.txt", header = TRUE, fill = TRUE) 

# EDA
dim(enriched_SPN4)
str(enriched_SPN4)

# There are 46122 genes
head(enriched_SPN4)

# Remove miRNAs, cTel, poorly annotated genes:
colnames(enriched_SPN4)
AllGenes_SPN4 <- enriched_SPN4 %>%
  filter(!str_detect(name, pattern = "cTel")) %>%
  filter(!str_detect(name, pattern = "21u")) %>%
  filter(!is.na(Stoeckius_twoCell_RPKM)) %>%
  filter(!is.na(spn4_enrichment))

# Now there are 14,563 genes
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
rownames(stringentMatrix) <- SPN4_wide_stringentCutoff2$name

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

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file3 <- paste("03_figures/", today, "_heatmap_strictly_changing.pdf", sep = "")

pdf(file3, height = 9, width = 6)
q
dev.off()

# What about clustering? How many clusters are there?
# Note - this isn't a figure in the paper

r <- pheatmap(stringentMatrix, 
              scale="row", 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              cutree_rows = 5,
              show_rownames = FALSE)

r

help(pheatmap)

# append clusters into the same data frame
cl = cutree(r$tree_row,k = 5)
cl
ann = data.frame(cl)
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
              cutree_rows = 5,
              show_rownames = FALSE,
              annotation_row=ann)

s

help(pheatmap)


#write.table(SPN4_wide_2cutoff_2FPKM, file = "03_figures/240731_clusters.txt", quote = FALSE, sep = "\t")
#help(write.table)


################# HEATMAP OF GENES WITH AT LEAST ONE BINDING SITE #################

# The above heatmap was a little too stringent because it required that there be changing dynamics between the enrichment IPs. What we really want to see instead is all mRNAs represented. So, the filter will be - is a transcript IP'd by at least one of these RBPs. Also, we're going to plot this so that there is no z-scaling. This will help us better visualize when a gene's enrichment is equal across all columns.

## Supplemental Figure 1

# Calculate enrichment sum (Sum of SPN4 == TRUE +  OMA1 == TRUE LIN41 == TRUE)

AllGenes_SPN4

AllGenes_SPN4$LIN41 <- as.logical(AllGenes_SPN4$LIN41)
SPN4_SUM <- cbind(AllGenes_SPN4, sum =  apply(AllGenes_SPN4[, 3:5], 1, sum))
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

#install.packages("UpSetR")
#library(UpSetR)

allMatrix <- AllGenes_SPN4[,c(3,4,5)]
rownames(allMatrix) <- AllGenes_SPN4$name

allNumericMatrix <- 1*allMatrix


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
length(AllGenes_SPN4$gene_ID)
write(AllGenes_SPN4$gene_ID, file = paste("04_output_data/", today, "_ALLlist.txt"))

# Gene IDs - bound by a single RBP independent of any other binding
OMA1_list <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE)
OMA1_list$gene_ID
length(OMA1_list$gene_ID)
write(OMA1_list$gene_ID, file = paste("04_output_data/", today, "_OMA1_list.txt"))

SPN4_list <- AllGenes_SPN4 %>%
  filter(SPN4 == TRUE)
SPN4_list$gene_ID
length(SPN4_list$gene_ID)
write(SPN4_list$gene_ID, file = paste("04_output_data/", today, "_SPN4_list.txt"))

LIN41_list <- AllGenes_SPN4 %>%
  filter(LIN41 == TRUE)
LIN41_list$gene_ID
length(LIN41_list$gene_ID)
write(LIN41_list$gene_ID, file = paste("04_output_data/", today, "_LIN41_list.txt"))


#Gene IDs - Single RBP binding only (no overlap with another RBP)
OMA1_only <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == FALSE)
OMA1_only$gene_ID
length(OMA1_only$gene_ID)
write(OMA1_only$gene_ID, file = paste("04_output_data/", today, "_OMA1_ONLY_list.txt"))

SPN4_only <- AllGenes_SPN4 %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == FALSE)
SPN4_only$gene_ID
length(SPN4_only$gene_ID)
write(SPN4_only$gene_ID, file = paste("04_output_data/", today, "_SPN4_ONLY_list.txt"))

LIN41_only <- AllGenes_SPN4 %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == TRUE)
LIN41_only$gene_ID
length(LIN41_only$gene_ID)
write(LIN41_only$gene_ID, file = paste("04_output_data/", today, "_LIN41_ONLY_list.txt"))

# Overlapping sets between two RBPs
LIN41_SPN4 <- AllGenes_SPN4 %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == TRUE)
LIN41_SPN4$gene_ID
length(LIN41_SPN4$gene_ID)
write(LIN41_SPN4$gene_ID, file = paste("04_output_data/", today, "_LIN41_and_SPN4_list.txt"))


OMA1_SPN4 <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == FALSE)
OMA1_SPN4$gene_ID
length(OMA1_SPN4$gene_ID)
write(OMA1_SPN4$gene_ID, file = paste("04_output_data/", today, "_OMA1_SPN4_list.txt"))

OMA1_LIN41 <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == TRUE)
OMA1_LIN41$gene_ID
length(OMA1_LIN41$gene_ID)
write(OMA1_LIN41$gene_ID, file = paste("04_output_data/", today, "_OMA1_LIN41_list.txt"))

OMA1_SPN4_LIN41 <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == TRUE)
OMA1_SPN4_LIN41$gene_ID
length(OMA1_SPN4_LIN41$gene_ID)
write(OMA1_SPN4_LIN41$gene_ID, file = paste("04_output_data/", today, "_OMA1_SPN4_LIN41_list.txt"))


################# HEATMAP OF SELECT GENES #################

## I want to make a heatmap of select genes that are assayed later in Figure 3

## No figure. Built into a Figure for Figure 3 in section below.

#Heatmap of key assayed genes:

assayedGenes <- c("chs-1", "lin-41", "cpg-2", "oma-2", "car-1", "mex-5", "C04B4.2", "nos-1", "egg-3", "Y37H2A.12", "ZK666.4", "pigv-1", "npr-35", "nasp-2", "R05H11.1", "cgh-1", "mcm-2", "vab-2", "Y19D2B.2")

length(assayedGenes)
length(intersect(AllGenes_SPN4$name, assayedGenes) )
# Not in the dataset: mex-5, Y37HA.12, pigv-1, mcm-2, nasp-2


filtered_SPN4 <- AllGenes_SPN4 %>%
  filter(AllGenes_SPN4$name %in% assayedGenes)

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
#today <- format(Sys.Date(),"%Y%m%d")
#file6 <- paste("03_figures/", today, "_heatmap_smFISH_selecteGenes.pdf", sep = "")

#pdf(file6, height = 9, width = 6)
#v
#dev.off()

########## Fold CHange V. SPN Enrich V. Express. Level #######

## Figure 3

#Heatmap of key assayed genes plotted to include the expression level as well as their expression levels as this seems to make a difference

# Subset
joined_filtered_matrix2 <- joined_filtered_matrix1 %>%
  arrange(desc(foldChange))

joined_filtered_matrix2

w <- pheatmap(joined_filtered_matrix2[,10], 
              scale="none", 
              color = colorRampPalette(c("darkgreen", "white", "orange"), space = "Lab")(100),
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              border_color = "white", 
              show_rownames = TRUE)

w

# Save the plot 
today <- format(Sys.Date(),"%Y%m%d")
file7 <- paste("03_figures/", today, "_heatmap_smFISH_over2lysate_PlusExp.pdf", sep = "")

pdf(file7, height = 9, width = 6)
w
dev.off()

######### SAVE SESSION INFO ###############

today <- format(Sys.Date(),"%Y%m%d")
file8 <- paste("04_output_data/", today, "_sessionInfo.txt", sep = "")


writeLines(capture.output(sessionInfo()), file8)

################################
##         END SCRIPT         ##
################################
