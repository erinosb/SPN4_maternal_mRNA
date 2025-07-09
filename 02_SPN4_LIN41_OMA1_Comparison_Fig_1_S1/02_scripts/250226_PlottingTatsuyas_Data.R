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

# HISTOGRAM

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



# Draw the HISTOGRAM and save as a .pdf file
ggsave(histo1, filename = "03_figures/240724_histogram.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")


############### X Y Scatter Facets #################



############### CORR Matrix ############################


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

pdf("03_figures/240724_corr_matrix_plots.pdf", height = 6, width = 7)
par(mfrow=c(1,1))
p

dev.off()

help(dist)

################# HEATMAP #################

## REVIEW:
## Create a matrix
# head(AllGenes_SPN4)
# str(AllGenes_SPN4)
#dim(AllGenes_SPN4)

## REVIEW:


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


pdf("03_figures/240724_heatmap_all.pdf", height = 9, width = 6)
q
dev.off()

# What about clustering?

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


write.table(SPN4_wide_2cutoff_2FPKM, file = "03_figures/240731_clusters.txt", quote = FALSE, sep = "\t")
help(write.table)


################# Make a HEATMAP of all GENES - no Z Score ##############

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

# Make the cutoffs
## Select for enrichments over -2 and more than Log2(2)
## Select for expression levels greater than e^2(4)
## Select for stdev > 0.6
SPN4_expressed <- SPN4_wide_stringentCutoff1 %>%
  filter(Log2_SPN.4_LYSATE_FPKM > 4)

SPN4_expressed

# Full heatmap of high-ish expression genes:
dim(SPN4_expressed)
head(SPN4_expressed)
summary(SPN4_expressed$spn4_enrichment)

SPN4_expressed

laxMatrix <- as.matrix(SPN4_expressed[,c(6,7,8)])
rownames(laxMatrix) <- SPN4_expressed$name

dim(laxMatrix)

t <- pheatmap(laxMatrix, 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              cutree_rows = 5,
              show_rownames = FALSE)

t

################# HEATMAP OF GENES WITH AT LEAST ONE BINDING SITE #################

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

today <- format(Sys.Date(), "%y%m%d")
filename1 <- paste("03_figures/", today, "_Heatmap_BoundTargets_noScale.pdf", sep = "")
pdf(filename1, height = 9, width = 6)
u
dev.off()



################# UPSET PLOT ##################

install.packages("UpSetR")
library(UpSetR)

allMatrix <- AllGenes_SPN4[,c(3,4,5)]
rownames(allMatrix) <- AllGenes_SPN4$name

allNumericMatrix <- 1*allMatrix


upset(allNumericMatrix, text.scale=2)

help(upset)
help(upset)
today <- format(Sys.Date(), "%y%m%d")
filename <- paste("03_figures/", today, "_upsetPlot.pdf", sep = "")

pdf(filename, height = 4, width = 8)
upset(allNumericMatrix, text.scale=2)
dev.off()


### SAVE a bunch of lists for GO Ontology Analysis #############


length(AllGenes_SPN4$gene_ID)
write(AllGenes_SPN4$gene_ID, file = "03_figures/250227_ALLlist.txt")

#Gene IDs - OMA1 only
OMA1_only <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == FALSE)
OMA1_only$gene_ID
length(OMA1_only$gene_ID)
write(OMA1_only$gene_ID, file = "03_figures/250226_OMA1list.txt")

SPN4_only <- AllGenes_SPN4 %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == FALSE)
SPN4_only$gene_ID
length(SPN4_only$gene_ID)
getwd()
write(SPN4_only$gene_ID, file = "03_figures/250226_SPN4list.txt")
help(write)

LIN41_only <- AllGenes_SPN4 %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == TRUE)
LIN41_only$gene_ID
length(LIN41_only$gene_ID)
write(LIN41_only$gene_ID, file = "03_figures/250226_LIN41list.txt")
help(write)

LIN41_SPN4 <- AllGenes_SPN4 %>%
  filter(OMA1 == FALSE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == TRUE)
LIN41_SPN4$gene_ID
length(LIN41_SPN4$gene_ID)
getwd()
write(LIN41_SPN4$gene_ID, file = "03_figures/250226_LIN4_and_SPN4_1list.txt")


OMA1_SPN4 <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == FALSE)
OMA1_SPN4$gene_ID
length(OMA1_SPN4$gene_ID)
getwd()
write(OMA1_SPN4$gene_ID, file = "03_figures/250227_OMA1_SPN4_1list.txt")

OMA1_LIN41 <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == FALSE) %>%
  filter(LIN41 == TRUE)
OMA1_LIN41$gene_ID
length(OMA1_LIN41$gene_ID)
getwd()
write(OMA1_LIN41$gene_ID, file = "03_figures/250227_OMA1_LIN41_1list.txt")

OMA1_SPN4_LIN41 <- AllGenes_SPN4 %>%
  filter(OMA1 == TRUE) %>%
  filter(SPN4 == TRUE) %>%
  filter(LIN41 == TRUE)
OMA1_SPN4_LIN41$gene_ID
length(OMA1_SPN4_LIN41$gene_ID)
getwd()
write(OMA1_SPN4_LIN41$gene_ID, file = "03_figures/250227_OMA1_SPN4_LIN41_1list.txt")


################# HEATMAP OF SELECT GENES #################

#Heatmap of key assayed genes:

assayedGenes <- c("chs-1", "lin-41", "cpg-2", "oma-2", "car-1", "mex-5", "C04B4.2", "nos-1", "egg-3", "Y37H2A.12", "ZK666.4", "pigv-1", "npr-35", "nasp-2", "R05H11.1", "cgh-1", "mcm-2", "vab-2", "Y19D2B.2")

length(assayedGenes)
length(intersect(AllGenes_SPN4$name, assayedGenes) )
# Not in the dataset: mex-5, Y37HA.12, pigv-1, mcm-2, nasp-2


filtered_SPN4 <- AllGenes_SPN4 %>%
  filter(AllGenes_SPN4$name %in% assayedGenes)

filtered_SPN4


# Input smFISH Fold Change data:
setwd("~/Library/CloudStorage/Dropbox/LABWORK/PROJECTS/EOP248_spn4_dataAnalysis/Parsing_Tatsuyas_data")

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

s <- pheatmap(joined_filtered_matrix2[,c(3,2,1)], 
              scale="none", 
              color = colorRampPalette(c("blue4", "white", "maroon2"), space = "Lab")(100),
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              border_color = "white", 
              show_rownames = TRUE)

s


pdf("03_figures/240724_heatmap_smFISH_selecteGenes.pdf", height = 9, width = 6)
s
dev.off()




########## Fold CHange V. SPN Enrich V. Express. Level #######





########### FILTERING BY EXPRESSION LEVEL ##################


##########################################

#Heatmap of key assayed genes:

# Subset
joined_filtered_matrix2 <- joined_filtered_matrix1 %>%
  arrange(desc(foldChange))

joined_filtered_matrix2

u <- pheatmap(joined_filtered_matrix2[,10], 
              scale="none", 
              color = colorRampPalette(c("darkgreen", "white", "orange"), space = "Lab")(100),
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              border_color = "white", 
              show_rownames = TRUE)

u

pdf("03_figures/heatmap_over2lysatePlusExp.pdf", height = 9, width = 6)
u
dev.off()

colnames(joined_filtered_matrix2)

sp2<-ggplot(joined_filtered_matrix2, aes(x=spn4_enrichment, y=foldChange, color=Log2_SPN.4_LYSATE_FPKM)) + 
  geom_point()
sp2
# Change the low and high colors
# Sequential color scheme

sp2+scale_color_gradientn(colours = rainbow(4))

