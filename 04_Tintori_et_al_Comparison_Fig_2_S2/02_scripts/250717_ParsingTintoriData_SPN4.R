
# Cluster Tintori data in a way that I like
# Filter Tintori data for SPN4, OMA1/2, or LIN-41 binding and determine "fate"

###########  LOAD PACKAGES  #####################

library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(vcd)

###########  LOAD PACKAGES  #####################

###########  READ IN THE DATA  #####################

# Set dir
getwd()
setwd("~/Library/CloudStorage/Dropbox/LABWORK/PROJECTS/EOP255_TintoriData")

# Import the metadata from Tintori et al.
getwd()
metadata_tintori <- read.table(file = "01_input/metadata/metadata_parsed_50306.txt", header = FALSE, fill = TRUE)
metadata_tintori

# Annotate & Resort the metadata
colnames(metadata_tintori) <- c("sample", "cellID", "stage", "rep")
head(metadata_tintori)

# Re-order by sample name
metadata_tintori <- metadata_tintori[order(metadata_tintori$sample),]
metadata_tintori

# Import the normalized RNA-seq data from Tintori et al.
tintori <- read.table(file = "01_input/TableS2_RPKMs.txt", header = TRUE)
dim(tintori)
head(tintori)

# Transform the RNA-seq data long-wise using tidyverse
tintori_long <- pivot_longer(tintori, cols = 2:220, names_to ="sample", values_to = "rpkm")
dim(tintori_long)
head(tintori_long)

# Join the Noramlized RNA-seq data to the metadata into tintori_merged
tintori_merged <- left_join(tintori_long, metadata_tintori)

# Remove tossed cellIDs
dim(tintori_merged)
head(tintori_merged)

# How many samples should be tossed?
length(which(metadata_tintori$cellID == "tossed"))
# 54

# How many entries in the merged set?
length(which(tintori_merged$cellID == "tossed"))
# 1,694682

# Remove "tossed" entries
tintori_retained_merged <- tintori_merged[which( !tintori_merged$cellID == "tossed"),]
dim(tintori_retained_merged)
head(tintori_retained_merged)

# Reduce the data structue down into a nested data frame:
tintori_nest_by_transcript <- tintori_retained_merged |>
  group_by(transcript) |>
  nest() 
tintori_nest_by_transcript

# Merge with a list of WormBase IDs:
# Import the list that I downloaded from simplemine. There is a screenshot in the same folder that shows what I selected
wormbaseIDs <- read.table(file = "01_input/download_from_simplemine/simplemine_results_WBGene_publicName.txt", header = TRUE, fill = TRUE)
head(wormbaseIDs)

# Merge the annotation genenames to the tintori nested data frame:          
tintori_nest_by_transcript <- left_join(tintori_nest_by_transcript, wormbaseIDs, join_by( transcript == Public_Name ), keep = TRUE) %>%
  select( c("transcript", "WormBase_Gene_ID", "Public_Name", "data"))
  
dim(tintori_nest_by_transcript)
head(tintori_nest_by_transcript)

#### Make a giant heatmap ####

# Filter for present genes greater than 5 summed across 
tintori_nest_present <- tintori_nest_by_transcript |>
  rowwise() |>
  mutate(total_rpkm = sum(data$rpkm)) |>
  filter(total_rpkm > 5)

# I have 14,776 genes that passed the filter
dim(tintori_nest_present)
head(tintori_nest_present)


# Filter for changing genes that have a variance over 10 
changing <- tintori_nest_present |>
  rowwise() |>
  mutate(var = 100*var(data$rpkm)/total_rpkm) |>
  filter(var > 10)  

unchanging <- tintori_nest_present |>
  rowwise() |>
  mutate(var = 100*var(data$rpkm)/total_rpkm) |>
  filter(var <= 10)  

# I have 12,783 changing genes
dim(changing)
head(changing)

# I have 1,993 unchanging genes
dim(unchanging)
head(unchanging)

# Unnest the changing datastructure and convert into a matrix for heatmap plotting
changing

# Unnest
changing_unnested <- changing |>
  select(transcript, data) |>
  unnest(data) |>
  select(transcript, cellID, rpkm, rep, stage) |>
  mutate( cell_combo = paste(cellID, rep, sep = "_"))
changing_unnested

# Here is a dataframe that has each column as a cell type + replicate and each row as a gene and each value as rpkm
changing_unnested_wide <- changing_unnested |>
  filter(stage != "16-cell") |>
  select(transcript, rpkm, cell_combo) |>
  pivot_wider(names_from = cell_combo, values_from = rpkm)
changing_unnested_wide

# Here is a dataframe that merges each column as a cell type (replicates merged), each row is a gene, and each value is the mean rpkm
changing_wide_df <- changing_unnested |>
  group_by(transcript, cellID) |>
  summarise(mean_rpkm = mean(rpkm)) |>
  pivot_wider(names_from = cellID, values_from = mean_rpkm)
changing_wide_df
dim(changing_wide_df)
str(changing_wide_df) #it's a tibble data frame

# change to a data frame
changing_wide_df <- as.data.frame(changing_wide_df)
changing_wide_df

# create rownames
rownames(changing_wide_df) <- changing_wide_df$transcript

# Remove transcript column
changing_wide_df <- changing_wide_df[ 2:dim(changing_wide_df)[2]]
changing_wide_df <- as.data.frame(changing_wide_df)
str(changing_wide_df)
dim(changing_wide_df)

# Change to a matrix:
changing_wide_mat <- as.matrix(changing_wide_df)

# Order the cell types by stage and anterior-> posterior like so:
cellType_sorted <- read.table(file = "cellTypes_sorted.txt", header = TRUE)
cellType_sorted

colnames(changing_wide_mat)

bestOrder <- c(23, 1, 24, 2, 7, 17, 25, 3, 5, 8, 10, 12, 16, 20, 26, 4, 6, 9, 11, 13, 14, 15, 18, 19, 21, 22, 27)

# Do I have 27 cell types?
length(bestOrder)
length(unique(bestOrder))
dim(cellType_sorted)[1]

# re-arrange changing_wide_mat columns by the best order 
changing_wide_mat <- changing_wide_mat[ , bestOrder]

changing_wide_mat

# changing_wide_mat is my matrix. Next steps will go on to make the heatmaps:




####################### This is a section where I tested ALL the combinations of heatmap clustering options ############
#pheatmap(changing_wide_mat, 
#         cluster_cols = FALSE,
#         scale = "row", 
#         clustering_distance_rows = "minkowski", 
#         clustering_method = "w", 
#         cutree_rows = 5)

#"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#euclidean - looks ok
#maximum - bad
#manhattan - looks ok not as good as euclidean
#canberra - good at sorting out the decreasing v. increasing v. cell-specific gene expression
#binary - oh wow, terrible
#minkowski - not bad at sorting out 

#pheatmap(changing_wide_mat, 
#         cluster_cols = FALSE,
#         scale = "row", 
#         clustering_distance_rows = "canberra", 
#         clustering_method = "median", 
#         cutree_rows = 5)

# Canberra looks best - Next step, vary teh method 

# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#wardD - looks good
#wardD2 - also looks good
#single - nope
# complete - prev - looks good
# average - nope 
# mcquitty - nope
# median - nope
# Centroid - nope

#
#pheatmap(changing_wide_mat, 
#         cluster_cols = FALSE,
#         scale = "row", 
#         clustering_distance_rows = "euclidean", 
#         clustering_method = "single", 
#         cutree_rows = 5)

# Euclidean also looks good
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#wardD - looks good
#wardD2 - looks good
#single- nope

# OF all these, the canberra complete looks the best at separataing out the "decay" group. Euclidean Complete is also very good.

# Here is my favorite heatmap: cannaberra + complete. this one is great for separating out the zygotically transcribed (global) and maternal decay (global) genes:
set.seed(2025)
heatmap_caco5 <- pheatmap(changing_wide_mat, 
         cluster_cols = FALSE,
         scale = "row", 
         show_rownames = FALSE,
         clustering_distance_rows = "canberra", 
         clustering_method = "complete", 
         cutree_rows = 5)

# I split it into five clusters. I want to annotate the cluster number onto the plot so we can ensure that the cluster numbers match. 
# I want what looks like the 2nd and 3rd clusters as the "maternal transcripts undergoing decay groups"

# obtain the cluster information and add it into the matrix
cluster_set = cutree(heatmap_caco5$tree_row, 5)
clustered_changing_wide_mat <- cbind(changing_wide_mat, 
                     cluster = cluster_set )

# add cluster data into its own data frame
ann = data.frame(cluster_set)
rownames(ann) = rownames(clustered_changing_wide_mat)

head(ann)
# Colors of the 
vcolors = plasma(7)[1:5]
vcolors[1]


my_colour = list(cluster_set = c( "1" = vcolors[1], "2" = vcolors[5], "3" = vcolors[3], "4" = vcolors[2], "5" = vcolors[4]))


# Make an annotated heatmap using canaberra + complete:
set.seed(2025)
heatmap_caco5_ann <- pheatmap(changing_wide_mat, 
         cluster_cols = FALSE,
         scale = "row", 
         annotation_colors = my_colour,
         show_rownames = FALSE,
         clustering_distance_rows = "canberra", 
         clustering_method = "complete", 
         cutree_rows = 5,
         annotation_row=ann)


# Save this plot: 
today <- format(Sys.Date(), "%y%m%d")
filename1 <- paste("03_output/", today, "_Heatmap_Tintori_changing_present_transcripts_clust_ann.pdf", sep = "")
pdf(filename1, height = 8, width = 6)

heatmap_caco5_ann
dev.off()

# OK, so the gene clusters I want are clusters 5 and 3: 


###################################################
# SAVE CLUSTER LISTS - PLOT LINEPLOTS BY CLUSTER LIST
###################################################

# filter out lists of genes from each cluster
head(clustered_changing_wide_mat)

# change to a data frame:
clustered_changing_wide_df <- as.data.frame(clustered_changing_wide_mat)
head(clustered_changing_wide_df)

# Set the cluster # as a factor:
clustered_changing_wide_df$cluster <- as.factor(clustered_changing_wide_df$cluster)

# How many in each cluster:
summary(clustered_changing_wide_df$cluster)

# subset clusters
cluster1 <- clustered_changing_wide_df |>
  filter(cluster == 1)
cluster2 <- clustered_changing_wide_df |>
  filter(cluster == 2)
cluster3 <- clustered_changing_wide_df |>
  filter(cluster == 3)
cluster4 <- clustered_changing_wide_df |>
  filter(cluster == 4)
cluster5 <- clustered_changing_wide_df |>
  filter(cluster == 5)

# Create a dataframe that is grouped by the cluster and calculate the mean values of each cell ID across the cluster:
mean_by_cluster <- clustered_changing_wide_df %>% 
  group_by(cluster) %>%
  summarise_all("mean")

# eda 
dim(mean_by_cluster)
meanbycluster <- as.data.frame(mean_by_cluster)
meanbycluster

# convert the meanbycluster into longform data to plot using ggplot2:
longer_means_by_cluster <- pivot_longer(meanbycluster, cols = 2:28, names_to ="cell",
             values_to = "intensity")

# Set cell ID as an ordered factor
longer_means_by_cluster$cell <- factor(longer_means_by_cluster$cell, levels = colnames(changing_wide_mat))

longer_means_by_cluster$cluster <- factor(longer_means_by_cluster$cluster, levels = c(1,4,3,5,2))

# create lineplots of the trend over cell type for each cluster
vcolors = plasma(7)[1:5]
vcolors
lineplot_clusters <- ggplot(longer_means_by_cluster, aes(x=cell, y=intensity, group=cluster, colour=cluster)) + 
  geom_line()+
  geom_point() +
  scale_color_manual(values = vcolors) +
  facet_grid(rows = vars(cluster))

lineplot_clusters
# Save theis plot to annotate the heatmap:
today <- format(Sys.Date(), "%y%m%d")
filename2 <- paste("03_output/", today, "_cluster_lineplots.pdf", sep = "")
pdf(filename2, width = 6, height = 6)

lineplot_clusters

dev.off()

#########################################################
# import the gene lists of SPN-4, OMA-1, and LIN-41 targets
#########################################################

# import the gene lists
OMA1_1 <- read.table(file = "01_input/tatsuyaLists/250226_OMA1list.txt")
SPN4_2 <- read.table(file = "01_input/tatsuyaLists/250226_SPN4list.txt")
LIN41_3 <- read.table(file = "01_input/tatsuyaLists/250226_LIN41list.txt")
OMA1LIN41_4 <- read.table(file = "01_input/tatsuyaLists/250227_OMA1_LIN41_1list.txt")
OMA1SPN4_5 <- read.table(file = "01_input/tatsuyaLists/250227_OMA1_SPN4_1list.txt")
LIN41SPN4_6 <- read.table(file = "01_input/tatsuyaLists/250226_LIN4_and_SPN4_1list.txt")
ALL_7 <- read.table(file = "01_input/tatsuyaLists/250227_OMA1_SPN4_LIN41_1list.txt")

# eda
str(OMA1_1)
str(SPN4_2)
str(LIN41_3)
str(OMA1LIN41_4)
str(OMA1SPN4_5)
str(LIN41SPN4_6)
str(ALL_7)

# annotate the dataframes:
OMA1_1 <- OMA1_1 %>% 
  mutate(IP = "OMA-1")

SPN4_2 <- SPN4_2 %>% 
  mutate(IP = "SPN-4")

LIN41_3 <- LIN41_3 %>% 
  mutate(IP = "LIN-41")

OMA1LIN41_4 <- OMA1LIN41_4 %>% 
  mutate(IP = "OMA-1_LIN-41")

OMA1SPN4_5 <- OMA1SPN4_5 %>% 
  mutate(IP = "OMA-1_SPN-4")

LIN41SPN4_6 <- LIN41SPN4_6 %>% 
  mutate(IP = "LIN-41_SPN-4")

ALL_7 <- ALL_7 %>% 
  mutate(IP = "OMA-1_SPN-4_LIN-41")
ALL_7

# merge the gene lists
IP_lookup <- rbind(OMA1_1, SPN4_2, LIN41_3, OMA1LIN41_4, OMA1SPN4_5, LIN41SPN4_6, ALL_7)
colnames(IP_lookup) <- c("WormBase_Gene_ID", "IP") 
dim(IP_lookup)
table(IP_lookup$IP)
head(IP_lookup)

# user wormbaseID dataframe to add transcript names to gene lists
head(wormbaseIDs)
IP_lookup_annotated <- left_join(IP_lookup, wormbaseIDs)
head(IP_lookup_annotated)
# merge gene lists with changing_wide_mat

# create an annotated heatmap
changing_wide_df <- rownames_to_column(as.data.frame(changing_wide_mat))
head(changing_wide_df)

changing_wide_IP <- left_join(changing_wide_df, IP_lookup_annotated, by = c("rowname"="Sequence_Name"))

changing_wide_df_plusIP <- column_to_rownames(changing_wide_IP[,1:30]) 
# Add the column changing_wide_IP$IP to the annotation table and re-map the heatmap 

# add cluster data into its own data frame
head(ann)
head(changing_wide_df_plusIP)
table(rownames(ann) == rownames(changing_wide_df_plusIP))

changing_wide_df_plusIP
ann2 <- cbind(ann, changing_wide_df_plusIP$IP)
colnames(ann2) <- c("cluster_set", "IP")

ann3 <- ann2

# Make an annotated heatmap using canaberra + complete:
set.seed(2025)
heatmap_caco5_ip <- pheatmap(changing_wide_mat, 
                              cluster_cols = FALSE,
                              scale = "row", 
                              show_rownames = FALSE,
                              annotation_colors = my_colour,
                              clustering_distance_rows = "canberra", 
                              clustering_method = "complete", 
                              cutree_rows = 5,
                              annotation_row=ann2)

dim(changing_wide_mat)
# Save the heatmap:
today <- format(Sys.Date(), "%y%m%d")
filename4 <- paste("03_output/", today, "_heatmap_clusters_and_IP.pdf", sep = "")
pdf(filename4, width = 6, height = 8)
heatmap_caco5_ip
    
dev.off()

# This didn't work. Had to use Export -> Export as Jpeg
filename5 <- paste("03_output/", today, "_heatmap_clusters_and_IP.jpg", sep = "")
jpeg(filename5, width = 6, height = 8)

heatmap_caco5_ip

dev.off()

help(pdf)

# Make lineplots:

ann3$cluster_set <- factor(ann3$cluster_set, levels = c(1, 4, 3, 5, 2))

ann3$cluster_set
ann3[which(is.na(ann3$IP)), 2] <- "none"
ann3$IP <- factor(ann3$IP, levels = c("OMA-1", "SPN-4", "LIN-41", "OMA-1_LIN-41", "OMA-1_SPN-4", "LIN-41_SPN-4", "OMA-1_SPN-4_LIN-41", "none"))
ann3
  
table(ann3)
ann3
str(ann3)

ann3tabs <- xtabs(~ cluster_set + IP, data = ann3)
mosaicplot1 <- mosaic(ann3tabs, gp = shading_max, 
                      split_horizontal = TRUE)

mosaic(ann3tabs, gp = shading_max, 
       split_horizontal = TRUE)

mosaicplot1
today <- format(Sys.Date(), "%y%m%d")
filename3 <- paste("03_output/", today, "_mosaicplot.pdf", sep = "")
filename3

#--> Note: this didn't look right. Had to save with Export -> Save As PDF
pdf(filename3, width = 8, height = 8)

mosaic(ann2tabs, gp = shading_max, 
       split_horizontal = TRUE)

dev.off()
dev.off()
############ CREATE LINEPLOTS group_by IP #############

# start with clustered_changing_wide_df and annotate with ann2 using left_join
head(clustered_changing_wide_df)
head(ann2)

clustered_changing_wide_ann2 <- left_join(rownames_to_column(clustered_changing_wide_df), rownames_to_column(ann3))
clustered_changing_wide_ann2 <- column_to_rownames(clustered_changing_wide_ann2) 
clustered_changing_wide_ann2 <- clustered_changing_wide_ann2 %>%
  select(!(cluster:cluster_set))


head(clustered_changing_wide_ann2)
# Create a dataframe that is grouped by the cluster and calculate the mean values of each cell ID across the cluster:

mean_by_IP <- clustered_changing_wide_ann2 %>% 
  group_by(IP) %>%
  summarise_all("mean")

# eda 
dim(mean_by_IP)
mean_by_IP
mean_by_IP <- as.data.frame(mean_by_IP)
mean_by_IP

# convert the meanbycluster into longform data to plot using ggplot2:
longer_means_by_IP <- pivot_longer(mean_by_IP, cols = 2:28, names_to = "cell",
                                        values_to = "intensity")

# Set cell ID as an ordered factor
longer_means_by_IP
longer_means_by_IP$cell <- factor(longer_means_by_IP$cell, levels = colnames(changing_wide_mat))

longer_means_by_IP

# create lineplots of the trend over cell type for each cluster
lineplot_IP <- ggplot(longer_means_by_IP, aes(x=cell, y=intensity, group=IP, colour=IP)) + 
  geom_line()+
  geom_point() +
  facet_grid(rows = vars(IP))

lineplot_IP
# Save theis plot to annotate the heatmap:
today <- format(Sys.Date(), "%y%m%d")
filename2 <- paste("03_output/", today, "_IP_lineplots.pdf", sep = "")
pdf(filename2, width = 6, height = 6)

lineplot_IP

dev.off()



########## DELETE BELOW HERE ###############


changing_bystage_df <- changing_unnested |>
  group_by(transcript, stage) |>
  summarise(mean_rpkm = mean(rpkm)) |>
  pivot_wider(names_from = stage, values_from = mean_rpkm)


changing_wide_df
dim(changing_wide_df)
changing_wide_df <- as.data.frame(changing_wide_df)
rownames(changing_wide_df) <- changing_wide_df$transcript
changing_wide_df <- changing_wide_df[ 2:dim(changing_wide_df)[2]]
changing_wide_df <- as.data.frame(changing_wide_df)

changing_wide_df
dim(changing_wide_df)

changing_wide_mat <- as.matrix(changing_wide_df)




changing_bystage_df
changing_wide_df
dim(changing_wide_df)
changing_wide_df <- as.data.frame(changing_wide_df)
rownames(changing_wide_df) <- changing_wide_df$transcript
changing_wide_df <- changing_wide_df[ 2:dim(changing_wide_df)[2]]
changing_wide_df <- as.data.frame(changing_wide_df)



##########



changing_unnested |>
  select(transcript, rpkm, cell_combo) |>
  pivot_wider(names_from = cell_combo, values_from = rpkm)


print(changing_unnested[1:1000,], n = 1000)
aagr2_only <- changing |>
  select(transcript, stage, data) |>
  unnest_wider(films, names_sep = "_") |>
  filter(transcript == "aagr-2")


{changing_unnested} |>
  dplyr::summarise(n = dplyr::n(), .by = c(transcript, cell_combo)) |>
  dplyr::filter(n > 1L) 

head(changing)

aagr2 <- changing_unnested |>
  filter (transcript == "aagr-2")
aagr2[1:165,]
print(aagr2, n = 165)

changing_unnested

pivot_wider(changing_unnested, names_from = changing_unnested$cellID,
            values_from = changing_unnested$rpkm) 


df <- tibble(
  x = 1:2,
  y = list(
    tibble(a = 1, b = 2),
    tibble(a = 3:4, b = 5:6)
  ),
  z = list(
    tibble(c = 1, d = 2),
    tibble(c = 3:4, d = 5:6)
  )
)
df[[2]]

df %>% unnest(c(y, z))

# Pheatmap




#### FILTER FOR THE DIFFERENT CATEGORIES:

# Import lists of sets:
SPN4 <- read.table(file = "", header = TRUE, fill = TRUE)
head(wormbaseIDs)
