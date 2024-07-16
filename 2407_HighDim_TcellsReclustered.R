#' R script for producing outputs in from Fabien's rewilded vs. SPF mice in T cells reclustered
#' Author: Line Wulff
#' Date (created): 23-08-18

#' #### ---- Initiate libraries ---- ####
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)
library(stringr)
library(uwot)
library(ggrastr)
library(groupdata2)
library(viridis)
library(scales)
library(tidyr)

#### ---- variables used throughout script ---- ####
projdir <- getwd()
RAID_dir <- "/Volumes/Promise RAID/Line/projects/2402_Fabien_HighDimFlow"
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# flow markers
marker_cols <- panel$marker[!panel$marker %in% c("L-D","CD45","Tetramer")]

#### Reclustering of T cells (non-downsampled) ####
# At resolution 8
# cl 5 are CD4 T cells 
# cl 6 are CD8 T cells

AllTcells <- readRDS("24_07_15_AllTcellsreclustered_SPFvsrewilded_NONdownsampled.rds")
#AllTcells <- readRDS("24_07_15_AllTcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds")
#AllTcells <- readRDS("24_07_16_AllTcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds")
CD4Tcells <- readRDS("24_07_15_CD4Tcellsreclustered_SPFvsrewilded_NONdownsampled.rds")
#AllTcells <- readRDS("24_07_16_CD4Tcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds")
#AllTcells <- readRDS("24_07_16_CD4Tcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds")
CD8Tcells <- readRDS("24_07_15_CD8Tcellsreclustered_SPFvsrewilded_NONdownsampled.rds")
#AllTcells <- readRDS("24_07_16_CD8Tcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds")
#AllTcells <- readRDS("24_07_16_CD8Tcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds")

nrow(AllTcells) #226.520 cells, re-wilded: 58966, SPF: 167554
nrow(CD4Tcells) #111.647 cells, re-wilded: 29903, SPF: 81744 
nrow(CD8Tcells) #114.873 cells, re-wilded: 29063, SPF: 85810 

#### --- All T cells ---- ####
#' #### ---- Running UMAP ---- ####
# possibility of downsampling
# downsample to equal amount of cells in each group/sample/what makes sense
# check group with fewest observations
table(AllTcells$sample)

# downsampled for overall analysis of CD45+
AllTcells <- downsample(AllTcells, cat_col = "sample")
# or non downsampled for T cell analysis
# just don't run above


# prepare data for umapr (matrix format required)
data_umap <- AllTcells[, marker_cols]
data_umap <- as.matrix(data_umap)
dups <- duplicated(data_umap)
data_umap <- data_umap[!dups, ]

# run umap, takes a long time with many cells, downsample first to test
umap_emb <- umap(data_umap, scale = "maxabs") 

# prepare umap embedding output data for plot
AllTcells_umap <- as.data.frame(umap_emb)
colnames(AllTcells_umap) <- c("UMAP_1", "UMAP_2")
head(AllTcells_umap);dim(AllTcells_umap)
AllTcells[,c("UMAP_1","UMAP_2")] <- AllTcells_umap


#### Plot each individual group as seperate contour on top of total cells
ggplot(AllTcells, aes(x = UMAP_1, y = UMAP_2, colour = mice))+ 
  geom_point_rast(size=0.1)+
  theme_classic()

ggplot(AllTcells, aes(x = UMAP_1, y = UMAP_2, colour=mice))+ 
  #geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme_classic()
ggsave(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_UMAP_Contour_rewildvsSPF_plot.pdf",sep=""), height = 4, width = 4)

# rewilded only
ggplot(AllTcells, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=AllTcells[AllTcells$mice=="re-wilded",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#F8766D")+
  theme_classic()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_UMAP_rewildedContour_plot.pdf",sep=""), height = 4, width = 4)

#SPF only
ggplot(AllTcells, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=AllTcells[AllTcells$mice=="SPF",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#00BFC4")+
  theme_classic()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_UMAP_SPFContour_plot.pdf",sep=""), height = 4, width = 4)


#### Colour by individual marker and save pdf versions
for (marker in marker_cols){
  plot_mark <- ggplot(AllTcells, aes(x = UMAP_1, y = UMAP_2, colour=data_umap[,marker]))+ 
    geom_point_rast(size=0.1)+
    scale_colour_viridis_c(option = "plasma")+
    theme_classic()+
    labs(colour=marker)+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_UMAP_",marker,"_plot.pdf", sep=""), height = 4, width = 4)
  print(plot_mark)
  dev.off()
}


#' #### Now do clustering ####
# Individual samples based on timestamps - "sampleID"
# FlowSOM variables
n_max <- 10 # set slightly higher than actual number as easier to combine some clusters after than having to few
seed <- 13

for (i in seq(3,n_max)){
  n_meta <- i
  # Compute the FlowSOM object
  fsom <- FlowSOM(input = as.matrix(data_umap),
                  scale = FALSE,
                  colsToUse = marker_cols,
                  seed = seed,
                  nClus = n_meta)
  # add to FC file
  AllTcells[paste("cluster", i, sep = "_")] <- GetMetaclusters(fsom)
  # plot and check matches
  clus_plot <- ggplot(AllTcells, aes(x=UMAP_1,y=UMAP_2,colour = AllTcells[,paste("cluster", i, sep = "_")]))+
    geom_point()+
    labs(colour=paste("cluster", i, sep = "_"))+
    theme_classic()+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  pdf(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_UMAP_","cluster", i,"_plot.pdf", sep=""), height = 4, width = 4)
  print(clus_plot)
  dev.off()
}

#' #### Plots of distributions ####
# run function from 2407_DistributionDF.R
cl10_dist <- perc_df_calc(AllTcells, "sample", "cluster_10")
cl10_dist$mice <- unlist(str_split(cl10_dist$sample,"_"))[seq(1,nrow(cl10_dist)*2,2)]
cl10_dist$cluster <- factor(cl10_dist$cluster, levels = 1:10)
# summarized table for errorbars
cl10_dist_av <- data_summary(cl10_dist, varname="percentage", 
             groupnames=c("mice", "cluster"))
cl10_dist_av$cluster <- factor(cl10_dist_av$cluster, levels = 1:10)
cl10_dist_av$mice <- factor(cl10_dist_av$mice, levels = c("SPF","re-wilded"))

## Plots
## stacked
ggplot(cl10_dist, aes(x=sample, y=percentage, fill=cluster))+ 
  geom_bar(position = "stack", stat="identity", colour="black")+
  theme_classic()+
  ylab("% of total T cells per sample")+
  theme(axis.text.x = element_text(angle=90))
ggsave(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_cl10distribution_stacked.pdf",sep=""), height = 4, width = 4)

## side by side
ggplot(cl10_dist_av, aes(x=mice, y=percentage, fill=cluster))+
  geom_bar(position = "dodge", stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=.2,
                position=position_dodge(.9))+
  facet_grid(.~cluster)+
  theme_classic()+
  xlab("")+ylab("% of total T cells per sample")+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_cl10distribution_sidebyside.pdf",sep=""), height = 4, width = 5)

#' #### Plots of markers - density per cluster ####
for (marker in marker_cols){
  plot_mark <- ggplot(data = AllTcells, aes(x = AllTcells[,marker], group = cluster_10, fill = cluster_10)) +
    geom_density(adjust = 1.5)+
    facet_wrap(~cluster_10)+
    theme_classic()+
    xlab(marker)+
    theme(axis.text.y = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"AllTcellsreclustered_cluster10_density_",marker,"_plot.pdf", sep=""), height = 6, width = 6)
  print(plot_mark)
  dev.off()
}

#' #### save objects as analysed for later analysis ####
saveRDS(AllTcells, file = paste(dato,"AllTcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds",sep = "_"))
saveRDS(AllTcells, file = paste(dato,"AllTcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds",sep = "_"))




#### --- CD4 T cells ---- ####
#' #### ---- Running UMAP ---- ####
# possibility of downsampling
# downsample to equal amount of cells in each group/sample/what makes sense
# check group with fewest observations
table(CD4Tcells$sample)

# downsampled for overall analysis of CD45+
CD4Tcells <- downsample(CD4Tcells, cat_col = "sample")
# or non downsampled for T cell analysis
# just don't run above


# prepare data for umapr (matrix format required)
data_umap <- CD4Tcells[, marker_cols]
data_umap <- as.matrix(data_umap)
dups <- duplicated(data_umap)
data_umap <- data_umap[!dups, ]

# run umap, takes a long time with many cells, downsample first to test
umap_emb <- umap(data_umap, scale = "maxabs") 

# prepare umap embedding output data for plot
CD4Tcells_umap <- as.data.frame(umap_emb)
colnames(CD4Tcells_umap) <- c("UMAP_1", "UMAP_2")
head(CD4Tcells_umap);dim(CD4Tcells_umap)
CD4Tcells[,c("UMAP_1","UMAP_2")] <- CD4Tcells_umap


#### Plot each individual group as seperate contour on top of total cells
ggplot(CD4Tcells, aes(x = UMAP_1, y = UMAP_2, colour = mice))+ 
  geom_point_rast(size=0.1)+
  theme_classic()

ggplot(CD4Tcells, aes(x = UMAP_1, y = UMAP_2, colour=mice))+ 
  #geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme_classic()
ggsave(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_UMAP_Contour_rewildvsSPF_plot.pdf",sep=""), height = 4, width = 4)

# rewilded only
ggplot(CD4Tcells, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=CD4Tcells[CD4Tcells$mice=="re-wilded",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#F8766D")+
  theme_classic()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_UMAP_rewildedContour_plot.pdf",sep=""), height = 4, width = 4)

#SPF only
ggplot(CD4Tcells, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=CD4Tcells[CD4Tcells$mice=="SPF",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#00BFC4")+
  theme_classic()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_UMAP_SPFContour_plot.pdf",sep=""), height = 4, width = 4)


#### Colour by individual marker and save pdf versions
for (marker in marker_cols){
  plot_mark <- ggplot(CD4Tcells, aes(x = UMAP_1, y = UMAP_2, colour=data_umap[,marker]))+ 
    geom_point_rast(size=0.1)+
    scale_colour_viridis_c(option = "plasma")+
    theme_classic()+
    labs(colour=marker)+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_UMAP_",marker,"_plot.pdf", sep=""), height = 4, width = 4)
  print(plot_mark)
  dev.off()
}


#' #### Now do clustering ####
# Individual samples based on timestamps - "sampleID"
# FlowSOM variables
n_max <- 10 # set slightly higher than actual number as easier to combine some clusters after than having to few
seed <- 13

for (i in seq(3,n_max)){
  n_meta <- i
  # Compute the FlowSOM object
  fsom <- FlowSOM(input = as.matrix(data_umap),
                  scale = FALSE,
                  colsToUse = marker_cols,
                  seed = seed,
                  nClus = n_meta)
  # add to FC file
  CD4Tcells[paste("cluster", i, sep = "_")] <- GetMetaclusters(fsom)
  # plot and check matches
  clus_plot <- ggplot(CD4Tcells, aes(x=UMAP_1,y=UMAP_2,colour = CD4Tcells[,paste("cluster", i, sep = "_")]))+
    geom_point()+
    labs(colour=paste("cluster", i, sep = "_"))+
    theme_classic()+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  pdf(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_UMAP_","cluster", i,"_plot.pdf", sep=""), height = 4, width = 4)
  print(clus_plot)
  dev.off()
}

#' #### Plots of distributions ####
# run function from 2407_DistributionDF.R
cl7_dist <- perc_df_calc(CD4Tcells, "sample", "cluster_7")
cl7_dist$mice <- unlist(str_split(cl7_dist$sample,"_"))[seq(1,nrow(cl7_dist)*2,2)]
cl7_dist$cluster <- factor(cl7_dist$cluster, levels = 1:10)
# summarized table for errorbars
cl7_dist_av <- data_summary(cl7_dist, varname="percentage", 
                             groupnames=c("mice", "cluster"))
cl7_dist_av$cluster <- factor(cl7_dist_av$cluster, levels = 1:10)
cl7_dist_av$mice <- factor(cl7_dist_av$mice, levels = c("SPF","re-wilded"))

## Plots
## stacked
ggplot(cl7_dist, aes(x=sample, y=percentage, fill=cluster))+ 
  geom_bar(position = "stack", stat="identity", colour="black")+
  theme_classic()+
  ylab("% of total T cells per sample")+
  theme(axis.text.x = element_text(angle=90))
ggsave(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_cl7distribution_stacked.pdf",sep=""), height = 4, width = 4)

## side by side
ggplot(cl7_dist_av, aes(x=mice, y=percentage, fill=cluster))+
  geom_bar(position = "dodge", stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=.2,
                position=position_dodge(.9))+
  facet_grid(.~cluster)+
  theme_classic()+
  xlab("")+ylab("% of total T cells per sample")+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_cl7distribution_sidebyside.pdf",sep=""), height = 4, width = 5)

#' #### Plots of markers - density per cluster ####
for (marker in marker_cols){
  plot_mark <- ggplot(data = CD4Tcells, aes(x = CD4Tcells[,marker], group = cluster_7, fill = cluster_10)) +
    geom_density(adjust = 1.5)+
    facet_wrap(~cluster_10)+
    theme_classic()+
    xlab(marker)+
    theme(axis.text.y = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"CD4Tcellsreclustered_cluster7_density_",marker,"_plot.pdf", sep=""), height = 6, width = 6)
  print(plot_mark)
  dev.off()
}

#' #### save objects as analysed for later analysis ####
saveRDS(CD4Tcells, file = paste(dato,"CD4Tcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds",sep = "_"))
saveRDS(CD4Tcells, file = paste(dato,"CD4Tcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds",sep = "_"))

#### --- CD8 T cells --- ####
#' #### ---- Running UMAP ---- ####
# possibility of downsampling
# downsample to equal amount of cells in each group/sample/what makes sense
# check group with fewest observations
table(CD8Tcells$sample)

# downsampled for overall analysis of CD45+
CD8Tcells <- downsample(CD8Tcells, cat_col = "sample")
CD8Tcells$mice <- factor(CD8Tcells$mice, levels = c("SPF","re-wilded"))
# or non downsampled for T cell analysis
# just don't run above


# prepare data for umapr (matrix format required)
data_umap <- CD8Tcells[, marker_cols]
data_umap <- as.matrix(data_umap)
dups <- duplicated(data_umap)
data_umap <- data_umap[!dups, ]

# run umap, takes a long time with many cells, downsample first to test
umap_emb <- umap(data_umap, scale = "maxabs") 

# prepare umap embedding output data for plot
CD8Tcells_umap <- as.data.frame(umap_emb)
colnames(CD8Tcells_umap) <- c("UMAP_1", "UMAP_2")
head(CD8Tcells_umap);dim(CD8Tcells_umap)
CD8Tcells[,c("UMAP_1","UMAP_2")] <- CD8Tcells_umap


#### Plot each individual group as seperate contour on top of total cells
ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2, colour = mice))+ 
  geom_point_rast(size=0.1)+
  theme_classic()

ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2, colour=mice))+ 
  #geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d()+
  xlim(c(-8,8))+ylim(c(-8,8))+
  theme_classic()
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_Contour_rewildvsSPF_plot.pdf",sep=""), height = 4, width = 4)

ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2))+
  geom_point(colour="darkgrey")+
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), alpha = 0.5, colour = "black",
                  bins = 5)+
  facet_grid(.~mice)+
  xlim(c(-8,8))+ylim(c(-8,8))+
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_classic()
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_Contour_rewildvsSPF_plot2blue.pdf",sep=""), height = 4, width = 8)


ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2, colour=mice))+
  geom_point(colour="lightgrey",alpha=0.3)+
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = mice), bins=5)+
  facet_grid(.~mice)+
  xlim(c(-8,8))+ylim(c(-8,8))+
  theme_classic()
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_Contour_rewildvsSPF_plot3sepcols.pdf",sep=""), height = 4, width = 8)


# rewilded only
ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=CD8Tcells[CD8Tcells$mice=="re-wilded",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#00BFC4")+
  theme_classic()+
  xlim(c(-8,8))+ylim(c(-8,8))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_rewildedContour_plot.pdf",sep=""), height = 4, width = 4)

#SPF only
ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=CD8Tcells[CD8Tcells$mice=="SPF",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#F8766D")+
  theme_classic()+
  xlim(c(-8,8))+ylim(c(-8,8))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_SPFContour_plot.pdf",sep=""), height = 4, width = 4)


#### Colour by individual marker and save pdf versions
for (marker in marker_cols){
  plot_mark <- ggplot(CD8Tcells, aes(x = UMAP_1, y = UMAP_2, colour=data_umap[,marker]))+ 
    geom_point_rast(size=0.1)+
    scale_colour_viridis_c(option = "plasma")+
    theme_classic()+
    labs(colour=marker)+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_",marker,"_plot.pdf", sep=""), height = 4, width = 4)
  print(plot_mark)
  dev.off()
}


#' #### Now do clustering ####
# Individual samples based on timestamps - "sampleID"
# FlowSOM variables
n_max <- 10 # set slightly higher than actual number as easier to combine some clusters after than having to few
seed <- 13

for (i in seq(3,n_max)){
  n_meta <- i
  # Compute the FlowSOM object
  fsom <- FlowSOM(input = as.matrix(data_umap),
                  scale = FALSE,
                  colsToUse = marker_cols,
                  seed = seed,
                  nClus = n_meta)
  # add to FC file
  CD8Tcells[paste("cluster", i, sep = "_")] <- GetMetaclusters(fsom)
  # plot and check matches
  clus_plot <- ggplot(CD8Tcells, aes(x=UMAP_1,y=UMAP_2,colour = CD8Tcells[,paste("cluster", i, sep = "_")]))+
    geom_point()+
    labs(colour=paste("cluster", i, sep = "_"))+
    theme_classic()+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  pdf(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_UMAP_","cluster", i,"_plot.pdf", sep=""), height = 4, width = 4)
  print(clus_plot)
  dev.off()
}

#' #### Plots of distributions ####
# run function from 2407_DistributionDF.R
cl10_dist <- perc_df_calc(CD8Tcells, "sample", "cluster_10")
cl10_dist$mice <- unlist(str_split(cl10_dist$sample,"_"))[seq(1,nrow(cl10_dist)*2,2)]
cl10_dist$cluster <- factor(cl10_dist$cluster, levels = 1:10)
# summarized table for errorbars
cl10_dist_av <- data_summary(cl10_dist, varname="percentage", 
                             groupnames=c("mice", "cluster"))
cl10_dist_av$cluster <- factor(cl10_dist_av$cluster, levels = 1:10)
cl10_dist_av$mice <- factor(cl10_dist_av$mice, levels = c("SPF","re-wilded"))

## Plots
## stacked
ggplot(cl10_dist, aes(x=sample, y=percentage, fill=cluster))+ 
  geom_bar(position = "stack", stat="identity", colour="black")+
  theme_classic()+
  ylab("% of total T cells per sample")+
  theme(axis.text.x = element_text(angle=90))
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_cl10distribution_stacked.pdf",sep=""), height = 4, width = 4)

## side by side
ggplot(cl10_dist_av, aes(x=mice, y=percentage, fill=cluster))+
  geom_bar(position = "dodge", stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=.2,
                position=position_dodge(.9))+
  facet_grid(.~cluster)+
  theme_classic()+
  xlab("")+ylab("% of total T cells per sample")+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_cl10distribution_sidebyside.pdf",sep=""), height = 4, width = 5)

#' #### Plots of markers - density per cluster ####
for (marker in marker_cols){
  plot_mark <- ggplot(data = CD8Tcells, aes(x = CD8Tcells[,marker], group = cluster_10, fill = cluster_10)) +
    geom_density(adjust = 1.5)+
    facet_wrap(~cluster_10)+
    theme_classic()+
    xlab(marker)+
    theme(axis.text.y = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"CD8Tcellsreclustered_cluster10_density_",marker,"_plot.pdf", sep=""), height = 6, width = 6)
  print(plot_mark)
  dev.off()
}

#' #### save objects as analysed for later analysis ####
saveRDS(CD8Tcells, file = paste(dato,"CD8Tcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds",sep = "_"))
saveRDS(CD8Tcells, file = paste(dato,"CD8Tcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds",sep = "_"))
