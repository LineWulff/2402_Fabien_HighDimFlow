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

AllTcells <- readRDS("24_07_10_AllTcellsreclustered_SPFvsrewilded_NONdownsampled.rds")
CD4Tcells <- readRDS("24_07_10_CD4Tcellsreclustered_SPFvsrewilded_NONdownsampled.rds")
CD8Tcells <- readRDS("24_07_10_CD8Tcellsreclustered_SPFvsrewilded_NONdownsampled.rds")

nrow(AllTcells) #349249 cells, re-wilded  - 83232, SPF - 266017 
nrow(CD4Tcells) #234376 cells, re-wilded - 54169, SPF - 180207
nrow(CD8Tcells) #114873 cells, re-wilded - 29063, SPF - 85810

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


#### Now do clustering ####
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

saveRDS(AllTcells, file = paste(dato,"AllTcellsreclustered_SPFvsrewilded_DOWNsampled_analysed.rds",sep = "_"))
saveRDS(AllTcells, file = paste(dato,"AllTcellsreclustered_SPFvsrewilded_NONdownsampled_analysed.rds",sep = "_"))


