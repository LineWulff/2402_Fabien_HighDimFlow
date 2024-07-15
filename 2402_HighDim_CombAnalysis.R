#' R script for producing outputs in from Fabien's rewilded vs. SPF mice
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

#' #### ---- Read in data, inspect and transform labels if necessary ---- ####
SPF_FC <- as.data.frame(flowCore::exprs(flowCore::read.FCS(
  filename = "/Volumes/Promise RAID/Line/projects/2402_Fabien_HighDimFlow/inputfcsfiles/SPF.fcs", 
  transformation = FALSE, truncate_max_range = FALSE)))

rewild_FC <- as.data.frame(flowCore::exprs(flowCore::read.FCS(
  filename = "/Volumes/Promise RAID/Line/projects/2402_Fabien_HighDimFlow/inputfcsfiles/re-wilded.fcs", 
  transformation = FALSE, truncate_max_range = FALSE)))

# add mice model information to the separate samples that are to be analysed together
SPF_FC$mice <- "SPF"
rewild_FC$mice <- "re-wilded"


# Individual samples based on timestamps - "sampleID"
# FlowSOM variables
n_meta <- 10 # set slightly higher than actual number as easier to combine some clusters after than having to few
seed <- 13
scaling <- FALSE

# Compute the FlowSOM object
fsom <- FlowSOM(input = as.matrix(SPF_FC[,c("SampleID","Time")]),
                scale = TRUE,
                colsToUse = c("SampleID","Time"),
                seed = seed,
                nClus = n_meta)
# add to FC file
SPF_FC$samples <- GetMetaclusters(fsom)
# plot and check matches
ggplot(SPF_FC,aes(x=Time,y=SampleID,colour=samples))+
  geom_point()+
  theme_classic()
# collapse 5+6
SPF_FC[SPF_FC$samples %in% c(5,6),]$samples <- 5


# Compute the FlowSOM object
fsom <- FlowSOM(input = as.matrix(rewild_FC[,c("SampleID","Time")]),
                scale = TRUE,
                colsToUse = c("SampleID","Time"),
                seed = seed,
                nClus = n_meta)
# add to FC file
rewild_FC$samples <- GetMetaclusters(fsom)
# plot and check matches
ggplot(rewild_FC,aes(x=Time,y=SampleID,colour=samples))+
  geom_point()+
  theme_classic()
# collapse 3+5
rewild_FC[rewild_FC$samples %in% c(3,5),]$samples <- 3

# now gather FC files in one data frame 
# check matchomg column names
colnames(SPF_FC) %in% colnames(rewild_FC)

data_FC <- rbind(SPF_FC,rewild_FC)

head(data_FC)
dim(data_FC)

## change colnames so they fit excel sheet EXACTLY and exchange with actual markers
#remoce FJcomp
colnames(data_FC) <- str_replace(colnames(data_FC),"FJComp-","")
#remove end -A
colnames(data_FC) <- str_replace(colnames(data_FC),"-A","")
#replace Alexa Fluor with AF
colnames(data_FC) <- str_replace(colnames(data_FC),"Alexa Fluor","AF")
#replace eFluor with EF
colnames(data_FC) <- str_replace(colnames(data_FC),"eFluor","EF")
#remove empty spaces - will make life easier for you downstream, programming languages hate empty spaces
colnames(data_FC) <- str_replace_all(colnames(data_FC)," ","")
# for samples to be unique per sample add mice model/condition before number
data_FC$samples <- paste(data_FC$mice, data_FC$samples, sep = "_")


# now open csv file of panel and overwrite column name with matching marker
panel <- read.csv(paste(projdir,"/input/FF_panel.csv",sep=""), row.names = 2)
# check non overlaps, should not be any of the markers/channels you used
colnames(data_FC)[!colnames(data_FC) %in% rownames(panel)]
# exchange markers for channels
colnames(data_FC)[colnames(data_FC) %in% rownames(panel)] <- panel[colnames(data_FC)[colnames(data_FC) %in% rownames(panel)],]
#remove empty spaces and / - will make life easier for you downstream, programming languages hate empty spaces
colnames(data_FC) <- str_replace_all(colnames(data_FC)," ","")
colnames(data_FC) <- str_replace_all(colnames(data_FC),"/","-")

# check worked
head(data_FC)

#make copy of data_FC to work with
data_FC_df <- data_FC

#' select protein marker columns to use for dim. reduc.
# markers should NOT include markers used in gating strategy 
#if short list simply list all exactly as they are named in the csv file, otherwise you can exclude as below
marker_cols <- panel$marker[!panel$marker %in% c("L-D","CD45","Tetramer")]

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)
asinh_scale <- 150
data_FC_df[, marker_cols] <- asinh(data_FC_df[, marker_cols] / asinh_scale)

summary(data_FC)

#' #### ---- Running UMAP ---- ####
# possibility of downsampling
# downsample to equal amount of cells in each group/sample/what makes sense
# check group with fewest observations
table(data_FC_df$samples)

# downsampled for overall analysis of CD45+
down_FC <- downsample(data_FC_df, cat_col = "samples")
# or non downsampled for T cell analysis
down_FC <- data_FC_df

#n_sub <- 6000
#set.seed(1234)
#ix <- sample(1:length(labels), n_sub)
ix <- rownames(down_FC)
## or choose all cells 
#ix <- rownames(data_FC_df)
#down_FC <- data_FC_df

# prepare data for umapr (matrix format required)
data_umap <- down_FC[ix, marker_cols]
data_umap <- as.matrix(data_umap)
dups <- duplicated(data_umap)
data_umap <- data_umap[!dups, ]

# run umap, takes a long time with many cells, downsample first to test
umap_emb <- umap(data_umap, scale = "maxabs") 

# prepare umap embedding output data for plot
data_plot_umap <- as.data.frame(umap_emb)
colnames(data_plot_umap) <- c("UMAP_1", "UMAP_2")
head(data_plot_umap);dim(data_plot_umap)

# add group labels and merged CI label
data_plot_umap[,"sample"] <- down_FC[ix,]$sample
data_plot_umap[,"mice"] <- down_FC[ix,]$mice

#### Plot each individual group as seperate contour on top of total cells
ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2, colour = mice))+ 
  geom_point_rast(size=0.1)+
  theme_classic()

ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2, colour=mice))+ 
  #geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme_classic()
ggsave(paste(RAID_dir,"/output/",dato,"_UMAP_Contour_rewildvsSPF_plot.pdf",sep=""), height = 4, width = 4)

# rewilded only
ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=data_plot_umap[data_plot_umap$mice=="re-wilded",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#F8766D")+
  theme_classic()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"_UMAP_rewildedContour_plot.pdf",sep=""), height = 4, width = 4)

#SPF only
ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=data_plot_umap[data_plot_umap$mice=="SPF",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="#00BFC4")+
  theme_classic()+
  xlim(c(-13,13))+ylim(c(-15,17))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(RAID_dir,"/output/",dato,"_UMAP_SPFContour_plot.pdf",sep=""), height = 4, width = 4)


#### Colour by individual marker and save pdf versions
for (marker in marker_cols){
  plot_mark <- ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2, colour=data_umap[,marker]))+ 
    geom_point_rast(size=0.1)+
    scale_colour_viridis_c(option = "plasma")+
    theme_classic()+
    labs(colour=marker)+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"_UMAP_",marker,"_plot.pdf", sep=""), height = 4, width = 4)
  print(plot_mark)
  dev.off()
}


#### Now do clustering ####
# Individual samples based on timestamps - "sampleID"
# FlowSOM variables
n_max <- 18 # set slightly higher than actual number as easier to combine some clusters after than having to few
seed <- 13

for (i in seq(6,n_max)){
  n_meta <- i
  # Compute the FlowSOM object
  fsom <- FlowSOM(input = as.matrix(data_umap),
                  scale = FALSE,
                  colsToUse = marker_cols,
                  seed = seed,
                  nClus = n_meta)
  # add to FC file
  data_plot_umap[paste("cluster", i, sep = "_")] <- GetMetaclusters(fsom)
  # plot and check matches
  clus_plot <- ggplot(data_plot_umap,aes(x=UMAP_1,y=UMAP_2,colour = data_plot_umap[,paste("cluster", i, sep = "_")]))+
    geom_point_rast()+
    labs(colour=paste("cluster", i, sep = "_"))+
    theme_classic()+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  pdf(paste(RAID_dir,"/output/",dato,"_UMAP_","cluster", i,"_plot.pdf", sep=""), height = 4, width = 4)
  print(clus_plot)
  dev.off()
}

## save clustering to rds object for later analaysis (indicate downsampled or full set for continued analysis of T cells)
saveRDS(data_plot_umap, file = paste(dato,"CD45+clustered_SPFvsrewilded_NONdownsampled.rds",sep = "_"))
# or downsampled
saveRDS(data_plot_umap, file = paste(dato,"CD45+clustered_SPFvsrewilded_Downsampled.rds",sep = "_"))


# Now histograms based of each marker based on the clusters
for (marker in marker_cols){
  plot_mark <- ggplot(data = data_plot_umap, aes(x = data_umap[,marker], group = cluster_10, fill = cluster_10)) +
    geom_density(adjust = 1.5)+
    facet_wrap(~cluster_10)+
    theme_classic()+
    xlab(marker)+
    theme(axis.text.y = element_blank(), axis.ticks = element_blank())
  if (str_detect(marker,"/")){marker <- str_replace(marker,"/","-")}
  pdf(paste(RAID_dir,"/output/",dato,"cluster10_density_",marker,"_plot.pdf", sep=""), height = 6, width = 6)
  print(plot_mark)
  dev.off()
}

# see slides based on these, for ID of populations

#### Reclustering of T cells (non-downsampled) ####
# At resolution 8
# cl 5 are CD4 T cells 
# cl 6 are CD8 T cells
AllTcells <- cbind(data_plot_umap[data_plot_umap$cluster_10 %in% c(7,8),],data_umap[data_plot_umap$cluster_10 %in% c(7,8),marker_cols])
CD4Tcells <- cbind(data_plot_umap[data_plot_umap$cluster_10 %in% c(7),],data_umap[data_plot_umap$cluster_10 %in% c(7),marker_cols])
CD8Tcells <- cbind(data_plot_umap[data_plot_umap$cluster_10 %in% c(8),],data_umap[data_plot_umap$cluster_10 %in% c(8),marker_cols])

saveRDS(AllTcells, file = paste(dato,"AllTcellsreclustered_SPFvsrewilded_NONdownsampled.rds",sep = "_"))
saveRDS(CD4Tcells, file = paste(dato,"CD4Tcellsreclustered_SPFvsrewilded_NONdownsampled.rds",sep = "_"))
saveRDS(CD8Tcells, file = paste(dato,"CD8Tcellsreclustered_SPFvsrewilded_NONdownsampled.rds",sep = "_"))

nrow(AllTcells) #226.520 cells
nrow(CD4Tcells) #111.647 cells
nrow(CD8Tcells) #114.873 cells



