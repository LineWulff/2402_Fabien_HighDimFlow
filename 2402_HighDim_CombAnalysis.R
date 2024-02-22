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

#### ---- variables used throughout script ---- ####
projdir <- getwd()
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

# check matchomg column names
#colnames(SPF_FC) %in% colnames(rewild_FC)

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



# samples based on timestamps - "sampleID"
# 17. Specify the FlowSOM variables
#SOM_x <- 10
#SOM_y <- 10
n_meta <- 8
seed <- 13
scaling <- FALSE

# 18. Compute the FlowSOM object
fsom <- FlowSOM(input = SPF_FC,
                scale = TRUE,
                colsToUse = c("SampleID","Time"),
                seed = seed,
                nClus = n_meta)

SPF_FC$samples <- GetMetaclusters(fsom)

SPF_FC <- as.data.frame(SPF_FC)
# plot and check matches
ggplot(SPF_FC,aes(x=Time,y=SampleID,colour=samples))+
  geom_point()+#scale_color_manual(values=c("red","blue","cyan","magenta"))+
  #geom_hline(yintercept=c(190000,90000))+
  #geom_vline(xintercept = c(222000,238000))
  theme_classic()

## correct names for protein expression
colnames(data_FC_df) <- sub("^FJComp-","",colnames(data_FC_df))
colnames(data_FC) <- sub("^FJComp-","",colnames(data_FC))
#correct laser channel to protein
channel_prot <- c("GATA-3","viability","CD4","T-bet","CD90.2","TCRb","FOXP3","RORyT","Lin","CD45")
names(channel_prot) <- colnames(data_FC_df)[5:14]

for (i in seq(5,length(channel_prot)+4)){
  colnames(data_FC_df)[i] <- channel_prot[colnames(data_FC_df)[i]]
}


#' select protein marker columns to use for dim. reduc.
# markers should NOT include markers used in gating strategy 
marker_cols <- c("FOXP3","RORyT","GATA-3","T-bet","CD4","TCRb")

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)
asinh_scale <- 150
data_FC_df[, marker_cols] <- asinh(data_FC_df[, marker_cols] / asinh_scale)

summary(data_FC)

#' #### ---- Running UMAP ---- ####
# possibility of downsampling
# downsample to equal amount of cells in each group
groupn <- c()
for (group in unique(data_FC_df$group)){
  groupn[group] <- nrow(data_FC_df[data_FC_df$group==group,])}
# check group with fewest observations
groupn

down_FC <- downsample(data_FC_df, cat_col = "group")

#n_sub <- 6000
#set.seed(1234)
#ix <- sample(1:length(labels), n_sub)
ix <- rownames(down_FC)

# prepare data for umapr (matrix format required)
data_umap <- down_FC[ix, marker_cols]
data_umap <- as.matrix(data_umap)
dups <- duplicated(data_umap)
data_umap <- data_umap[!dups, ]

umap_emb <- umap(data_umap)

# prepare umap embedding output data for plot
data_plot_umap <- as.data.frame(umap_emb)
colnames(data_plot_umap) <- c("UMAP_1", "UMAP_2")
head(data_plot_umap);dim(data_plot_umap)

# add group labels and merged CI label
data_plot_umap[,"sample"] <- down_FC[ix,]$sample
data_plot_umap[,"group"] <- down_FC[ix,]$group
data_plot_umap$mergedCI <- down_FC$group
data_plot_umap[startsWith(data_plot_umap$mergedCI,"CI"),]$mergedCI <- "CI"

#### Plot each individual group as seperate contour on top of total cells
group_col <- c("CI8"="lightslateblue","CI7-1"="darkslategray4","GF"="gray60")
ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d(data=data_plot_umap[data_plot_umap$group=="GF",],
                  aes(x = UMAP_1, y = UMAP_2, colour=group))+
  scale_color_manual(values = group_col)+
  theme_classic()+
  xlim(c(-10,12))+ylim(c(-13,11))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"lymphoid","UMAP_GF_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)

ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d(data=data_plot_umap[data_plot_umap$group=="CI7-1",],
                  aes(x = UMAP_1, y = UMAP_2, colour=group))+
  scale_color_manual(values = group_col)+
  theme_classic()+
  xlim(c(-10,12))+ylim(c(-13,11))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"lymphoid","UMAP_CI7-1_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)

ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d(data=data_plot_umap[data_plot_umap$group=="CI8",],
                  aes(x = UMAP_1, y = UMAP_2, colour=group))+
  scale_color_manual(values = group_col)+
  theme_classic()+
  xlim(c(-10,12))+ylim(c(-13,11))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"lymphoid","UMAP_CI8_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)


ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=data_plot_umap[data_plot_umap$mergedCI=="CI",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="skyblue3")+
  theme_classic()+
  xlim(c(-10,12))+ylim(c(-13,11))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"lymphoid","UMAP_CImerged_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)

#### Colour by individual marker and save pdf versions
for (marker in marker_cols){
  plot_mark <- ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2, colour=data_umap[,marker]))+ 
    geom_point_rast(size=1)+
    scale_colour_viridis_c(option = "plasma")+
    theme_classic()+
    labs(colour=marker)+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  pdf(paste(dato,"lymphoid","UMAP",marker,"plot.pdf", sep="_"), height = 4, width = 4)
  print(plot_mark)
  dev.off()
}


