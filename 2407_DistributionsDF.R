## Function to calculate distributions of cell populations between different samples
##

perc_df_calc <- function(input_df, sample_col, perc_col){
  perc_df <- data.frame()
  for (samp in unique(input_df[,sample_col])){
    samp_dist <- table(input_df[input_df[,sample_col]==samp, perc_col])/nrow(input_df[input_df[,sample_col]==samp,])*100
    perc_df <- rbind(perc_df, samp_dist)
  }
  rownames(perc_df) <- unique(input_df[,sample_col])
  colnames(perc_df) <- names(table(input_df[,perc_col]))
  perc_df$sample <- rownames(perc_df)
  perc_df <- gather(perc_df, cluster, percentage, 1:length(table(input_df[,perc_col])))
  return(perc_df)
}

##
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
