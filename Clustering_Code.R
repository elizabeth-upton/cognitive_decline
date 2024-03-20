

library(reshape2)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(dplyr)
library(purrr)
library(stringr)
library(cluster)
library(factoextra)
library(ggplot2)
library(lme4)
library(mclust)
library(tidyr)
library(mixtools)
require(ggplot2)
require(kml)
require(traj)
require(lcmm)
require(lmerTest)
require(plyr)
require(tidyr)
require(psych)
require(fpc)
require(mclust)
require(rcompanion)

#This is an example of how we did the K-means clustering to get clusters distinguishing the 
#gender effect on subtest trajectories (see figure 4)
#Load Stacked Results In

#setwd("/Users/samuelliu/Downloads/Results_Sept8_Final/")
names <- list.files(pattern = ".rds") 
df <- purrr::map(names, readRDS)

names <- str_replace(names, "results_","") 
names <- str_replace(names, "__B1000_startseed1.rds","")
names <- str_replace(names, "_B1000_startseed1.rds","")


fitted.values <- matrix(nrow = 336)

# Put all data into one large dataframe, saving stacked mean for each age/gender/education level combo
for(j in 1:17){
  
  df_one <- df[[j]]
  new_list<- list()
  
  for (i in 1:1000){
    new_list[[i]] <- df_one[[i]]$PlotData
  }
  
  wide.preds <- do.call(rbind, new_list) %>% 
    group_by(age,Female, CollegeOrLess, AdvancedDeg) %>% 
    dplyr::summarize( ST.mean = mean(ST)) %>% 
    dplyr::mutate(Sex = ifelse(Female == 0, "M", "W")) 
  
  fitted.values <- cbind(fitted.values,wide.preds$ST.mean)
}

#Rename Variables for Ease of Understanding
fitted.values <- fitted.values[,-1]
fitted.values <- cbind(fitted.values, wide.preds[,c(1,3,4,6)], drop = F)
fitted.values <- fitted.values %>%
  mutate(Gender = Sex) %>%
  dplyr::select(-Sex)

colnames(fitted.values)[1:17] <- names
fitted.values <- as.data.frame(fitted.values) %>%
  mutate(gonogoten = -gonogoten, trailA = -trailA, trailB = -trailB)

#Pivoting data to longer for k-means clustering program
fitted.valueslong=fitted.values%>%pivot_longer(cols=(colnames(fitted.values)[1:17]),names_to='subtest.name',values_to='mean.score')
fitted.valueslong$education=''

fitted.valueslong$education <- ifelse(fitted.valueslong$CollegeOrLess == 0 & fitted.valueslong$AdvancedDeg == 0, "HS",
                                      ifelse(fitted.valueslong$AdvancedDeg == 1, "Advanced Degree",
                                             ifelse(fitted.valueslong$CollegeOrLess == 1, "College", fitted.valueslong$education)))



#Deriving the difference in stacked mean scores between college men and women, and adding it to dataframe
meanscore=(fitted.valueslong %>% 
             filter(education=='College') %>% 
             filter(Gender=='M'))$mean.score-(fitted.valueslong %>% filter(education=='College')%>% filter(Gender=='W'))$mean.score

fitted.valuesclus=fitted.valueslong %>% 
  filter(education=='College') %>%
  filter(Gender=='W')

fitted.valuesclus$mean.score=meanscore

fitted.valuesclus=fitted.valuesclus %>%
  pivot_wider(id_cols=subtest.name,names_from=age,values_from=mean.score)



dat_cld <- cld(fitted.valuesclus, timeInData = 2:57,idAll=fitted.valuesclus$subtest.name)   # function cld() builds a ClusterLongData object for clustering

set.seed(12345)
#choosing Euclidean distance
(option2 <- parALGO(distanceName="euclidean"))
#running the k-means clustering, allowing for 
kml(dat_cld, nbRedraw = 200, nbClusters = 2:7, toPlot = "none",parAlgo = option2)  # function kml() runs k-means with various starts (and clusters)

#noting different metrics for cluster number selection. Combination of metrics as
#well as comparing cluster analyses across different number of clusters
#results in us choosing 5 clusters.
plotAllCriterion(dat_cld, criterion=CRITERION_NAMES[c(1:8)], standardized=TRUE)

## add clusters to data frame
fitted.valuesclus$clust3kml <- getClusters(dat_cld, 5)
table(fitted.valuesclus$clust3kml)
fitted.valuesclus[,c(1,58)]
datlong=fitted.valueslong %>%
  filter(education=='College') %>%
  filter(Gender=='W')
datlong$mean.score=meanscore

