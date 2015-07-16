#Lab 6 practice

#set working directory
setwd("~/Desktop/SISMID11/Lab-06")

#load saved data
load("~/Desktop/SISMID11/Lab-06/STAT.RData") #that didn't work

load("~/Desktop/SISMID11/Lab-04/STAT.RData")

#load packages
library("phyloseq")
packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

library("plyr")
packageVersion("plyr")

theme_set(theme_bw())

#explore dataset a little
sample_names(phy)


fecal.phy <- subset_samples(phy, Location == 'fecal')

#Compute Jensen-Shannon Divergence distance of the fecal samples
jsdDist <- distance(phy, method = "jsd")

sample_names(fecal.phy)
rownames(phenotypes)
rownames(phenotypes) = paste("fecal", rownames(phenotypes), sep="_")
phenotypes = phenotypes[sample_names(fecal.phy),]

#set seed because you want it to be reproducable


#