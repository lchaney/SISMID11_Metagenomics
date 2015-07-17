#Lab 9 R Code

library(phyloseq)
load("STAT.RData")

load("~/Desktop/SISMID11/Lab-04/STAT.RData")

library(ade4)
library(vegan)

##compute JSD distance matrix
normalizeSample = function(x){x / sum(x)}
physeq.norm = transformSampleCounts(phy, normalizeSample)
jsdnormphy <- distance(physeq.norm, method = "jsd")

#use adonis function of the ade4 package to
#test for difference between the levels of Location

#my try didn't work
#adonis( jsdnormphy ~ Location, data=physeq.norm, perm = 999)

samp.data <- data.frame(sample_data(physeq.norm))
adonis(jsdnormphy ~ Location, data=samp.data, permutations = 9999)

#output is similar to a regular anova table
#the p-value is low showing that in all the permutation we didn't see any___ - so really significant

adonis(jsdnormphy ~ Treatment, data=samp.data, permutations = 9999)


adonis(jsdnormphy ~ Treatment*Location, data=samp.data, permutations = 9999)
#not right maybe?

adonis(jsdnormphy ~ Treatment, strata = samp.data$Location, data=samp.data, permutations = 9999)
#this specificies that permutations only look within each strata
#it takes both groups (cecal and fecal) and centers them

jsd.pco = dudi.pco(cailliez(jsdnormphy), scannf=F, nf=2)


jsd.wca <- wca(jsd.pco, samp.data$Location)

s.class(jsd.wca$li, samp.data$Location)
s.class(jsd.wca$li, samp.data$Treatment)

#now to see where the differences are
dist.mat <-  as.matrix(jsdnormphy)
phoc.test = function(x){
  spair = (samp.data$Treatment %in% x)
  phoc - adonis(as.dist(dist.mat[spair, spair]) ~ Treatment[spair],
                data=samp.data, strata = samp.data$Location[spair])
  c(x, phoc$aov.tab[,6][1])
}
combn(levels(Treatment), 2, phoc.test)
#didn't work -- error did not find Treatment
