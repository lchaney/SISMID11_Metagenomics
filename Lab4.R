setwd("~/Desktop/SISMID11/Lab-04")

#load data
otufile = "seqs_otu_table.txt"
mapfile = "obesity_mapping.txt"
trefile = "seqs_rep_set.tre.txt"

lab4data <- import_qiime(otufile, mapfile, trefile)

phenfile = "stat_phenotype.txt"
phendata <- read.table()



####to save time we have been given the data!
load("~/Desktop/SISMID11/Lab-04/STAT.RData")

#normalize the data
phy.norm <- transformSampleCounts(phy, function(x) x/sum(x))

sample_sums(phy.norm) #everything now adds up to 1

phylum.phy <- tax_glom(phy.norm, taxrank = "Phylum")
order.phy <- tax_glom(phy.norm, taxrank = "Order")

taxon = otu_table(phylum.phy)[1,]
test.res <- wilcox.test(c(taxon) ~ sample_data(phylum.phy)$Location)
test.res
names(test.res)
test.res$statistic
test.res$p.value

#testing difference of one taxa between location (fecal vs cecal)
#one fault in this is these are not complete independent samples a rank sum test would be more appropriate
#for pairs -- each mouse had a fecal and cecal sample
Location = sample_data(phylum.phy)$Location
test.taxon = function(taxon){
  test.res = wilcox.test(c(taxon) ~ Location)
  medians = tapply(c(taxon), Location, median)
  ns = tapply(c(taxon), Location, length)
  auc = test.res$statistic/(ns[1]*ns[2])
  c(medians, U=test.res$statistic, AUC=auc, p.value=test.res$p.value)
}

test.taxon(taxon)

#run the test for all taxa
phylum.location.test = apply(otu_table(phylum.phy), 1, test.taxon) #R is really slow when you do for loops, apply gets around that
phylum.location.test
phylum.location.test = t(phylum.location.test)
phylum.location.test

#add the bonferroni correction
phylum.location.test = cbind()