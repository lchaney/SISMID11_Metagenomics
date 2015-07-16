#this is code for Lab 2
setwd("~/Desktop/SISMID11/Lab-02-01-phyloseq-intro")


NTaxaSim = 10
NSampSim = 8
randomNBValues = rnbinom(NTaxaSim*NSampSim, 0.15, 0.01)
randomNBValues
otumat = matrix(randomNBValues, nrow = NTaxaSim, ncol = NSampSim)
otumat


rownames(otumat)
colnames(otumat)
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat

taxmat = matrix(sample(letters[1:10], 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat[, "Domain"] <- "Bacteria"
taxmat[, "Phylum"] <- sample(x = c("Firmicutes", "Bacteroidetes", "Proteobacteria"),
                             size = nrow(taxmat), replace = TRUE)
class(otumat)
class(taxmat)

otumat
taxmat

#not sure here
install.packages("devtools")
library("devtools")
install_github("phyloseq", "joey711")


library("phyloseq")
packageVersion("phyloseq")


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX
physeq = phyloseq(OTU, TAX)
physeq


plot_bar(physeq, fill = "Family")

plot_heatmap(physeq)

sampledata = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = rgeom(nsamples(physeq), 0.2),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata

sample_names(sampledata) %in% sample_names(physeq)

random_tree = ape::rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1

physeq2 = phyloseq(OTU, TAX, sampledata, random_tree)
physeq2

identical(physeq1, physeq2)


plot_tree(physeq1, 
          color="Location",
          label.tips="taxa_names",
          size = "Abundance",
          justify = "left",
          ladderize="left",
          plot.margin=0.3)
          
plot_heatmap(physeq1, taxa.label="Phylum")


biomfile = "rich_sparse_otu_table.biom"
treefile = "biom-tree.phy"
import_biom(biomfile, treefile, parseFunction=parse_taxonomy_greengenes)
		
		#source("http://bioconductor.org/biocLite.R")
		#biocLite("BiocUpgrade")
		
biomfile = "otu_table.biom"
treefile = "rep_set.tre"
ps0 = import_biom(biomfile, treefile, parseFunction=parse_taxonomy_default)
ps0 <- merge_phyloseq(ps0, import_qiime_sample_data("Fasting_Map.txt"))
ps0

plot_tree(ps0, color = "Treatment", ladderize = "left", justify = "left")

OTU = import_usearch_uc("closed_map.uc")
SD = import_qiime_sample_data("Fasting_Map.txt")
tree = read_tree_greengenes("13_8_97_otus_unannotated.tree")
closedps = phyloseq(OTU, SD, tree)
closedps

plot_tree(closedps, color = "Treatment", ladderize = "left", justify = "left")

OTU = import_usearch_uc("denovo_map.uc")
dim(OTU)


library("data.table")
taxDT = fread("13_8_97_otu_taxonomy.txt", 
                        sep="\t", header=FALSE, colClasses="character") 
setnames(taxDT, c("OTU", "taxstring"))
setkey(taxDT, "OTU")
closedtaxDT = taxDT[taxa_names(closedps), ]
closedtaxlist = lapply(
  lapply(closedtaxDT$taxstring,
         function(x){strsplit(x, split="; ", fixed=TRUE)[[1]]}
         ), parse_taxonomy_greengenes)
names(closedtaxlist) <- closedtaxDT$OTU
closedtax = build_tax_table(closedtaxlist)

closedps <- merge_phyloseq(closedps, closedtax)
closedps

plot_tree(closedps, color = "Treatment", label.tips = "Phylum",
          ladderize = "left", justify = "left")
          
          
otufile = system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
mapfile = system.file("extdata", "master_map.txt", package="phyloseq")
trefile = system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
rs_file = system.file("extdata", "qiime500-refseq.fasta", package="phyloseq")
qiimedata = import_qiime(otufile, mapfile, trefile, rs_file)
qiimedata