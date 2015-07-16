data("soilrep")
ntaxa(soilrep)
nsamples(soilrep)
sample_names(soilrep)
taxa_names(soilrep)[1:20] #truncated

rank_names(soilrep) #doesn't work because there is no tax table with this data
sample_variables(soilrep)


otu_table(soilrep)[1:5, 1:5] #truncated length and width
sample_data(soilrep)[1:5,] #truncated length
tax_table(soilrep)[1:5, 1:4] #not available
phy_tree(soilrep) #not available


dt = data.table(psmelt(soilrep))
dt

ggplot(dt[Abundance > 0L], aes(Treatment, Abundance)) + 
  geom_violin() +
  geom_point(size = 3,
             alpha = 0.5,
             position = position_jitter(width = 0.1, height = 0.01)) +
  scale_y_sqrt()


cpsf = filter_taxa(soilrep, function(x) mean(x) > 1, TRUE)

mcps = merge_samples(soilrep, "Treatment")
sample_data(mcps)
sample_names(mcps)

otu_table(mcps)[, 1:5] #truncated width


get_taxa_unique(closedps, "Phylum") #no phylum in this data
taxg = tax_glom(closedps, "Phylum") #no phylum in this data
taxg #no phylum in this data