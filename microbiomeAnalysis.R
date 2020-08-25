# microbiome community diversity - Sadeq, Mills, Beckerman

# libraries
library(tidyverse)
library(vegan)
library(gplots)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('phyloseq')
library(phyloseq)

# import data at various scales (CHOOSE) ----

trim_data<- read_delim("./RTLGenomics_Guts/AnalysisB/TrimmedTaxa.genus.counts.txt",
                          delim = "\t")
# trim_data <- read_delim("./RTLGenomics_Guts/AnalysisB/TrimmedTaxa.family.counts.txt",
#                          delim = "\t")
# trim_data <- read_delim("./RTLGenomics_Guts/AnalysisB/TrimmedTaxa.order.counts.txt",
#                         delim = "\t")
# trim_data <- read_delim("./RTLGenomics_Guts/AnalysisB/TrimmedTaxa.class.counts.txt",
#                         delim = "\t")

glimpse(trim_data)
head(trim_data)

# prepare data for physeq  ----
an_scale <- "genus"
an_scale_labels <- c("Kingdom","Phylum","Class","Order","Family","Genus")

# label column 1 via taxonomy
names(trim_data)[1] <- an_scale
trim_data$genus

# check
head(trim_data)

# create sensible taxa names from scale
tmp_names <- trim_data %>% 
  separate(an_scale, sep = " ; ",
           into = an_scale_labels,
           fill = 'right') %>% 
  unite(taxa, an_scale_labels, sep = "_")

tmp_names$taxa

# generate physeq components -----

otus <- tmp_names %>% select(taxa)
otu_mat <- select(tmp_names, -taxa) %>% as.matrix
rownames(otu_mat) <- otus$taxa
colnames(otu_mat) <- c("Control_AB", "Cu_AB", "CuJuJu_AB", "JuJu_AB",
                       "Control", "Cu", "Cu_JuJu", "JuJu")
otu_mat

tax_mat <- otus %>% separate(taxa,into = c("Kingdom","Phylum","Class","Order","Family","Genus"),
                             sep = "_") %>% 
  as.matrix

rownames(tax_mat) <- otus$taxa

explanatory1 <- expand.grid(
  stress = factor(c("Control", "Cu", "Cu-Predation", "Predation"), levels = c("Control", "Cu", "Predation", "Cu-Predation")),
  antibiotic = factor(c("Antibiotic", "Control"), levels = c("Control", "Antibiotic"))) %>% 
  mutate(stress_ant = paste(stress, antibiotic, sep = "_"))

rownames(explanatory1) <- colnames(otu_mat)
explanatory1

# combine physeq components ------
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)                  
SAM = sample_data(explanatory1)

guts  <-  phyloseq(OTU, TAX, SAM)
guts

# Further Cleaning -----

guts2 <- prune_taxa(taxa_sums(guts) > 0, guts) %>% 
  subset_taxa(.,Class != "No Hit" & Class !="Unclassified" &
                Class != "Unknown")

# heat map ------

# 200 or more counts
guts3 <- prune_taxa(taxa_sums(guts) > 200, guts) %>% 
  subset_taxa(.,Class != "No Hit" & Class !="Unclassified" &
                Class != "Unknown")

# 200 bar - extra
plot_bar(guts3, fill = "Class", x = "stress", 
         facet_grid = ~antibiotic)

# hclust defaults for clustering: nice col-treatment labels
pheatmap(otu_table(guts3), 
         labels_col = guts3@sam_data$stress_ant)

# NMDS ordination using jaccard (or Bray dist mat) ----
# order y-axis by family and x by stress ant
# plot_heatmap(guts3, "NMDS", "jaccard", "stress_ant", "Genus")

bray_pcoa <-  ordinate(guts2, "NMDS", "bray")

plot_ordination(guts2, bray_pcoa, "samples", color="antibiotic")+
  geom_point(size = 5)+ xlim(-1.2,1.2)+ ylim(-0.75, 0.75)+
  geom_text(label = explanatory1$stress, size = 3, vjust = -1.2, hjust = -0.5)+
  geom_polygon(aes(fill=antibiotic), alpha= 0.3)+
  scale_fill_manual(values = c(Control = "grey", Antibiotic = "orange"))+
  annotate(geom = "text", x = 0.7, y= 0.25, label = "Point A")+
  annotate(geom = "text", x = 0.1, y= 0.05, label = "Point B")+
  annotate(geom = "text", x = -1, y= -0.62, label = "Point C-1,-2,-3 ->", angle = 90)+
  geom_segment(aes(x = 0.59, y = -0.74, xend = -0.096, yend = 0.058), 
               col = 'grey50', linetype = "dotted",
               arrow = arrow(length = unit(0.02, "npc")))+
  theme_bw(base_size = 15)

