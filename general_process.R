# Install necessary packages if not already installed
# BiocManager::install("MicrobiotaProcess")

# Load required libraries
library(ggplot2)
library(MicrobiotaProcess)
library(patchwork)
library(dada2)
library(phyloseq)
library(ggpubr)
library(tibble)
library(dplyr)
library(microViz)
library(vegan)
library(pairwiseAdonis)
library(rstatix)
library(data.table)
library(tidyr)
library(microbiome)
library(reshape2)
library(ComplexHeatmap)
library(gridExtra)
library(scales)
library(vcd)
library(stats)
library(car)
library(Hmisc)
library(magrittr)
library(microbiomeMarker)
library(microbiomeutilities)
library(jeevanuDB)
library(ANCOMBC)
library(tidyverse)
library(nlme)
library(compositions)
library(SpiecEasi)
library(seqtime)
library(MicEco)
library(iNEXT)
library(decontam)

# Read in data necessary for creating phyloseq object
otu_mat <- read_excel("seqALL.xlsx")
tax_mat <- read_excel("tax_All_filter.xlsx")
samples_df <- read_excel("metadata_human.xlsx")

# Transform OTU and taxonomy matrices to use the first column as row names
otu_mat <- otu_mat %>% tibble::column_to_rownames("otu")
tax_mat <- tax_mat %>% tibble::column_to_rownames("otu")
samples_df <- samples_df %>% tibble::column_to_rownames("sample")

# Convert data to matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Create phyloseq object components
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

# Create phyloseq object
PS <- phyloseq(OTU, TAX, samples)

# Check if any taxa have zero counts and summarize
any(taxa_sums(PS) == 0)
sum(taxa_sums(PS) == 0)
summarize_phyloseq(PS)

# Count ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "phylum"], exclude = NULL)

# Count reads for off-target eukaryotes and archaea
by((otu_table(PS)), tax_table(PS)[, "phylum"], sum)

# Filtering
# Remove samples with zero reads
PS <- prune_samples(sample_sums(PS) > 0, PS)
summarize_phyloseq(PS)

# Remove samples with less than 100 reads
PS <- prune_samples(sample_sums(PS) >= 100, PS)
summarize_phyloseq(PS)

# Filter out unclassified or unwanted taxa
PS <- subset_taxa(PS, !is.na(phylum) & !phylum %in% c("", "uncharacterized", "NA", "p__Vertebrata"))

# Save and reload the phyloseq object
saveRDS(PS, "PS.Rds")
PS <- readRDS("PS.Rds")

# General check for richness across groups
plot_richness(PS, x= "Blastocystis.Illumina", color = "Blastocystis.Illumina", measures = c("Observed","Simpson", "Shannon")) +
  geom_jitter(alpha= 0.005) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90))

# Create a prevalence data frame to examine low prevalent taxa
Prevdf <- apply(X = otu_table(PS), MARGIN = 1, FUN = function(x){sum(x > 0)})
Prevdf <- data.frame(Prevalence = Prevdf, TotalAbundance = taxa_sums(PS), tax_table(PS))

# Visualize prevalence vs abundance
ggplot(Prevdf, aes(TotalAbundance, Prevalence / nsamples(PS), color = phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +
  xlab("Log 10 Total Reads") +
  ylab("Prevalence [Prop. of Samples]") +
  theme_bw() +
  facet_wrap(~phylum) +
  theme(legend.position="none")

# Remove low prevalent ASVs (appearing in less than 10% of samples with more than 5 reads)
wh0 <- genefilter_sample(PS, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS))
PS <- prune_taxa(wh0, PS)

# Normalize the data
set.seed(1024)
rareres <- get_rarecurve(obj=PS, chunks=400)

# Create rarefaction curve
p_rare <- ggrarecurve(obj=rareres, indexNames=c("Observe", "Shannon", "ACE")) +
  theme(legend.spacing.y=unit(0.01,"cm"), legend.text=element_text(size=4))

# Barplot of relative abundance by groups
classtaxa <- get_taxadf(obj=PS.Norm, taxlevel=5)
pclass <- ggbartax(obj=classtaxa, facetNames="Burden_Illumina", topn=20) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
pclass

# Venn diagram of common taxa
vennlist <- get_vennlist(obj=PS, factorNames="Burden_Illumina")
library(VennDiagram)
vennp <- venn.diagram(vennlist, height=5, width=5, filename=NULL, fill=c("#00AED7", "#FD9347"), 
                      cat.col=c("#00AED7", "#FD9347"), alpha = 0.85, cex = 1.2, cat.cex = 1.3)

grid::grid.draw(vennp)

# Alpha diversity analysis
alphaobj <- get_alphaindex(PS)
p_alpha <- ggbox(alphaobj, geom="violin", factorNames="Burden_Illumina") + theme(strip.background = element_rect(colour=NA, fill="grey"))
p_alpha

# Beta diversity analysis with PCoA
pcoares <- get_pcoa(obj=PS.Norm, distmethod="bray", method ="hellinger")
pcoaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE, factorNames=c("Burden_Illumina"), ellipse=TRUE) +
  scale_color_manual(values=c("#FD9347","#00AED7",  "red","green")) +
  scale_fill_manual(values=c("#FD9347","#00AED7", "red","green"))
pcoaplot1

# PERMANOVA (Permutational Multivariate Analysis of Variance)
distme <- get_dist(PS.Norm, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(PS.Norm), check.names=FALSE)
set.seed(1024)
adores <- adonis2(distme ~ Burden_Illumina, data=sampleda, permutation=9999)
adores

# Differential abundance analysis using DESeq2
dds <- phyloseq_to_deseq2(PS.Norm, ~ Burden_Illumina)
geoMeans <- apply(counts(dds), 1, function(x) exp(sum(log(x[x > 0])) / length(x)))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds)

# Compare Infected vs Non-Infected
res_Infected_vs_NonInfected <- results(dds, contrast = c("Burden_Illumina", "Positive", "Negative"))

# Merge results with taxonomy table
tax_table_df <- as.data.frame(tax_table(PS.Norm))
res_df_Infected_vs_NonInfected <- as.data.frame(res_Infected_vs_NonInfected)
res_df_Infected_vs_NonInfected$OTU <- rownames(res_df_Infected_vs_NonInfected)
res_tax_Infected_vs_NonInfected <- merge(res_df_Infected_vs_NonInfected, tax_table_df, by.x = "OTU", by.y = "row.names")

# Volcano plot
EnhancedVolcano(res_tax_Infected_vs_NonInfected, lab = as.character(res_tax_Infected_vs_NonInfected$Genus), 
                x = 'log2FoldChange', y = 'pvalue', title = "", subtitle = "Volcano plot", 
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3.0, labSize = 3.0)

