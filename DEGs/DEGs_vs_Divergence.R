# 5.2.2023

# Cavity nesting gene expression vs. divergence time

# load in libraries
library(ggplot2)
library(ggrepel)

# read in data
DEGs_vs_Divergence <-read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/Divergence_DEGs_logfoldchange0.5.csv")

# plot
ggplot(DEGs_vs_Divergence, aes(x=Phylogeny.Divergence,y=DEGs)) + geom_point() + geom_text_repel(aes(label = Family))
DEGs_vs_Divergence.plot <-ggplot(DEGs_vs_Divergence, aes(x=Phylogeny.Divergence,y=DEGs)) + geom_point() + geom_text_repel(aes(label = Family)) + theme_classic()
DEGs_vs_Divergence.plot

# export as ppt figure
library (officer)
library (rvg)
editable_graph <- dml(ggobj = DEGs_vs_Divergence.plot)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "DEGs_vs_Divergence.plot.pptx")

# correlate divergence with gene expression
cor.test(DEGs_vs_Divergence$Phylogeny.Divergence, DEGs_vs_Divergence$DEGs, method = "pearson") 
#t = 2.5681, df = 3, p-value = 0.08263, r2 = 0.8290565 
