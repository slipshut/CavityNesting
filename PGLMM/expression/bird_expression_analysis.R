library(ape)
library(tidyverse)
library(phytools)
library(MCMCglmm)

### This script conducts PGLM analysis for expression of a single gene,
### with the gene's numerical index speciffied as a command-line 
### argument. Written by Mark Hibbins

args = commandArgs(trailingOnly=TRUE)

### Load datasets ###
exp <- read.csv("normalizedCounts_Species.csv")
trees <- ape::read.nexus("10sp_CNtree_1000.nex")
agg <- read.csv("10sp_Agg_avg.csv")

### Covariance matrix for consensus tree ###

#Consensus tree needs to be read in from newick 
#to circumvent issue with internal node labels
consensus_tree <- ape::read.tree("bird_consensus.newick")
inv.phylo <- inverseA(consensus_tree, nodes = "TIPS", scale = TRUE)

### Data tidying ###

exp_long <- exp %>%
  pivot_longer(!ZebraFinchProteinID:GeneDescription, #long format 
  names_to = "species", values_to = "expression") %>%
    mutate(log_expression = log(expression), #log expression
           sex = ifelse(grepl("M", species), "M", "F"), #sex
           species_scientific = case_when(grepl("BB", species) ~ "Sialia_sialis", #Species scientific names
					  grepl("BS", species) ~ "Hirundo_rustica",
					  grepl("CW", species) ~ "Thryothorus_ludovicianus",
					  grepl("ET", species) ~ "Passer_montanus",
					  grepl("HS", species) ~ "Passer_domesticus",
					  grepl("HW", species) ~ "Troglodytes_aedon",
					  grepl("PW", species) ~ "Protonotaria_citrea",
					  grepl("RO", species) ~ "Turdus_migratorius",
					  grepl("KYTS", species) ~ "Tachycineta_bicolor",
					  grepl("YW", species) ~ "Dendroica_petechia"),
	   nesting_strat = case_when(grepl("BB", species) ~ "ob_cavity", #Nesting strategy
				     grepl("BS", species) ~ "open",
				     grepl("CW", species) ~ "fac_cavity",
				     grepl("ET", species) ~ "ob_cavity",
				     grepl("HS", species) ~ "fac_cavity",
				     grepl("HW", species) ~ "ob_cavity",
				     grepl("PW", species) ~ "ob_cavity",
				     grepl("RO", species) ~ "open",
				     grepl("KYTS", species) ~ "ob_cavity",
				     grepl("YW", species) ~ "open"))		

exp_long <- as.data.frame(exp_long)

### Add attack rate data ###

nrows <- dim(exp_long)[1]
attrate <- rep(NA, nrows)

for (i in 1:nrows) {
    species <- (exp_long$species_scientific)[i]
    sex <- (exp_long$sex)[i]  
    attrate_index <- which(agg$Species.Scientific == species)
    sex_index <- which(agg$Sex == sex)
    index <- intersect(attrate_index, sex_index)
    attrate[i] = (agg$Attack.Rate)[index]
}

exp_long <- cbind(exp_long, attrate)
print(head(exp_long))
	   				 
## Function for phylogenetic analysis on a gene ### 

expression_phylo <- function(exp_df, phylo, gene_index) {

  gene_names <- unique(exp_df$GeneName) #list of gene names
  gene <- gene_names[gene_index] #get subset based on slice
  exp_df_subset <- subset(exp_df, GeneName == gene) #subset expression data

  exp_model <- MCMCglmm(log_expression ~ nesting_strat*sex*attrate,
 			random = ~species_scientific,
			family = "gaussian",
			ginverse=list(species_scientific = phylo$Ainv),
			data = exp_df_subset, nitt = 100000, verbose = FALSE)
  return(list(unique(exp_df_subset$ZebraFinchProteinID), gene, exp_model))

}

#Run model

results <- expression_phylo(exp_long, inv.phylo, as.numeric(args[1]))
model <- results[[3]]
model_check <- heidel.diag(model$Sol)

#Save results
sink(paste("bird_expression_agg_model_gene_", args[1], ".txt", sep = ""))
print(paste("Protein ID: ", results[[1]]))
print(paste("Gene name: ", results[[2]]))
print(model_check)
print(summary(model))
sink()


