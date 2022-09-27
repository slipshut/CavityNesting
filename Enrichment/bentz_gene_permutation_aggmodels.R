remove(list=ls())
library(tidyverse)

### This script does permutation analyses to test for enrichment ###
### of phyloDEGs in lists of candidate aggression genes. 
### Written by Mark Hibbins ### 

# Read in phyloDEG and all-genes datasets

phyloDEG <- read.csv(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/results/bird_expression_agg_models_fdr_results.csv")
normalizedCounts_Species <- read.csv(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/datasets/expression/normalizedCounts_Species.csv")

# Read in and clean candidate gene dataset

Bentz_TRES_genes <- read.csv(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/datasets/cand_genes/Bentz_TRES_genes.csv", header=FALSE)

Bentz_TRES_genes[1,1] = "Symbol"
colnames(Bentz_TRES_genes) <- Bentz_TRES_genes[1,]
Bentz_TRES_genes <- Bentz_TRES_genes[-1,]

# Get candidate gene list for each set

HYPO_Day0_DEG_set <- unique(subset(Bentz_TRES_genes$Symbol,
                            Bentz_TRES_genes$Set == "HYPO_Day0_DEG"))
HYPO_Day2_DEG_set <- unique(subset(Bentz_TRES_genes$Symbol,
                                   Bentz_TRES_genes$Set == "HYPO_Day2_DEG"))
VMT_WGCNA_Blue_set <- unique(subset(Bentz_TRES_genes$Symbol,
                                   Bentz_TRES_genes$Set == "VMT_WGCNA_Blue"))
VMT_WGCNA_Brown_set <- unique(subset(Bentz_TRES_genes$Symbol,
                                    Bentz_TRES_genes$Set == "VMT_WGCNA_Brown"))
HYPO_WGCNA_Green_set <- unique(subset(Bentz_TRES_genes$Symbol,
                                    Bentz_TRES_genes$Set == "HYPO_WGCNA_Green"))
VMT_Day0_DEG_set <- unique(subset(Bentz_TRES_genes$Symbol,
                                   Bentz_TRES_genes$Set == "VMT_Day0_DEG"))
VMT_Day2_DEG_set <- unique(subset(Bentz_TRES_genes$Symbol,
                                  Bentz_TRES_genes$Set == "VMT_Day2_DEG"))
pooled_set <- unique(Bentz_TRES_genes$Symbol)

### Get list of phyloDEGs associated with nesting strat 

is_nestingstrat_gene <- function(gene) {
  
  if (grepl("nesting_strat", gene[4])) {
    return(gene[2])
  } 
  else {
    return(NA)
  }
}

obs_genes <- apply(phyloDEG, 1, is_nestingstrat_gene)
obs_genes <- unique(obs_genes[!is.na(obs_genes)])
obs_genes <- unique(unlist(lapply(obs_genes, toupper)))

## Get total pool of bird genes 

all_bird_genes <- normalizedCounts_Species$GeneName

## Function to calculate gene set overlap 

gene_set_overlap <- function(query, cand_set) {
  
  n_overlapping <- length(intersect(query, cand_set))
  overlap <- n_overlapping/length(query)
  
  return(overlap)
  
}

## Function to calculate overlap for a random draw of genes 

random_set_overlap <- function(pool, n_genes, cand_set) {
  
  random_set <- sample(pool, n_genes, replace = FALSE)
  random_overlap <- gene_set_overlap(random_set, cand_set)
  
  return(random_overlap)
}

## Generate sets of random overlaps 

HYPO_Day0_DEG_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      HYPO_Day0_DEG_set))
HYPO_Day2_DEG_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      HYPO_Day2_DEG_set))
HYPO_WGCNA_Green_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      HYPO_WGCNA_Green_set))
VMT_Day0_DEG_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      VMT_Day0_DEG_set))
VMT_Day2_DEG_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      VMT_Day2_DEG_set))
VMT_WGCNA_Blue_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      VMT_WGCNA_Blue_set))
VMT_WGCNA_Brown_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      VMT_WGCNA_Brown_set))
pooled_random_overlaps <- replicate(10000,
                                   random_set_overlap(all_bird_genes,
                                                      length(obs_genes),
                                                      pooled_set))

#Observed overlaps for each dataset
HYPO_Day0_DEG_obs_overlap <- gene_set_overlap(obs_genes, HYPO_Day0_DEG_set)
HYPO_Day2_DEG_obs_overlap <- gene_set_overlap(obs_genes, HYPO_Day2_DEG_set)
HYPO_WGCNA_Green_obs_overlap <- gene_set_overlap(obs_genes, HYPO_WGCNA_Green_set)
VMT_Day0_DEG_obs_overlap <- gene_set_overlap(obs_genes, VMT_Day0_DEG_set)
VMT_Day2_DEG_obs_overlap <- gene_set_overlap(obs_genes, VMT_Day2_DEG_set)
VMT_WGCNA_Blue_obs_overlap <- gene_set_overlap(obs_genes, VMT_WGCNA_Blue_set)
VMT_WGCNA_Brown_obs_overlap <- gene_set_overlap(obs_genes, VMT_WGCNA_Brown_set)
pooled_obs_overlap <- gene_set_overlap(obs_genes, pooled_set)

#Plot results 

overlaps <- c(HYPO_Day0_DEG_random_overlaps, HYPO_Day2_DEG_random_overlaps,
              HYPO_WGCNA_Green_random_overlaps, VMT_Day0_DEG_random_overlaps,
              VMT_Day2_DEG_random_overlaps, VMT_WGCNA_Blue_random_overlaps,
              VMT_WGCNA_Brown_random_overlaps, pooled_random_overlaps)

datasets <- rep(c("HYPO_Day0_DEG", "HYPO_Day2_DEG", "HYPO_WGCNA_Green",
                  "VMT_Day0_DEG", "VMT_Day2_DEG", "VMT_WGNCA_Blue",
                  "VMT_WGCNA_Brown", "pooled"), each = 10000)

random_df <- as.data.frame(cbind(as.numeric(overlaps), datasets))

observed <- c(HYPO_Day0_DEG_obs_overlap, HYPO_Day2_DEG_obs_overlap,
              HYPO_WGCNA_Green_obs_overlap, VMT_Day0_DEG_obs_overlap,
              VMT_Day2_DEG_obs_overlap, VMT_WGCNA_Blue_obs_overlap,
              VMT_WGCNA_Brown_obs_overlap, pooled_obs_overlap)

obs_dataset <- c("HYPO_Day0_DEG", "HYPO_Day2_DEG", "HYPO_WGCNA_Green",
                 "VMT_Day0_DEG", "VMT_Day2_DEG", "VMT_WGNCA_Blue",
                 "VMT_WGCNA_Brown", "pooled")

obs_df <- as.data.frame(cbind(as.numeric(observed), datasets = obs_dataset))



enrichment_plot <- ggplot() + 
  geom_density(data = random_df, aes(x = overlaps), size=1) +
  geom_vline(data = obs_df, mapping = aes(xintercept = observed),
             linetype = "dashed", size = 1) + 
  facet_wrap(~datasets, scales = "free") + 
  theme_bw(base_size = 20)

enrichment_plot

#Report p-values

get_enrichment_pval <- function(obs_overlap, random_overlaps) {
  
  enrichment_rank = 0
  
  for (i in 1:length(random_overlaps)) {
    if (obs_overlap >= random_overlaps[i]) {
      enrichment_rank <- enrichment_rank + 1
    }
    else {
      next 
    }
  }
  
  p = 1 - 2*abs(0.5 - (enrichment_rank/10000))
  
  return(p)
}

HYPO_Day0_DEG_pval <- get_enrichment_pval(HYPO_Day0_DEG_obs_overlap,
                                          HYPO_Day0_DEG_random_overlaps)
HYPO_Day2_DEG_pval <- get_enrichment_pval(HYPO_Day2_DEG_obs_overlap,
                                          HYPO_Day2_DEG_random_overlaps)
HYPO_WGCNA_Green_pval <- get_enrichment_pval(HYPO_WGCNA_Green_obs_overlap,
                                             HYPO_WGCNA_Green_random_overlaps)
VMT_Day0_DEG_pval <- get_enrichment_pval(VMT_Day0_DEG_obs_overlap,
                                          VMT_Day0_DEG_random_overlaps)
VMT_Day2_DEG_pval <- get_enrichment_pval(VMT_Day2_DEG_obs_overlap,
                                         VMT_Day2_DEG_random_overlaps)
VMT_WGCNA_Blue_pval <- get_enrichment_pval(VMT_WGCNA_Blue_obs_overlap,
                                         VMT_WGCNA_Blue_random_overlaps)
VMT_WGCNA_Brown_pval <- get_enrichment_pval(VMT_WGCNA_Brown_obs_overlap,
                                           VMT_WGCNA_Brown_random_overlaps)
pooled_pval <- get_enrichment_pval(pooled_obs_overlap,
                                   pooled_random_overlaps)





