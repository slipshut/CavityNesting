remove(list=ls())
library(tidyverse)
library(readxl)

### This script does permutation analyses to test for enrichment ###
### of phyloDEGs in lists of candidate aggression genes. 
### Written by Mark Hibbins ### 

# Read in phyloDEG and all-genes datasets

phyloDEG <- read.csv(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/results/bird_expression_agg_models_fdr_results.csv")
normalizedCounts_Species <- read.csv(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/datasets/expression/normalizedCounts_Species.csv")

# Read in candidate gene datasets using a function I found at 
# https://www.geeksforgeeks.org/how-to-read-a-xlsx-file-with-multiple-sheets-in-r/

read_multisheet <- function(fpath) {
  
  sheets <- readxl::excel_sheets(fpath)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fpath, sheet = x))
  df <- lapply(tibble, as.data.frame)
  names(df) <- sheets
  
  return(df) 
}

cand_df <- read_multisheet(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/datasets/cand_genes/CandidateAggGeneLists.xlsx")

rittschof <- cand_df$`Rittschof Table S4`
zhangjames <- cand_df$`Zhang-James Table S4 `
filby <- cand_df$Filby

### Get and clean aggression candidate gene list

colnames(rittschof) <- rittschof[1,]
rittschof <- rittschof[-1,]
rittschof_genes <- rittschof$`Mouse symbol`
rittschof_genes <- rittschof_genes[as.numeric(rittschof[,7]) < 0.05]
rittschof_genes <- unique(unlist(lapply(rittschof_genes, toupper)))


filby_genes <- filby$Gene[!is.na(filby$Gene)]
filby_genes <- unique(unlist(lapply(filby_genes, toupper)))

zhangjames_genes <- unique(zhangjames$Genes)

cand_genes <- unique(c(rittschof_genes, filby_genes, zhangjames_genes))


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
filby_random_overlaps <- replicate(10000,
                             random_set_overlap(all_bird_genes,
                                                length(obs_genes),
                                                filby_genes))

rittschof_random_overlaps <- replicate(10000,
                             random_set_overlap(all_bird_genes,
                                                length(obs_genes),
                                                rittschof_genes))

zhangjames_random_overlaps <- replicate(10000,
                             random_set_overlap(all_bird_genes,
                                                length(obs_genes),
                                                zhangjames_genes))


pooled_random_overlaps <- replicate(10000,
                             random_set_overlap(all_bird_genes,
                                             length(obs_genes),
                                             cand_genes))

#Observed overlaps for each dataset
filby_obs_overlap <- gene_set_overlap(obs_genes, filby_genes)
rittschof_obs_overlap <- gene_set_overlap(obs_genes, rittschof_genes)
zhangjames_obs_overlap <- gene_set_overlap(obs_genes, zhangjames_genes)
pooled_obs_overlap <- gene_set_overlap(obs_genes, cand_genes)

#Plot results 

overlaps <- c(filby_random_overlaps, rittschof_random_overlaps,
              zhangjames_random_overlaps, pooled_random_overlaps)

datasets <- rep(c("filby", "rittschof", "zhangjames", "pooled"), each = 10000)
random_df <- as.data.frame(cbind(as.numeric(overlaps), datasets))

observed <- c(filby_obs_overlap, rittschof_obs_overlap, 
              zhangjames_obs_overlap, pooled_obs_overlap)
obs_dataset <- c("filby", "rittschof", "zhangjames", "pooled")
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

filby_pval <- get_enrichment_pval(filby_obs_overlap, filby_random_overlaps)

rittschof_pval <- get_enrichment_pval(rittschof_obs_overlap, 
                                      rittschof_random_overlaps)
zhangjames_pval <- get_enrichment_pval(zhangjames_obs_overlap,
                                       zhangjames_random_overlaps)
pooled_pval <- get_enrichment_pval(pooled_obs_overlap,
                                   pooled_random_overlaps)





