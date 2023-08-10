library(ape)
library(tidyverse)
library(phytools)
library(MCMCglmm)

eigengenes_all_ind <- read.csv(
   "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/wgcna/eigengenes_all_ind.csv")

bird_trees <- ape::read.nexus(
  "C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/datasets/10sp_CNtree_1000.nex"
)

agg <- read.csv("C:/Users/18126/OneDrive - University of Toronto/Projects/bird_expression/datasets/10sp_Agg_avg.csv")

### Get majority-rules consensus tree ### 

consensus_tree <- ape::consensus(bird_trees, p = 0.5, rooted = TRUE)
consensus_tree <- phytools::consensus.edges(bird_trees, consensus.tree = consensus_tree)
consensus_tree$node.label <- NULL
plot(consensus_tree)

### Phylogenetic covariance matrix ###
inv.phylo <- inverseA(consensus_tree, nodes = "TIPS", scale = TRUE)

### Create tidy dataframe ### 

eigenvalues_long <- eigengenes_all_ind %>%
  pivot_longer(!X, #long format 
               names_to = "module", values_to = "eigenvalue") %>%
  mutate(sex = ifelse(grepl("M", X), "M", "F"), #sex
         species_scientific = case_when(grepl("BB", X) ~ "Sialia_sialis", #Species scientific names
                                        grepl("BS", X) ~ "Hirundo_rustica",
                                        grepl("CW", X) ~ "Thryothorus_ludovicianus",
                                        grepl("ET", X) ~ "Passer_montanus",
                                        grepl("HS", X) ~ "Passer_domesticus",
                                        grepl("HW", X) ~ "Troglodytes_aedon",
                                        grepl("PW", X) ~ "Protonotaria_citrea",
                                        grepl("RO", X) ~ "Turdus_migratorius",
                                        grepl("KYTS", X) ~ "Tachycineta_bicolor",
                                        grepl("YW", X) ~ "Dendroica_petechia"),
         nesting_strat = case_when(grepl("BB", X) ~ "ob_cavity", #Nesting strategy
                                   grepl("BS", X) ~ "other",
                                   grepl("CW", X) ~ "other",
                                   grepl("ET", X) ~ "ob_cavity",
                                   grepl("HS", X) ~ "other",
                                   grepl("HW", X) ~ "ob_cavity",
                                   grepl("PW", X) ~ "ob_cavity",
                                   grepl("RO", X) ~ "other",
                                   grepl("KYTS", X) ~ "ob_cavity",
                                   grepl("YW", X) ~ "other"))

eigenvalues_long <- as.data.frame(eigenvalues_long[,2:6])

nrows <- dim(eigenvalues_long)[1]
attrate <- rep(NA, nrows)

for (i in 1:nrows) {
  species <- (eigenvalues_long$species_scientific)[i]
  sex <- (eigenvalues_long$sex)[i]  
  attrate_index <- which(agg$Species.Scientific == species)
  sex_index <- which(agg$Sex == sex)
  index <- intersect(attrate_index, sex_index)
  attrate[i] = (agg$Attack.Rate)[index]
}

eigenvalues_long <- cbind(eigenvalues_long, attrate)

#Want to fit a separate model for each module - split by module
module_eigenvals <- split(eigenvalues_long, eigenvalues_long$module)

### Function for phylogenetic analysis on an eigengene

eigen_phylo <- function(eigen_df, phylo) {
  
  eigen_model <- MCMCglmm(eigenvalue ~ nesting_strat*sex*attrate,
                        random = ~species_scientific,
                        family = "gaussian",
                        ginverse=list(species_scientific = phylo$Ainv),
                        data = eigen_df, nitt = 100000, verbose = FALSE)
  return(list(unique(eigen_df$module), eigen_model))
  
}

### Fit models 

model_modules <- c()
model_terms <- c()
model_coefs <- c()
model_pvals <- c()

eigen_test <- summary(eigen_phylo(module_eigenvals[[1]], inv.phylo)[[2]])$solutions

for (i in 1:length(module_eigenvals)) {
  
  eigen_model <- eigen_phylo(module_eigenvals[[i]], inv.phylo)
  eigen_model_results <- eigen_model[[2]]
  eigen_model_check <- heidel.diag(eigen_model_results$Sol)
  coefs <- summary(eigen_model_results)$solutions
  terms <- rownames(coefs)
  post_means <- unname(coefs[,1])
  pvals <- unname(coefs[,5])
  module_name <- eigen_model[[1]]
  
  for (j in 1:length(post_means)) {
    model_modules <- append(model_modules, module_name)
    model_terms <- append(model_terms, terms[j])
    model_coefs <- append(model_coefs, post_means[j])
    model_pvals <- append(model_pvals, pvals[j])
  }
  
  print(paste("Finished writing ", as.character(i), 
              " out of ", as.character(length(module_eigenvals))))

}

model_results_df <- as.data.frame(cbind(model_modules, model_terms, model_coefs, model_pvals))
model_results_df$p_adj <- p.adjust(model_results_df$model_pvals, method = "BH")
write.csv(model_results_df, "wgcna_pglm_bothsexes_2levels_results.csv")
