iteremove(list=ls())
library(ape)
library(tidyverse)
library(phytools)
library(MCMCglmm)

### This script does PGLM analysis for aggression, testosterone ###
### and distance datasets. Written by Mark Hibbins ### 

### Load datasets ### 

agg <- read.csv(
  "10sp_Agg_MountSex.csv")

tes <- read.csv(
  "10sp_logT.csv")

bird_trees <- ape::read.nexus(
  "10sp_CNtree_1000.nex")

### Get majority-rules consensus tree ### 

consensus_tree <- ape::consensus(bird_trees, p = 0.5, rooted = TRUE)
consensus_tree <- phytools::consensus.edges(bird_trees, consensus.tree = consensus_tree)
consensus_tree$node.label <- NULL

plotTree(consensus_tree) #Plot the tree


### Zero-inflated binomial GLMM for attack rate - 2 level model ### 

#Set up dataframe

agg <- agg[c("Species.Scientific", "Nest.Type.Specific",
             "Nest.Type", "Sex", "Physical.Contact",
             "Distance..raw.epochs.", "Mount.Sex")]

# If running 2-level model (Obligate vs. Non-Obligate) instead of 3-level model (Obligate vs. Facultative vs. Open):
agg$Nest.Type.Combined <- ifelse(agg$Nest.Type.Specific == "Obligate.Cavity",
                                                           "Obligate.Cavity",
                                                           "Non_Ob")

agg$Species.Scientific <- as.factor(agg$Species.Scientific)

agg["No.Contact"] <- 60 - agg["Physical.Contact"]


#Inverse matrix for phylogeny as a random effect

inv.phylo <- inverseA(consensus_tree, nodes = "TIPS", scale = TRUE)

#Model priors; fixed random effect and residual variance

prior0 <- list(R = list(V = diag(2), nu = 0.002, fix = 2),
               G = list(G1 = list(V = 1, nu = 0.002)))

#Model fit

agg_model_2levels<- MCMCglmm(c(Physical.Contact, No.Contact) ~ trait - 1 + 
                                at.level(trait, 1): Nest.Type.Combined + 
                                at.level(trait, 1):Sex + 
                                at.level(trait, 1):Nest.Type.Combined:Sex,
                              rcov = ~idh(trait):units,
                              random = ~Species.Scientific,
                              family = "zibinomial",
                              prior = prior0,
                              ginverse=list(Species.Scientific = inv.phylo$Ainv),
                              data = agg, nitt = 2000000, verbose = TRUE)

#Check for convergence 

agg_check <- heidel.diag(agg_model_2levels$Sol)
print(agg_check)

#Save model results to file

sink("agg_model_results_2levels.txt")
print(summary(agg_model_2levels))
sink()


### Zero-inflated binomial GLMM for attack rate - 3 level model ### 

#Set up dataframe

agg <- agg[c("Species.Scientific", "Nest.Type.Specific",
             "Nest.Type", "Sex", "Physical.Contact",
             "Distance..raw.epochs.", "Mount.Sex")]

agg$Species.Scientific <- as.factor(agg$Species.Scientific)

agg["No.Contact"] <- 60 - agg["Physical.Contact"]


#Inverse matrix for phylogeny as a random effect

inv.phylo <- inverseA(consensus_tree, nodes = "TIPS", scale = TRUE)

#Model priors; fixed random effect and residual variance

prior0 <- list(R = list(V = diag(2), nu = 0.002, fix = 2),
               G = list(G1 = list(V = 1, nu = 0.002)))

#Model fit

agg_model_3levels <- MCMCglmm(c(Physical.Contact, No.Contact) ~ trait - 1 + 
                        at.level(trait, 1):Nest.Type.Specific + 
                        at.level(trait, 1):Sex + 
                        at.level(trait, 1):Nest.Type.Specific:Sex,
                      rcov = ~idh(trait):units,
                      random = ~Species.Scientific,
                      family = "zibinomial",
                      prior = prior0,
                      ginverse=list(Species.Scientific = inv.phylo$Ainv),
                      data = agg, nitt = 2000000, verbose = TRUE)

#Check for convergence 

agg_check <- heidel.diag(agg_model_3levels$Sol)
print(agg_check)

#Save model results to file

sink("agg_model_results_3levels.txt")
print(summary(agg_model_3levels))
sink()



### Gaussian GLMM for distance  - 2 levels ### 

#Model fit

# If running 2-level model (Obligate vs. Non-Obligate) instead of 3-level model (Obligate vs. Facultative vs. Open):
agg$Nest.Type.Combined <- ifelse(agg$Nest.Type.Specific == "Obligate.Cavity",
                                                           "Obligate.Cavity",
                                                           "Non_Ob")

dis_model <- MCMCglmm(Distance..raw.epochs.~Nest.Type.Combined*Sex,
                      random = ~Species.Scientific,
                      family = "gaussian",
                      ginverse=list(Species.Scientific = inv.phylo$Ainv),
                      data = agg, nitt = 100000, verbose = TRUE)

#Check convergence

dis_check <- heidel.diag(dis_model$Sol)
print(dis_check)

#Save results to file

sink("dis_model_results_2level.txt")
print(summary(dis_model))
sink()


### Gaussian GLMM for distance  - 3 levels ### 

#Model fit

dis_model_3levels <- MCMCglmm(Distance..raw.epochs.~Nest.Type.Specific*Sex,
                      random = ~Species.Scientific,
                      family = "gaussian",
                      ginverse=list(Species.Scientific = inv.phylo$Ainv),
                      data = agg, nitt = 100000, verbose = TRUE)

#Check convergence

dis_check <- heidel.diag(dis_model_3levels$Sol)
print(dis_check)

#Save results to file

sink("dis_model_results_3levels.txt")
print(summary(dis_model_3levels))
sink()


### Gaussian GLMM for testosterone - 2 level model ### 

#Set up dataframe

tes <- tes[c("Species.Scientific", "Nest.Type.Specific",
             "Sex", "logTestosterone")]

# If running 2-level model (Obligate vs. Non-Obligate) instead of 3-level model (Obligate vs. Facultative vs. Open):
agg$Nest.Type.Combined <- ifelse(agg$Nest.Type.Specific == "Obligate.Cavity",
                                                           "Obligate.Cavity",
                                                           "Non_Ob")

tes$Species.Scientific <- as.factor(tes$Species.Scientific)

#Model fit

tes_model <- MCMCglmm(logTestosterone~Nest.Type.Combined*Sex,
                      random = ~Species.Scientific,
                      family = "gaussian",
                      ginverse=list(Species.Scientific = inv.phylo$Ainv),
                      data = tes, nitt = 100000, verbose = TRUE)

#Check convergence

tes_check <- heidel.diag(tes_model$Sol)
print(tes_check)

#Save results to file

sink("tes_model_results_2levels.txt")
print(summary(tes_model))
sink()


### Gaussian GLMM for testosterone - 3 level model ### 

#Set up dataframe

tes <- tes[c("Species.Scientific", "Nest.Type.Specific",
             "Sex", "logTestosterone")]

tes$Species.Scientific <- as.factor(tes$Species.Scientific)

#Model fit

tes_model_3levels <- MCMCglmm(logTestosterone~Nest.Type.Specific*Sex,
                      random = ~Species.Scientific,
                      family = "gaussian",
                      ginverse=list(Species.Scientific = inv.phylo$Ainv),
                      data = tes, nitt = 100000, verbose = TRUE)

#Check convergence

tes_check <- heidel.diag(tes_model_3levels$Sol)
print(tes_check)

#Save results to file

sink("tes_model_results_3levels.txt")
print(summary(tes_model))
sink()
