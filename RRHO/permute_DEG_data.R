remove(list=ls())
library(dplyr)

#Load in DEG datasets

sparrow_results <- read.csv(
  "C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/CN_Species_DEGs_Sparrows.csv", header = TRUE)

swallow_results <- read.csv(
  "C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/CN_Species_DEGs_Swallows.csv", header = TRUE)

thrush_results <- read.csv(
  "C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/CN_Species_DEGs_Thrushes.csv", header = TRUE)

warbler_results <- read.csv(
  "C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/CN_Species_DEGs_Warblers.csv", header = TRUE)

wren_results <- read.csv(
  "C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/CN_Species_DEGs_Wrens.csv", header = TRUE)


permute_log2foldchange <- function(df) {
  
  perm <- sample(1:nrow(df))
  
  df[, 4:6] <- df[perm, 4:6]
  
  return(df)
}

sparrow_permuted <- permute_log2foldchange(sparrow_results)
swallow_permuted <- permute_log2foldchange(swallow_results)
thrush_permuted <- permute_log2foldchange(thrush_results)
warbler_permuted <- permute_log2foldchange(warbler_results)
wren_permuted <- permute_log2foldchange(wren_results)

write.csv(sparrow_permuted, "sparrow_permuted_DEGs.csv", row.names = FALSE)
write.csv(swallow_permuted, "swallow_permuted_DEGs.csv", row.names = FALSE)
write.csv(thrush_permuted, "thrush_permuted_DEGs.csv", row.names = FALSE)
write.csv(warbler_permuted, "warbler_permuted_DEGs.csv", row.names = FALSE)
write.csv(wren_permuted, "wren_permuted_DEGs.csv", row.names = FALSE)