remove(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)

### Function to read in excel file with multiple tabs

multiplesheets <- function(fname) {

  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)

  # assigning names to data frames
  names(data_frame) <- sheets

  # print data frame
  return(data_frame)
}

#List of gene counts

RRHO.Permuted.Both.Sexes.Down <- read.csv("C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/RRHO Permuted Both Sexes Down.csv")
RRHO.Permuted.Both.Sexes.Up <- read.csv("C:/Users/mhibb/OneDrive - University of Toronto/Projects/bird_expression/datasets/RRHO/RRHO Permuted Both Sexes Up.csv")

RRHO_pivot_list <- list(RRHO.Permuted.Both.Sexes.Down,
                        RRHO.Permuted.Both.Sexes.Up)

#Count number of genes for each category

count_genes_RRHO <- function(counts_df) {

  gene_counts <- table(counts_df[,2])[2:10]
  gene_counts_df <- as.data.frame(matrix(gene_counts))

  return(gene_counts_df)
}

count_genes_RRHO_list <- lapply(RRHO_pivot_list, count_genes_RRHO)
num_genes <- dplyr::bind_rows(count_genes_RRHO_list)
num_comparisons <- rep(c(2, 3, 4, 5, 6, 7, 8, 9, 10), times = 2)
direction <- rep(c("Lower", "Higher"), each = 9)

counts_df <- as.data.frame(cbind(num_genes, num_comparisons, direction))
names(counts_df)[names(counts_df) == 'V1'] <- 'num_genes'

#Plot

counts_plot <- ggplot(counts_df, aes(x = num_comparisons, y = num_genes)) +
                geom_line(aes(linetype = direction), size = 1.5) +
                theme_minimal(base_size = 20) +
                labs(x = "Number of RRHO comparisons", y = "Number of shared genes",
                     linetype = "Sign", shape = "Sign") 

#Get list of genes shared by at least 5 comparisons

bothsexes_up <- subset(RRHO_pivot_list[[1]], gene_count >= 5)
bothsexes_up <- bothsexes_up[2:nrow(bothsexes_up), 1]
write.csv(bothsexes_up, "RRHO_shared_bothsexes_up.csv")

bothsexes_down <- subset(RRHO_pivot_list[[2]], gene_count >= 5)
bothsexes_down <- bothsexes_down[2:nrow(bothsexes_down), 1]
write.csv(bothsexes_down, "RRHO_shared_bothsexes_down.csv")

females_up <- subset(RRHO_pivot_list[[3]], gene_count >= 5)
females_up <- females_up[2:nrow(females_up), 1]
write.csv(females_up, "RRHO_shared_females_up.csv")

females_down <- subset(RRHO_pivot_list[[4]], gene_count >= 5)
females_down <- females_down[2:nrow(females_down), 1]
write.csv(females_down, "RRHO_shared_females_down.csv")

males_up <- subset(RRHO_pivot_list[[5]], gene_count >= 5)
males_up <- males_up[2:nrow(males_up), 1]
write.csv(males_up, "RRHO_shared_males_up.csv")

males_down <- subset(RRHO_pivot_list[[6]], gene_count >= 5)
males_down <- males_down[2:nrow(males_down), 1]
write.csv(males_down, "RRHO_shared_males_down.csv")
