# CavityNesting_Analysis

Analysis scripts for Lipshutz, Hibbins, et al. 2023. A brief description of the purpose 
of each script is outlined here. Scripts were run in Python v3 and R v4.

## Behavior_Testosterone

**CN_behavior_analyses.R** - Examines and visuzalizes attack rate and distance from mount by 
nest strategy and sex, from simulated territorial intrusions.

**CN_testosterone_analyses.R** - Examines and visuzalizes plasma testosterone levels by 
nest strategy and sex, from enzyme immunoassay.

**Cavity.Agg.Collect.Allyears.csv** - Input file for behavior and testosterone analyses

## DEGs

**DEGs_vs_Divergence.R** - Visualizes and correlates number of differentially expressed genes between 
species pairs by degree of evolutionary divergence, for each of five families

**Divergence_DEGs.csv** - Input file for DEGs vs. Divergence

## Enrichment
**bentz_gene_permutation_aggmodels.R** - Tests for significant overlap between our PhyloDEGs
associated with nesting strategy and a list of bird behavioral candidate genes generated by 
Bentz et al. 2021.

**gene_permutation_notbirds_aggmodel.R** - Tests for significant overlap between our PhyloDEGs
associated with nesting strategy and lists of behavioral candidate genes from several non-bird
studies. 

## PGLM
**aggression_phyloglm.R** - Fits phylogenetic generalized linear mixed models (PGLMs) for the 
relationship of attack rate with nesting strategy, testosterone, and physical distance. 

### WGCNA
**bird_wgcna_pglm.R** - Fits PGLMs with three nesting strategies to WGCNA modules inferred with
data from both sexes.

**bird_wgcna_pglm_F.R** - Fits PGLMs with three nesting strategies to WGCNA modules inferred for
females.

**bird_wgcna_pglm_M.R** - Fits PGLMs with three nesting strategies to WGCNA modules inferred for 
males. 

**wgcna_2level_pglm.R** - Fits PGLMs with two nesting strategies to WGCNA modules inferred with
data from both sexes.

**wgcna_pglm_2levels_F.R** - Fits PGLMs with two nesting strategies to WGCNA modules inferred for
females.

**wgcna_pglm_2levels_M.R** - Fits PGLMs with two nesting strategies to WGCNA modules inferred for 
males. 

### Expression

**bird_expression_2levels_model_fdr_correction.py** - Applies false discovery rate correction to 
expression PGLMs with two nesting strategies.

**bird_expression_agg_model_fdr_correction.py** - Applies false discovery rate correction to 
expression PGLMs with three nesting strategies. 

**bird_expression_analyis_2levels.R** - Fits a PGLM with two nesting strategies to log-expression
of an individual gene. 

**bird_expression_analysis_attrate.R** - Fits a PGLM with three nesting strategies to log-expression
of an individual gene. 

**get_convergent_agg_genes.R** - Gets putatively convergent genes significantly associated with 
aggression in three-strategy models. 

**get_convergent_nestingstrat_2level_genes.R** - Gets putatively convergent genes significantly 
associated with nesting strategy in two-strategy models.

**get_convergent_nestingstrat_genes.R** - Gets putatively convergent genes significantly
associated with nesting strategy in three-strategy models. 

**job_scripts/** - SLURM scripts for running gene expression PGLM analyses on high-performance
computing clusters. 

**summarize_bird_expression_2level_results.py** - Counts genes associated with each model term 
in two-strategy PGLMs, including the sign of the inferred coefficient. 

**summarize_bird_expression_agg_results.py** - Counts genes associated with each model term
in three-strategy PGLMs, including the sign of the inferred coefficient. 

## RRHO 

**CavNesting_RankDivergence.Rmd** - Conducts RRHO analysis for both sexes combined. 

**CavNesting_RankDivergence_Females.Rmd** - Conducts RRHO analysis for female expression. 

**CavNesting_RankDivergence_Males.Rmd** - Conducts RRHO analysis for male expression. 

**count_RRHO_shared.R** - Counts highly concordantly expressed genes shared across multiple
RRHO family comparisons. 

## WGCNA

**CN_raw_counts.csv** - Expression data input file for WGCNA, for both sexes

**CN_WGCNA_traits.csv** - Trait data input file for WGCNA, for both sexes

**WGCNA_cavity_Nesting_allind_cytoscape.R** - Runs weighted gene co-expression network analysis
for gene expression using data from both sexes. Includes output for cytoscape network visualization.

**WGCNA_cavity_Nesting_male.R** - Runs weighted gene co-expression network analysis
for gene expression using data from males only.

**WGCNA_cavity_nesting_female.R** - Runs weighted gene co-expression network analysis
for gene expression using data from females only.

**WGCNA_downstream_analyses.R** - Visualizes WGCNA by nest strategy, aggression, and sex. For samples 
with behavioral assays followed by immediate collection, includes latency analysis.

**WGCNA_latency.csv** - Input file for WGCNA latency analysis, including module eigengenes.






