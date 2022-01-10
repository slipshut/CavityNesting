import sys

### This script parses the raw output from bird_expression_analysis.R
### and writes the results to a single, more easily readable line.
### Written by Mark Hibbins

model_output = []
protein_ID = []
gene_name = []
gene_num = sys.argv[2]

with open(sys.argv[1], "r") as results_file:

	model_results = False 

	for line in results_file:

		splitline = line.replace('"', '').strip().split()

		if len(splitline) > 2:

			if splitline[1] == 'Protein':
				protein_ID.append(splitline[-1])

			if splitline[1] == 'Gene':
				gene_name.append(splitline[-1])

		
			if splitline[-1] == 'pMCMC':
				model_results = True

			if model_results == True:
				model_output.append(splitline)
			else:
				continue

model_output = model_output[1:-2]
protein_ID = protein_ID[0]
gene_name = gene_name[0]

for line in model_output:

	if line[5] == '<1e-04':
		pval = 0.0001
	else:
		pval = line[5]

	print(protein_ID, gene_name, gene_num, line[0], line[1], line[2], line[3], line[4], pval)
