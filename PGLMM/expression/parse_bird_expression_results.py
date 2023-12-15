import sys

### This script parses the raw output from bird_expression_analysis.R
### and writes the results to a single, more easily readable line.
### Written by Mark Hibbins

model_output = []
protein_ID = []
gene_name = []
gene_num = sys.argv[2]
term = []
coef = []
pval = []

with open(sys.argv[1], "r") as results_file:

	append_coef = False
	append_pval = False 

	for line in results_file:
								
		splitline = line.replace('"', '').strip().split()
		#print(splitline)	

		if len(splitline) > 2:

			if splitline[1] == 'Protein':
				protein_ID.append(splitline[-1])

			if splitline[1] == 'Gene':
				gene_name.append(splitline[-1])

		
		if len(splitline) == 1 and splitline[0] == 'pMCMC':
			append_pval = True
			append_coef = False
		elif len(splitline) > 1 and splitline[0] == 'Location':
			append_coef = True
			append_pval = False
		elif len(splitline) == 1 and splitline[0] == "---":
			append_coef = False
			append_pval = False

		if append_coef == True and len(splitline) > 1:
			term.append(splitline[0])
			coef.append(splitline[1])
		elif append_pval == True and len(splitline) > 1:
			if splitline[1] == "<":
				pval.append(splitline[2])
			else:
				pval.append(splitline[1])
		else:
			continue





model_output = model_output[1:-2]
protein_ID = protein_ID[0]
gene_name = gene_name[0]
term = term[2:]
coef = coef[2:]

for i in range(len(pval)):

	if pval[i] == '1e-04':
		pval[i] = 0.0001

	print(protein_ID, gene_name, gene_num, term[i], coef[i], pval[i])
