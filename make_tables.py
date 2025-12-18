# Usage: python make_tables.py
# Purpose: count GOSlim terms in a representative dataset 
# Flags:
#    -i, --input    : input best blast hits file, tab delimited, three columns 
#    -d, --database : mapping file for locus/gene to GO Slim category 


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = "input")
parser.add_argument('-d', '--database', dest = "database", default = './TAIR_2021/ATH_GO_GOSLIM_2021.txt')
args = parser.parse_args()
SLIM_DB = args.database
HITS = args.input


import pandas as pd

# Read in best blast hits and subset to include the subject of each hit, replace with column name that matches the GOSlim database 
best_hits = pd.read_table(HITS)
best_hits_S = best_hits[['Subject']]
best_hits_S.columns = ['ObjectName']

# Read in the GOSlim database and add names to each column 
GOSlim_DB = pd.read_table(SLIM_DB, sep='\t', comment='!')
GOSlim_DB.columns = ['locus','TAIRAccession', 'ObjectName', 'RelationshipType','GOTerm','GOID','TAIRKeywordID','Aspect','GOSlimTerm','EvidenceCode','EvidenceDescription','EvidenceWidth','Ref','Annotator','Date']

# Generate a new df with the locus information 
best_hits_S_locus = pd.DataFrame(best_hits_S['ObjectName'].str.replace('\.[0-9]','', regex = True))
best_hits_S_locus.columns = ['locus']

#Subset the GOSlim database to include only the ObjectName and corresponding GOSlimTerm
GOSlim_DB_OGOSlim_gene = GOSlim_DB[['ObjectName', 'GOSlimTerm']]
GOSlim_DB_GOSlim_locus = GOSlim_DB[['locus', 'GOSlimTerm']]

#Merge the GOSlim database and best hits files
GOSlim_merge_gene = best_hits_S.merge(GOSlim_DB_OGOSlim_gene,'inner', on='ObjectName')  
GOSlim_merge_locus = best_hits_S_locus.merge(GOSlim_DB_GOSlim_locus, 'inner', on='locus')

# Get raw and frequency counts at gene-level only 
rawcounts_gene = GOSlim_merge_gene['GOSlimTerm'].value_counts().sort_index()
rawcounts_gene.to_csv(HITS + '.rawcounts-gene.tsv',sep='\t')

freqcounts_gene = GOSlim_merge_gene['GOSlimTerm'].value_counts(normalize=True).sort_index()
freqcounts_gene.to_csv(HITS + '.freqcounts-gene.tsv',sep='\t')


# Get raw and frequency counts at locus-level 
rawcounts_locus = GOSlim_merge_locus['GOSlimTerm'].value_counts().sort_index()
rawcounts_locus.to_csv(HITS + '.rawcounts-locus.tsv',sep='\t')

freqcounts_locus = GOSlim_merge_locus['GOSlimTerm'].value_counts(normalize=True).sort_index()
freqcounts_locus.to_csv(HITS + '.freqcounts-locus.tsv',sep='\t')


# Below is code added by E Yaklich to create an output file with geneID and GO term mapping
# read in best hits table again
best_hits = pd.read_table(HITS)

#create a copy when subsetting to avoid SettingWithCopyWarning
best_hits_subset = best_hits[['Query', 'Subject']].copy()


# just keep transcript ID
best_hits_subset['Query'] = best_hits_subset['Query'].str.extract(r'([^\::]+)')


# rename columns for clarity
best_hits_subset.columns = ['gene_name', 'TAIR_ID']

# clean the 'TAIR_ID' column by removing version numbers (e.g., AT2G07760.1 -> AT2G07760)....i.e. going from transcript to locus. This will remove any number of digits or letters, ex ID.t18 or ID.UORF1
# this is needed to map to the GOSlim database
best_hits_subset['TAIR_ID'] = best_hits_subset['TAIR_ID'].str.replace(r'\..*', '', regex=True)

# subset the GOSlim database to include the TAIRAccession (Arabidopsis TAIR ID), the associated GO term and GOID, and Evidence Code
GOSlim_DB_subset = GOSlim_DB[['locus', 'GOID', 'GOTerm','EvidenceCode']].copy()  # make a copy to avoid SettingWithCopyWarning

# Check example sample cleaned data from both files
print("Sample Best Hits (TAIR IDs cleaned):")
print(best_hits_subset.head())

print("Sample GOSlim DB (TAIRAccession cleaned):")
print(GOSlim_DB_subset.head())

# merge based on cleaned-up 'TAIR_ID' and 'locus' (from GOSlim DB)
GOSlim_merge_gene = best_hits_subset.merge(GOSlim_DB_subset, 'inner', left_on='TAIR_ID', right_on='locus')

# check the first few rows of the merged result to ensure there are matches
print("Merged Data (first few rows):")
print(GOSlim_merge_gene.head())


# drop duplicate rows
GOSlim_merge_gene = GOSlim_merge_gene.drop_duplicates()

# save to file: keeping gene name, TAIR ID, and GOSlim terms, Go Term, Evidence Code
output_file = HITS + '.genes_to_GO_terms.tsv'
GOSlim_merge_gene[['gene_name', 'TAIR_ID', 'GOID', 'GOTerm','EvidenceCode']].to_csv(output_file, sep='\t', index=False)

print(f"Results saved to {output_file}")












