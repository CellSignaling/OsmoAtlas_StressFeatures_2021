#!/usr/bin/env bash

# Script for obtaining the motifs enriched in Hog1 upregulated dependent genes over independent genes

# Download S.cerevisiae annotation in gff
wget ftp://ftp.ensembl.org/pub/release-89/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.89.gff3.gz -P data/original_data

# Unzip
gunzip data/original_data/Saccharomyces_cerevisiae.R64-1-1.89.gff3.gz

# Convert to BED

gff2bed < data/original_data/Saccharomyces_cerevisiae.R64-1-1.89.gff3 > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89.bed

# Keep only genes
awk '{if ($4 ~ /gene/) print $0}' data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89.bed > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes.bed

# Remove last column since it is not necessary and could cause troubles
awk -F'\t' 'BEGIN { OFS = FS }; NF { NF -= 1 }; 1' data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes.bed > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_9col.bed

# Remove : from gene:YXXXXX since it causes trouble with René program
cat data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_9col.bed | sed "s/gene://g" > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_9col_fixed.bed

# Add chr to the chromosome column and Mito to ChrM
awk '{print "chr"$0}' data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_9col_fixed.bed | sed 's/chrMito/chrM/g' > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_9col_fixed_format.bed

# Create upstream flanks
bedtools flank -i data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_9col_fixed_format.bed -g data/original_data/sacCer3_chromsizes.tsv -l 400 -r 0 -s > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_upstream.bed

# Divide in stress groups
for i in data/derived_data/*toGrep.txt; do
	echo $i
	f=$(basename "$i" _toGrep.txt)
	grepVF.pl -c 4 -a $i -b data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_upstream.bed > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_upstream_${f}_toHOMER.bed
done


# Dependent using independent as background
findMotifsGenome.pl data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_upstream_upregulated_Hog1_dependent_toHOMER.bed  \
sacCer3 results/S_cerevisiae/HOMER/ -size given \
-bg data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_genes_upstream_upregulated_Hog1_independent_toHOMER.bed


# Parse results Hog1 dep vs ind as background
## Upstream
files=$(find results/S_cerevisiae/HOMER/homerResults | grep -v "RV" | grep -v "similar" | grep ".motif$" | sort -V)

# Make full motif file
cat $files > results/S_cerevisiae/HOMER/homerMotifs_upstream_upregulated_dep_vs_ind_bg_full.motifs