#!/bin/sh

# Prune by pairwise LD
# do not consider --maf 0.05 (population too small)
plink --vcf $HOME/Projects/DHR/data/Genetics/Illumina_VCF/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only.vcf \
	--geno 0.05 --indep-pairwise 50 5 0.2 \
		--out $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_pruned

plink --vcf $HOME/Projects/DHR/data/Genetics/Illumina_VCF/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only.vcf \
	--extract $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_pruned.prune.in \
		--make-bed --out $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_pruned_input_PCA

# Do PCA
plink --bfile $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_pruned_input_PCA \
	--pca --out $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_pruned_PCA
