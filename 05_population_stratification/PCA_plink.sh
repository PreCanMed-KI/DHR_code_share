#!/bin/sh

# check mismatching call betweem 1000 genomes and my samples
plink --vcf $HOME/Projects/DHR/data/Genetics/Illumina_VCF/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only.vcf \
	--bmerge $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR \
		--merge-mode 6 \
			--out $HOME/Projects/DHR/out/SNP_PCA/merge_mismatches

# merge and remove mismatching SNPs
plink --vcf $HOME/Projects/DHR/data/Genetics/Illumina_VCF/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only.vcf \
	--exclude $HOME/Projects/DHR/out/SNP_PCA/merge_mismatches.missnp \
		--make-bed --out $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_tmp

plink --bfile $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR \
	--exclude $HOME/Projects/DHR/out/SNP_PCA/merge_mismatches.missnp \
		--make-bed --out $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_tmp
 
 
plink --bfile $HOME/Projects/DHR/out/SNP_PCA/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only_tmp \
	 --bmerge $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_tmp \
		 --make-bed --out $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged

# Prune by pairwise LD
plink --bfile $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged \
	--geno 0.05 --maf 0.05 \
		--indep-pairwise 50 5 0.2 \
			--out $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged_pruned

plink --bfile $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged \
	--extract $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged_pruned.prune.in \
		--make-bed --out $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged_input_PCA

# Do PCA
plink --bfile $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged_input_PCA\
	--pca --out $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR_merged_pruned_PCA
