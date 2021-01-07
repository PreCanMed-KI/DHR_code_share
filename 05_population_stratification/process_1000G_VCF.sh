#!/bin/sh

# Download 1000 genomes VCF files and PED file
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;

for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;


# Download GRCh37 reference genome
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;
gunzip human_g1k_v37.fasta.gz ;


# Convert the VCF files to PLINK format
for chr in {1..22}; do
    plink --vcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr --split-x b37 no-fail \
    --extract $HOME/Projects/DHR/data/Genetics/Illumina_VCF/DHR_ILLUMINA_PLUS_96IID_chr1-22_snps-only.snplist \
    --make-bed --out $HOME/Projects/DHR/out/SNP_PCA/ALL.chr"${chr}".phase3_1000G_intersect_DHR ;
done

# Merge into single file
ls $HOME/Projects/DHR/out/SNP_PCA/*.bim > $HOME/Projects/DHR/out/SNP_PCA/merge_list.txt ;
sed -i 's/.bim//g' $HOME/Projects/DHR/out/SNP_PCA/merge_list.txt ;
plink --merge-list $HOME/Projects/DHR/out/SNP_PCA/merge_list.txt --out $HOME/Projects/DHR/out/SNP_PCA/ALL_1000G_intersect_DHR ;
