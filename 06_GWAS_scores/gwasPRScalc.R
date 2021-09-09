### Author: Tojo James, Karolinska Institutet
#Unix CMD: /usr/bin/R --slave --no-restore --vanilla --file=gwasPRScalc.R --args GCST004988

setwd('/home/tojo.james/Genetics/Imputation/')
library(data.table)
library(rbgen)
library(BSgenome)
library(rtracklayer)
library(BiocInstaller)
#biocLite("BSgenome.Hsapiens.UCSC.hg19")

load("GWAS_all_diseaseTrait.RData")


args <- commandArgs(TRUE)
#Study_Accession <- "GCST004988"
Study_Accession <- as.character(args[1])

Disease_PRS_chrom_sample<- list()
Disease_PRS_chrom_sample[["PRS_file_version"]] <-Study_Accession

diseaselabel <- as.character(SelectedGWASTraitSNP$TraitSNPCount[SelectedGWASTraitSNP$TraitSNPCount$Accession == Study_Accession,][["Trait"]])
diseaselabel <- paste(Study_Accession,diseaselabel,sep = " ") 
diseaselabel <- gsub(" ", "_", diseaselabel, fixed = TRUE)

PRS_chr_merged <- list()
for (i in 1:22){
  chrN <- i
  if (chrN < 10){
    chr0N <-chrN
  } else{
    chr0N <-chrN 
  }
  
  
  

  inputFileName <-paste0("DHR_1000G_SISu_imputed_chr",chrN,"_v.bgen")
  getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19",chr=chr){
    g <- getBSgenome(genome, masked=FALSE)
    seqlengths(g)[as.numeric(chr)]
  }
  
  
  chrstart = 1
  chrend = getChrLength(chr=chrN)
  print(c(chrstart,chrend))
  
  
  
  
  ranges = data.frame(chromosome = as.character(chr0N), start = as.numeric(chrstart), end = as.numeric(chrend))
  geneticdata = bgen.load(inputFileName, ranges)
  print(paste0("Read completed for genotype file ",inputFileName))
  print(length(row.names(geneticdata$variants)))
  

  

  selected_chrGWAS <-SelectedGWASTraitSNP[[Study_Accession]][SelectedGWASTraitSNP[[Study_Accession]]$seqnames == paste0("chr",chrN),]
  #selected_variant_alleles <-subset(geneticdata$variants, chromosome == chr0N & position == as.numeric(selected_chrGWAS$start))
  geneticvariant_PRS <-merge(geneticdata$variants,selected_chrGWAS ,by.x = "position",by.y = "start",all.y = TRUE)
  
  head(geneticvariant_PRS)
  missinggenotypecount <-nrow(geneticvariant_PRS[is.na(geneticvariant_PRS$rsid),])
  geneticvariant_PRS <- unique(geneticvariant_PRS[!is.na(geneticvariant_PRS$rsid),][c('chromosome','position','rsid','allele0','allele1','risksnp','OR_BETA','X95_Text')])
  geneticvariant_PRS <- geneticvariant_PRS[!duplicated(geneticvariant_PRS$position),]
  matchedgenotypecount <-nrow(geneticvariant_PRS)
  
  increase_betavar <-geneticvariant_PRS[grepl("increase",geneticvariant_PRS[,"X95_Text"]),]
  
  decrease_betavar <-geneticvariant_PRS[grepl("decrease",geneticvariant_PRS[,"X95_Text"]),]
  decrease_betavar$OR_BETA <- as.numeric(decrease_betavar$OR_BETA)*(-1)
  
  odds_var <-geneticvariant_PRS[!grepl(c("crease"),geneticvariant_PRS[,"X95_Text"]),]
  odds_var$OR_BETA <-log(as.numeric(odds_var$OR_BETA))
  
  geneticvariant_PRS <- rbind(increase_betavar,decrease_betavar,odds_var)
  
  
  
  geneticvariant_PRS$effect_allele <- sapply(strsplit(as.character(geneticvariant_PRS$risksnp), "-", fixed=T), "[", 2)
  #geneticvariant_PRS$rsid <- sapply(strsplit(as.character(geneticvariant_PRS$risksnp), "-", fixed=T), "[", 1)
  
  
  geneticvariant_PRS_A1 <- geneticvariant_PRS[as.character(geneticvariant_PRS$allele0) == as.character(geneticvariant_PRS$effect_allele),]
  geneticvariant_PRS_A2 <- geneticvariant_PRS[as.character(geneticvariant_PRS$allele1) == as.character(geneticvariant_PRS$effect_allele),]
  dim(geneticvariant_PRS_A1);dim(geneticvariant_PRS_A2)
  
  geneticdata_g0 <-data.frame(geneticdata$data[as.character(geneticvariant_PRS_A1$rsid),,"g=0"])
  geneticdata_g0_1 <-data.frame(geneticdata$data[as.character(geneticvariant_PRS_A1$rsid),,"g=1"])
  geneticdata_g2 <-data.frame(geneticdata$data[as.character(geneticvariant_PRS_A2$rsid),,"g=2"])
  geneticdata_g2_1 <-data.frame(geneticdata$data[as.character(geneticvariant_PRS_A2$rsid),,"g=1"])
  colnames(geneticdata_g0)<-gsub("X.anonymous_|\\.", "",colnames(geneticdata_g0))
  colnames(geneticdata_g0_1)<-gsub("X.anonymous_|\\.", "",colnames(geneticdata_g0_1))
  colnames(geneticdata_g2)<-gsub("X.anonymous_|\\.", "",colnames(geneticdata_g2))
  colnames(geneticdata_g2_1)<-gsub("X.anonymous_|\\.", "",colnames(geneticdata_g2_1))
  geneticdata_g0$rsid <- row.names(geneticdata_g0)
  geneticdata_g0_1$rsid <- row.names(geneticdata_g0_1)
  geneticdata_g2$rsid <- row.names(geneticdata_g2)
  geneticdata_g2_1$rsid <- row.names(geneticdata_g2_1)
  
  
  c(as.character(geneticvariant_PRS_A1$rsid),geneticdata_g0$rsid)
  
  
  
  genticvariant_A1_g0 <-merge(geneticvariant_PRS_A1,geneticdata_g0,by= "rsid")
  genticvariant_A1_g0_1 <-merge(geneticvariant_PRS_A1,geneticdata_g0_1,by = c("rsid"))
  
  genticvariant_A1_g2 <-merge(geneticvariant_PRS_A2,geneticdata_g2,by = c("rsid"))
  genticvariant_A1_g2_1 <-merge(geneticvariant_PRS_A2,geneticdata_g2_1,by = c("rsid"))
  
  dosage_g0 <- 2*genticvariant_A1_g0[-(1:9)] + genticvariant_A1_g0_1[-(1:9)]
  dosage_g2 <- 2*genticvariant_A1_g2[-(1:9)] + genticvariant_A1_g2_1[-(1:9)]
  
  SNP_PRS_g0 <- genticvariant_A1_g0[["OR_BETA"]]*dosage_g0
  SNP_PRS_g2 <- genticvariant_A1_g2[["OR_BETA"]]*dosage_g2
  
  sum_A1_g0<-colSums(SNP_PRS_g0)
  sum_A1_g2<-colSums(SNP_PRS_g2)
  
  PRS_chr <- sum_A1_g0 + sum_A1_g2
  
  PRS_chr_merged[[i]] <-data.frame(PRS_chr)
  colnames(PRS_chr_merged[[i]]) <-as.character(paste0("chr",chrN))
  
  Disease_PRS_chrom_sample[[paste0("PRS_Sample_chr", i)]]<-as.data.frame(PRS_chr)
  Disease_PRS_chrom_sample[[paste0("SampleGenotypeSize_chr", i)]]<-as.numeric(length(row.names(geneticdata$variants)))
  Disease_PRS_chrom_sample[[paste0("PRSVariantSize_chr", i)]]<-as.numeric(nrow(selected_chrGWAS))
  Disease_PRS_chrom_sample[[paste0("MatchingGenotypeCount_chr", i)]]<-as.numeric(matchedgenotypecount)
  Disease_PRS_chrom_sample[[paste0("MissingGenotypeCount_chr", i)]]<-as.numeric(missinggenotypecount)
  print(paste0("Completed chr_",i))
}
  
PRS_Disease_sample <- data.frame(PRS_chr_merged)
PRS_Disease_sample$PR_Score<-rowSums(PRS_Disease_sample)
Disease_PRS_chrom_sample[["Result"]] <-PRS_Disease_sample
write.csv(Disease_PRS_chrom_sample[["Result"]],file = paste0(diseaselabel,"_GWAS_PRS_dosage_filtered.csv"),row.names = TRUE)
rm(list=setdiff(ls(),c("Disease_PRS_chrom_sample","outdatafile","diseaselabel")))
save(Disease_PRS_chrom_sample, file=paste0(diseaselabel,"_GWAS_PRS_filtered.csv")) 
 
print("######### Results written into .RData ########")
print("######### COMPLETED ########")


