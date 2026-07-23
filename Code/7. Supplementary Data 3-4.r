########################################
## Supplementary Data 3
########################################
dl=list()



dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"

dir_fh="G:/My Drive/DB/Firehose/LUAD/" # The location where the TCGA(obtained from Firehose) data is stored
dir_tcga="G:/My Drive/DB/TCGA/mRNA_protein_coding/" # TCGA (obtained from GDC) dataset directory
dir_cmap="G:/My Drive/DB/CMap/" # Dataset's directory The dataset directory was obtained from Connectivity Map (CMap) Genomics.



## MetMap Sample for MP EDA
load(file=paste0(dir_data,"0. all sample information(DepMap).Rdata")) # variable name: sinfo
dl=c(dl , list(`MetMap EDA cell line`=sort(sinfo$Sample)))


## MetMap Sample for Metastasis Siganture extraction
dl[["MetMap signature cell line"]]=sort(sinfo$Sample[sinfo$`Cancer Type`=='Lung' & sinfo$`meta type to all5`!='Weakly Metastatic (Low Confidence)'])


## TCGA patient samples collected from Firehose for validation analysis
load(paste0(dir_data,"0.1. sample information(TCGA-Firehose).Rdata"))
psinfo=sinfo
load(paste0(dir_fh,"sinfo_merged.Rdata"))
psinfo$`Neoadjuvant Treatment`=sinfo$`patient.history_of_neoadjuvant_treatment`[match(psinfo$Sample,sinfo$patient.bcr_patient_barcode)]
sinfo=psinfo

cols=c('Sample','Age','Metastasis Ocurrence','Metastasis Free Survival','Recurrence Ocurrence','Relapse Free Survival','Vital Status','Overall Survival','Pathologic Stage','Race','Gender','Tissue')
sinfo=sinfo[,cols]
sinfo=sinfo[apply(sinfo, 1, function(v) sum(is.na(v))!=(length(cols)-1)),]
dl[['TCGA-Firehose validation samples']]=sort(sinfo$Sample)


## TCGA tissue samples collected from GDC for matched-normal DEG test
load(paste0(dir_tcga,"clinical/LUAD.Rdata")) # sample information; variable name: sinfo
load(paste0(dir_tcga,"tpm/LUAD.Rdata")) # TPM expression; variable name: tpm
mids=intersect(sinfo$patient[grep("Tumor",sinfo$definition)] , sinfo$patient[grep("Normal",sinfo$definition)])
csinfo=sinfo[sinfo$patient %in% mids & grepl("Normal",sinfo$definition),]
cid=sort(csinfo[!duplicated(csinfo$patient),"barcode"])
tsinfo=sinfo[sinfo$patient %in% mids & grepl("Tumor",sinfo$definition),]
tid=sort(tsinfo[!duplicated(tsinfo$patient),"barcode"])

dl[['TCGA-GDC matched-normal samples']]=sort(c(cid,tid))


## DepMap cell-line for Drug sensitivity analysis
load(paste0(dir_data,"5. Predicted Metastasis risk score with random gene set.Rdata"))
dl[['DepMap drug-sensitivity cell lines']]=sort(rownames(mdf))


## CMap Signature ID for Drug repurposing
load(paste0(dir_cmap,"level5(gold)_compound_perturbation(2021-11-23).Rdata"))
load(paste0(dir_data,"6. 1st Drug repurposing res (based drug signature).Rdata")) # cdf

sinfo=sinfo[sinfo$Pooling=="tumor",]
tinfo=tinfo[tinfo$CC_q75>=0.8 & tinfo$Replicate_self_rank_q25<=0.05 & !is.na(tinfo$CC_q75) & !is.na(tinfo$Replicate_self_rank_q25) & tinfo$Replicate_no>=3 & tinfo$QC_pass>=1,] # remove low quality; Replicate_self_rank_q25=pct_self_rank_q25 
tinfo=tinfo[tinfo$Perturbation_type=="trt_cp",] # only perturbation
tinfo=tinfo[tinfo$Cell_Line %in% sinfo$Cell_Line,]

pm=pm[,colnames(pm) %in% tinfo$Signature_id]

dl[['CMap drug-repurposing signatures']]=sort(intersect(tinfo$Signature_id[tinfo$Perturbagen_cmap %in% cdf$drug] , colnames(pm)))


## Supplementary Data 3
gmt = vapply(names(dl), function(x) { paste(c(x, "NA", dl[[x]]), collapse = "\t")}, FUN.VALUE = character(1))
# writeLines(gmt, con = "G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/Supplementary Data/Supplementary Data 3.gmt")



########################################
## Supplementary Data 4
########################################
dir_dep='G:/My Drive/DB/DepMap/'

## DepMap cell-line for Drug sensitivity analysis
load(paste0(dir_dep,"sinfo_24Q4.Rdata")) # sinfo
load(paste0(dir_data,"5. Predicted Metastasis risk score with random gene set.Rdata")) # mdf

sdf=data.frame(`Cell line ID`=rownames(mdf) , `Tumor Lineage`=sinfo$OncotreeLineage[match(rownames(mdf) , sinfo$ModelID)] , `Primary Disease`=sinfo$OncotreePrimaryDisease[match(rownames(mdf) , sinfo$ModelID)] , `Disease Subtype`=sinfo$OncotreeSubtype[match(rownames(mdf) , sinfo$ModelID)], check.names=F)
sdf=sdf[order(sdf$`Cell line ID`),]
# write.csv(sdf, file="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/Supplementary Data/Supplementary Data 4.csv", row.names=F, quote=F)