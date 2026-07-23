########################################
## Directory
########################################
dir_tcga="G:/My Drive/DB/TCGA/mRNA_protein_coding/" # The location where the TCGA(obtained from GDC) data is stored
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"



########################################
## Library
########################################
library(dplyr)
library(stringr)



########################################
## Load Data
########################################

## sample information
load(paste0(dir_tcga,"clinical/pan_cancer.Rdata")) # sinfo
gsinfo=sinfo

## Read count
load(paste0(dir_tcga,"read_count/pan_cancer.Rdata")) # rcm

## TPM
load(paste0(dir_tcga,"tpm/pan_cancer.Rdata")) # tpm



########################################
## Dataset Preview
########################################
tpm[1:5,1:5]
#          TCGA-RW-A680-01A-11R-A35K-07 TCGA-S7-A7WO-01A-11R-A35L-07
# TSPAN6                      2.3522218                    4.5646769
# TNMD                        0.2107629                    0.1826923
# DPM1                        4.3254008                    6.5985126
# SCYL3                       1.2570710                    2.4935192
# C1orf112                    0.3207735                    0.8487982
#          TCGA-W2-A7HH-01A-11R-A35L-07 TCGA-SR-A6MX-06A-11R-A35L-07
# TSPAN6                      2.9797144                    2.6657566
# TNMD                        0.8020691                    0.1256511
# DPM1                        6.2033793                    6.1109040
# SCYL3                       2.8687460                    2.8512393
# C1orf112                    1.3679860                    1.2668765
#          TCGA-TT-A6YO-01A-11R-A35L-07
# TSPAN6                      4.4141761
# TNMD                        0.0000000
# DPM1                        7.2540649
# SCYL3                       1.0533199
# C1orf112                    0.3523066

rcm[1:5,1:5]
#          TCGA-RW-A680-01A-11R-A35K-07 TCGA-S7-A7WO-01A-11R-A35L-07
# TSPAN6                            401                         1548
# TNMD                                5                            3
# DPM1                              495                         1743
# SCYL3                             206                          480
# C1orf112                           32                           72
#          TCGA-W2-A7HH-01A-11R-A35L-07 TCGA-SR-A6MX-06A-11R-A35L-07
# TSPAN6                            427                          361
# TNMD                               15                            2
# DPM1                             1199                         1224
# SCYL3                             593                          637
# C1orf112                          129                          125
#          TCGA-TT-A6YO-01A-11R-A35L-07
# TSPAN6                           1619
# TNMD                                0
# DPM1                             3215
# SCYL3                             130
# C1orf112                           29



########################################
## Sample information data
########################################

## only Lung & Tumor
nid=gsinfo$barcode[grep("Normal|New",gsinfo$sample_type)]
gsinfo=gsinfo[!(gsinfo$barcode %in% nid) , ]
sinfo=setNames(gsinfo[,c("patient","barcode","name")] , c("Sample","Tissue","Cancer Type"))


## Age
sinfo$Age=gsinfo$age_at_index


## recurrence
rdf=summarise(group_by(gsinfo,patient),
        `Recurrence Ocurrence`=ifelse(sum(paper_local.recurrence=="Locoregional Recurrence",na.rm=T)>0 , "Locoregional Recurrence" ,
                ifelse(sum(paper_Recurrence=="NO",na.rm=T)+sum(paper_recurred_progressed=="NO",na.rm=T)+sum(paper_distant.recurrence!="Distant Metastasis",na.rm=T)+sum(!(paper_distant_metastasis_pathologic_spread %in% c("M1","M1a")),na.rm=T) + sum(paper_Pathologic.Spread.Distant.Metastasis!="M1",na.rm=T) + sum(paper_Clinical.Spread..Distant.Metastases..M.!="M1",na.rm=T) + sum(paper_AJCC.metastasis.category!="M1",na.rm=T) + sum(paper_ajcc_metastasis_pathologic_pm!="M1",na.rm=T)==0,"No Recurrence",NA)),
        `Relapse Free Survival`=ifelse(sum(paper_local.recurrence=="Locoregional Recurrence",na.rm=T)>0, min(as.numeric(paper_local.recurrence.days),na.rm=T),max(days_to_last_follow_up,na.rm=T)))

sinfo$`Recurrence Ocurrence`=rdf$`Recurrence Ocurrence`[match(sinfo$Sample,rdf$patient)]
sinfo$`Relapse Free Survival`=rdf$`Relapse Free Survival`[match(sinfo$Sample,rdf$patient)]
sinfo$`Relapse Free Survival`=ifelse(!is.finite(sinfo$`Relapse Free Survival`),NA,sinfo$`Relapse Free Survival`)


# metastasis
mdf=summarise(group_by(gsinfo,patient),
        `Metastasis Ocurrence`=ifelse(sum(paper_distant.recurrence=="Distant Metastasis",na.rm=T)+sum(paper_Pathologic.Spread.Distant.Metastasis=="M1")+sum(paper_AJCC.metastasis.category=="M1",na.rm=T)+sum(paper_Clinical.Spread..Distant.Metastases..M.=="M1",na.rm=T) + sum(paper_ajcc_metastasis_pathologic_pm=="M1",na.rm=T)>0 , "Distant Metastasis" ,
                ifelse(sum(paper_Recurrence=="NO",na.rm=T)+sum(paper_recurred_progressed=="NO",na.rm=T)+sum(paper_local.recurrence!="Locoregional Recurrence",na.rm=T)==0,"No Metastasis",NA)),
        `Metastasis Free Survival`=ifelse(sum(paper_distant.recurrence=="Distant Metastasis",na.rm=T)>0, min(as.numeric(paper_distant.recurrence.days),na.rm=T),max(days_to_last_follow_up,na.rm=T)))

sinfo$`Metastasis Ocurrence`=mdf$`Metastasis Ocurrence`[match(sinfo$Sample,mdf$patient)]
sinfo$`Metastasis Free Survival`=mdf$`Metastasis Free Survival`[match(sinfo$Sample,mdf$patient)]
sinfo$`Metastasis Free Survival`=ifelse(!is.finite(sinfo$`Metastasis Free Survival`),NA,sinfo$`Metastasis Free Survival`)

msdf=summarise(group_by(gsinfo, patient), `Metastasis Site`=paste(unique(paper_distant.recurrence.site[!is.na(paper_distant.recurrence.site)]),collapse="&"))
sinfo$`Metastasis Site`=msdf$`Metastasis Site`[match(sinfo$Sample,msdf$patient)]


## survival yes/no
sinfo$`Vital Status`=str_to_title(gsinfo$vital_status)
sdays=gsinfo$days_to_last_follow_up
ddays=gsinfo$days_to_death
sinfo$`Overall Survival`=as.numeric(ifelse(sinfo$`Vital Status`=="Dead", ddays , ifelse(sinfo$`Vital Status`=="Alive" , sdays , NA)))


## stage
sinfo$`AJCC Stage`=gsub("Stage |A$|B$|C$","",gsinfo$ajcc_pathologic_stage)


## TNM stage
sinfo$`TNM-M Stage`=ifelse(as.character(gsinfo$paper_M.stage) %in% c("M0","M1"),as.character(gsinfo$paper_M.stage),NA)
gsinfo$paper_T.stage=gsub("a$|b$","",gsinfo$paper_T.stage)
sinfo$`TNM-T Stage`=ifelse(as.character(gsinfo$paper_T.stage) %in% c("T0","T1","T3","T4"),as.character(gsinfo$paper_T.stage),NA)

## Race
sinfo$Race=gsinfo$race
sinfo$Race[sinfo$Race=="not reported"]=NA


## Gender
sinfo$Gender=str_to_title(gsinfo$gender)


## first collected sample
dim(sinfo) # 10524 16
sum(duplicated(sinfo$Sample)) # 162

sdf=summarise(group_by(sinfo,Sample),
            "Vital Status"=ifelse(sum(`Vital Status`=="Dead")>0,"Dead","Alive") ,
            "Overall Survival"=ifelse(sum(`Vital Status`=="Dead")>0,min(`Overall Survival`[`Vital Status`=="Dead"]),max(`Overall Survival`[`Vital Status`=="Alive"])),
            "AJCC Stage"=ifelse(sum(!is.na(`AJCC Stage`))==0,NA,sort(`AJCC Stage`[!is.na(`AJCC Stage`)],decreasing=T)[1]),
            "TNM-M Stage"=ifelse(sum(!is.na(`TNM-M Stage`))==0,NA,sort(`TNM-M Stage`[!is.na(`TNM-M Stage`)],decreasing=T)[1]),
            "TNM-T Stage"=ifelse(sum(!is.na(`TNM-T Stage`))==0,NA,sort(`TNM-T Stage`[!is.na(`TNM-T Stage`)],decreasing=T)[1])) 

dsams=unique(sinfo$Sample[duplicated(sinfo$Sample)])
for(sam in dsams) { 
    esinfo=gsinfo[gsinfo$patient==sam,]
    if(sum(!is.na(esinfo$days_to_collection))==0) {bcode=sort(esinfo$barcode)[1]} else {bcode=sort(setdiff(esinfo$barcode[esinfo$days_to_collection==min(esinfo$days_to_collection,na.rm=T)],NA))[1]}
    sinfo=rbind(sinfo[!(sinfo$Sample %in% sam),] , sinfo[sinfo$Tissue==bcode,])
}

dim(sinfo) # 10362 16
sum(duplicated(sinfo$Sample)) # 0

sinfo=left_join(sinfo[,setdiff(colnames(sinfo),setdiff(colnames(sdf),c("Sample")))] , sdf)



########################################
## expression data
########################################

## read count
rcm=rcm[,match(sinfo$Tissue,colnames(rcm))]
colnames(rcm)=sinfo$Sample
dugs=unique(rownames(rcm)[duplicated(rownames(rcm))])
for(dug in dugs) {
    ercm=rcm[rownames(rcm)==dug,]    
    rcm=rbind(rcm[rownames(rcm)!=dug,] , matrix(apply(ercm,2,mean), nrow=1, dimnames=list(dug , colnames(ercm))))    
}


## tpm
tpm=tpm[,match(sinfo$Tissue,colnames(tpm))]
colnames(tpm)=sinfo$Sample
dugs=unique(rownames(tpm)[duplicated(rownames(tpm))])
for(dug in dugs) {
    etpm=tpm[rownames(tpm)==dug,]    
    tpm=rbind(tpm[rownames(tpm)!=dug,] , matrix(apply(etpm,2,mean), nrow=1, dimnames=list(dug , colnames(etpm))))    
}



########################################
## Save
########################################

## common sample
csam=Reduce(intersect , list(sinfo$Sample,colnames(rcm),colnames(tpm)))
sinfo=sinfo[match(csam,sinfo$Sample),]
rcm=rcm[,match(csam,colnames(rcm))]
tpm=tpm[,match(csam,colnames(tpm))]


## save
# save(sinfo, file=paste0(dir_data,"0.2. sample information(TCGA-GDC).Rdata"))
# save(rcm, file=paste0(dir_data,"0.2. read count(TCGA-GDC).Rdata"))
# save(tpm, file=paste0(dir_data,"0.2. tpm(TCGA-GDC).Rdata"))