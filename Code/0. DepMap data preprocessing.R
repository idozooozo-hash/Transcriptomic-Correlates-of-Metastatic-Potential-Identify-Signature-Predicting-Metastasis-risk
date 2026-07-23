########################################
## Directory
########################################
dir_dep="G:/My Drive/DB/DepMap/" # The location where the DepMap data is stored
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"



########################################
## Library
########################################
library(stringr)
library(data.table)
library(dplyr)



########################################
## Load Data
########################################
load(paste0(dir_dep,"231212_metmap500.Rdata")) # load(paste0(dir_dep,"231212_metmap500.Rdata")) # R variable name: mts; source: DepMap MetMap_500_Metastatic_Potential.csv, 2020 release
load(paste0(dir_dep,"tpm_24Q4.Rdata")) # R variable name: tpm; source: DepMap OmicsExpressionProteinCodingGenesTPMLogp1.csv, 24Q4 release
load(paste0(dir_dep,"read_count_24Q4.Rdata")) # R variable name: rcm; source: DepMap OmicsExpressionGenesExpectedCountProfile.csv, 24Q4 release
load(paste0(dir_dep,"sinfo_24Q4.Rdata")) # R variable name: sinfo; source: DepMap 'Model.csv', 'OmicsDefaultModelProfiles.csv', 'OmicsDefaultModelConditionProfiles.csv' and 'ModelCondition.csv', 24Q4 release



########################################
## Data preview
########################################
head(mts)
#                  all5       bone     brain      kidney      liver         lung
# ACH-000824 -0.3990943 -0.4150925 -3.588466 -2.70126910 -2.8143728 -1.960035800
# ACH-001307  0.7004353  0.3598184 -3.087720  0.34062310 -1.0269101 -0.354831550
# ACH-000756 -3.6238341 -3.6238341 -3.623834 -3.62383410 -3.6238341 -3.623834100
# ACH-000681  0.7230241  0.4708827 -1.318449  0.06815637 -0.2603131 -0.251146500
# ACH-000444  1.4668136 -0.3881296  1.144734 -0.32549992  1.1296722 -0.008543288
# ACH-000142 -3.9475245 -3.9475245 -3.947525 -3.94752450 -3.9475245 -3.947524500


tpm[1:5,1:5]
#          ACH-001113 ACH-001289 ACH-001339 ACH-001979 ACH-002438
# TSPAN6     4.331992  4.5674238   3.150560   4.240314   4.032101
# TNMD       0.000000  0.5849625   0.000000   0.000000   0.000000
# DPM1       7.364660  7.1066415   7.379118   5.681168   6.674687
# SCYL3      2.792855  2.5434959   2.333424   2.063503   2.117695
# C1orf112   4.471187  3.5046204   4.228049   1.641546   3.003602


rcm[1:5,1:5]
#          ACH-001113 ACH-001289 ACH-001339 ACH-001979 ACH-002438
# TSPAN6      2383.00    2529.00    1552.00    3032.00    2616.00
# TNMD           0.00      13.00       0.00       0.00       0.00
# DPM1        5332.80    3978.00    8303.00    2157.00    4012.90
# SCYL3        961.53     555.83     849.46     727.48     618.33
# C1orf112    1518.50     832.17    2511.50     252.52     798.67


head(sinfo[,c("ProfileID","ModelID","Datatype","PatientID","CellLineName","OncotreeLineage","OncotreeSubtype","OncotreePrimaryDisease","Age","Sex","PrimaryOrMetastasis","SampleCollectionSite")])
#   ProfileID    ModelID Datatype PatientID CellLineName      OncotreeLineage
# 1 PR-01r7OM ACH-000957      rna PT-FXUrcz       LS 180                Bowel
# 2 PR-02XmLG ACH-002785      rna PT-6sPicj  NCC-LMS1-C1          Soft Tissue
# 3 PR-045poV ACH-003273      rna PT-44yhk7  CHRF-288-11              Myeloid
# 4 PR-04VvBz ACH-001289      wes PT-773uN4   COG-AR-359            CNS/Brain
# 5 PR-08baLC ACH-000237      wgs PT-wwKJYr       JHOM-1 Ovary/Fallopian Tube
# 6 PR-09gmEI ACH-000520      rna PT-tv8Ku2          59M Ovary/Fallopian Tube
#                    OncotreeSubtype    OncotreePrimaryDisease Age    Sex
# 1             Colon Adenocarcinoma Colorectal Adenocarcinoma  58 Female
# 2                   Leiomyosarcoma            Leiomyosarcoma  42 Female
# 3  Acute Megakaryoblastic Leukemia    Acute Myeloid Leukemia   1   Male
# 4 Atypical Teratoid/Rhabdoid Tumor           Embryonal Tumor  NA   Male
# 5          Mucinous Ovarian Cancer  Ovarian Epithelial Tumor  88 Female
# 6 High-Grade Serous Ovarian Cancer  Ovarian Epithelial Tumor  65 Female
#   PrimaryOrMetastasis   SampleCollectionSite
# 1             Primary                  Colon
# 2          Metastatic                   bone
# 3             Primary                       
# 4                     central_nervous_system
# 5             Primary                  ovary
# 6          Metastatic                ascites



########################################
## sample information Preprocessing
########################################
sinfo=sinfo[,c("ModelID","PatientID","CellLineName","OncotreeLineage","OncotreeSubtype","OncotreePrimaryDisease","PatientSubtypeFeatures","Age","AgeCategory","Sex","PatientRace","PrimaryOrMetastasis","SampleCollectionSite","Stage")]
colnames(sinfo)=c("Sample","Patient","Cell Line","Cancer Type","Cancer Subtype","Primary Disease","Subtype Features","Age","Age Category","Sex","Race","Cell Line Origin","Sample Collection Site","Stage")

sinfo$`Cancer Type`=gsub("/"," or ",sinfo$`Cancer Type`)
sinfo$`Cancer Subtype`=gsub("/"," or ",sinfo$`Cancer Subtype`)
sinfo$`Sample Collection Site`=str_to_title(gsub("_"," ",sinfo$`Sample Collection Site`))
sinfo$Race=str_to_title(sinfo$Race)



########################################
## MetMap500 Preprocessing
########################################
mts=cbind(data.frame(Sample=rownames(mts)),setNames(mts,paste0("meta to ",colnames(mts))))
mts_binary=mts[,grep("meta",colnames(mts))]
mts_binary=ifelse(mts_binary<=(-4),"Non Metastatic",ifelse(mts_binary>=(-2),"Metastatic","Weakly Metastatic (Low Confidence)"))
colnames(mts_binary)=gsub("meta to ","meta type to ",colnames(mts_binary))

mts=cbind(mts , mts_binary)
sinfo=left_join(mts,sinfo)
sinfo=distinct(sinfo)



########################################
## select common samples
########################################
csams=Reduce(intersect , list(sinfo$Sample,mts$Sample,colnames(tpm),colnames(rcm))) # common sample ID

tpm=tpm[,csams]
rcm=rcm[match(rownames(tpm),rownames(rcm)),csams]

dim(tpm) # 19193 481
dim(rcm) # 19193 481



########################################
## Save
########################################
# save(tpm, file=paste0(dir_data,"0. tpm(DepMap).Rdata"))
# save(rcm, file=paste0(dir_data,"0. read count(DepMap).Rdata"))
# save(sinfo, file=paste0(dir_data,"0. all sample information(DepMap).Rdata"))

sinfo=sinfo[match(csams,sinfo$Sample),]
# save(sinfo, file=paste0(dir_data,"0. filtered sample information(DepMap).Rdata"))
