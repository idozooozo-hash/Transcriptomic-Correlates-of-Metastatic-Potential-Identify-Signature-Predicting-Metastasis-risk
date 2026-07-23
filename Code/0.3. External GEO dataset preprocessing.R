########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"



########################################
## library
########################################
library(GEOquery)
library(GSVA)
library(dplyr)



########################################
## load data
########################################
load(paste0(dir_data,"2. Metastasis biomarker.Rdata")) # variable name: fdegs



########################################
## LUAD (GSE31210)
########################################

## get data
dir.create(paste0(dir_data,"GSE31210/"))
setwd(paste0(dir_data,"GSE31210/"))
gse = getGEO("GSE31210", GSEMatrix=TRUE, destdir=paste0(dir_data,"GSE31210/"), getGPL=T)
gse = gse[[1]]

em=exprs(gse)
sinfo=pData(phenoData(gse)) 


## annotation
adf=fread(file=paste0(dir_data,"1.14.15. GSE31210/GSE31210_family.xml/GPL570-tbl-1.txt"))
adf=adf[,c("V1","V11")]
colnames(adf)=c("probe","symbol")
rownames(em)=adf$symbol[match(rownames(em),adf$probe)]
rownames(em)=gsub("///.*","",rownames(em))
colnames(sinfo)[colnames(sinfo)=="geo_accession"]="Sample"


## processing duplicated gene
dugs=unique(rownames(em)[duplicated(rownames(em))])
for(gene in dugs) {
    print(which(dugs %in% gene))

    mexps=apply(em[rownames(em)==gene,] , 2 , mean )
    em=em[!(rownames(em) %in% gene),]
    em=rbind(em , t(as.data.frame(mexps)))
    rownames(em)[nrow(em)]=gene
}


## commona samples
em=em[,sinfo$Sample]


## Preprocessing
# (1) stage
colnames(sinfo)[colnames(sinfo)=="pathological stage:ch1"]="Pathologic Stage"
sinfo$`Pathologic Stage`=gsub("A$|B$","",sinfo$`Pathologic Stage`)
# (2) relapse
colnames(sinfo)[colnames(sinfo)=="relapse:ch1"]="Recurrence Occurrence"
sinfo$`Recurrence Occurrence`[sinfo$`Recurrence Occurrence`=="relapsed" & !is.na(sinfo$`Recurrence Occurrence`)]="Recurrence"
sinfo$`Recurrence Occurrence`[sinfo$`Recurrence Occurrence`=="not relapsed" & !is.na(sinfo$`Recurrence Occurrence`)]="No Recurrence"
colnames(sinfo)[colnames(sinfo)=="days before relapse/censor:ch1"]="Relapse Free Survival"
sinfo$`Relapse Free Survival`=as.numeric(sinfo$`Relapse Free Survival`)/365
# (3) Overall survival
colnames(sinfo)[colnames(sinfo)=="days before death/censor:ch1"]="Overall Survival"
sinfo$`Overall Survival`=as.numeric(sinfo$`Overall Survival`)/365
colnames(sinfo)[colnames(sinfo)=="death:ch1"]="Vital Status"
sinfo$`Vital Status`[sinfo$`Vital Status`=="dead" & !is.na(sinfo$`Vital Status`)]="Dead"
sinfo$`Vital Status`[sinfo$`Vital Status`=="alive" & !is.na(sinfo$`Vital Status`)]="Alive"
# (4) Age
colnames(sinfo)[colnames(sinfo)=="age (years):ch1"]="Age"
# (5) Gender
sinfo$`Sex`=toupper(sinfo$`gender:ch1`)


## Store
sinfol=list()
eml=list()
sinfol[["LUAD (GSE31210)"]]=sinfo
eml[["LUAD (GSE31210)"]]=em



########################################
## NSCLC (GSE50081)
########################################

## get data
dir.create(paste0(dir_data,"GSE50081/"))
setwd(paste0(dir_data,"GSE50081/"))
gse = getGEO("GSE50081", GSEMatrix=TRUE, destdir=paste0(dir_data,"1.14.21. GSE50081/"), getGPL=T)
gse = gse[[1]]

em=exprs(gse)
sinfo=pData(phenoData(gse)) 


## annotation
adf=fread(file=paste0(dir_data,"GSE50081/GSE50081_family.xml/GPL570-tbl-1.txt"))
adf=adf[,c("V1","V11")]
colnames(adf)=c("probe","symbol")
rownames(em)=adf$symbol[match(rownames(em),adf$probe)]
rownames(em)=gsub(" ///.* ","",rownames(em))
em=em[rownames(em)!="",]
colnames(sinfo)[colnames(sinfo)=="geo_accession"]="Sample"


## processing duplicated gene
dugs=unique(rownames(em)[duplicated(rownames(em))])
for(gene in dugs) {
    print(which(dugs %in% gene))

    mexps=apply(em[rownames(em)==gene,] , 2 , mean )
    em=em[!(rownames(em) %in% gene),]
    em=rbind(em , t(as.data.frame(mexps)))
    rownames(em)[nrow(em)]=gene
}


## common samples
em=em[,sinfo$Sample]


## Preprocessing
# (1) stage
colnames(sinfo)[colnames(sinfo)=="Stage:ch1"]="Stage"
sinfo$`Stage`=gsub("A$|B$","",sinfo$`Stage`)
stages=setNames(c("1","2","3","4"),c("I","II","III","IV"))
sinfo$Stage=names(stages)[match(sinfo$Stage,stages)]
# (2) relapse
colnames(sinfo)[colnames(sinfo)=="recurrence:ch1"]="Recurrence Occurrence"
sinfo$`Recurrence Occurrence`[sinfo$`Recurrence Occurrence`=="Y" & !is.na(sinfo$`Recurrence Occurrence`)]="Recurrence"
sinfo$`Recurrence Occurrence`[sinfo$`Recurrence Occurrence`=="N" & !is.na(sinfo$`Recurrence Occurrence`)]="No Recurrence"
colnames(sinfo)[colnames(sinfo)=="disease-free survival time:ch1"]="Relapse Free Survival"
sinfo$`Relapse Free Survival`=ifelse(sinfo$`Recurrence Occurrence`=="Recurrence",as.numeric(sinfo$`Relapse Free Survival`),apply(data.frame(as.numeric(sinfo$`Relapse Free Survival`),as.numeric(sinfo$`survival time:ch1`)),1,max))
# (3) Overall survival
colnames(sinfo)[colnames(sinfo)=="status:ch1"]="Vital Status"
sinfo$`Vital Status`=str_to_title(sinfo$`Vital Status`)
colnames(sinfo)[colnames(sinfo)=="survival time:ch1"]="Overall Survival"
sinfo$`Overall Survival`=ifelse(sinfo$`Vital Status`=="Dead",as.numeric(sinfo$`Overall Survival`),apply(data.frame(as.numeric(sinfo$`Overall Survival`),as.numeric(sinfo$`Relapse Free Survival`)),1,max))
# (4) Sex
sinfo$Sex=ifelse(sinfo$`Sex:ch1`=='F','FEMALE',ifelse(sinfo$`Sex:ch1`=='M','MALE',NA))
# (5) Age
sinfo$Age=as.numeric(sinfo$`age:ch1`)


## Store
sinfol[["NSCLC (GSE50081)"]]=sinfo
eml[["NSCLC (GSE50081)"]]=em



########################################
## NSCLC (GSE42127)
########################################

## get data
setwd(paste0(dir_data,"GSE42127/"))
gse = getGEO("GSE42127", GSEMatrix=TRUE, destdir=paste0(dir_data,"GSE42127/"), getGPL=T)
gse = gse[[1]]

em=exprs(gse)
sinfo=pData(phenoData(gse)) 


## annotation
adf=fread(file=paste0(dir_data,"GSE42127/GSE42127_family.xml/GPL6884-tbl-1.txt"))
adf=adf[,c("V1","V6")]
colnames(adf)=c("probe","symbol")
rownames(em)=paste0(adf$symbol[match(rownames(em),adf$probe)]," (",rownames(em),")")

colnames(sinfo)[colnames(sinfo)=="geo_accession"]="Sample"


## processing duplicated gene
rownames(em)=sapply(str_split(rownames(em)," \\(| /// "), function(g) g[1])
dugs=unique(rownames(em)[duplicated(rownames(em))])

for(gene in dugs) {
    print(which(dugs %in% gene))

    mexps=apply(em[rownames(em)==gene,] , 2 , mean )
    em=em[!(rownames(em) %in% gene),]
    em=rbind(em , t(as.data.frame(mexps)))
    rownames(em)[nrow(em)]=gene
}


## Preprocessing
# (1) Survival
sinfo$`Vital Status`=ifelse(sinfo$`survival status:ch1`=='D','Dead',ifelse(sinfo$`survival status:ch1`=='A','Alive',NA))
sinfo$`Overall Survival`=as.numeric(sinfo$`overall survival months:ch1`)/12 # month to year
# (2) Age
sinfo$Age=as.numeric(sinfo$`age at surgery:ch1`)
# (3) Sex
sinfo$Sex=ifelse(sinfo$`gender:ch1`=='F','FEMALE',ifelse(sinfo$`gender:ch1`=='M','MALE',NA))
# (4) Stage
sinfo$`Pathologic Stage`=gsub("A$|B$","",sinfo$`final.pat.stage:ch1`)
sinfo$`Pathologic Stage`[sinfo$`Pathologic Stage`=='unknown']=NA


## Save
sinfol[["NSCLC (GSE42127)"]]=sinfo
eml[["NSCLC (GSE42127)"]]=em



########################################
## NSCLC (GSE37745)
########################################

dir.create(paste0(dir_data,"GSE37745/"))

## get data
setwd(paste0(dir_data,"GSE37745/"))
gse = getGEO("GSE37745", GSEMatrix=TRUE, destdir=paste0(dir_data,"GSE37745/"), getGPL=T)
gse = gse[[1]]

em=exprs(gse)
sinfo=pData(phenoData(gse)) 


## annotation
adf=fread(file=paste0(dir_data,"GSE37745/GSE37745_family.xml/GPL570-tbl-1.txt"))
adf=adf[,c("V1","V11")]
colnames(adf)=c("probe","symbol")
rownames(em)=adf$symbol[match(rownames(em),adf$probe)]
rownames(em)=gsub(" ///.* ","",rownames(em))
em=em[rownames(em)!="",]

colnames(sinfo)[colnames(sinfo)=="geo_accession"]="Sample"


## processing duplicated gene
dugs=unique(rownames(em)[duplicated(rownames(em))])
for(gene in dugs) {
    print(which(dugs %in% gene))
    mexps=apply(em[rownames(em)==gene,] , 2 , mean )
    em=em[!(rownames(em) %in% gene),]
    em=rbind(em , t(as.data.frame(mexps)))
    rownames(em)[nrow(em)]=gene
}


## common samples
em=em[,sinfo$Sample]


## Preprocessing
# (1) drug
colnames(sinfo)[colnames(sinfo)=="adjuvant treatment:ch1"]="Treatment Status"
sinfo$`Treatment Status`[sinfo$`Treatment Status`=="yes"]="Treatment"
sinfo$`Treatment Status`[sinfo$`Treatment Status`=="no"]="No Treatment"
# (2) stage
colnames(sinfo)[colnames(sinfo)=="tumor stage:ch1"]="Stage"
sinfo$`Stage`=gsub("a$|b$","",sinfo$`Stage`)
stages=setNames(c("1","2","3","4"),c("I","II","III","IV"))
sinfo$Stage=names(stages)[match(sinfo$Stage,stages)]
# (3) relapse
colnames(sinfo)[colnames(sinfo)=="recurrence:ch1"]="Recurrence Occurrence"
sinfo$`Recurrence Occurrence`[sinfo$`Recurrence Occurrence`=="yes" & !is.na(sinfo$`Recurrence Occurrence`)]="Recurrence"
sinfo$`Recurrence Occurrence`[sinfo$`Recurrence Occurrence`=="no" & !is.na(sinfo$`Recurrence Occurrence`)]="No Recurrence"
colnames(sinfo)[colnames(sinfo)=="days to recurrence / to last visit:ch1"]="Relapse Free Survival"
sinfo$`Relapse Free Survival`=as.numeric(sinfo$`Relapse Free Survival`)/365
# (4) Overall survival
colnames(sinfo)[colnames(sinfo)=="dead:ch1"]="Vital Status"
sinfo$`Vital Status`[sinfo$`Vital Status`=="yes" & !is.na(sinfo$`Vital Status`)]="Dead"
sinfo$`Vital Status`[sinfo$`Vital Status`=="no" & !is.na(sinfo$`Vital Status`)]="Alive"
colnames(sinfo)[colnames(sinfo)=="days to determined death status:ch1"]="Overall Survival"
sinfo$`Overall Survival`=ifelse(sinfo$`Vital Status`=="Dead",as.numeric(sinfo$`Overall Survival`)/365,sinfo$`Relapse Free Survival`)
# (5) Sex
sinfo$Sex=toupper(sinfo$`gender:ch1`)
# (6) Age
sinfo$Age=as.numeric(sinfo$`age:ch1`)


## Save
sinfol[["NSCLC (GSE37745)"]]=sinfo
eml[["NSCLC (GSE37745)"]]=em



########################################
## Metastasis Risk calculation
########################################

## Run ssGSEA
sinfol=lapply(names(eml), function(set) {
    em=eml[[set]]
    sinfo=sinfol[[set]]

    param=ssgseaParam(em, list('Metastasis Risk'=fdegs))
    mdf=cbind(data.frame(Sample=colnames(em)),t(as.data.frame(gsva(param,verbose=T))))
    sinfo=left_join(sinfo,mdf)
    sinfo
})
names(sinfol)=names(eml)



########################################
## Save
########################################
# save(eml, sinfol, file=paste0(dir_data,"0.3. external dataset(GEO).Rdata"))