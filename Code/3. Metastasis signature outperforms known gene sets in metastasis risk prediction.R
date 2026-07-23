########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/figure/"

dir_fh="G:/My Drive/DB/Firehose/LUAD/" # The location where the TCGA(obtained from Firehose) data is stored
dir_tcga="G:/My Drive/DB/TCGA/mRNA_protein_coding/" # TCGA (obtained from GDC) dataset directory

## Additional gene sets directory
dir_msig="G:/My Drive/DB/MSigDB/homosapiens/" # MSigDB datasets directory
dir_csea="G:/My Drive/DB/cancerSEA/" # cancerSEA datasets directory
dir_ctd="G:/My Drive/DB/ctd-Comparative Toxicgenomics Database/" # CTD datasets directory
dir_creeds="G:/My Drive/DB/CREEDS/" # CREEDS dataset directory
dir_emtome="G:/My Drive/DB/EMTome/" # EMTome dataset directory



########################################
## Library
########################################
library(data.table)
library(stringr)
library(colorRamp2)
library(ComplexHeatmap)
library(GSVA)
library(survival)
library(survminer)
library(venndir)



########################################
## Load Data
########################################
load(paste0(dir_data,"2. Metastasis biomarker.Rdata")) # variable name: fdegs
load(paste0(dir_data,"2. DEG analysis res.Rdata")) # variable name: ddf

## TCGA data
tpml=c()
load(paste0(dir_data,"0.1. tpm(TCGA-Firehose).Rdata"))
tpml[["Firehose"]]=tpm
load(paste0(dir_data,"0.2. tpm(TCGA-GDC).Rdata"))
tpml[["GDC"]]=tpm

sinfol=c()
load(paste0(dir_data,"0.1. sample information(TCGA-Firehose).Rdata"))
sinfol[["Firehose"]]=sinfo
load(paste0(dir_data,"0.2. sample information(TCGA-GDC).Rdata"))
sinfol[["GDC"]]=sinfo



########################################
## Additional Gene sets
########################################

# (1) EMTome Gene sets
emtome=qusage::read.gmt(paste0(dir_emtome,"EMTome_signatures.gmt"))
names(emtome)=paste0("EMT (",names(emtome),")")
names(emtome)=gsub("_et_al."," et al. ",names(emtome))
names(emtome)=paste0(names(emtome)," | EMTome")
emtome=lapply(emtome, function(g) gsub("\\ $","",g))
# (2) CTD EMT-related Gene sets
cdf=fread(paste0(dir_ctd,"CTD_genes_diseases.csv/CTD_genes_diseases.csv"), check.names=F)
cdf=cdf[grep("Metastasis|Invasiveness",cdf$V3),]
colnames(cdf)[1]="Gene"
ctd=split(cdf$Gene,cdf$V3)
names(ctd)=paste0(names(ctd)," | CTD")
# (3) MSigDB-Hallmark 'EMT' Gene set
load(paste0(dir_msig,"hm_2024.1.Rdata")) # variable name: hm
hm=list(`EMT | Hallmark`=hm[["EPITHELIAL MESENCHYMAL TRANSITION"]])
# (4)  CancerSEA EMT-related Gene sets
list.files(dir_csea)
# [1] "Angiogenesis-LUAD-Guillaumet-Adkins A. Genome Biol. 2017 (PDX_LUAD).csv"
# [2] "Angiogenesis-LUAD-Kim KT. Genome Biol. 2015 (PDX).csv"                  
# [3] "EMT-LUAD-Guillaumet-Adkins A. Genome Biol. 2017 (PDX_LUAD).csv"         
# [4] "EMT-LUAD-Kim KT. Genome Biol. 2015 (PDX).csv"                           
# [5] "Invasion-LUAD-Guillaumet-Adkins A. Genome Biol. 2017 (PDX_LUAD).csv"    
# [6] "Invasion-LUAD-Kim KT. Genome Biol. 2015 (PDX).csv"                      
# [7] "Metastasis-LUAD-Kim KT. Genome Biol. 2015 (PDX).csv"  
fns=list.files(dir_csea) 
sdfl=lapply(fns, function(fn) read.csv(paste0(dir_csea,fn)))
names(sdfl)=gsub(".csv","",fns)
names(sdfl)=gsub("LUAD-","LUAD (",names(sdfl))
names(sdfl)=paste0(gsub("\\. [0-9].*","., ",names(sdfl)),regmatches(names(sdfl), regexpr("(?<=\\. )[0-9]+", names(sdfl), perl = TRUE)),")")
names(sdfl)=paste(names(sdfl),"|","cancerSEA")
csea=lapply(sdfl, function(sdf) sdf$Symbol[sdf$FDR==0 & sdf$Correlation>0])

## Total Gene set list
gsl=c(list(`Metastasis Signature`=fdegs),ctd,hm,csea,emtome) 
gsl=lapply(gsl, function(gs) unique(gs[gs %in% ddf$Gene]))
gsl=gsl[sapply(gsl,length)!=0]


length(gsl)-1 # additional gene set No.: 96 (except metastasis biomarker set)


## Save
# save(gsl, file=paste0(dir_data,"3. Additional Gene set list.Rdata"))



########################################
## Overlap between the metastasis signature and previously published gene sets
########################################

## Over-representation analysis (ORA)
odf=data.frame(Reference=NA, `Gene Sets`=names(gsl),pvalue=NA,padj=NA,OR=NA,Overlap=NA, check.names=F)
odf$Reference=gsub(".*\\| ","",odf$`Gene Sets`)
tgs=unique(ddf$Gene) # total protein-coding gene (length: 17080)

for(set in names(gsl)) {
    sgs=gsl[[set]]
    
    a=length(unique(intersect(sgs,fdegs))) # both
    b=length(setdiff(sgs,fdegs)) # only sgs
    c=length(setdiff(fdegs,sgs)) # only fdegs
    d=length(setdiff(tgs,c(fdegs,sgs))) # not any

    mtx = matrix(c(a, b, c, d), nrow=2)
    fres=fisher.test(mtx, alternative="greater")
    odf$pvalue[odf$`Gene Sets`==set]=fres$p.value
    odf$OR[odf$`Gene Sets`==set]=fres$estimate # Odds ratio

    odf$Overlap[odf$`Gene Sets`==set]=a
}
odf$padj=p.adjust(odf$pvalue,method="BH")
# save(odf, file=paste0(dir_data,"3. ORA with additional gene sets.Rdata"))


## Viz: UpSet plot
# 1. In this UpSet plot, each gene is assigned to a single, mutually exclusive intersection group.
# 2. To avoid displaying an excessive number of combinations, only the full intersection and intersections missing at most one set are visualized.
cuts=sapply(unique(odf$Reference), function(ref) {
    cut=sort(odf$Overlap[odf$Reference==ref],decreasing=T)[1:3]
    min(cut[!is.na(cut)])
})

vodf=do.call(rbind,
    lapply(names(cuts), function(ref) {
        vdf=odf[odf$Reference==ref,]
        vdf[vdf$Overlap>=cuts[ref],]
    })
)
vodf=vodf[c(1,setdiff(order(vodf$Reference),1)),] # sorting with reference

gs_viz=vodf$`Gene Sets`
groups=vodf$Reference
vgsl=gsl[gs_viz]

ccols_dir=setNames(c("#ffdb3a","#000dff"),c("Up","Down"))
fcols_dir=setNames(c("#ffdf50","#403dff"),c("Up","Down"))

dir="Up"

if(length(setdiff(gs_viz,names(vgsl)))>0) {for(gs in setdiff(gs_viz,names(vgsl))) {vgsl[[gs]]=character(0)}}
tgenes=unique(unlist(vgsl))
if(length(tgenes)==0) {
    tiff(filename=paste0(dir_fig,".tif"), width=10, height=10, units='cm',res=300)
    plot.new()
    text(x=0.5, y=0.5, labels="no DEGs", cex=3)
    dev.off() 
    next
}
vdf=data.frame(lapply(vgsl, function(genes) ifelse(tgenes %in% genes,1,0)), check.names=F)
rownames(vdf)=tgenes
tgenes = setNames(rownames(vdf) , apply(vdf,1, function(n) paste0(colnames(vdf)[n==1],collapse="&")))
if(sum(sapply(vgsl, length)!=0)==1) {
    isize=length(unlist(vgsl))
    ssize=sapply(vgsl, length) # set size
    gmtx=matrix(ifelse(gs_viz %in% names(vgsl)[sapply(vgsl, length)!=0],1,0),ncol=1,nrow=length(ssize))
    colnames(gmtx)=names(vgsl)[sapply(vgsl, length)!=0]
    rownames(gmtx)=gs_viz
    hmtx=matrix(1:length(vgsl), ncol=1, nrow=nrow(gmtx))
    colnames(hmtx)=colnames(gmtx)
    rownames(hmtx)=rownames(gmtx)
} else {
    combis=unique(names(tgenes))
    isize=sapply(combis, function(combi) sum(names(tgenes)==combi)) # intersection size
    ssize=sapply(vgsl, length) # set size

    gmtx=as.matrix(data.frame(lapply(combis, function(combi) ifelse(gs_viz %in% unlist(str_split(combi,"&")) , 1:length(gs_viz) ,0)))) 
    colnames(gmtx)=combis
    rownames(gmtx)=gs_viz

    gmtx=gmtx[,order(colnames(gmtx),decreasing=T)] # sorting
    gmtx=gmtx[,order(str_count(colnames(gmtx),"&"),decreasing=T)]
    gmtx=gmtx[gs_viz,]
    gmtx=gmtx[,grepl("Signature",colnames(gmtx))]

    if("numeric" %in% class(gmtx)) { gmtx=matrix(gmtx,ncol=1,nrow=length(gmtx), dimnames=list(names(gmtx),combis[grep("&",combis)])) }
    isize=isize[match(colnames(gmtx),names(isize))]

    hmtx=matrix(rep(1:nrow(gmtx),ncol(gmtx)), nrow=nrow(gmtx), ncol=ncol(gmtx))
    colnames(hmtx)=colnames(gmtx)
    rownames(hmtx)=rownames(gmtx)
}
ssize=ssize[match(rownames(hmtx),names(ssize))]

hfils=setNames(c('#ffa5a5','#ffdcbf','#edffc3','#eccfff','#b1e6ff'),c("Metastasis Signature","cancerSEA","CTD","EMTome","Hallmark"))
pfils=setNames(c('#ff3d3d','#ffa153','#c6ff43','#ca7aff','#22b9ff'),c("Metastasis Signature","cancerSEA","CTD","EMTome","Hallmark"))
pcols=setNames(c('#ff0000','#ff963f','#aee729','#c164ff','#1269ff'),c("Metastasis Signature","cancerSEA","CTD","EMTome","Hallmark"))

fill.ht=colorRamp2(1:length(gs_viz) , hfils[match(groups,names(hfils))])
fill.bp=colorRamp2(1:length(gs_viz) , rev(pfils[match(groups,names(pfils))]))
color.bp=colorRamp2(1:length(gs_viz) , rev(pcols[match(groups,names(pcols))]))

color.pt=colorRamp2(c(1:length(gs_viz),0) , c(rep("#000000",length(gs_viz)),"#a4a4a4"))
fill.pt=colorRamp2(c(1:length(gs_viz),0) , c(rep("#3f3f3f",length(gs_viz)),"#d4d4d4"))

# top annotation barplot
tagap=ifelse(max(isize)>=100,100,50)
tan.ht=HeatmapAnnotation( `Intersection\nSize`=anno_barplot( isize, height=unit(3.5, "cm"), gp=gpar(fill=fcols_dir[dir], col=ccols_dir[dir], lwd=2), border=FALSE, axis=FALSE, add_numbers=TRUE, numbers_gp=gpar(fontsize=21, col="#000000"), numbers_rot=0, numbers_offset=unit(0.3, "cm") ),
    annotation_name_side="left", annotation_name_gp=gpar(fontsize=18, col="#000000"), annotation_name_rot=90, annotation_name_offset=c(`Intersection\nSize`="0.3cm") )

# row annotation barplot
ragap=ifelse(max(ssize)>=100,100,50)
ran.ht=rowAnnotation(`DEGs\nCount` = anno_barplot(ssize, height=unit(4.3,"cm"), width = unit(2.3, "cm"), gp=gpar(fill=rev(attr(fill.bp, "colors")), col=rev(attr(color.bp, "colors")), lwd=2), border=F, axis=F, axis_param=list(gp=gpar(fontsize=13), at=seq(0,max(ssize),ragap), labels=seq(0,max(ssize),ragap)), annotation_name_gp=gpar(fontsize=18, col="#000000"), annotation_name_rot=90,
    add_numbers=T, numbers_gp=gpar(fontsize=21, col='#000000'), numbers_rot=0, numbers_offset=unit(0.3,"cm")),
                    annotation_name_offset=c(`DEGs count`="0.3cm"), annotation_name_gp=gpar(fontsize=18, col="#000000"), annotation_name_rot=0)

# heatmap
hp=Heatmap(hmtx, show_heatmap_legend=F, col=fill.ht, border=F, cluster_columns=F, cluster_rows=F, show_row_names=T, show_column_names=F, cluster_column_slices=F, show_row_dend=F, column_names_rot=85, row_names_side = "left", row_labels=gsub(" \\|.*","",rownames(hmtx)),
            width=ncol(gmtx)*unit(1.5,"cm"), height=nrow(gmtx)*unit(1.5,"cm"), rect_gp=gpar(col="#80808000", lwd=0, lty=1), row_names_gp=gpar(fontsize=22), column_names_gp=gpar(fontsize=25),
            cell_fun = function(j, i, x, y, width, height, fill) {
                current_value= gmtx[i, j]
                # grid.point
                point_color = color.pt(current_value)
                point_fill = fill.pt(current_value)
                grid.points(x, y, pch = 21, size = unit(1, "cm"), gp = gpar(col = point_color, fill=point_fill, lwd=2.5))
                # grid.semgent
                threshold = 0.5
                if (current_value > threshold) {
                    col_values = gmtx[, j]
                    significant_rows = which(col_values > threshold)
                    if (length(significant_rows) > 1) {
                        y_coords = (nrow(gmtx) - significant_rows + 0.5) / nrow(gmtx)
                        y_top = unit(min(y_coords), "npc")
                        y_bottom = unit(max(y_coords), "npc")

                        if (i == significant_rows[which.max(y_coords)]) {
                            grid.segments(x, y_top, x, y_bottom, gp = gpar(col = "#3f3f3f", lwd = 3))
                        }
                    }
                }
            },
            top_annotation=tan.ht,
            right_annotation=ran.ht
        )
tiff(filename=paste0(dir_fig,"Main Figure 3(A).tif"), width=(27+1.5*ncol(hmtx)), height=32, units = 'cm',res=300)
draw(hp, merge_legend=F, padding=unit(c(0, 5, 0, 0), "cm"), gap=unit(0.05,"cm")) 
dev.off()

# Legend
lg=Legend(labels=names(pfils)[c(1,setdiff(order(names(pfils)),1))], legend_gp=gpar(fill=pfils[c(1,setdiff(order(names(pfils)),1))]), title="Batch", border=pcols[c(1,setdiff(order(names(pcols)),1))], labels_gp=gpar(col='#000000', fontsize=12), grid_height=unit(0.5,"cm"), grid_width=unit(0.5,"cm"))
tiff(paste0(dir_fig,"Main Figure 3(A)_legend.tiff"), width = 5, height = 12, units="cm", res=300) 
draw(lg)
dev.off()


## Supplementary Data 1-2
sdf=data.frame(`Metastasis Signature`=sort(vgsl[["Metastasis Signature"]]), check.names=F)
# write.csv(sdf, file="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/Supplementary Data/Supplementary Data 1.csv", row.names=F)
sgsl=sapply(vgsl[names(vgsl)!='Metastasis Signature'], function(gs) {
    gs=c(sort(gs) , rep("",(max(sapply(vgsl,length))-length(gs))))
    data.frame(gs)
})
sdf=do.call(cbind , sgsl)
colnames(sdf)=setdiff(names(vgsl),'Metastasis Signature')
gs = lapply(sdf, function(x) {
  x = trimws(as.character(x))
  unique(x[!is.na(x) & x != ""])
})
gs = vapply(names(gs), function(x) { paste(c(x, "NA", gs[[x]]), collapse = "\t") }, FUN.VALUE = character(1))
# writeLines(gs, con = "G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/Supplementary Data/Supplementary Data 2.gmt", useBytes = TRUE)



########################################
## Metastasis Risk Prediction in Patient(TCGA) dataset
########################################

ggsl=c(gsl["Metastasis Signature"],list(`Known metastasis signature`=gsl[["Metastasis Signature"]][gsl[["Metastasis Signature"]] %in% unlist(gsl[2:length(gsl)])]),gsl[2:length(gsl)])
names(ggsl)[duplicated(names(ggsl))]=gsub(")"," v2)",names(ggsl)[duplicated(names(ggsl))])

mdfl=c()
effective_size=c()
for(db in names(tpml)) {
    tpm=tpml[[db]]
    param=ssgseaParam(tpm, ggsl)

     effective_size[[db]]=data.frame(
        Database=db,
        Signature=names(ggsl),
        Effective_gene_count=lengths(ggsl),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )
}

## Effective size
esdf=do.call(rbind,effective_size)
rownames(esdf)=1:nrow(esdf)

# save(mdfl, esdf, file=paste0(dir_data,"3. Predicted Metastasis risk.Rdata"))



########################################
## Add Treatment info
########################################
load(paste0(dir_fh,"sinfo_merged.Rdata"))
sinfol[['Firehose']]$`Neoadjuvant Treatment`=sinfo$`patient.history_of_neoadjuvant_treatment`[match(sinfol[["Firehose"]]$Sample,sinfo$patient.bcr_patient_barcode)]

load(paste0(dir_tcga,"clinical/pan_cancer.Rdata"))
sinfol[['GDC']]$`Neoadjuvant Treatment`=sinfo$`paper_History.of.neoadjuvant.treatment`[match(sinfol[["GDC"]]$Sample,sub("-[^-]+$", "", sinfo$bcr_patient_barcode))]
sinfol[['GDC']]$`Radiation Therapy`=sinfo$`paper_radiation.treatment.adjuvant`[match(sinfol[["GDC"]]$Sample,sub("-[^-]+$", "", sinfo$bcr_patient_barcode))]
sinfol[['GDC']]$`Radiation Therapy`[sinfol[['GDC']]$`Radiation Therapy`=='[Unknown]' & !is.na(sinfol[['GDC']]$`Radiation Therapy`)]=NA


## No. of Samples
table(sinfol[['Firehose']]$`Radiation Therapy`, useNA='ifany')
#   No  Yes <NA> 
#   14   14  487 
table(sinfol[['Firehose']]$`Neoadjuvant Treatment`, useNA='ifany')
#  no yes 
# 512   3 

table(sinfol[['GDC']]$`Radiation Therapy`, useNA='ifany')
#    NO   YES  <NA> 
#   136    57 10169 
table(sinfol[['GDC']]$`Neoadjuvant Treatment`, useNA='ifany')
#   No  Yes <NA> 
#  389   10 9963 



########################################
## Continuous survival analysis Cases
########################################
atypes=c("mfs","rfs","os")
occurences=setNames(c("Metastasis Ocurrence","Recurrence Ocurrence","Vital Status"), atypes)
free_survivals=setNames(c("Metastasis Free Survival","Relapse Free Survival","Overall Survival"), atypes)
nons=setNames(c("No Metastasis","No Recurrence","Alive"), atypes)
adjust_template=c("Age","Gender","Race","Stage","Cancer Type", "Neoadjuvant Treatment")
categorical_template=c("Gender","Race","Stage","Cancer Type", "Neoadjuvant Treatment")

quote_term=function(x) paste0("`",x,"`")



########################################
## Construct administratively censored outcome
########################################
make_endpoint_data=function(db, atype, horizon_year) {
    sinfo=sinfol[[db]]
    mdf=mdfl[[db]]

    dat=merge(sinfo, mdf, by="Sample", all=FALSE, sort=FALSE)
    if(!"Stage" %in% colnames(dat)) {
        if("Pathologic Stage" %in% colnames(dat)) {
            colnames(dat)[colnames(dat)=="Pathologic Stage"]="Stage"
        } else if("AJCC Stage" %in% colnames(dat)) {
            colnames(dat)[colnames(dat)=="AJCC Stage"]="Stage"
        }
    }

    event_col=occurences[atype]
    time_col=free_survivals[atype]

    dat=dat[!is.na(dat[[event_col]]) & !is.na(dat[[time_col]]), , drop=F]
    dat$day=as.numeric(dat[[time_col]])
    dat$event_original=ifelse(dat[[event_col]]==nons[atype], 0, 1)

    dcut = 365.25 * horizon_year

    dat$time=pmin(dat$day,dcut)
    dat$event=as.integer(dat$day<=dcut & dat$event_original==1)

    dat=dat[!is.na(dat$time) & dat$time>=0 & !is.na(dat$event), , drop=F]
    dat
}



########################################
## Common complete-case model data
########################################
prepare_analysis_data=function(dat,score_vars,adjust_vars) {
    score_vars=intersect(score_vars,colnames(dat))
    adjust_vars=intersect(adjust_vars,colnames(dat))

    repeat {
        vars=c("time","event",score_vars,adjust_vars)
        x=dat[complete.cases(dat[,vars,drop=FALSE]), vars, drop=F]
        if(nrow(x)==0) { return(list(data=x,adjust_vars=character(0))) }
        invalid=adjust_vars[vapply(adjust_vars, function(v) length(unique(x[[v]]))<2, logical(1))]
        if(length(invalid)==0) break
        adjust_vars=setdiff(adjust_vars,invalid)
    }
    categorical_vars=intersect(categorical_template, adjust_vars)

    for(v in categorical_vars) { x[[v]]=droplevels(factor(x[[v]])) }
    list(data=x,adjust_vars=adjust_vars)
}



########################################
## Extract C-index and 95% CI
########################################
extract_cindex=function(cfit) {
    est=as.numeric(coef(cfit))[1]
    vc=as.matrix(vcov(cfit))
    se=sqrt(vc[1,1])

    data.frame(C_index=est, C_low=max(0,est-1.96*se), C_high=min(1,est+1.96*se), C_SE=se, check.names=F)
}



########################################
## HR per SD and C-index for one signature
########################################
continuous_signature_analysis=function(dat, score_name, adjust_vars=adjust_template) {
    prep=prepare_analysis_data(dat=dat, score_vars=score_name, adjust_vars=adjust_vars)
    x=prep$data
    adjust_vars=prep$adjust_vars

    if(nrow(x)==0 || length(unique(x$event))<2) return(NULL)
    if(sum(x$event)<3 || sum(x$event==0)<3) return(NULL)
    score_sd=sd(x[[score_name]],na.rm=TRUE)

    if(is.na(score_sd) || score_sd==0) return(NULL)
    x$ScoreSD=as.numeric(scale(x[[score_name]]))

    rhs=c("ScoreSD", quote_term(adjust_vars))

    cox_formula=as.formula( paste( "Surv(time,event) ~", paste(rhs,collapse=" + ")))

    fit=tryCatch(coxph(cox_formula, data=x, ties="efron", x=T, y=T, model=T), error=function(e) NULL)
    if(is.null(fit)) return(NULL)

    csum=summary(fit)
    if(!"ScoreSD" %in% rownames(csum$coefficients)) return(NULL)

    ## Adjusted HR per one-SD increase
    hr=csum$coefficients["ScoreSD","exp(coef)"]
    hr_p=csum$coefficients["ScoreSD","Pr(>|z|)"]
    hr_low=csum$conf.int["ScoreSD","lower .95"]
    hr_high=csum$conf.int["ScoreSD","upper .95"]

    ## Primary C-index: score only
    marker_cfit=concordance(Surv(time,event)~ScoreSD, data=x, reverse=TRUE, timewt="n")
    marker_c=extract_cindex(marker_cfit)

    ## Secondary C-index: clinical covariates plus score
    model_cfit=concordance(fit, timewt="n")
    model_c=extract_cindex(model_cfit)

    data.frame(Signature=score_name, N=nrow(x), Events=sum(x$event), HR_per_SD=hr, HR_low=hr_low, HR_high=hr_high, HR_P=hr_p, C_score=marker_c$C_index, C_score_low=marker_c$C_low, C_score_high=marker_c$C_high, C_model=model_c$C_index, C_model_low=model_c$C_low, C_model_high=model_c$C_high, stringsAsFactors=FALSE, check.names=FALSE
    )
}



########################################
## Run continuous analyses
########################################
horizons=c(1,5,10)
continuous_results=list()
result_index=1

for(db in names(mdfl)) {
    score_names=setdiff(colnames(mdfl[[db]]),"Sample")

    for(atype in atypes) {
        for(horizon in horizons) {
            dat=make_endpoint_data(db=db, atype=atype, horizon_year=horizon)

            for(score_name in score_names) {
                res=continuous_signature_analysis(dat=dat, score_name=score_name)

                if(is.null(res)) next
                res$Database=db
                res$Endpoint=atype
                res$Horizon_year=horizon
                continuous_results[[result_index]]=res
                result_index=result_index+1
            }
        }
    }
}
continuous_results=do.call(rbind,continuous_results)

hr_group=interaction(continuous_results$Database, continuous_results$Endpoint, continuous_results$Horizon_year, drop=T)
continuous_results$HR_FDR=ave(continuous_results$HR_P, hr_group, FUN=function(x) p.adjust(x,method="BH"))
continuous_results=merge(continuous_results, esdf, by.x=c("Database","Signature"), by.y=c("Database","Signature"), all.x=T, sort=F)
# save(continuous_results, file=paste0(dir_data, "3. Continuous HR and C-index results.Rdata"))
# write.csv(continuous_results, paste0(dir_data,"3. Continuous HR and C-index results.csv"), row.names=F)



########################################
## Formal paired Delta C-index test
########################################

## Define function
delta_cindex_test=function(dat, reference_score="Metastasis Signature", comparator_score, adjusted=FALSE, adjust_vars=adjust_template) {
    prep=prepare_analysis_data(dat=dat, score_vars=c(reference_score,comparator_score), adjust_vars=adjust_vars)
    x=prep$data
    adjust_vars=prep$adjust_vars

    if(nrow(x)==0 || length(unique(x$event))<2) return(NULL)
    if(sum(x$event)<3 || sum(x$event==0)<3) return(NULL)

    if(sd(x[[reference_score]])==0) return(NULL)
    if(sd(x[[comparator_score]])==0) return(NULL)

    x$ReferenceSD=as.numeric(scale(x[[reference_score]]))
    x$ComparatorSD=as.numeric(scale(x[[comparator_score]]))

    if(!adjusted) {
        cfit=concordance(Surv(time,event)~ReferenceSD+ComparatorSD, data=x, reverse=TRUE, timewt="n")
    } else {
        rhs_ref=c("ReferenceSD", quote_term(adjust_vars))
        rhs_cmp=c("ComparatorSD", quote_term(adjust_vars))
        fit_ref=coxph(as.formula(paste("Surv(time,event) ~", paste(rhs_ref,collapse=" + "))), data=x,  ties="efron", x=TRUE, y=TRUE)
        fit_cmp=coxph(as.formula(paste("Surv(time,event) ~", paste(rhs_cmp,collapse=" + ") )), data=x, ties="efron", x=T, y=T)
        cfit=concordance(fit_ref, fit_cmp, timewt="n")
    }

    cvalues=as.numeric(coef(cfit))
    vc=as.matrix(vcov(cfit))

    contrast=c(1,-1)
    delta=as.numeric(contrast%*%cvalues)
    delta_se=sqrt(as.numeric(contrast%*%vc%*%contrast))
    z=delta/delta_se
    pvalue=2*pnorm(abs(z),lower.tail=FALSE)

    data.frame(Reference=reference_score, Comparator=comparator_score, Model=ifelse(adjusted, "Clinical covariates + score", "Score only"), N=nrow(x), Events=sum(x$event), C_reference=cvalues[1], C_comparator=cvalues[2], Delta_C=delta, Delta_C_low=delta-1.96*delta_se, Delta_C_high=delta+1.96*delta_se, Delta_C_SE=delta_se, Z=z, P_value=pvalue, stringsAsFactors=F, check.names=F)
}


## Run function
delta_results=list()
result_index=1

for(db in names(mdfl)) {
    score_names=setdiff(colnames(mdfl[[db]]),c("Sample","Metastasis Signature"))

    for(atype in atypes) {
        for(horizon in horizons) {
            dat=make_endpoint_data(db=db,atype=atype,horizon_year=horizon)

            for(comparator in score_names) {
                res1=delta_cindex_test(dat=dat,comparator_score=comparator,adjusted=FALSE)

                if(!is.null(res1)) {
                    res1$Database=db
                    res1$Endpoint=atype
                    res1$Horizon_year=horizon
                    delta_results[[result_index]]=res1
                    result_index=result_index+1
                }

                res2=delta_cindex_test(dat=dat,comparator_score=comparator,adjusted=TRUE)

                if(!is.null(res2)) {
                    res2$Database=db
                    res2$Endpoint=atype
                    res2$Horizon_year=horizon
                    delta_results[[result_index]]=res2
                    result_index=result_index+1
                }
            }
        }
    }
}
delta_results=do.call(rbind,delta_results)

delta_group=interaction(delta_results$Database, delta_results$Endpoint, delta_results$Horizon_year, delta_results$Model, drop=T)
delta_results$Adjusted_P_value=ave( delta_results$P_value, delta_group, FUN=function(x) p.adjust(x,method="BH") )

# save(delta_results, file=paste0(dir_data, "3. Formal Delta C-index comparison.Rdata"))
# write.csv(delta_results, paste0(dir_data, "3. Formal Delta C-index comparison.csv"), row.names=F)



########################################
## Size-matched 13-gene subset controls
########################################

## C-index from a continuous score
cindex_from_score=function(time,event,score) {
    dat=data.frame(time=time,event=event,score=as.numeric(score))
    dat=dat[complete.cases(dat),,drop=F]

    if(nrow(dat)<2) return(NA_real_)
    if(length(unique(dat$event))<2) return(NA_real_)
    if(sum(dat$event)<3 || sum(dat$event==0)<3) return(NA_real_)

    score_sd=sd(dat$score)
    if(is.na(score_sd) || score_sd==0) return(NA_real_)

    dat$score=as.numeric(scale(dat$score))

    cfit=tryCatch(concordance(Surv(time,event)~score, data=dat, reverse=T, timewt="n"), error=function(e) NULL)
    if(is.null(cfit)) return(NA_real_)

    as.numeric(coef(cfit))[1]
}


## Reproducible seed from text
seed_from_text=function(text, base_seed=123) {
    values=utf8ToInt(text)
    seed=base_seed+sum(values*seq_along(values))
    seed=as.integer(seed%%2147483646)+1
    seed
}


## Generate unique size-matched subsets
generate_size_matched_sets=function(genes,target_size,B,seed=123) {
    genes=sort(unique(genes))
    n_gene=length(genes)

    if(n_gene<target_size) { return(list( sets=NULL, mode="Fewer genes than target size", total_possible=0 )) }
    if(n_gene==target_size) { return(list( sets=NULL, mode="Already size-matched", total_possible=1 )) }
    total_possible=choose(n_gene,target_size)

    # Enumerate every possible subset when the number is small
    if(is.finite(total_possible) && total_possible<=B && total_possible<=5000) {
        sets=combn(genes,target_size,simplify=FALSE)
        return(list( sets=sets, mode="All possible subsets", total_possible=total_possible ))
    }

    # Otherwise generate B unique random subsets
    set.seed(seed)

    sets=list()
    set_keys=character(0)
    iteration=0
    max_iteration=B*100

    while(length(sets)<B && iteration<max_iteration) {
        iteration=iteration+1

        candidate=sort(sample(genes,size=target_size,replace=FALSE))
        candidate_key=paste(candidate,collapse="\001")
        if(!candidate_key %in% set_keys) { sets[[length(sets)+1]]=candidate ; set_keys=c(set_keys,candidate_key) }
    }

    list( sets=sets, mode="Unique random subsets", total_possible=total_possible )
}


## Prepare all outcome cases for one database
make_size_control_cases=function(db,tpm_db,atypes,horizons,reference_score="Metastasis Signature") {
    case_list=list()
    case_log=list()
    case_index=1

    for(atype in atypes) {
        for(horizon in horizons) {
            dat=make_endpoint_data(db=db,atype=atype,horizon_year=horizon)
            status="Included"
            reason=NA_character_

            if(!reference_score %in% colnames(dat)) { status="Skipped" ; reason="Reference score was not found" }

            if(status=="Included") {
                common_samples=intersect(dat$Sample,colnames(tpm_db))
                if(length(common_samples)==0) { status="Skipped" ; reason="No common expression and survival samples" }
            }

            if(status=="Included") {
                dat=dat[match(common_samples,dat$Sample),,drop=F]
                dat=dat[!duplicated(dat$Sample),,drop=F]
                dat=dat[ !is.na(dat$time) & !is.na(dat$event) & !is.na(dat[[reference_score]]), , drop=F ]

                if(nrow(dat)==0) {
                    status="Skipped"
                    reason="No complete survival observations"
                } else if(length(unique(dat$event))<2) {
                    status="Skipped"
                    reason="Only one event status"
                } else if(sum(dat$event)<3) {
                    status="Skipped"
                    reason="Fewer than three events"
                } else if(sum(dat$event==0)<3) {
                    status="Skipped"
                    reason="Fewer than three non-events"
                }
            }

            if(status=="Included") {
                sample_index=match(dat$Sample,colnames(tpm_db))
                if(any(is.na(sample_index))) { status="Skipped" ; reason="Expression sample matching failed" }
            }

            if(status=="Included") {
                target_c=cindex_from_score(time=dat$time,event=dat$event,score=dat[[reference_score]])
                if(is.na(target_c)) { status="Skipped" ; reason="Reference C-index could not be estimated" }
            }

            case_name=paste(db,atype,horizon,sep="__")

            if(status=="Included") { case_list[[case_name]] = list(Database=db, Endpoint=atype, Horizon_year=horizon, data=dat, sample_index=sample_index, Target_C_index=target_c) }

            case_log[[case_index]] = data.frame( Database=db, Endpoint=atype, Horizon_year=horizon, Status=status, Reason=reason, N=if(exists("dat")) nrow(dat) else NA_integer_, Events=if(exists("dat") && "event" %in% colnames(dat)) sum(dat$event,na.rm=T) else NA_integer_, stringsAsFactors=F, check.names=F )
            case_index=case_index+1
        }
    }
    list(cases=case_list,log=do.call(rbind,case_log))
}


## Evaluate one collection of control sets
evaluate_control_sets_all_cases=function(tpm_db,case_list, gene_set_list, parent_signature, control_type, parent_gene_count, sampling_mode, total_possible, batch_size=100, label="") {
    if(length(gene_set_list)==0 || length(case_list)==0) return(NULL)
    batches=split(seq_along(gene_set_list), ceiling(seq_along(gene_set_list)/batch_size))

    output=list()
    output_index=1

    cat("\n",label,"\n",sep="")
    cat("Total gene sets:",length(gene_set_list),"\n")
    cat("Total survival cases:",length(case_list),"\n")
    cat("Total GSVA batches:",length(batches),"\n")

    pb=txtProgressBar(min=0,max=length(batches),style=3)
    start_time=Sys.time()

    for(batch_no in seq_along(batches)) {
        batch=batches[[batch_no]]
        sets_batch=gene_set_list[batch]

        param=ssgseaParam(exprData=tpm_db,geneSets=sets_batch,minSize=1,maxSize=Inf,alpha=0.25,normalize=T)

        score_batch=tryCatch(
            gsva(param,verbose=F),
            error=function(e) { message("\nGSVA error: ",conditionMessage(e)) ; NULL }
        )

        if(is.null(score_batch)) {
            setTxtProgressBar(pb,batch_no)
            next
        }

        score_batch=as.matrix(score_batch)

        # Ensure the score columns follow the TPM sample order
        common_score_samples=intersect(colnames(tpm_db),colnames(score_batch))
        score_batch=score_batch[,match(common_score_samples,colnames(score_batch)),drop=F]

        for(case_name in names(case_list)) {
            case=case_list[[case_name]]
            case_score_index=match(case$data$Sample,colnames(score_batch))

            if(any(is.na(case_score_index))) next

            parent_full_c=NA_real_

            if(parent_signature %in% colnames(case$data)) {
                parent_full_c=cindex_from_score(time=case$data$time,event=case$data$event,score=case$data[[parent_signature]])
            }

            cvalues=vapply(
                seq_len(nrow(score_batch)),
                function(i) {
                    cindex_from_score(time=case$data$time,event=case$data$event,score=score_batch[i,case_score_index])
                },
                numeric(1)
            )

            output[[output_index]] = data.frame(Database=case$Database, Endpoint=case$Endpoint, Horizon_year=case$Horizon_year, Control_type=control_type, Parent_signature=parent_signature, Parent_gene_count=parent_gene_count, Sampling_mode=sampling_mode, Total_possible_subsets=total_possible, Set=rownames(score_batch), C_index=cvalues, Target_C_index=case$Target_C_index, Parent_full_C_index=parent_full_c, N=nrow(case$data), Events=sum(case$data$event), stringsAsFactors=F, check.names=F)
            output_index=output_index+1
        }
        setTxtProgressBar(pb,batch_no)
    }

    close(pb)

    elapsed=as.numeric(difftime(Sys.time(), start_time,units="mins"))
    cat("\nCompleted in",round(elapsed,2),"minutes\n")

    if(length(output)==0) return(NULL)

    do.call(rbind,output)
}


## Run every database, endpoint and horizon
run_size_matched_controls=function(dbs, atypes, horizons, comparators, ggsl, tpml, mdfl, target_genes, gene_universe, B=1000, batch_size=100, base_seed=123, reference_score="Metastasis Signature", checkpoint_file=NULL) {
    result_parts=list()
    result_index=1
    skip_parts=list()
    skip_index=1
    case_log_parts=list()
    case_log_index=1

    comparators=unique(comparators)
    comparators=setdiff(comparators,reference_score)
    comparators=intersect(comparators,names(ggsl))
    target_size=length(unique(target_genes))

    for(db in dbs) {
        cat("\n========================================\n")
        cat("Database:",db,"\n")
        cat("========================================\n")

        if(!db %in% names(tpml) || !db %in% names(mdfl)) {
            skip_parts[[skip_index]] = data.frame( Database=db, Parent_signature=NA_character_, Control_type=NA_character_, Effective_gene_count=NA_integer_, Reason="Database was not found in tpml or mdfl", stringsAsFactors=F )

            skip_index=skip_index+1
            next
        }

        tpm_db=as.matrix(tpml[[db]])

        db_samples=intersect(colnames(tpm_db),unique(mdfl[[db]]$Sample))

        tpm_db=tpm_db[,db_samples,drop=F]

        case_object=make_size_control_cases(db=db,tpm_db=tpm_db,atypes=atypes,horizons=horizons,reference_score=reference_score)

        case_list=case_object$cases
        case_log_parts[[case_log_index]]=case_object$log
        case_log_index=case_log_index+1

        if(length(case_list)==0) {
            skip_parts[[skip_index]] = data.frame( Database=db, Parent_signature=NA_character_, Control_type=NA_character_, Effective_gene_count=NA_integer_, Reason="No valid survival cases", stringsAsFactors=F )

            skip_index=skip_index+1
            next
        }

        # Comparator-derived 13-gene subsets
        for(comparator in comparators) {
            comparator_genes = intersect(unique(ggsl[[comparator]]),rownames(tpm_db))

            n_comparator=length(comparator_genes)

            if(n_comparator<target_size) {
                skip_parts[[skip_index]] = data.frame( Database=db, Parent_signature=comparator, Control_type="13-gene subset from comparator", Effective_gene_count=n_comparator, Reason=paste0( "Fewer than ", target_size, " effective genes" ), stringsAsFactors=F )

                skip_index=skip_index+1
                next
            }

            if(n_comparator==target_size) {
                skip_parts[[skip_index]] = data.frame( Database=db, Parent_signature=comparator, Control_type="13-gene subset from comparator", Effective_gene_count=n_comparator, Reason="Already size-matched; direct C-index comparison is used", stringsAsFactors=F )

                skip_index=skip_index+1
                next
            }

            comparator_seed=seed_from_text(paste(db,comparator,sep="__"),base_seed=base_seed)
            generated=generate_size_matched_sets(genes=comparator_genes,target_size=target_size,B=B,seed=comparator_seed)
            random_sets=generated$sets

            if(is.null(random_sets) || length(random_sets)==0) next

            names(random_sets)=paste0(make.names(comparator),"__subset_",sprintf("%04d",seq_along(random_sets)))

            temp = evaluate_control_sets_all_cases( tpm_db=tpm_db, case_list=case_list, gene_set_list=random_sets, parent_signature=comparator, control_type="13-gene subset from comparator", parent_gene_count=n_comparator, sampling_mode=generated$mode, total_possible=generated$total_possible, batch_size=batch_size, label=paste0(db," | ",comparator) )
            if(!is.null(temp)) { result_parts[[result_index]]=temp ; result_index=result_index+1 }
        }

        # Background random 13-gene sets
        background_genes = intersect(unique(gene_universe),rownames(tpm_db))
        background_genes = setdiff(background_genes,unique(target_genes))

        n_background=length(background_genes)

        if(n_background>=target_size) {
            background_seed=seed_from_text( paste(db,"Random expressed genes",sep="__"), base_seed=base_seed )

            generated = generate_size_matched_sets(genes=background_genes,target_size=target_size,B=B,seed=background_seed)

            background_sets=generated$sets

            if(!is.null(background_sets) && length(background_sets)>0) {
                names(background_sets) = paste0("Background_13gene_",sprintf("%04d",seq_along(background_sets)))

                temp = evaluate_control_sets_all_cases( tpm_db=tpm_db, case_list=case_list, gene_set_list=background_sets, parent_signature="Random expressed genes", control_type="Random background 13-gene set", parent_gene_count=n_background, sampling_mode=generated$mode, total_possible=generated$total_possible, batch_size=batch_size, label=paste0(db," | Random background genes") )

                if(!is.null(temp)) {
                    result_parts[[result_index]]=temp
                    result_index=result_index+1
                }
            }
        } else {
            skip_parts[[skip_index]] = data.frame( Database=db, Parent_signature="Random expressed genes", Control_type="Random background 13-gene set", Effective_gene_count=n_background, Reason=paste0( "Fewer than ", target_size, " background genes" ), stringsAsFactors=F )

            skip_index=skip_index+1
        }

        # Checkpoint after each database
        if(!is.null(checkpoint_file)) {
            checkpoint_results=if(length(result_parts)>0) do.call(rbind,result_parts) else NULL
            checkpoint_skip=if(length(skip_parts)>0) do.call(rbind,skip_parts) else NULL
            checkpoint_cases=if(length(case_log_parts)>0) do.call(rbind,case_log_parts) else NULL

            saveRDS(list(results=checkpoint_results,skip_log=checkpoint_skip,case_status=checkpoint_cases),checkpoint_file,compress=F)
            cat("Checkpoint saved:",checkpoint_file,"\n")
        }

        rm(tpm_db)
        gc()
    }

    results=if(length(result_parts)>0) do.call(rbind,result_parts) else NULL
    skip_log=if(length(skip_parts)>0) do.call(rbind,skip_parts) else NULL
    case_status=if(length(case_log_parts)>0) do.call(rbind,case_log_parts) else NULL

    list(results=results,skip_log=skip_log,case_status=case_status)
}


## Summarize empirical distributions
summarize_size_control_results=function(size_control_results) {
    if(is.null(size_control_results) || nrow(size_control_results)==0) return(NULL)

    split_key=interaction(size_control_results$Database,size_control_results$Endpoint,size_control_results$Horizon_year,size_control_results$Control_type,size_control_results$Parent_signature,drop=T,lex.order=T)

    summary_result=do.call( rbind,
        lapply(
            split(size_control_results,split_key),
            function(x) {
                valid_c=x$C_index[!is.na(x$C_index)]
                target_c=x$Target_C_index[1]

                if(length(valid_c)==0) {
                    return(data.frame(Database=x$Database[1], Endpoint=x$Endpoint[1], Horizon_year=x$Horizon_year[1], Control_type=x$Control_type[1], Parent_signature=x$Parent_signature[1], Parent_gene_count=x$Parent_gene_count[1], Sampling_mode=x$Sampling_mode[1], Total_possible_subsets=x$Total_possible_subsets[1], Iterations=0, Mean_C=NA_real_, Median_C=NA_real_, SD_C=NA_real_, C_2.5=NA_real_, C_97.5=NA_real_, Target_C=target_c, Parent_full_C=x$Parent_full_C_index[1], Target_minus_Median_C=NA_real_, Target_percentile=NA_real_, Empirical_P=NA_real_,N=x$N[1], Events=x$Events[1], stringsAsFactors=F, check.names=F))
                }

                empirical_p=(1+sum(valid_c>=target_c))/(1+length(valid_c))

                data.frame(Database=x$Database[1], Endpoint=x$Endpoint[1], Horizon_year=x$Horizon_year[1], Control_type=x$Control_type[1], Parent_signature=x$Parent_signature[1], Parent_gene_count=x$Parent_gene_count[1], Sampling_mode=x$Sampling_mode[1], Total_possible_subsets=x$Total_possible_subsets[1], Iterations=length(valid_c),Mean_C=mean(valid_c), Median_C=median(valid_c), SD_C=sd(valid_c), C_2.5=unname(quantile(valid_c,0.025)), C_97.5=unname(quantile(valid_c,0.975)), Target_C=target_c, Parent_full_C=x$Parent_full_C_index[1], Target_minus_Median_C=target_c-median(valid_c), Target_percentile=100*mean(valid_c<=target_c), Empirical_P=empirical_p, N=x$N[1], Events=x$Events[1], stringsAsFactors=F, check.names=F)
            }
        )
    )

    adjustment_group=interaction(summary_result$Database, summary_result$Endpoint, summary_result$Horizon_year, summary_result$Control_type, drop=T)
    summary_result$Adjusted_P=ave(summary_result$Empirical_P, adjustment_group,
        FUN=function(x) {
            adjusted=rep(NA_real_,length(x))
            valid=!is.na(x)
            adjusted[valid]=p.adjust(x[valid],method="BH")
            adjusted
        }
    )

    summary_result
}


## Run all cases
# Gene universe must have been defined before ddf is overwritten
analysis_gene_universe=unique(ddf$Gene)

stopifnot(exists("analysis_gene_universe"))
stopifnot(exists("ggsl"))
stopifnot(exists("vodf"))

# Selected signatures shown in the comparison figure
primary_comparators=unique(c( "Known metastasis signature", vodf$`Gene Sets`))
primary_comparators=setdiff( primary_comparators, "Metastasis Signature" )
primary_comparators=intersect(primary_comparators,names(ggsl))

# All available databases
size_control_dbs=intersect(names(tpml),names(mdfl))

size_control_run = run_size_matched_controls( dbs=size_control_dbs, atypes=atypes, horizons=horizons, comparators=primary_comparators, ggsl=ggsl, tpml=tpml, mdfl=mdfl, target_genes=fdegs, gene_universe=analysis_gene_universe, B=1000, batch_size=100, base_seed=123, reference_score="Metastasis Signature", checkpoint_file=paste0( dir_data, "3. Size-matched control checkpoint.rds" ) )

size_control_results=size_control_run$results
size_control_skip_log=size_control_run$skip_log
size_control_case_status=size_control_run$case_status

size_control_summary=summarize_size_control_results( size_control_results )


## Save
# save(size_control_results, size_control_summary, size_control_skip_log, size_control_case_status, file=paste0( dir_data, "3. Size-matched 13-gene control results.Rdata"), compress="gzip" )
# fwrite(size_control_summary, paste0( dir_data, "3. Size-matched 13-gene control summary.csv"))
# fwrite(size_control_case_status, paste0( dir_data, "3. Size-matched 13-gene control case status.csv"))
# if(!is.null(size_control_skip_log)) { fwrite( size_control_skip_log, paste0( dir_data, "3. Size-matched 13-gene control skip log.csv" ) ) }

# Iteration-level CSV can be very large.
# Rdata already contains the complete iteration-level result.
# fwrite(size_control_results, paste0(dir_data, "3. Size-matched 13-gene control iteration results.csv"))



########################################
## Visualize Results
########################################

## Viz setting
display_signature=function(x) {
    gsub(" \\| "," from ",x)
}

safe_filename=function(x) {
    x=gsub("[\\\\/:*?\"<>|]","_",x)
    x=gsub("[[:space:]]+"," ",x)
    trimws(x)
}

format_p=function(x) {
    ifelse(
        is.na(x),
        "NA",
        ifelse(
            x<0.001,
            formatC(x,format="e",digits=2),
            formatC(x,format="f",digits=3)
        )
    )
}

get_ctype=function(db) {
    if(db=="Firehose") return("Lung Cancer")
    if(db=="GDC") return("Pan Cancer")
    db
}

endpoint_labels=setNames(
    c("Metastasis-free survival","Relapse-free survival","Overall survival"),
    c("mfs","rfs","os")
)

endpoint_short=setNames(
    c("MFS","RFS","OS"),
    c("mfs","rfs","os")
)

plot_signatures=unique(c(
    "Metastasis Signature",
    "Known metastasis signature",
    vodf$`Gene Sets`
))

plot_signatures=intersect(
    plot_signatures,
    unique(continuous_results$Signature)
)


## Output directories
dir_hr=paste0(dir_fig,"Continuous analysis/HR per SD/")
dir_cindex=paste0(dir_fig,"Continuous analysis/C-index/")
dir_delta=paste0(dir_fig,"Continuous analysis/Delta C-index/")
dir_size=paste0(dir_fig,"Continuous analysis/Size-matched control/")
dir_km_table=paste0(dir_fig,"Continuous analysis/KM exploratory with table/")
dir_km_notable=paste0(dir_fig,"Continuous analysis/KM exploratory no table/")

dir.create(dir_hr,recursive=T,showWarnings=F)
dir.create(dir_cindex,recursive=T,showWarnings=F)
dir.create(dir_delta,recursive=T,showWarnings=F)
dir.create(dir_size,recursive=T,showWarnings=F)
dir.create(dir_km_table,recursive=T,showWarnings=F)
dir.create(dir_km_notable,recursive=T,showWarnings=F)


## Gene set order (plotting signature)
plot_signatures=unique(c('Metastasis Signature',unique(with(subset(continuous_results, Database=="Firehose" & Endpoint=="mfs" & Horizon_year==1 & Signature %in% plot_signatures & !is.na(HR_per_SD)), Signature[order(HR_per_SD,decreasing=T)]))))



########################################
## HR per SD forest plots
########################################
for(db in unique(continuous_results$Database)) {
    for(atype in unique(continuous_results$Endpoint)) {
        for(horizon in sort(unique(continuous_results$Horizon_year))) {

            hdf=continuous_results[ continuous_results$Database==db & continuous_results$Endpoint==atype & continuous_results$Horizon_year==horizon & continuous_results$Signature %in% plot_signatures, , drop=F ]
            hdf=hdf[complete.cases(hdf[,c("HR_per_SD","HR_low","HR_high")]),,drop=F]

            if(nrow(hdf)==0) next

            current_order=plot_signatures[ plot_signatures %in% hdf$Signature]
            hdf$Feature=display_signature(hdf$Signature)
            hdf$Feature=factor(hdf$Feature,levels=rev(display_signature(current_order)))

            hdf$Significant=!is.na(hdf$HR_P) & hdf$HR_P<0.05
            hdf$HR_label=sprintf("%.3f",hdf$HR_per_SD)

            gp=ggplot(hdf, aes(x=Feature,y=HR_per_SD)) +
                geom_errorbar( aes(ymin=HR_low,ymax=HR_high), width=0.2, linewidth=1.2) +
                geom_point( aes(fill=Significant), size=5, shape=21, color="black", stroke=1) +
                geom_hline( yintercept=1, linetype="dashed", color="red", linewidth=1.2 ) +
                geom_text( aes(label=HR_label), vjust=-0.75, size=5 ) +
                scale_fill_manual( values=c( `FALSE`="#3e3e3e", `TRUE`="#ffdf50" ), guide="none" ) +
                scale_y_log10( expand=expansion(mult=c(0.08,0.22)) ) +
                coord_flip(clip="off") +
                labs(title=paste0( get_ctype(db), " - ", horizon, "-year ", endpoint_short[atype] ), x="Gene signature", y="Hazard ratio per SD" ) +
                theme_bw() +
                theme(
                    plot.title=element_text(size=10,face="bold"),
                    axis.title=element_text(size=18,color="black"),
                    axis.text.x=element_text(size=14,color="black"),
                    axis.text.y=element_text(size=15,color="black"),
                    panel.grid=element_blank(),
                    plot.margin=margin(10,35,10,10)
                )

            filename=paste0(dir_hr,safe_filename(paste(db,toupper(atype),paste0(horizon,"year"),sep="_")),".tiff")
            ggsave(filename,plot=gp,width=22,height=22,units="cm",dpi=300,limitsize=F,compression="lzw")
        }
    }
}



########################################
## Score-only Harrell's C-index plots
########################################
for(db in unique(continuous_results$Database)) {
    for(atype in unique(continuous_results$Endpoint)) {
        for(horizon in sort(unique(continuous_results$Horizon_year))) {

            cdf=continuous_results[continuous_results$Database==db & continuous_results$Endpoint==atype &
                continuous_results$Horizon_year==horizon & continuous_results$Signature %in% plot_signatures,,drop=F]
            cdf=cdf[complete.cases(cdf[,c("C_score","C_score_low","C_score_high")]),,drop=F]
            if(nrow(cdf)==0) next

            current_order=plot_signatures[plot_signatures %in% cdf$Signature]
            cdf$Feature=display_signature(cdf$Signature)
            cdf$Feature=factor(cdf$Feature,levels=rev(display_signature(current_order)))
            cdf$C_label=sprintf("%.3f",cdf$C_score)

            cdf$Top3=FALSE
            top_index=order(cdf$C_score,decreasing=T)[seq_len(min(3,nrow(cdf)))]
            cdf$Top3[top_index]=TRUE

            c_range=range(c(cdf$C_score_low,cdf$C_score_high,0.5),na.rm=T)
            c_pad=max(0.02,diff(c_range)*0.12)
            lower_limit=max(0,c_range[1]-c_pad)
            upper_limit=min(1,c_range[2]+c_pad)

            gp=ggplot(cdf,aes(x=Feature,y=C_score)) +
                geom_errorbar(aes(ymin=C_score_low,ymax=C_score_high),width=0.2,linewidth=1.2,color="black") +
                geom_point(aes(fill=Top3),size=5,shape=21,color="black",stroke=1) +
                geom_hline(yintercept=0.5,linetype="dashed",color="red",linewidth=1.2) +
                geom_text(aes(label=C_label),vjust=-0.75,size=5) +
                scale_fill_manual(values=c(`FALSE`="#BDBDBD",`TRUE`="#FF72BB"),guide="none") +
                scale_y_continuous(limits=c(lower_limit,upper_limit),expand=expansion(mult=c(0.08,0.22))) +
                coord_flip(clip="off") +
                labs(title=paste0(get_ctype(db)," - ",horizon,"-year ",endpoint_short[atype]),
                    x="Gene signature",y="Harrell's C-index") +
                theme_bw() +
                theme(plot.title=element_text(size=10,face="bold"),
                    axis.title=element_text(size=18,color="black"),
                    axis.text.x=element_text(size=14,color="black"),
                    axis.text.y=element_text(size=15,color="black"),
                    panel.grid=element_blank(),plot.margin=margin(10,35,10,10))

            filename=paste0(dir_cindex,safe_filename(paste(db,toupper(atype),
                paste0(horizon,"year"),"Score_only",sep="_")),".tiff")

            ggsave(filename,plot=gp,width=22,height=22,units="cm",dpi=300,
                limitsize=F,compression="lzw")
        }
    }
}



########################################
## Adjusted-model Harrell's C-index plots
########################################
for(db in unique(continuous_results$Database)) {
    for(atype in unique(continuous_results$Endpoint)) {
        for(horizon in sort(unique(continuous_results$Horizon_year))) {

            cdf=continuous_results[continuous_results$Database==db & continuous_results$Endpoint==atype &
                continuous_results$Horizon_year==horizon & continuous_results$Signature %in% plot_signatures,,drop=F]
            cdf=cdf[complete.cases(cdf[,c("C_model","C_model_low","C_model_high")]),,drop=F]
            if(nrow(cdf)==0) next

            current_order=plot_signatures[plot_signatures %in% cdf$Signature]
            cdf$Feature=display_signature(cdf$Signature)
            cdf$Feature=factor(cdf$Feature,levels=rev(display_signature(current_order)))
            cdf$C_label=sprintf("%.3f",cdf$C_model)

            cdf$Top3=FALSE
            top_index=order(cdf$C_model,decreasing=T)[seq_len(min(3,nrow(cdf)))]
            cdf$Top3[top_index]=TRUE

            c_range=range(c(cdf$C_model_low,cdf$C_model_high,0.5),na.rm=T)
            c_pad=max(0.02,diff(c_range)*0.12)
            lower_limit=max(0,c_range[1]-c_pad)
            upper_limit=min(1,c_range[2]+c_pad)

            gp=ggplot(cdf,aes(x=Feature,y=C_model)) +
                geom_errorbar(aes(ymin=C_model_low,ymax=C_model_high),width=0.2,linewidth=1.2,color="black") +
                geom_point(aes(fill=Top3),size=5,shape=21,color="black",stroke=1) +
                geom_hline(yintercept=0.5,linetype="dashed",color="red",linewidth=1.2) +
                geom_text(aes(label=C_label),vjust=-0.75,size=5) +
                scale_fill_manual(values=c(`FALSE`="#BDBDBD",`TRUE`="#FF72BB"),guide="none") +
                scale_y_continuous(limits=c(lower_limit,upper_limit),expand=expansion(mult=c(0.08,0.22))) +
                coord_flip(clip="off") +
                labs(title=paste0(get_ctype(db)," - ",horizon,"-year ",endpoint_short[atype]),
                    x="Gene signature",y="Adjusted-model Harrell's C-index") +
                theme_bw() +
                theme(plot.title=element_text(size=10,face="bold"),
                    axis.title=element_text(size=18,color="black"),
                    axis.text.x=element_text(size=14,color="black"),
                    axis.text.y=element_text(size=15,color="black"),
                    panel.grid=element_blank(),plot.margin=margin(10,35,10,10))

            filename=paste0(dir_cindex,safe_filename(paste(db,toupper(atype),
                paste0(horizon,"year"),"Clinical_factor_addition",sep="_")),".tiff")

            ggsave(filename,plot=gp,width=22,height=22,units="cm",dpi=300,
                limitsize=F,compression="lzw")
        }
    }
}



########################################
## Formal paired Delta C-index plots
########################################
delta_models=c("Score only","Clinical covariates + score")
delta_range=range(c(ddf_plot$Delta_C_low,ddf_plot$Delta_C_high,0),na.rm=T)
delta_pad=max(0.01,diff(delta_range)*0.12)
lower_limit=delta_range[1]-delta_pad
upper_limit=delta_range[2]+delta_pad
delta_breaks=sort(unique(c(delta_range[1],0,delta_range[2])))

for(model_name in delta_models) {
    for(db in unique(delta_results$Database)) {
        for(atype in unique(delta_results$Endpoint)) {
            for(horizon in sort(unique(delta_results$Horizon_year))) {

                ddf_plot=delta_results[delta_results$Database==db & delta_results$Endpoint==atype &
                    delta_results$Horizon_year==horizon & delta_results$Model==model_name &
                    delta_results$Comparator %in% setdiff(plot_signatures,"Metastasis Signature"),,drop=F]

                ddf_plot=ddf_plot[complete.cases(ddf_plot[,c("Delta_C","Delta_C_low","Delta_C_high")]),,drop=F]
                if(nrow(ddf_plot)==0) next

                current_order=setdiff(plot_signatures,"Metastasis Signature")
                current_order=current_order[current_order %in% ddf_plot$Comparator]

                ddf_plot$Feature=display_signature(ddf_plot$Comparator)
                ddf_plot$Feature=factor(ddf_plot$Feature,levels=rev(display_signature(current_order)))
                ddf_plot$Point_group=ifelse(ddf_plot$Delta_C<0,"Comparator higher",
                    ifelse(ddf_plot$Delta_C>0,"Metastasis Signature higher","No difference"))
                ddf_plot$Delta_label=sprintf("%.3f",ddf_plot$Delta_C)

                delta_range=range(c(ddf_plot$Delta_C_low,ddf_plot$Delta_C_high,0),na.rm=T)
                delta_pad=max(0.01,diff(delta_range)*0.12)
                lower_limit=delta_range[1]-delta_pad
                upper_limit=delta_range[2]+delta_pad

                model_file=ifelse(model_name=="Score only","Score_only","Clinical_adjusted")

                gp=ggplot(ddf_plot,aes(x=Feature,y=Delta_C)) +
                    geom_errorbar(aes(ymin=Delta_C_low,ymax=Delta_C_high),
                        width=0.2,linewidth=1.2,color="black") +
                    geom_point(aes(fill=Point_group),size=5,shape=21,color="black",stroke=1) +
                    geom_hline(yintercept=0,linetype="dashed",color="red",linewidth=1.2) +
                    geom_text(aes(label=Delta_label),vjust=-0.75,size=5) +
                    scale_fill_manual(values=c("Comparator higher"="#BDBDBD",
                        "Metastasis Signature higher"="#FFA5C4"),guide="none") +
                    scale_y_continuous(limits=c(lower_limit,upper_limit),
                        breaks=round(delta_breaks,1),
                        labels=scales::label_number(accuracy=0.1),
                        expand=expansion(mult=c(0.08,0.22))) +
                    coord_flip(clip="off") +
                    labs(title=paste0(get_ctype(db)," - ",horizon,"-year ",endpoint_short[atype]),
                        x="Comparator signature",
                        y=expression(Delta*"C-index (Metastasis Signature - comparator)")) +
                    theme_bw() +
                    theme(plot.title=element_text(size=10,face="bold"),
                        axis.title=element_text(size=9,color="black"),
                        axis.text.x=element_text(size=12,color="black"),
                        axis.text.y=element_text(size=15,color="black"),
                        panel.grid=element_blank(),plot.margin=margin(10,35,10,10))

                filename=paste0(dir_delta,safe_filename(paste(db,toupper(atype), paste0(horizon,"year"),model_file,sep="_")),".tiff")
                ggsave(filename,plot=gp,width=22,height=22,units="cm",dpi=300, limitsize=F,compression="lzw")
            }
        }
    }
}



########################################
## Size-matched control plots
########################################
endpoint_short=setNames(c("MFS","RFS","OS"),c("mfs","rfs","os"))

get_ctype=function(db) if(db=="Firehose") "Lung Cancer" else if(db=="GDC") "Pan Cancer" else db
safe_filename=function(x) trimws(gsub("[[:space:]]+"," ",gsub("[\\\\/:*?\"<>|]","_",x)))
signature_key=function(x) tolower(gsub("[^a-zA-Z0-9]","",x))

canonical_signature_key=function(x) {
    k=signature_key(x)
    k[grepl("known",k) & grepl("metastasis",k) & grepl("signature",k)]="knownmetastasissignature"
    k[grepl("random",k) & (grepl("express",k) | grepl("gene",k))]="randomexpressedgenes"
    k[k=="metastasissignature"]="metastasissignature"
    k
}

find_signature_name=function(x,key,fallback) {
    idx=match(key,canonical_signature_key(x))
    if(is.na(idx)) fallback else as.character(x[idx])
}

if(!exists("short_signature")) short_signature=function(x) gsub(" \\| "," from ",x)

plot_size_control=function(db,atype,horizon) {
    all_names=unique(c(plot_signatures,size_control_summary$Parent_signature,continuous_results$Signature))
    metastasis_sig=find_signature_name(all_names,"metastasissignature","Metastasis Signature")
    known_sig=find_signature_name(all_names,"knownmetastasissignature","Known metastasis signature")
    random_sig=find_signature_name(all_names,"randomexpressedgenes","Random expressed genes")

    special_keys=c("metastasissignature","knownmetastasissignature","randomexpressedgenes")
    remaining=unique(plot_signatures)
    remaining=remaining[!canonical_signature_key(remaining) %in% special_keys]
    ordered_signatures=c(metastasis_sig,known_sig,random_sig,remaining)
    ordered_signatures=ordered_signatures[!duplicated(canonical_signature_key(ordered_signatures))]
    allowed_keys=canonical_signature_key(ordered_signatures)

    sdf=size_control_summary[
        size_control_summary$Database==db &
        size_control_summary$Endpoint==atype &
        size_control_summary$Horizon_year==horizon,,drop=F]
    sdf$Signature_key=canonical_signature_key(sdf$Parent_signature)
    sdf=sdf[sdf$Signature_key %in% allowed_keys,,drop=F]
    sdf=sdf[!duplicated(sdf$Signature_key),,drop=F]

    odf=continuous_results[
        continuous_results$Database==db &
        continuous_results$Endpoint==atype &
        continuous_results$Horizon_year==horizon,,drop=F]
    odf$Signature_key=canonical_signature_key(odf$Signature)
    odf=odf[odf$Signature_key %in% allowed_keys & !is.na(odf$C_score),,drop=F]
    odf=odf[!duplicated(odf$Signature_key),,drop=F]

    plot_df=data.frame(
        Parent_signature=ordered_signatures,
        Signature_key=canonical_signature_key(ordered_signatures),
        stringsAsFactors=F
    )

    sdf_idx=match(plot_df$Signature_key,sdf$Signature_key)
    odf_idx=match(plot_df$Signature_key,odf$Signature_key)

    plot_df$Median_C=sdf$Median_C[sdf_idx]
    plot_df$Adjusted_P=sdf$Adjusted_P[sdf_idx]
    plot_df$Original_C=odf$C_score[odf_idx]
    plot_df$Original_low=odf$C_score_low[odf_idx]
    plot_df$Original_high=odf$C_score_high[odf_idx]

    metastasis_idx=match("metastasissignature",plot_df$Signature_key)
    if(!is.na(metastasis_idx) && !is.finite(plot_df$Original_C[metastasis_idx])) {
        target_values=unique(sdf$Target_C[is.finite(sdf$Target_C)])
        if(length(target_values)>0) plot_df$Original_C[metastasis_idx]=target_values[1]
    }

    plot_df=plot_df[is.finite(plot_df$Original_C) | is.finite(plot_df$Median_C),,drop=F]
    if(nrow(plot_df)==0) return(NULL)

    display_keys=unique(allowed_keys[allowed_keys %in% plot_df$Signature_key])
    plot_df$Plot_position=length(display_keys)-match(plot_df$Signature_key,display_keys)+1
    plot_df$Parent_label=short_signature(plot_df$Parent_signature)
    axis_df=plot_df[order(plot_df$Plot_position),c("Plot_position","Parent_label"),drop=F]

    original_points=plot_df[is.finite(plot_df$Original_C),,drop=F]
    matched_points=plot_df[is.finite(plot_df$Median_C),,drop=F]
    original_ci=plot_df[is.finite(plot_df$Original_low) & is.finite(plot_df$Original_high),,drop=F]
    matched_points$Significant=!is.na(matched_points$Adjusted_P) & matched_points$Adjusted_P<0.05

    all_c=c(plot_df$Original_C,plot_df$Original_low,plot_df$Original_high,plot_df$Median_C,0.5)
    all_c=all_c[is.finite(all_c)]
    if(length(all_c)==0) return(NULL)

    data_range=max(all_c)-min(all_c)
    if(data_range==0) data_range=0.1

    axis_min=max(0,floor((min(all_c)-max(0.025,data_range*0.15))*100)/100)
    axis_max=min(1,ceiling((max(all_c)+max(0.035,data_range*0.20))*100)/100)

    if(axis_max-axis_min<0.12) {
        axis_mid=(axis_min+axis_max)/2
        axis_min=max(0,axis_mid-0.06)
        axis_max=min(1,axis_mid+0.06)
    }

    axis_breaks=sort(unique(c(axis_min,0.5,axis_max)))
    signature_ticks=element_line(color=ifelse(horizon==1,"black","white"))

    gp=ggplot()

    if(nrow(original_ci)>0) {
        gp=gp+geom_errorbar(
            data=original_ci,
            aes(x=Plot_position,ymin=Original_low,ymax=Original_high),
            inherit.aes=F,width=0.2,linewidth=1.2,color="black")
    }

    if(nrow(original_points)>0) {
        gp=gp+geom_point(
            data=original_points,
            aes(x=Plot_position,y=Original_C),
            inherit.aes=F,shape=21,size=5,color="black",fill="#A8CBFF",stroke=1)
    }

    if(nrow(matched_points)>0) {
        gp=gp+geom_point(
            data=matched_points,
            aes(x=Plot_position,y=Median_C,fill=Significant),
            inherit.aes=F,shape=23,size=5,color="black",stroke=1)
    }

    gp=gp+
        geom_hline(yintercept=0.5,linetype="dashed",color="red",linewidth=1.2)+
        scale_x_continuous(
            limits=c(0.5,length(display_keys)+0.5),
            breaks=axis_df$Plot_position,
            labels=axis_df$Parent_label,
            expand=c(0,0))+
        scale_fill_manual(values=c(`FALSE`="#ffdf50",`TRUE`="#3e3e3e"),guide="none")+
        scale_y_continuous(
            limits=c(axis_min,axis_max),
            breaks=axis_breaks,
            labels=function(x) sprintf("%.1f",x),
            expand=expansion(mult=c(0.05,0.10)))+
        coord_flip(clip="off")+
        labs(
            title=paste0(get_ctype(db)," - ",horizon,"-year ",endpoint_short[atype]),
            x="Signature",y="Harrell's C-index")+
        theme_bw()+
        theme(
            plot.title=element_text(size=10,face="bold"),
            axis.title=element_text(size=18,color="black"),
            axis.text.x=element_text(size=13,color="black"),
            axis.text.y=element_text(size=15,color="black"),
            axis.ticks.y=signature_ticks,
            legend.position="none",
            panel.grid=element_blank(),
            plot.margin=margin(10,35,10,10))

    filename=paste0(dir_size, safe_filename(paste(db,toupper(atype),paste0(horizon,"year"),sep="_")),".tiff")
    ggsave(filename,plot=gp,width=22,height=22,units="cm",dpi=300,limitsize=F,compression="lzw")
    gp
}

for(db in unique(size_control_summary$Database)) {
    for(atype in unique(size_control_summary$Endpoint)) {
        for(horizon in sort(unique(size_control_summary$Horizon_year))) {
            gp=plot_size_control(db,atype,horizon)
            if(!is.null(gp)) print(gp)
        }
    }
}




########################################
## Exploratory top/bottom-quartile KM analysis
########################################

## Function
fit_quartile_survival=function(dat, score_name, adjust_vars=adjust_template) {
    if(!score_name %in% colnames(dat)) return(NULL)

    qdat=dat[complete.cases(dat[,c("time","event","day",score_name),drop=F]),,drop=F]
    if(nrow(qdat)==0) return(NULL)

    cuts=quantile(qdat[[score_name]], probs=c(0.25,0.75), na.rm=T, names=F)
    if(any(is.na(cuts)) || cuts[1]>=cuts[2]) return(NULL)

    qdat$Level=NA_character_
    qdat$Level[qdat[[score_name]]<=cuts[1]]="Low"
    qdat$Level[qdat[[score_name]]>=cuts[2]]="High"
    qdat=qdat[!is.na(qdat$Level),,drop=F]
    qdat$Level=factor(qdat$Level, levels=c("Low","High"))

    level_table=table(qdat$Level)
    if(length(level_table)<2 || any(level_table<3)) return(NULL)
    if(length(unique(qdat$event))<2) return(NULL)
    if(sum(qdat$event)<3 || sum(qdat$event==0)<3) return(NULL)

    sres=survfit(Surv(time,event)~Level, data=qdat)
    logrank=survdiff(Surv(time,event)~Level, data=qdat)
    pv=pchisq(logrank$chisq, df=length(logrank$n)-1, lower.tail=F)

    ## Adjusted top-versus-bottom HR
    adjust_vars=intersect(adjust_vars, colnames(qdat))

    repeat {
        cox_vars=c("Level",adjust_vars)
        vars=c("time","event",cox_vars)
        cox_dat=qdat[complete.cases(qdat[,vars,drop=F]),vars,drop=F]

        if(nrow(cox_dat)==0) break

        invalid=adjust_vars[vapply(
            adjust_vars,
            function(v) length(unique(cox_dat[[v]]))<2,
            logical(1)
        )]

        if(length(invalid)==0) break
        adjust_vars=setdiff(adjust_vars,invalid)
    }

    hr=NA_real_
    hr_low=NA_real_
    hr_high=NA_real_
    cox=NULL

    if(exists("cox_dat") && nrow(cox_dat)>0) {
        categorical_vars=intersect(categorical_template,adjust_vars)

        for(v in categorical_vars) {
            cox_dat[[v]]=droplevels(factor(cox_dat[[v]]))
        }

        cox_dat$Level=relevel(factor(cox_dat$Level),ref="Low")

        if(length(unique(cox_dat$Level))==2 &&
           length(unique(cox_dat$event))==2 &&
           sum(cox_dat$event)>=3) {

            rhs=c("Level",quote_term(adjust_vars))
            cox_formula=as.formula(
                paste("Surv(time,event) ~",paste(rhs,collapse=" + "))
            )

            cox=tryCatch(
                coxph(cox_formula,data=cox_dat,ties="efron"),
                error=function(e) NULL
            )

            if(!is.null(cox)) {
                csum=summary(cox)

                if("LevelHigh" %in% rownames(csum$coefficients)) {
                    hr=csum$coefficients["LevelHigh","exp(coef)"]
                    hr_low=csum$conf.int["LevelHigh","lower .95"]
                    hr_high=csum$conf.int["LevelHigh","upper .95"]
                }
            }
        }
    }

    list( data=qdat, survfit=sres, logrank_p=pv, cox=cox, hr=hr, hr_low=hr_low, hr_high=hr_high, cut_low=cuts[1], cut_high=cuts[2] )
}


## Viz
km_cols=setNames(c("#403dff","#ff3d3d"),c("Low","High"))

for(db in names(mdfl)) {
    for(atype in atypes) {
        for(horizon in horizons) {
            dat=make_endpoint_data(db=db,atype=atype,horizon_year=horizon)
            score_names=plot_signatures[plot_signatures %in% colnames(dat)]

            for(score_name in score_names) {
                kres=fit_quartile_survival(dat=dat,score_name=score_name)
                if(is.null(kres)) next

                qdat=kres$data
                sres=kres$survfit
                pv_label=format_p(kres$logrank_p)

                if(is.na(kres$hr)) {
                    hr_label="NA"
                } else {
                    hr_label=ifelse(kres$hr>100,'HR > 100',round(kres$hr,3))
                }

                dcut=365.25*horizon

                if(horizon==1) {
                    xlab="Months"
                    xscale_value=365.25/12
                    break_value=365.25/12
                } else {
                    xlab="Years"
                    xscale_value="d_y"
                    break_value=365.25
                }

                ymin=max(0,floor(min(sres$surv,na.rm=T)*10)/10-0.05)
                annotation_y=ymin+0.12*(1-ymin)

                plot_title=paste0( get_ctype(db)," | ", display_signature(score_name)," | ", horizon,"-year ",endpoint_short[atype] )

                gplot=ggsurvplot(
                    sres, data=qdat, conf.int=T, risk.table=T,
                    legend.labs=c("Low","High"),
                    legend.title="Metastasis risk level",
                    palette=unname(km_cols[c("Low","High")]),
                    title=plot_title, xlab=xlab, xscale=xscale_value,
                    break.time.by=break_value, risk.table.height=0.3,
                    surv.median.line="hv", size=0.65, fontsize=4,
                    censor=T, censor.shape="+", censor.size=5,
                    ylim=c(ymin,1), xlim=c(0,dcut)
                )

                gplot$plot=gplot$plot +
                    labs(y=paste0(endpoint_labels[atype]," probability")) +
                    annotate( "text", label=paste0( "Log-rank P = ",pv_label, "\nAdjusted HR = ",hr_label), x=0.03*dcut, y=annotation_y, hjust=0, color="black", size=5, fontface="italic" ) +
                    theme(
                        plot.title=element_text(size=13,face="bold",hjust=0.5),
                        axis.title=element_text(size=13),
                        axis.text=element_text(size=12.5),
                        panel.grid=element_blank()
                    )

                gplot$table=gplot$table +
                    labs(title="Number at risk",y="Risk level") +
                    theme(
                        plot.title=element_text(size=11,face="italic"),
                        axis.title=element_text(size=11),
                        axis.text=element_text(size=10)
                    )

                file_base=safe_filename(paste(db,toupper(atype),paste0(horizon,"year"),score_name,sep="_"))

                # With risk table
                tiff( filename=paste0(dir_km_table,file_base,".tiff"), width=15, height=17, units="cm", res=300, compression="lzw" )
                print(gplot)
                dev.off()

                # Without risk table
                ggsave( paste0(dir_km_notable,file_base,".tiff"), plot=gplot$plot, width=14, height=12, units="cm", dpi=300, limitsize=F, compression="lzw" )
            }
        }
    }
}



########################################
## Expression levels of the Metastasis signature in Patient samples
########################################

## Load TCGA (GDC) LUAD origina dataset
load(paste0(dir_tcga,"clinical/LUAD.Rdata")) # sample information; variable name: sinfo
load(paste0(dir_tcga,"tpm/LUAD.Rdata")) # TPM expression; variable name: tpm


## DEG analysis in patient samples
fcut=log2(2)
pcut=0.05

mids=intersect(sinfo$patient[grep("Tumor",sinfo$definition)] , sinfo$patient[grep("Normal",sinfo$definition)])

sinfo=sinfo[order(sinfo$days_to_last_follow_up),]
sinfo=sinfo[order(sinfo$days_to_collection),]

csinfo=sinfo[sinfo$patient %in% mids & grepl("Normal",sinfo$definition),]
cid=sort(csinfo[!duplicated(csinfo$patient),"barcode"])

tsinfo=sinfo[sinfo$patient %in% mids & grepl("Tumor",sinfo$definition),]
tid=sort(tsinfo[!duplicated(tsinfo$patient),"barcode"])

ddf=data.frame(Gene=rownames(tpm),stat=NA,pvalue=NA,padj=NA,deg=NA)
for(i in 1:nrow(tpm)) {
      print(i)
      gene=rownames(tpm)[i]

      wres=wilcox.test(tpm[i,tid],tpm[i,cid],paired=T,alternative="greater")
      ddf$stat[ddf$Gene==gene]=wres$stat
      ddf$pvalue[ddf$Gene==gene]=wres$p.value
}
ddf$padj=p.adjust(ddf$pvalue,method="BH")
ddf$deg=ifelse(!is.na(ddf$padj) & ddf$padj<pcut,T,F)
ddf=ddf[order(ddf$pvalue),]
ddf=ddf[order(ddf$padj),]
# save(ddf, file=paste0(dir_data,"3. LUAD Patient sample DEG analysis res.Rdata"))


## Metastasis biomarker
ddf[ddf$Gene %in% fdegs,]
#           Gene stat       pvalue         padj   deg
# 17888  SLC12A8 1630 1.032706e-09 2.607092e-08  TRUE
# 963       NGEF 1565 2.017366e-08 2.889870e-07  TRUE
# 12615    NPAS2 1528 9.812349e-08 1.114578e-06  TRUE
# 17463 SERPINB5 1468 1.077503e-06 9.199238e-06  TRUE
# 14132    OPLAH 1281 5.000512e-04 2.291309e-03  TRUE
# 6573     LPIN3 1121 2.009833e-02 6.103157e-02 FALSE
# 15823     DNER  892 3.902281e-01 7.794804e-01 FALSE
# 16903   SFT2D1  821 6.067976e-01 1.000000e+00 FALSE
# 2071      NRP1  765 7.594560e-01 1.000000e+00 FALSE
# 999   TNFRSF1A  751 7.918756e-01 1.000000e+00 FALSE
# 5946      RAC2  490 9.976994e-01 1.000000e+00 FALSE
# 15875  PRKAR1B  369 9.999185e-01 1.000000e+00 FALSE
# 6067     KLK10  157 1.000000e+00 1.000000e+00 FALSE

## Viz by boxplot
degs=intersect(ddf$Gene[ddf$deg==T],fdegs)

tiff(filename=paste0(dir_fig, "Supplementary Figure 3(G).tif"),  width= 27, height = 27 , units = 'cm', res=300)
par(mfrow=c(3,length(degs)),mar=c(2,2.2,2.7,0.5),mgp=c(1.1,0.3,0),cex.lab=1,cex.axis=1.7,cex.main=2,tck=-0.02,las=1,bty="l")
for(g in c(degs,sort(setdiff(fdegs,degs)))) {
      ne=tpm[g,cid] # normal expression
      de=tpm[g,tid] # disease expression
      names(ne)=substr(names(ne),1,12)
      names(de)=substr(names(de),1,12)
      de=de[match(names(ne),names(de))] # disease expression value ordering by TCGA
      tt=wilcox.test(de, ne, alternative="greater", paired=T)
      pv=ddf$padj[ddf$Gene==g]
      df=cbind(data.frame(normal=log2(ne+1)), data.frame(cancer=log2(de+1)))
      ymin=0
      ymax=4
      ylim=c(ymin,ymax)
      boxplot(df$normal,df$cancer, xlab="", ylab="", ylim=ylim, outline=F, col="white", border="white", boxwex=0.6, frame=T, xaxt="n")
      grid(NA, NULL, lty=3, lwd=1, col="darkgrey")
      par(new=T)
      boxplot(main=g, df$normal,df$cancer, xlab="", ylab="", ylim=ylim, outline=F, col=NA, border="black", boxwex=0.6, frame=T, xaxt="n")
      for(i in 1:ncol(df)) {stripchart(df[,i], at=i, cex=1.5, lwd=0.3, pch = 21, col=c("#4dc7ff","#ff478d")[i], bg=c("#82d7ff","#ff81b1")[i], method="stack", vertical=TRUE, add = TRUE)}
      segments(x0=rep(1,nrow(df)), x1=rep(2,nrow(df)), y0=df$normal, y1=df$cancer, col=rgb(169/255,169/255,169/255,alpha=0.7), lwd=1.1)
      if(pv < 0.05) {
      exp = ceiling(log10(pv))-1
      numord = round(pv* (10^-exp), digits=2)
      text(x=1.5, y=ymax*9.42/10, bquote(.(numord)~"x"~10^.(exp)), cex=2.1, adj=0.5,xpd=T)
      } else {
      pv=round(pv,2)
      text(x=1.5, y=ymax*9.42/10, bquote(.(pv) ), cex=2.1, adj=0.5,xpd=T)
      } 
      segments(x0=1,x1=2,y0=ymax*8.7/10,y1=ymax*8.7/10,col='black', xpd=T)
      segments(x0=c(1,2),x1=c(1,2),y0=c(ymax*8.5/10,ymax*8.5/10),y1=c(ymax*8.7/10,ymax*8.7/10),col='black', xpd=T)
}
dev.off()



########################################
## overlap among the Gene sets
########################################

## Gene list
mlist=c(gsl["Metastasis Signature"] , list(`Known Signature`=unique(unlist(gsl[names(gsl) %in% setdiff(gs_viz,"Metastasis Signature")])) , `Normal vs. Tumor DEGs`=ddf$Gene[ddf$padj < pcut & !is.na(ddf$padj)]))


## Viz
tiff(filename=paste0(dir_fig, "Main Figure 3(G).tif"),  width= 20, height = 20 , units = 'cm', res=300)
par(mfrow=c(3,length(degs)),mar=c(2,2.2,2.7,0.5),mgp=c(1.1,0.3,0),cex.lab=1,cex.axis=1.7,cex.main=2,tck=-0.02,las=1,bty="l")
vg=venndir(mlist, overlap_type="overlap", proportional=F, plot_style="gg", show_segments=F, font_cex=c(1.2,1.2),label_style="lite",set_colors=c('#FF5D5D','#EDFFC3','#FF81B1'))
dev.off()