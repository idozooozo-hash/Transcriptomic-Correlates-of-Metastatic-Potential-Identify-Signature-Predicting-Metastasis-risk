########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/figure/"
dir_msigdb="G:/My Drive/DB/MSigDB/homosapiens/" # MSigDB gene set directory

## Additional gene sets directory
dir_csea="G:/My Drive/DB/cancerSEA/" # cancerSEA datasets directory



########################################
## Library
########################################
library(ComplexHeatmap)
library(car)
library(edgeR)
library(data.table)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(dplyr)
library(colorRamp2)
library(huge)
library(igraph)
library(RCy3)
library(GSVA)
library(FSA)
library(stringr)

## Python library
library(reticulate)
use_python( "C:/Users/idozo/anaconda3/envs/metastatsis/python.exe", required = TRUE )
umap = import("umap") # conda install conda-forge::umap-learn
umap = umap$UMAP



########################################
## Load Data
########################################
load(file=paste0(dir_data,"0. tpm(DepMap).Rdata")) # variable name: tpm
load(file=paste0(dir_data,"0. read count(DepMap).Rdata")) # variable name: rcm
load(file=paste0(dir_data,"0. filtered sample information(DepMap).Rdata")) # variable name: sinfo

load(paste0(dir_msigdb,"hm_2024.1.Rdata")) # load MSigDB '2024.1 version' Hallmark gene set; variable name: hm
class(hm) # list
unique(lapply(hm,class)) # character
length(hm) # 50



########################################
## UMAP Visualization of transcriptomic profiles
########################################

## input
sinfo=sinfo[match(colnames(tpm),sinfo$Sample),]


## Dimension Reduction
random_state=123
set.seed(random_state)
ures=umap(n_components = as.integer(2), random_state=as.integer(random_state))$fit_transform(t(tpm))
rownames(ures)=colnames(tpm)
colnames(ures)=c("UMAP 1","UMAP 2")
# save(ures, file=paste0(dir_data,"2. Dimension Reduction res.Rdata"))


## Colored by Cancer Type
pcols=c('Esophagus or Stomach'='#f08120ff', 'Lung'='#ff311aff', 'Thyroid'='#F1948A', 'CNS or Brain'='#a8372dff', 'Bladder or Urinary Tract'='#F4D03F', 'Ovary or Fallopian Tube'='#85C1E9', 'Bowel'='#8E44AD', 'Uterus'='#0e99f7ff', 'Liver'='#ec991bff', 'Pancreas'='#00ff6aff', 'Breast'='#ff76e8ff', 'Head and Neck'='#1bb15aff', 'Soft Tissue'='#13cea5ff', 'Skin'='#ff4784ff', 'Pleura'='#b5c2baff', 'Bone'='#c38dd8ff', 'Prostate'='#6e6e6eff', 'Kidney'='#dd7e75ff', 'Peripheral Nervous System'='#F1948A', 'Biliary Tract'='#135079ff', 'Fibroblast'='#000000ff', 'Ampulla of Vater'='#82E5AA')

tiff(filename=paste0(dir_fig,"Supplementary Figure 2(A).tif"), width=10, height=10, units = 'cm',res=300)
par(mfcol=c(1,1), mar=c(2.5,2.5,1.5,0.5), mgp=c(1.2,0.1,0), tck=-0.005)
plot(ures, main="UMAP (Colored by Cancer Type)", bg=pcols[sinfo[match(rownames(ures),sinfo$Sample),"Cancer Type"]], col='#b9b9b9ff', pch=21, lwd=0.7, cex=0.9, cex.main=0.9)
dev.off()

lg=Legend(labels=sort(names(pcols)), legend_gp = gpar(fill=pcols[order(names(pcols))]),border="#b9b9b9ff")
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(A)_legend.tiff"), width=10, height=10, units = 'cm',res=300)
draw(packLegend(list=list(lg))) 
dev.off() 


## Colored by Metastatic potential category
pcols=setNames(c('#ff3d3d','#d4d4d4','#403dff') , unique(sinfo$`meta type to all5`))

tiff(filename=paste0(dir_fig,"Supplementary Figure 2(B).tif"), width=10, height=10, units = 'cm',res=300)
par(mfcol=c(1,1), mar=c(2.5,2.5,1.5,0.5), mgp=c(1.2,0.1,0), tck=-0.005)
plot(ures, main="UMAP (Colored by Metastatic Potential Type)", bg=pcols[sinfo[match(rownames(ures),sinfo$Sample),"meta type to all5"]], col='#b9b9b9ff', pch=21, lwd=0.7, cex=0.9, cex.main=0.9)
dev.off()

lg=Legend(labels=sort(names(pcols)), legend_gp = gpar(fill=pcols[order(names(pcols))]),border="#b9b9b9ff")
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(B)_legend.tiff"), width=10, height=10, units = 'cm',res=300)
draw(packLegend(list=list(lg))) 
dev.off()



########################################
## Metastasis Potential by Subtype in Lung cancer
########################################
psinfo=sinfo
psinfo=psinfo[psinfo$`Cancer Type`=='Lung',]
psinfo=psinfo[psinfo$`meta type to all5`!='Weakly Metastatic (Low Confidence)',]


## Subtype composition table
subtypes = c(
    "Lung Adenocarcinoma" = "LUAD",
    "Giant Cell Carcinoma of the Lung" = "Giant Cell Carcinoma",
    "Large Cell Lung Carcinoma" = "Large Cell Carcinoma",
    "Lung Squamous Cell Carcinoma" = "LUSC",
    "Non-Small Cell Lung Cancer" = "NSCLC",
    "Small Cell Lung Cancer" = "SCLC",
    "Lung Adenosquamous Carcinoma" = "Adenosquamous Carcinoma",
    "Lung Carcinoid" = "Lung Carcinoid",
    "Mucoepidermoid Carcinoma of the Lung" = "Mucoepidermoid Carcinoma",
    "SMARCA4-deficient undifferentiated tumor" = "SMARCA4 DUT"
)
psinfo$`Cancer Subtype`=subtypes[match(psinfo$`Cancer Subtype`,names(subtypes))]
tb=table(psinfo[,c("Cancer Subtype","meta type to all5")])
tb
#                           meta type to all5
# Cancer Subtype             Metastatic Non Metastatic
#   Giant Cell Carcinoma              1              0
#   Large Cell Carcinoma              7              0
#   LUAD                             28              2
#   LUSC                             12              0
#   Mucoepidermoid Carcinoma          1              0
#   NSCLC                             3              0
#   SCLC                              2              1
#   SMARCA4 DUT                       1              0


## Fisher exact test by subtype
set.seed(456)
tb = tb[rowSums(tb) > 0, , drop = FALSE]
fres = fisher.test(tb, simulate.p.value = TRUE, B = 1000000)
fres
# data:  tb
# p-value = 0.4554
# alternative hypothesis: two.sided


## Anova test by subtype
ldf=psinfo[,c("meta to all5","meta type to all5","Cancer Subtype")] # data frame for linear modeling
ldf=ldf[!is.na(ldf[,"Cancer Subtype"]) & ldf[,"Cancer Subtype"]!="Unknown" & ldf[,"Cancer Subtype"]!="",]
ldf[,"Cancer Subtype"]=factor(ldf[,"Cancer Subtype"])
lm=lm(ldf[,"meta to all5"]~ldf[,"Cancer Subtype"]) # linear modeling
ares=Anova(lm, type=2)
pv=ares[1,"Pr(>F)"] 
pv # 0.03421726


## Kruskal wallis test by subtype
kres=kruskal.test(`meta to all5`~`Cancer Subtype`, data=ldf)
pv=kres$p.value 
pv # 0.04651811
# Kruskal wallis test - <post-hoc test>
ntab=table(ldf$`Cancer Subtype`)
keep_subtype=names(ntab[ntab>=3])
ldf2=ldf[ldf$`Cancer Subtype` %in% keep_subtype,]
ldf2$`Cancer Subtype`=droplevels(ldf2$`Cancer Subtype`)
kruskal.test(`meta to all5`~`Cancer Subtype`,data=ldf2)
dres2=dunnTest(`meta to all5`~`Cancer Subtype`, data=ldf2, method="holm")$res
dres2
#                      Comparison           Z     P.unadj      P.adj
# 1   SCLC - Large Cell Carcinoma -0.83992446 0.400950742 1.00000000
# 2                  SCLC - NSCLC  0.66254135 0.507624344 1.00000000
# 3  Large Cell Carcinoma - NSCLC  1.62385396 0.104406981 0.83525585
# 4                   SCLC - LUAD  1.05830052 0.289918454 1.00000000
# 5   Large Cell Carcinoma - LUAD  2.90752663 0.003642993 0.03642993
# 6                  NSCLC - LUAD  0.16492995 0.868999117 1.00000000
# 7                   SCLC - LUSC  0.93475464 0.349914762 1.00000000
# 8   Large Cell Carcinoma - LUSC  2.48737855 0.012868838 0.11581955
# 9                  NSCLC - LUSC  0.09669876 0.922965628 0.92296563
# 10                  LUAD - LUSC -0.10964608 0.912690062 1.00000000

## Save
# save(fres, lm, ares, kres, dres2, file=paste0(dir_data,"2. Metastatic Potential Differences Across Lung Cancer Subtypes.RData"))



########################################
## Lung sample filtering
########################################

sinfo=sinfo[sinfo$`Cancer Type`=='Lung',]
nrow(sinfo) # No. of samples: 93



########################################
## Differentially expressed gene analysis
########################################

## DEG analysis
fcut=log2(2)
pcut=0.05

cid=sinfo$Sample[sinfo$`meta type to all5`=="Non Metastatic"]
tid=sinfo$Sample[sinfo$`meta type to all5`=="Metastatic"]
dg = rbind(data.frame(sample = cid, group = "control"), data.frame(sample = tid, group = "treat"))
counts = rcm[, dg$sample]
group = factor(dg$group)
dge = DGEList(counts=counts, group=group)
keep = filterByExpr(dge)
dge = dge[keep,,keep.lib.sizes=FALSE]
dge = calcNormFactors(dge)
design = model.matrix(~group) 
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
lrt = glmLRT(fit)
ddf = topTags(lrt, n=Inf)$table
colnames(ddf) = c("log2FoldChange","logCPM","LR","pvalue","padj")
ddf$Gene = gsub("_[0-9]+$", "", rownames(ddf))
# sort
ddf=ddf[order(abs(ddf$log2FoldChange), decreasing=T),]
ddf=ddf[order(ddf$pvalue),]
ddf=ddf[order(ddf$padj),]
# DEG
ddf$deg=(abs(ddf$log2FoldChange)>fcut & ddf$padj<pcut & !is.na(ddf$log2FoldChange) & !is.na(ddf$padj))
# save(ddf, file=paste0(dir_data,"2. DEG analysis res.Rdata"))


## Reference set
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


## Viz
fcut=log2(2)
pcut=0.05
tcols=setNames(c("#ff0000","#000dff","#000000","#d4d4d4"),c("Up","Down","DEG","Gene"))
pcols=setNames(c("#FF3F3F","#3E5BFF","#000000","#d4d4d4"),c("Up","Down","DEG","Gene"))
pfils=setNames(c("#FF6060","#548AFF","#000000","#d4d4d4"),c("Up","Down","DEG","Gene"))

pdf=data.frame(gene=ddf$Gene, logfc=ddf$log2FoldChange, logp=-log10(ddf$padj), pv=ddf$pvalue, deg=ddf$deg)
pdf$stat = sign(pdf$logfc) * (-log10(pdf$pv))
pdf=pdf[order(pdf$logp * sign(pdf$logfc), decreasing=T),]
pdf$rank=1:nrow(pdf)
pdf$stat=ifelse(abs(pdf$stat)>=7, 7*sign(pdf$stat), pdf$stat)

xcut = sort(c(max(pdf$rank[pdf$deg == TRUE & pdf$stat > 0], na.rm=TRUE), min(pdf$rank[pdf$deg == TRUE & pdf$stat < 0], na.rm=TRUE)))

up_n = sum(pdf$logfc > fcut & pdf$logp > log10(1/pcut), na.rm=TRUE)
down_n = sum(pdf$logfc < -fcut & pdf$logp > log10(1/pcut), na.rm=TRUE)

fpdf=pdf[pdf$deg==T,]
fpdf=fpdf[c(which(fpdf$gene %in% unlist(csea))[1:10],which(order(fpdf$rank,decreasing=T) %in% c(1:10))) , ] # Only Up/Down Top 10 
fpdf=fpdf[!is.na(fpdf$gene),]

gp=ggplot(pdf, aes(x=rank, y = stat)) +
      theme_bw() +
      theme(plot.title = element_blank(),
            axis.title.x=element_text(color="black", size=18, margin=margin(t=10)),
            axis.title.y=element_text(color="black", size=18, margin=margin(r=10)),
            axis.text.x=element_text(color="black", size=18),
            axis.text.y=element_text(color="black", size=18),
            axis.ticks.y = element_line(linewidth = 0.5, color="black"), 
            axis.ticks.x = element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_line(color="black", linewidth=0.5, linetype=3),
            panel.grid.minor.y=element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            plot.margin = unit(c(0.2,1,0.4,1), "cm")) +
      labs(x="Rank", y="sign(log2FC) × -log10(p-value)") +
      scale_y_continuous( breaks=sort(unique(c(seq(-6, 6, by=2), -1, 1))), expand=expansion(mult=c(0.05, 0.05)) ) +
      scale_x_continuous( breaks=xcut, expand=expansion(mult=c(0.05, 0.05)) ) +
      geom_point(data=pdf[pdf$logp<=log10(1/pcut),], shape = 21, col = "darkgrey", fill = pcols["Gene"], alpha = 0.5, size = 3,  stroke = 0.2) +  
      geom_point(data=pdf[pdf$logfc>fcut & pdf$deg==T,], shape = 21, col = pcols["Up"], fill = pfils["Up"], alpha = 1, size = 3,  stroke = 0.2) +
      geom_point(data=pdf[pdf$logfc<(fcut*-1) & pdf$deg==T,], shape = 21, col = pcols["Down"], fill = pfils["Down"], alpha = 1, size = 3,  stroke = 0.2) +
      coord_cartesian(ylim=c(-7,7), clip="off") +
      geom_hline(yintercept=c(-log10(pcut), log10(pcut)), linetype='solid', color='black', linewidth=0.6) +
      geom_vline( xintercept=xcut, linetype="solid", color="black", linewidth=0.6 ) +
      geom_text_repel(data=fpdf , aes(label = gene), size=5.5, color = c(rep(tcols["Up"],sum(fpdf$stat>0)),rep(tcols["Down"],sum(fpdf$stat<0))), segment.color=c(rep(pcols["Up"],sum(fpdf$stat>0)),rep(pcols["Down"],sum(fpdf$stat<0))), segment.size=0.4, box.padding=0.55, force=3.5, max.overlaps=Inf) +
      annotate("text", x=Inf, y=6.7, label=paste0("Up DEG: ", up_n), color=tcols["Up"], size=6.5, fontface="bold", hjust=1) +
      annotate("text", x=Inf, y=5.7, label=paste0("Down DEG: ", down_n), color=tcols["Down"], size=6.5, fontface="bold", hjust=1)
# ggsave(filename=paste0(dir_fig, "Main Figure 2(A).tiff"), plot = gp, width=16, height=14, units = "cm", dpi = 300, limitsize = TRUE)



########################################
## Pathway analysis
########################################

if(!("stat" %in% colnames(ddf))) {ddf$stat = sign(ddf$log2FoldChange)*-log10(ddf$pvalue)} 
st=setNames(ddf$stat, ddf$Gene)
st=sort(st[!is.na(st)], decreasing=T)  

stt=Sys.time()
gdf=as.data.frame(fgseaMultilevel(pathways=hm, stats=st, minSize=10, maxSize=500, nPermSimple=100000))
gdf=gdf[order(abs(gdf$NES), decreasing=T),]
gdf=gdf[order(gdf$pval),]
gdf=gdf[order(gdf$padj),]
gdf=gdf[!duplicated(gdf$pathway),]
gdf$score=-log10(gdf$padj)
ett=Sys.time()
ett-stt # 29.1841 secs
# save(gdf, file=paste0(dir_data,"2. Pathway analysis res.Rdata")) 


## Viz (Up-regulation Pathway)
pcut=0.001
pdf=filter(gdf,NES>0)[1:15,]

dir=rev(sign(pdf$NES))
score=rev(pdf$score*sign(pdf$NES))
path=rev(pdf$pathway)
pcol="#ff3d3d"
xmax=ceiling(max(abs(pdf$score),na.rm=T))
xrange=c(xmax*(-1),xmax)
xax=seq(0,xrange[2],ifelse(xmax<15,2,4))
tpos=xrange[2]*0.02

tiff(filename=paste0(dir_fig,"Main Figure 2(B)_Up.tif"), width=20, height=10, units="cm", res=300)
par(plt=c(0.1,0.1,0.9,0.9), mar=c(2.5,19.5,2,0.5), mgp=c(1.5,0.2,0), tck=-0.01) 
bp=barplot(score, xlim=xrange, horiz=T, xaxt='n', yaxt='n', xlab=bquote(-log[10]~("FDR")), names.arg=NA, main="TOP 15 Up-Regulated Pathways", width=1.7, space=0.3, border=NA, col=pcol)
axis(1, at=xax, labels=xax, cex.axis=0.85, las=1)
abline(v=0, lty=1)
segments(xax, (0-(bp[2]-bp[1])), xax, (max(bp)+bp[length(bp)-1]), col="#00000031", lty=3)
text(x=((-1)*tpos), y=bp[dir>0]+0.01, labels=path[dir>0], col="black", xpd=T, cex=0.72, adj=1)
dev.off()


## Viz (Down-regulation Pathway)
pcut=0.001
pdf=filter(gdf,NES<0)[1:15,]

score=rev(pdf$score)
path=rev(pdf$pathway)
pcol="#a4a4a4"
xmax=ceiling(max(abs(pdf$score),na.rm=T))
xrange=c(xmax*(-1),xmax)
xax=seq(0,xrange[2],ifelse(xmax<15,2,4))
tpos=xrange[2]*0.02

tiff(filename=paste0(dir_fig,"Main Figure 2(B)_Down.tif"), width=20, height=10, units="cm", res=300)
par(plt=c(0.1,0.1,0.9,0.9), mar=c(2.5,19.5,2,0.5), mgp=c(1.5,0.2,0), tck=-0.01) 
bp=barplot(score, xlim=xrange, horiz=T, xaxt='n', yaxt='n', xlab=bquote(-log[10]~italic("(FDR)")), names.arg=NA, main="Down-Regulated Pathways", width=1.7, space=0.3, border=NA, col=pcol)
axis(1, at=xax, labels=xax, cex.axis=0.85, las=1)
abline(v=0, lty=1)
segments(xax, (0-(bp[2]-bp[1])), xax, (max(bp)+bp[length(bp)-1]), col="#00000031", lty=3)
text(x=((-1)*tpos), y=bp[dir>0]+0.01, labels=path[dir>0], col="black", xpd=T, cex=0.72, adj=1)
dev.off() 



########################################
## Activation of EMT-related pathways
########################################
sets=c("EPITHELIAL MESENCHYMAL TRANSITION","TGF BETA SIGNALING","ANGIOGENESIS") # EMT related gene sets

for(set in sets) {
    # data
    gp=plotEnrichment(hm[[set]], st, ticksSize=0.7)
    pdf=gp$data
    # input
    xseq=seq(0,length(st),1)
    xrange=c(0,max(xseq))
    yseq=unique(c(seq(floor(min(pdf$ES)*10)*0.1,0,0.1),seq(0,ceiling(max(pdf$ES)*10)*0.1,0.1)))
    yrange=c(min(yseq),max(yseq))
    pv=ifelse(round(gdf$padj[gdf$pathway==set],3)>=0.0005,round(gdf$padj[gdf$pathway==set],3),formatC(gdf$padj[gdf$pathway==set], format = "e", digits = 4))
    # color
      color_lim=quantile(st, c(0.02,0.98), na.rm=T)
      st_color=pmax(pmin(st,color_lim[2]),color_lim[1])
    pcols=colorRamp2(c(color_lim[1],0,color_lim[2]),c('#000dff','#FFFFFF','#ff0000'))(st_color)

    # draw
    dir_fout=paste0(dir_fig,"Main Figure 2(C)_",set)
    tiff(filename=paste0(dir_fout,"Enrichment Plot (",set,").tif"), width=13, height=10, units="cm", res=300)
    par(plt=c(0.1,0.1,0.9,0.9), mar=c(5,2,2,0.5), mgp=c(1.5,0.2,0), tck=-0.01) 
    # plot
    plot(pdf$rank, pdf$ES, typ='l', col='#ffc400', xlim=xrange, ylim=yrange, main=set, lwd=4, xaxt='n', yaxt='n', xlab=NA,ylab=NA,frame.plot=F)
    segments(x0=1:max(xrange), x1=1:max(xrange), y0=min(yrange)-abs(max(yrange)-min(yrange))*0.12, y1=min(yrange)-abs(max(yrange)-min(yrange))*0.12*2, col=pcols, xpd=T)
    segments(x0=pdf$rank, x1=pdf$rank, y0=min(yrange)-abs(max(yrange)-min(yrange))*0.001, y1=min(yrange)-abs(max(yrange)-min(yrange))*0.12, col='black', xpd=T, lwd=1)
    # border
    segments(x0=min(xrange),x1=max(xrange),y0=yseq,y1=yseq,lty=2,col='grey',lwd=1.2)
    segments(x0=c(1,max(xrange),1,1,1,1),x1=c(1,max(xrange),max(xrange),max(xrange),max(xrange),max(xrange)),y0=c(min(yrange)-abs(max(yrange)-min(yrange))*0.12*2,min(yrange)-abs(max(yrange)-min(yrange))*0.12*2,min(yrange)-abs(max(yrange)-min(yrange))*0.12*2,max(yrange),min(yrange)-abs(max(yrange)-min(yrange))*0.12,min(yrange)),y1=c(max(yrange),max(yrange),min(yrange)-abs(max(yrange)-min(yrange))*0.12*2,max(yrange),min(yrange)-abs(max(yrange)-min(yrange))*0.12,min(yrange)),col='black', xpd=T, lwd=1.2)
    segments(x0=min(xrange),x1=max(xrange),y0=0,y1=0,col='#d92847',lty=2,lwd=1.2)
    # annotation text
    text(labels=c(paste0("FDR: ",pv,"\n","NES: ",round(gdf$NES[gdf$pathway==set],3))), x=max(xrange)-(max(xrange)-min(xrange))*0.05, y=max(yrange)-(max(yrange)-min(yrange))*0.15, adj=1, cex=1.3)
    # under plot
    text(labels="Metastatic", x=c(min(xrange)+(max(xrange)-min(xrange))*0.01), y=min(yrange)-abs(max(yrange)-min(yrange))*0.17, cex=1.2, adj=0, xpd=T)
    text(labels="Non Metastatic", x=c(max(xrange)-(max(xrange)-min(xrange))*0.01), y=min(yrange)-abs(max(yrange)-min(yrange))*0.17, cex=1.2, adj=1, xpd=T)
    # axis
    text(labels=c("Enrichment Score"),x=min(xrange)-abs(max(xrange)-min(xrange))*0.09, y=mean(c(max(yrange),min(yrange)-abs(max(yrange)-min(yrange))*0.12*2)),srt=90, xpd=T, cex=1.25)
    text(labels=c("Gene Rank"), x=mean(xrange), y=min(yrange)-abs(max(yrange)-min(yrange))*0.4 , srt=0, xpd=T, cex=1.25)
    axis(side=1, at=c(seq(min(xseq),max(xseq),5000),max(xseq)), labels=c(seq(min(xseq),max(xseq),5000),max(xseq)), cex.axis=0.8, las=1, mgp=c(3,2.38,2.35))
    text(labels=yseq, x=min(xrange)-abs(max(xrange)-min(xrange))*0.005, y=yseq, cex=0.8, adj=1, xpd=T)
    dev.off()
}



########################################
## Gene-gene interaction network inference
########################################

## Input data
degs=ddf$Gene[ddf$log2FoldChange>0 & ddf$deg==T] # up DEG
tid=sinfo$Sample[sinfo$`meta type to all5`=='Metastatic'] # Metastatic sample ID
exp=t(tpm[degs,tid]) # filter expression
dim(exp) # 55 96


## Network inference
# (1) Select optimal lambda
set.seed(123)
lams = seq(0.1, 1, length.out = 10)
out = huge(exp, method = "glasso", lambda = lams)
sel = huge.select(out, criterion = "stars")  # EBIC, stars, ric...
# (2) gene-gene interaction
sel$opt.lambda # 0.1
fit = huge(exp, method = "glasso", lambda = sel$opt.lambda)
adj = as.matrix(fit$path[[1]])
rownames(adj) = colnames(adj) = colnames(exp)
# (3) Edge list
edf = which(adj != 0, arr.ind = TRUE)
edf = edf[edf[, 1] < edf[, 2], ]  # remove duplicate
ggi = data.frame(from = rownames(adj)[edf[, 1]], to = rownames(adj)[edf[, 2]])
# (4) as network
graph = graph_from_data_frame(ggi, directed = T)
graph = igraph::simplify(graph)


## Cluster analysis
clusters=cluster_walktrap(graph) # non-random
length(communities(clusters)) # Group number: 5
cluster_sizes = sizes(clusters) # cluster size: 5
V(graph)$cluster=membership(clusters)[match(names(V(graph)) , names(membership(clusters)))]


## Centrality calculation
clusters_to_keep = which(cluster_sizes >= 10) # Set the member-count cutoff for PUBLIC group removal
ndf=data.frame()
for(c in clusters_to_keep) {
  nodes=V(graph)$name[(membership(clusters) ==c)] 
  egraph=delete_vertices(graph, setdiff(V(graph)$name,nodes))
  degree = degree(egraph)
  closeness = closeness(egraph, mode="all")
  betweenness = betweenness(egraph)
  eigenvector = eigen_centrality(egraph)$vector
  aut=hits_scores(egraph)$authority
  hub = hits_scores(egraph)$hub
  edf=data.frame(Gene=nodes,cluster=c,degree=degree,closeness=closeness,aut=aut,betweenness=betweenness,eigenvector=eigenvector, hub=hub)
  ndf=rbind(ndf,edf)
}
ndf=group_by(ndf, cluster)
ndf = mutate(ndf, z_degree = scale(degree), z_betweenness = scale(betweenness), z_eigenvector = scale(eigenvector), z_hub = scale(hub))
ndf = mutate(ndf, nscore = z_degree + z_betweenness + z_eigenvector + z_hub)


## Save
# save(graph, clusters, ndf, cluster_sizes, clusters_to_keep, file=paste0(dir_data,"2. Gene interaction Network analysis res.Rdata"))


## Viz - Main Figure 2(D)
nodes_to_keep=V(graph)$name[V(graph)$cluster %in% names(cluster_sizes[cluster_sizes>=10])]
hubs=unlist(lapply(unique(ndf$cluster), function(m) {
    endf=ndf[ndf$cluster==m,]
    cuts=sort(unique(endf$closeness),decreasing=T)[1]
    endf$Gene[endf$closeness>=min(cuts)]
  })
)
ggraph=graph
V(ggraph)$hub=ifelse(V(ggraph)$name %in% hubs,1,0)
V(ggraph)$label=ifelse(V(ggraph)$name %in% hubs,V(ggraph)$name,"")
comm = igraph::make_clusters(ggraph, membership = as.integer(V(ggraph)$cluster))
memberships=membership(comm)
E(ggraph)$external=apply(ends(ggraph, E(ggraph)), 1, function(int) {
  cluster1=V(ggraph)$cluster[V(ggraph)$name==int[1]]
  cluster2=V(ggraph)$cluster[V(ggraph)$name==int[2]]
  ifelse(cluster1==cluster2,cluster1,"external")
})
comm = igraph::make_clusters(ggraph, membership = as.integer(V(ggraph)$cluster))
dedges = delete_edges(ggraph, E(ggraph)[crossing(comm, ggraph)])
ggraph = subgraph_from_edges(ggraph, E(ggraph)[!crossing(comm, ggraph)], delete.vertices = FALSE)
V(ggraph)$hub=ifelse(V(ggraph)$name %in% hubs,1,0)
createNetworkFromIgraph(ggraph, title = "Main Figure 2(D)", collection = "Main Figure 2(D)")



########################################
## Metastatsis biomarker final selecion
########################################
fdegs=unlist(lapply(unique(ndf$cluster), function(m) {
    endf=ndf[ndf$cluster==m,]
    cuts=sort(unique(endf$closeness),decreasing=T)[1]
    endf$Gene[endf$closeness>=min(cuts)]
  })
)
# save(fdegs, file=paste0(dir_data,"2. Metastasis biomarker.Rdata"))



########################################
## Centrality (Closeness) Rank
########################################

pcols = colorRampPalette(c('#d92847','#221331'))

tiff(filename=paste0(dir_fig,"Main Figure 2(E).tif"), width=17, height=15, units="cm", res=300)
par(mfcol=c(2,3), plt=c(0.1,0.1,0.9,0.9), mar=c(3.5,1.5,3.5,1.5), mgp=c(1.5,0.2,0), tck=-0.01)  
for(c in clusters_to_keep) {
  gdf=ndf[ndf$cluster==c,]
  gdf=gdf[order(gdf$closeness,decreasing=T),][1:7,]
  xrange=c(0,max(gdf$closeness))
  xseq=seq(0,max(gdf$closeness),0.01)
  ranks=rank(1-gdf$closeness, ties.method = "min")
  gpcols=pcols(max(ranks)) 
  bp=barplot(rev(gdf$closeness), xlim=xrange, horiz=T, xaxt='n', yaxt='n', xlab="Closeness", names.arg=NA, main=paste("Cluster",c), width=1.7, space=0.3, border=NA, col=rev(gpcols[ranks]), cex.main=1.5)
  text(x=0.001, y=bp, labels=rev(gdf$Gene), col="white", xpd=T, cex=1.5, adj=0)
  axis(1, at=xseq, labels=xseq, cex.axis=0.8, las=1)
}
dev.off()



########################################
## 13 Genes' expression pattern
########################################
nid=sinfo$Sample[sinfo$`meta type to all5`=="Non Metastatic"]
mid=sinfo$Sample[sinfo$`meta type to all5`=="Metastatic"]

tiff(filename=paste0(dir_fig, "Supplementary Figure 2(C).tif"),  width= 27, height = 15 , units = 'cm', res=300)
par(mfrow=c(2,7),mar=c(2.2,2.2,5,0.5),mgp=c(1.1,0.3,0),cex.lab=1,cex.axis=1.7,cex.main=2,tck=-0.02,las=1,bty="l")
for(g in sort(fdegs)) {
      ne=tpm[g,nid] # normal expression
      me=tpm[g,mid] # disease expression

      pv=ddf$padj[ddf$Gene==g]
      bdf=data.frame(Group=c(rep("Non metastatic",length(ne)),rep("Metastatic",length(me))) , Expression=c(ne,me) , check.names=F)
      ymin=0
      ymax=max(bdf$Expression)*1.15
      ylim=c(ymin,ymax)
      boxplot(bdf$Expression[bdf$Group=="Non metastatic"] , bdf$Expression[bdf$Group=="Metastatic"],
            xlab="", ylab="", ylim=ylim, outline=F, col="white", border="white", boxwex=0.6, frame=T, xaxt="n")
      grid(NA, NULL, lty=3, lwd=1, col="darkgrey")
      par(new=T)
      boxplot(main=g, bdf$Expression[bdf$Group=="Non metastatic"] , bdf$Expression[bdf$Group=="Metastatic"],
            xlab="", ylab="", ylim=ylim, outline=F, col=NA, border="black", boxwex=0.6, frame=T, xaxt="n")
      for(i in 1:length(unique(bdf$Group))) {stripchart(bdf$Expression[bdf$Group==unique(bdf$Group)[i]], at=i, cex=1.5, lwd=0.3, pch = 21, col=c("#000dff","#ff0000")[i], bg=c("#403dff","#ff3d3d")[i], method="jitter", vertical=TRUE, add = TRUE)}
      if(pv < 0.05) {
      exp = ceiling(log10(pv))-1
      numord = round(pv* (10^-exp), digits=2)
      typos=ymax+(ymax-ymin)*0.02
      text(x=1.5, y=typos, bquote(.(numord)~"x"~10^.(exp)), cex=1.7, adj=0.5,xpd=T)
      } else {
      pv=round(pv,2)
      text(x=1.5, y=typos, bquote(.(pv) ), cex=2.1, adj=0.5,xpd=T)
      } 
      sypos=ymax-(ymax-ymin)*0.06
      sypos2=ymax-(ymax-ymin)*0.08
      segments(x0=1,x1=2,y0=sypos,y1=sypos,col='black', xpd=T)
      segments(x0=c(1,2),x1=c(1,2),y0=c(sypos2,sypos2),y1=c(sypos,sypos),col='black', xpd=T)
}
dev.off()



########################################
## Metastasis risk score calculation in each metastatic type group
########################################

## Metastasis risk score calculation
param=ssgseaParam(tpm, list(`Risk_score`=fdegs))
mdf=cbind(data.frame(Sample=colnames(tpm)),t(as.data.frame(gsva(param,verbose=T))))
# save(mdf, file=paste0(dir_data,"2. cancer cell line metastasis risk score.Rdata"))

mdf=left_join(sinfo,mdf)
nrisk=mdf$Risk_score[mdf$`meta type to all5`=="Non Metastatic"]
mrisk=mdf$Risk_score[mdf$`meta type to all5`=="Metastatic"]
wres=wilcox.test(mrisk,nrisk , alternative="greater")
pv=wres$p.value


## Viz
# col=c("#000dff","#ff0000"), bg=c("#403dff","#ff3d3d")
pcols=c("#403dff","#ff3d3d")

tiff(filename=paste0(dir_fig, "Supplementary Figure 2(D).tif"),  width= 8, height = 12 , units = 'cm', res=300)
par(mfrow=c(1,1),mar=c(2.2,3.2,5,0.5),mgp=c(1.1,0.3,0),cex.lab=1,cex.axis=1.7,cex.main=2,tck=-0.02,las=1,bty="l")
bdf=data.frame(Group=c(rep("Non metastatic",length(nrisk)),rep("Metastatic",length(mrisk))) , Risk_score=c(nrisk,mrisk) , check.names=F)
ymin=round(min(bdf$Risk_score),3)
ymax=max(bdf$Risk_score)*1.15
ylim=c(ymin,ymax)
boxplot(bdf$Risk_score[bdf$Group=="Non metastatic"] , bdf$Risk_score[bdf$Group=="Metastatic"], xlab="", ylab="", ylim=ylim, outline=F, col="white", border="white", boxwex=0.6, frame=T, xaxt="n")
grid(NA, NULL, lty=3, lwd=1, col="darkgrey")
par(new=T)
boxplot(main="Risk score", bdf$Risk_score[bdf$Group=="Non metastatic"] , bdf$Risk_score[bdf$Group=="Metastatic"], xlab="", ylab="", ylim=ylim, outline=F, col=NA, border="black", boxwex=0.6, frame=T, xaxt="n")
for(i in 1:length(unique(bdf$Group))) {stripchart(bdf$Risk_score[bdf$Group==unique(bdf$Group)[i]], at=i, cex=1.5, lwd=0.3, pch = 21, col=c("#000dff","#ff0000")[i], bg=c("#403dff","#ff3d3d")[i], method="jitter", vertical=TRUE, add = TRUE)}
if(pv < 0.05) {
    exp = ceiling(log10(pv))-1
    numord = round(pv* (10^-exp), digits=2)
    typos=ymax+(ymax-ymin)*0.01
    text(x=1.5, y=typos, bquote(.(numord)~"x"~10^.(exp)), cex=1.7, adj=0.5,xpd=T)
} else {
    pv=round(pv,2)
    text(x=1.5, y=typos, bquote(.(pv) ), cex=2.1, adj=0.5,xpd=T)
}
sypos=ymax-(ymax-ymin)*0.06
sypos2=ymax-(ymax-ymin)*0.08  
segments(x0=1,x1=2,y0=sypos,y1=sypos,col='black', xpd=T)
segments(x0=c(1,2),x1=c(1,2),y0=c(sypos2,sypos2),y1=c(sypos,sypos),col='black', xpd=T)
dev.off()





########################################
## Gene-gene interaction network inference
########################################

## Input data
degs=ddf$Gene[ddf$log2FoldChange>0 & ddf$deg==T] # up DEG
tid=sinfo$Sample[sinfo$`meta type to all5`=='Metastatic'] # Metastatic sample ID
exp=t(tpm[degs,tid]) # filter expression
dim(exp) # 55 96


## Network inference
# (1) Select optimal lambda
set.seed(123)
lams = seq(0.1, 1, length.out = 10)
out = huge(exp, method = "glasso", lambda = lams)
sel = huge.select(out, criterion = "stars")  # EBIC, stars, ric...
# (2) gene-gene interaction
sel$opt.lambda # 0.1
fit = huge(exp, method = "glasso", lambda = sel$opt.lambda)
adj = as.matrix(fit$path[[1]])
rownames(adj) = colnames(adj) = colnames(exp)
# (3) Edge list
edf = which(adj != 0, arr.ind = TRUE)
edf = edf[edf[, 1] < edf[, 2], ]  # remove duplicate
ggi = data.frame(from = rownames(adj)[edf[, 1]], to = rownames(adj)[edf[, 2]])
# (4) as network
graph = graph_from_data_frame(ggi, directed = T)
graph = igraph::simplify(graph)


## Cluster analysis
clusters=cluster_walktrap(graph) # non-random
length(communities(clusters)) # Group number: 5
cluster_sizes = sizes(clusters) # cluster size: 5
V(graph)$cluster=membership(clusters)[match(names(V(graph)) , names(membership(clusters)))]


## Centrality calculation
clusters_to_keep = which(cluster_sizes >= 10) # Set the member-count cutoff for PUBLIC group removal
ndf=data.frame()
for(c in clusters_to_keep) {
  nodes=V(graph)$name[(membership(clusters) ==c)] 
  egraph=delete_vertices(graph, setdiff(V(graph)$name,nodes))
  degree = degree(egraph)
  closeness = closeness(egraph, mode="all")
  betweenness = betweenness(egraph)
  eigenvector = eigen_centrality(egraph)$vector
  aut=hits_scores(egraph)$authority
  hub = hits_scores(egraph)$hub
  edf=data.frame(Gene=nodes,cluster=c,degree=degree,closeness=closeness,aut=aut,betweenness=betweenness,eigenvector=eigenvector, hub=hub)
  ndf=rbind(ndf,edf)
}
ndf=group_by(ndf, cluster)
ndf = mutate(ndf, z_degree = scale(degree), z_betweenness = scale(betweenness), z_eigenvector = scale(eigenvector), z_hub = scale(hub))
ndf = mutate(ndf, nscore = z_degree + z_betweenness + z_eigenvector + z_hub)
ndf = as.data.frame(ndf)


## Save
# save(graph, clusters, ndf, cluster_sizes, clusters_to_keep, file=paste0(dir_data,"2. Gene interaction Network analysis res.Rdata"))


## Viz - Main Figure 2(D)
nodes_to_keep=V(graph)$name[V(graph)$cluster %in% names(cluster_sizes[cluster_sizes>=10])]
hubs=unlist(lapply(unique(ndf$cluster), function(m) {
    endf=ndf[ndf$cluster==m,]
    cuts=sort(unique(endf$closeness),decreasing=T)[1]
    endf$Gene[endf$closeness>=min(cuts)]
  })
)
ggraph=graph
V(ggraph)$hub=ifelse(V(ggraph)$name %in% hubs,1,0)
V(ggraph)$label=ifelse(V(ggraph)$name %in% hubs,V(ggraph)$name,"")
comm = igraph::make_clusters(ggraph, membership = as.integer(V(ggraph)$cluster))
memberships=membership(comm)
E(ggraph)$external=apply(ends(ggraph, E(ggraph)), 1, function(int) {
  cluster1=V(ggraph)$cluster[V(ggraph)$name==int[1]]
  cluster2=V(ggraph)$cluster[V(ggraph)$name==int[2]]
  ifelse(cluster1==cluster2,cluster1,"external")
})
comm = igraph::make_clusters(ggraph, membership = as.integer(V(ggraph)$cluster))
dedges = delete_edges(ggraph, E(ggraph)[crossing(comm, ggraph)])
ggraph = subgraph_from_edges(ggraph, E(ggraph)[!crossing(comm, ggraph)], delete.vertices = FALSE)
V(ggraph)$hub=ifelse(V(ggraph)$name %in% hubs,1,0)
createNetworkFromIgraph(ggraph, title = "Main Figure 4(D)", collection = "Main Figure 4(D)")



########################################
## Metastatsis biomarker final selecion
########################################
fdegs=unlist(lapply(unique(ndf$cluster), function(m) {
    endf=ndf[ndf$cluster==m,]
    cuts=sort(unique(endf$closeness),decreasing=T)[1]
    endf$Gene[endf$closeness>=min(cuts)]
  })
)
# save(fdegs, file=paste0(dir_data,"2. Metastasis biomarker.Rdata"))



########################################
## Centrality Rank
########################################

pcols = colorRampPalette(c('#d92847','#221331'))

tiff(filename=paste0(dir_fig,"Main Figure 2(E).tif"), width=17, height=15, units="cm", res=300)
par(mfcol=c(2,3), plt=c(0.1,0.1,0.9,0.9), mar=c(3.5,1.5,3.5,1.5), mgp=c(1.5,0.2,0), tck=-0.01)  
for(c in clusters_to_keep) {
  gdf=ndf[ndf$cluster==c,]
  gdf=gdf[order(gdf$closeness,decreasing=T),][1:7,]
  xrange=c(0,max(gdf$closeness))
  xseq=seq(0,max(gdf$closeness),0.01)
  ranks=rank(1-gdf$closeness, ties.method = "min")
  gpcols=pcols(max(ranks)) 
  bp=barplot(rev(gdf$closeness), xlim=xrange, horiz=T, xaxt='n', yaxt='n', xlab="Closeness", names.arg=NA, main=paste("Cluster",c), width=1.7, space=0.3, border=NA, col=rev(gpcols[ranks]), cex.main=1.5)
  text(x=0.001, y=bp, labels=rev(gdf$Gene), col="white", xpd=T, cex=1.5, adj=0)
  axis(1, at=xseq, labels=xseq, cex.axis=0.8, las=1)
}
dev.off()



########################################
## 13 Genes' expression pattern
########################################
nid=sinfo$Sample[sinfo$`meta type to all5`=="Non Metastatic"]
mid=sinfo$Sample[sinfo$`meta type to all5`=="Metastatic"]

tiff(filename=paste0(dir_fig, "Supplementary Figure 2(C).tif"),  width= 27, height = 15 , units = 'cm', res=300)
par(mfrow=c(2,7),mar=c(2.2,2.2,5,0.5),mgp=c(1.1,0.3,0),cex.lab=1,cex.axis=1.7,cex.main=2,tck=-0.02,las=1,bty="l")
for(g in sort(fdegs)) {
      ne=tpm[g,nid] # normal expression
      me=tpm[g,mid] # disease expression

      pv=ddf$padj[ddf$Gene==g]
      bdf=data.frame(Group=c(rep("Non metastatic",length(ne)),rep("Metastatic",length(me))) , Expression=c(ne,me) , check.names=F)
      ymin=0
      ymax=max(bdf$Expression)*1.15
      ylim=c(ymin,ymax)
      boxplot(bdf$Expression[bdf$Group=="Non metastatic"] , bdf$Expression[bdf$Group=="Metastatic"],
            xlab="", ylab="", ylim=ylim, outline=F, col="white", border="white", boxwex=0.6, frame=T, xaxt="n")
      grid(NA, NULL, lty=3, lwd=1, col="darkgrey")
      par(new=T)
      boxplot(main=g, bdf$Expression[bdf$Group=="Non metastatic"] , bdf$Expression[bdf$Group=="Metastatic"],
            xlab="", ylab="", ylim=ylim, outline=F, col=NA, border="black", boxwex=0.6, frame=T, xaxt="n")
      for(i in 1:length(unique(bdf$Group))) {stripchart(bdf$Expression[bdf$Group==unique(bdf$Group)[i]], at=i, cex=1.5, lwd=0.3, pch = 21, col=c("#000dff","#ff0000")[i], bg=c("#403dff","#ff3d3d")[i], method="jitter", vertical=TRUE, add = TRUE)}
      if(pv < 0.05) {
      exp = ceiling(log10(pv))-1
      numord = round(pv* (10^-exp), digits=2)
      typos=ymax+(ymax-ymin)*0.02
      text(x=1.5, y=typos, bquote(.(numord)~"x"~10^.(exp)), cex=1.7, adj=0.5,xpd=T)
      } else {
      pv=round(pv,2)
      text(x=1.5, y=typos, bquote(.(pv) ), cex=2.1, adj=0.5,xpd=T)
      } 
      sypos=ymax-(ymax-ymin)*0.06
      sypos2=ymax-(ymax-ymin)*0.08
      segments(x0=1,x1=2,y0=sypos,y1=sypos,col='black', xpd=T)
      segments(x0=c(1,2),x1=c(1,2),y0=c(sypos2,sypos2),y1=c(sypos,sypos),col='black', xpd=T)
}
dev.off()



########################################
## Metastasis risk score calculation in each metastatic type group
########################################

## Metastasis risk score calculation
param=ssgseaParam(tpm, list(`Risk_score`=fdegs))
mdf=cbind(data.frame(Sample=colnames(tpm)),t(as.data.frame(gsva(param,verbose=T))))
# save(mdf, file=paste0(dir_data,"2. cancer cell line metastasis risk score.Rdata"))

mdf=left_join(sinfo,mdf)
nrisk=mdf$Risk_score[mdf$`meta type to all5`=="Non Metastatic"]
mrisk=mdf$Risk_score[mdf$`meta type to all5`=="Metastatic"]
wres=wilcox.test(mrisk,nrisk , alternative="greater")
pv=wres$p.value


## Viz
# col=c("#000dff","#ff0000"), bg=c("#403dff","#ff3d3d")
pcols=c("#403dff","#ff3d3d")

tiff(filename=paste0(dir_fig, "Supplementary Figure 2(D).tif"),  width= 8, height = 12 , units = 'cm', res=300)
par(mfrow=c(1,1),mar=c(2.2,3.2,5,0.5),mgp=c(1.1,0.3,0),cex.lab=1,cex.axis=1.7,cex.main=2,tck=-0.02,las=1,bty="l")
bdf=data.frame(Group=c(rep("Non metastatic",length(nrisk)),rep("Metastatic",length(mrisk))) , Risk_score=c(nrisk,mrisk) , check.names=F)
ymin=round(min(bdf$Risk_score),3)
ymax=max(bdf$Risk_score)*1.15
ylim=c(ymin,ymax)
boxplot(bdf$Risk_score[bdf$Group=="Non metastatic"] , bdf$Risk_score[bdf$Group=="Metastatic"], xlab="", ylab="", ylim=ylim, outline=F, col="white", border="white", boxwex=0.6, frame=T, xaxt="n")
grid(NA, NULL, lty=3, lwd=1, col="darkgrey")
par(new=T)
boxplot(main="Risk score", bdf$Risk_score[bdf$Group=="Non metastatic"] , bdf$Risk_score[bdf$Group=="Metastatic"], xlab="", ylab="", ylim=ylim, outline=F, col=NA, border="black", boxwex=0.6, frame=T, xaxt="n")
for(i in 1:length(unique(bdf$Group))) {stripchart(bdf$Risk_score[bdf$Group==unique(bdf$Group)[i]], at=i, cex=1.5, lwd=0.3, pch = 21, col=c("#000dff","#ff0000")[i], bg=c("#403dff","#ff3d3d")[i], method="jitter", vertical=TRUE, add = TRUE)}
if(pv < 0.05) {
    exp = ceiling(log10(pv))-1
    numord = round(pv* (10^-exp), digits=2)
    typos=ymax+(ymax-ymin)*0.01
    text(x=1.5, y=typos, bquote(.(numord)~"x"~10^.(exp)), cex=1.7, adj=0.5,xpd=T)
} else {
    pv=round(pv,2)
    text(x=1.5, y=typos, bquote(.(pv) ), cex=2.1, adj=0.5,xpd=T)
}
sypos=ymax-(ymax-ymin)*0.06
sypos2=ymax-(ymax-ymin)*0.08  
segments(x0=1,x1=2,y0=sypos,y1=sypos,col='black', xpd=T)
segments(x0=c(1,2),x1=c(1,2),y0=c(sypos2,sypos2),y1=c(sypos,sypos),col='black', xpd=T)
dev.off()



########################################
## 1st Evaluation: Leave-one out (LOO)
########################################

ddfl=c() # DEG analysis res
graphl=c()
clusterl=c()
ndfl=c()
cluster_sizel=c()
clusters_to_keepl=c()
fdegl=c()

cids=sinfo$Sample[sinfo$`meta type to all5`=="Non Metastatic"]
fcut=log2(2)
pcut=0.05
set.seed(123)

for(eid in cids) {
  # (1) DEG analysis
  cid=setdiff(cids,eid)
  tid=sinfo$Sample[sinfo$`meta type to all5`=="Metastatic"]
  dg = rbind(data.frame(sample = cid, group = "control"), data.frame(sample = tid, group = "treat"))
  counts = rcm[, dg$sample]
  group = factor(dg$group)
  dge = DGEList(counts=counts, group=group)
  keep = filterByExpr(dge)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  dge = calcNormFactors(dge)
  design = model.matrix(~group) 
  dge = estimateDisp(dge, design)
  fit = glmFit(dge, design)
  lrt = glmLRT(fit)
  vddf = topTags(lrt, n=Inf)$table
  colnames(vddf) = c("log2FoldChange","logCPM","LR","pvalue","padj")
  vddf$Gene = gsub("_[0-9]+$", "", rownames(vddf))
  # sort
  vddf=vddf[order(abs(vddf$log2FoldChange), decreasing=T),]
  vddf=vddf[order(vddf$pvalue),]
  vddf=vddf[order(vddf$padj),]
  # DEG
  vddf$deg=(abs(vddf$log2FoldChange)>fcut & vddf$padj<pcut & !is.na(vddf$log2FoldChange) & !is.na(vddf$padj))
  ddfl[[paste0("Excluding ",eid)]]=vddf

  # (2) Gene interaction network construction
  degs=vddf$Gene[vddf$log2FoldChange>0 & vddf$deg==T] # up DEG
  exp=t(tpm[degs,tid]) # filter expression
  # (2-1) Select optimal lambda
  lams = seq(0.1, 1, length.out = 10)
  out = huge(exp, method = "glasso", lambda = lams)
  sel = huge.select(out, criterion = "stars")  # EBIC, stars, ric...
  # (2-2) gene-gene interaction
  sel$opt.lambda # 0.1
  fit = huge(exp, method = "glasso", lambda = sel$opt.lambda)
  adj = as.matrix(fit$path[[1]])
  rownames(adj) = colnames(adj) = colnames(exp)
  # (2-3) Edge list
  edf = which(adj != 0, arr.ind = TRUE)
  edf = edf[edf[, 1] < edf[, 2], ]  # remove duplicate
  ggi = data.frame(from = rownames(adj)[edf[, 1]], to = rownames(adj)[edf[, 2]])
  # (2-4) as network
  vggi=ggi[ggi$from %in% vddf$Gene[abs(vddf$log2FoldChange)>fcut & !is.na(vddf$log2FoldChange) & !is.na(vddf$padj) & vddf$padj<pcut] & ggi$to %in% vddf$Gene[abs(vddf$log2FoldChange)>fcut & !is.na(vddf$log2FoldChange) & !is.na(vddf$padj) & vddf$padj<pcut],]
  vgraph = graph_from_data_frame(vggi, directed = T)
  vgraph = igraph::simplify(vgraph)
  clusters=cluster_walktrap(vgraph) # non-random
  length(communities(clusters)) 
  cluster_sizes = sizes(clusters) 
  V(vgraph)$cluster=membership(clusters)[match(names(V(vgraph)) , names(membership(clusters)))]

  # (3) Centrality calculation
  clusters_to_keep = which(cluster_sizes>= 10) # Set the member-count cutoff for PUBLIC group removal
  vndf=data.frame()
  for(c in clusters_to_keep) {
    nodes=V(vgraph)$name[(membership(clusters) ==c)] 
    evgraph=delete_vertices(vgraph, setdiff(V(vgraph)$name,nodes))
    degree = degree(evgraph)
    closeness = closeness(evgraph, mode="all")
    betweenness = betweenness(evgraph)
    eigenvector = eigen_centrality(evgraph)$vector
    aut=hits_scores(evgraph)$authority
    hub = hits_scores(evgraph)$hub
    edf=data.frame(Gene=nodes,cluster=c,degree=degree,closeness=closeness,aut=aut,betweenness=betweenness,eigenvector=eigenvector, hub=hub)
    vndf=rbind(vndf,edf)
  }
  vndf=group_by(vndf, cluster)
  vndf = mutate(vndf, z_degree = scale(degree), z_betweenness = scale(betweenness), z_eigenvector = scale(eigenvector), z_hub = scale(hub))
  vndf = mutate(vndf, nscore = z_degree + z_betweenness + z_eigenvector + z_hub)

  # (4) Hub gene selection
  hubs=unlist(lapply(unique(vndf$cluster), function(m) {
      evndf=vndf[vndf$cluster==m,]
      cuts=sort(unique(evndf$closeness),decreasing=T)[1]
      evndf$Gene[evndf$closeness>=min(cuts)]
    })
  )
  gvgraph=vgraph
  V(gvgraph)$hub=ifelse(V(gvgraph)$name %in% hubs,1,0)
  V(gvgraph)$label=ifelse(V(gvgraph)$name %in% hubs,V(gvgraph)$name,"")
  comm = igraph::make_clusters(gvgraph, membership = as.integer(V(gvgraph)$cluster))
  memberships=membership(comm)
  E(gvgraph)$external=apply(ends(gvgraph, E(gvgraph)), 1, function(int) {
    cluster1=V(gvgraph)$cluster[V(gvgraph)$name==int[1]]
    cluster2=V(gvgraph)$cluster[V(gvgraph)$name==int[2]]
    ifelse(cluster1==cluster2,cluster1,"external")
  })
  comm = igraph::make_clusters(gvgraph, membership = as.integer(V(gvgraph)$cluster))
  dedges = delete_edges(gvgraph, E(gvgraph)[crossing(comm, gvgraph)])
  gvgraph = subgraph_from_edges(gvgraph, E(gvgraph)[!crossing(comm, gvgraph)], delete.vertices = FALSE)
  V(gvgraph)$hub=ifelse(V(gvgraph)$name %in% hubs,1,0)

  graphl[[paste0("Excluding ",eid)]]=vgraph
  clusterl[[paste0("Excluding ",eid)]]=clusters
  ndfl[[paste0("Excluding ",eid)]]=vndf
  cluster_sizel[[paste0("Excluding ",eid)]]=cluster_sizes
  clusters_to_keepl[[paste0("Excluding ",eid)]]=clusters_to_keep

  # (5) Metastatsis biomarker final selecion
  vfdegs=unlist(lapply(unique(vndf$cluster), function(m) {
      evndf=vndf[vndf$cluster==m,]
      cuts=sort(unique(evndf$closeness),decreasing=T)[1]
      evndf$Gene[evndf$closeness>=min(cuts)]
    })
  )
  fdegl[[paste0("Excluding ",eid)]]=vfdegs
}

lapply(fdegl, function(gs) intersect(gs,fdegs))
# $`Excluding ACH-000800`
# [1] "OPLAH" "LPIN3"

# $`Excluding ACH-000343`
# [1] "NPAS2" "OPLAH"

# $`Excluding ACH-000021`
# character(0)

sort(table(unlist(fdegl)))
#   CD109   CDCP1  CHST15  COMMD7    CYBA FAM120A    GALM    HRH1   LPIN3     MET 
#       1       1       1       1       1       1       1       1       1       1 
# MFSD14B  NAP1L5     NHS   NPAS2  PLXNA1  SLC6A6  TCIRG1   TOR4A  WASHC5    WWP2 
#       1       1       1       1       1       1       1       1       1       1 
#   LAMP1   OPLAH 
#       2       2 


## Save
# save(ddfl, graphl, clusterl, ndfl, cluster_sizel, clusters_to_keepl, fdegl, file=paste0(dir_data,"2. 1st Gene selection performance evaulation (LOO).Rdata"))


## Performance Viz
fm = sapply(ddfl, function(x) { x$log2FoldChange[match(sort(fdegs), x$Gene)] })
fm = cbind(ddf$log2FoldChange[match(sort(fdegs),ddf$Gene)] , fm)
rownames(fm) = sort(fdegs)
colnames(fm)[1]='Metastasis Signature'

pm = sapply(ddfl, function(x) { x$padj[match(sort(fdegs), x$Gene)] })
pm = cbind(ddf$padj[match(sort(fdegs),ddf$Gene)] , pm)
rownames(pm) = sort(fdegs)

tm = matrix("", nrow = nrow(pm), ncol = ncol(pm))
rownames(tm) = rownames(pm)
colnames(tm) = colnames(pm)
tm[pm < 0.05] = "*"
tm[pm < 0.01] = "**"

dcolor='#403dff'
ucolor="#ff3d3d"
pcols=colorRamp2(c(-2,0,max(fm)), c(dcolor,'#FFFFFF',ucolor))
csplit=ifelse(!grepl("Excluding",colnames(fm)),"Metastasis Signature","Excluding One Sample")
csplit=factor(csplit, levels=c("Metastasis Signature","Excluding One Sample"))
hp=Heatmap(fm, name="hp", show_heatmap_legend=F, col=pcols, border=F, cluster_columns=F, cluster_rows=F, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_row_dend=F, column_names_rot=85, 
            row_gap=unit(0.3,"cm"), column_gap=unit(0.3,'cm'), column_title=NULL, width=ncol(fm)*unit(1.1,"cm"), height=nrow(fm)*unit(1,"cm"), rect_gp = gpar(col = "#7A7A7A", lwd = 2, lty=1),
            column_split=csplit,
            column_title_gp = gpar(fontsize = 26, fontface = "bold"),
            column_title_rot = 65,
            row_title_gp = gpar(fontsize = 28, col = c("#FF3E58", "#378EFF"), fontface = "bold"),
            row_names_gp = gpar(fontsize=23),
            column_names_gp = gpar(fontsize=22),
            cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(pm,2)[i, j], x, y, gp = gpar(fontsize=13, col="#000000")) })
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(E).tiff"), width=20, height=25, units = 'cm',res=300)
draw(hp, merge_legend=F, heatmap_legend_side="left", annotation_legend_side='left', padding=unit(c(0, 0, 0, 0), "cm"), gap=unit(0.05,"cm"))
dev.off()


## Legend
pcols=colorRamp2(c(-2,0,max(fm)), c(dcolor,'#FFFFFF',ucolor))
lg_hp=Legend(at=c(-2,0,max(fm)), labels=c(-2,0,round(max(fm),1)), col_fun = pcols, title = expression(Log[2] * FC), border='#7A7A7A', title_position="topcenter")
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(E)_legend.tiff"), width=7, height=15, units = 'cm',res=300)
draw(packLegend(list=list(lg_hp))) 
dev.off() 


## Corelation with FCs
sapply(seq_len(ncol(fm)), function(i) {
  cor(ddf$log2FoldChange[match(rownames(fm),ddf$Gene)], fm[, i], method = "pearson", use = "complete.obs")
})
# 0.9521472 0.9862175 0.9105194


for(i in setdiff(seq_len(ncol(fm)),1)) {
        x=ddf$log2FoldChange[match(rownames(fm),ddf$Gene)]
        y=fm[, i]
        cres=cor.test(x,y, method="pearson", use = "complete.obs")
        ly=lm(y~x)
        r=round(cres$estimate,3)
        pv=cres$p.value
        pv=sprintf("%.2e", pv)

        xlim=c(min(x),max(x))
        ylim=c(min(y),max(y))
        xpos=min(x)+(max(x)-min(x))*0.01
        ypos=max(y)-(max(y)-min(y))*0.06

        tiff(filename=paste0(dir_fig,"Supplementary Figure 2(F)_",colnames(fm)[i],".tif"), width=12, height=13, units = 'cm',res=300)
        par(mar=c(4,4,3,0.25),  mgp=c(2.2,0.25,0), cex.axis=1.5, cex.lab=1.2, cex.main=1.4, tck=-0.01)
        plot(x=x, y=y, xlim=xlim, ylim=ylim, main=colnames(fm)[i], xlab=colnames(fm)[1], ylab=colnames(fm)[i], lwd=1.9, cex=2.2, col="#FF3E58", bg="#FF3E5880", pch=21)
        abline(ly, lty=3, col='#d92847', lwd=5.2)
        text(x=xpos, y=ypos, paste0("r: ",r,"\n","p-value: ",pv), cex=1.7, adj=0, col='#000000')
        dev.off()
}



########################################
## 2nd Evaluation: Centrality
########################################
load((paste0(dir_data,"2. Gene interaction Network analysis res.Rdata")))
centralities=c("closeness","degree","betweenness","eigenvector","HITS hub")

fdegl=lapply(centralities, function(centrality) {
  unlist(lapply(unique(ndf$cluster), function(m) {
      endf=ndf[ndf$cluster==m,]
      cuts=sort(unique(endf[,centrality]),decreasing=T)[1]
      endf$Gene[endf[,centrality]>=min(cuts)]
    })
  )
})
names(fdegl)=centralities


## Save
# save(fdegl, file=paste0(dir_data,"2. 2nd Gene selection performance evaulation (Centrality).Rdata"))


## Viz
genes=sort(unique(unlist(fdegl)))
cm=do.call(cbind , lapply(centralities, function(centrality) genes %in% fdegl[[centrality]]))
rownames(cm)=genes
colnames(cm)=str_to_title(centralities)
cm=ifelse(cm==F,0,1)
original=ifelse(rownames(cm) %in% fdegs,"Metastasis Signature","Non Metastasis Signature")
original=factor(original, levels=sort(unique(original)))
csplit=c("Closeness",rep("Others",(ncol(cm)-1)))
csplit=factor(csplit, levels=c("Closeness","Others"))

pcols=colorRamp2(c(0,1), c('#221331','#d92847'))
hp=Heatmap(cm, name="hp", show_heatmap_legend=F, col=pcols, border=F, cluster_columns=F, cluster_rows=F, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_row_dend=F, column_names_rot=85,  
            row_split=original, row_gap=unit(0.3,"cm"), column_gap=unit(0.3,'cm'), column_split=csplit,
            width=ncol(cm)*unit(1.1,"cm"), height=nrow(cm)*unit(1,"cm"), rect_gp = gpar(col = "#7A7A7A", lwd = 2, lty=1),
            column_title=NULL, column_title_gp = gpar(fontsize = 27, fontface = "bold"), column_names_side='bottom', 
            row_title_rot = 90, row_title=gsub(" Signature","\nSignature",sort(unique(original))), row_title_gp = gpar(fontsize = 28, col = c("#ff3d3d", "darkgrey"), fontface = "bold"),
            row_names_gp = gpar(fontsize=23),
            column_names_gp = gpar(fontsize=23))
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(G).tiff"), width=20, height=30, units = 'cm',res=300)
draw(hp, merge_legend=F, heatmap_legend_side="left", annotation_legend_side='left', padding=unit(c(0, 0, 0, 0), "cm"), gap=unit(0.05,"cm"))
dev.off()

lg=Legend(labels=sort(c("Hub","Non Hub")), legend_gp = gpar(fill=c('#d92847','#221331')),border="#b9b9b9ff")
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(G)_legend.tiff"), width=4, height=4, units = 'cm',res=300)
draw(packLegend(list=list(lg))) 
dev.off()



########################################
## 3rd Evaluation: Clustering method
########################################
load((paste0(dir_data,"2. Gene interaction Network analysis res.Rdata")))

ndfl=c()
graphl=c()
clusterl=c()
cluster_sizel=c()
clusters_to_keepl=c()


## Leiden
clusters=cluster_leiden(graph, objective_function="modularity", resolution=1) # non-random
length(communities(clusters)) # Group number: 6
cluster_sizes = sizes(clusters) 
V(graph)$cluster=membership(clusters)[match(names(V(graph)) , names(membership(clusters)))]
clusters_to_keep = which(cluster_sizes >= 10) # Set the member-count cutoff for PUBLIC group removal
ndf=data.frame()

for(c in clusters_to_keep) {
  nodes=V(graph)$name[(membership(clusters) ==c)] 
  egraph=delete_vertices(graph, setdiff(V(graph)$name,nodes))
  degree = degree(egraph)
  closeness = closeness(egraph, mode="all")
  betweenness = betweenness(egraph)
  eigenvector = eigen_centrality(egraph)$vector
  aut=hits_scores(egraph)$authority
  hub = hits_scores(egraph)$hub
  edf=data.frame(Gene=nodes,cluster=c,degree=degree,closeness=closeness,aut=aut,betweenness=betweenness,eigenvector=eigenvector, hub=hub)
  ndf=rbind(ndf,edf)
}
ndf=group_by(ndf, cluster)
ndf = mutate(ndf, z_degree = scale(degree), z_betweenness = scale(betweenness), z_eigenvector = scale(eigenvector), z_hub = scale(hub))
ndf = mutate(ndf, nscore = z_degree + z_betweenness + z_eigenvector + z_hub)
ndf = as.data.frame(ndf)

ndfl[["Leiden"]]=ndf
graphl[['Leiden']]=graph
clusterl[['Leiden']]=clusters
cluster_sizel[['Leiden']]=cluster_sizes
clusters_to_keepl[['Leiden']]=clusters_to_keep


## Louvain
clusters=cluster_louvain(graph) # non-random
length(communities(clusters)) # Group number: 6
cluster_sizes = sizes(clusters) 
V(graph)$cluster=membership(clusters)[match(names(V(graph)) , names(membership(clusters)))]
clusters_to_keep = which(cluster_sizes >= 10) # Set the member-count cutoff for PUBLIC group removal
ndf=data.frame()
for(c in clusters_to_keep) {
  nodes=V(graph)$name[(membership(clusters) ==c)] 
  egraph=delete_vertices(graph, setdiff(V(graph)$name,nodes))
  degree = degree(egraph)
  closeness = closeness(egraph, mode="all")
  betweenness = betweenness(egraph)
  eigenvector = eigen_centrality(egraph)$vector
  aut=hits_scores(egraph)$authority
  hub = hits_scores(egraph)$hub
  edf=data.frame(Gene=nodes,cluster=c,degree=degree,closeness=closeness,aut=aut,betweenness=betweenness,eigenvector=eigenvector, hub=hub)
  ndf=rbind(ndf,edf)
}
ndf=group_by(ndf, cluster)
ndf = mutate(ndf, z_degree = scale(degree), z_betweenness = scale(betweenness), z_eigenvector = scale(eigenvector), z_hub = scale(hub))
ndf = mutate(ndf, nscore = z_degree + z_betweenness + z_eigenvector + z_hub)
ndf = as.data.frame(ndf)

ndfl[["Louvain"]]=ndf
graphl[['Louvain']]=graph
clusterl[['Louvain']]=clusters
cluster_sizel[['Louvain']]=cluster_sizes
clusters_to_keepl[['Louvain']]=clusters_to_keep


## Save
# save(ndfl, graphl, clusterl, cluster_sizel, clusters_to_keepl, file=paste0(dir_data,"2. 3rd Gene selection performance evaulation (Clustering).Rdata"))


fdegl=lapply(names(ndfl), function(clusetring_method) {
  ndf=ndfl[[clusetring_method]]
  unlist(lapply(unique(ndf$cluster), function(m) {
      endf=ndf[ndf$cluster==m,]
      cuts=sort(unique(endf$closeness),decreasing=T)[1]
      endf$Gene[endf$closeness>=min(cuts)]
    })
  )
})
names(fdegl)=names(ndfl)

genes=sort(unique(c(unlist(fdegl),fdegs)))
cm=do.call(cbind , lapply(names(ndfl), function(clusetring_method) genes %in% fdegl[[clusetring_method]]))
rownames(cm)=genes
cm=cm[match(sort(fdegs),rownames(cm)),]
colnames(cm)=str_to_title(names(ndfl))
cm=cm[,c("Leiden","Louvain")]
cm=ifelse(cm==F,0,1)

pcols=colorRamp2(c(0,1), c('#221331','#d92847'))
hp=Heatmap(cm, name="hp", show_heatmap_legend=F, col=pcols, border=F, cluster_columns=F, cluster_rows=F, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_row_dend=F, column_names_rot=85,  
            width=ncol(cm)*unit(1.1,"cm"), height=nrow(cm)*unit(1,"cm"), rect_gp = gpar(col = "#7A7A7A", lwd = 2, lty=1),
            column_title_gp = gpar(fontsize = 27, fontface = "bold"), column_names_side='bottom', 
            row_names_gp = gpar(fontsize=23),
            column_names_gp = gpar(fontsize=23))
tiff(filename=paste0(dir_fig,"Supplementary Figure 2(H).tiff"), width=20, height=40, units = 'cm',res=300)
draw(hp, merge_legend=F, heatmap_legend_side="left", annotation_legend_side='left', padding=unit(c(10, 0, 0, 0), "cm"), gap=unit(0.05,"cm"))
dev.off()
