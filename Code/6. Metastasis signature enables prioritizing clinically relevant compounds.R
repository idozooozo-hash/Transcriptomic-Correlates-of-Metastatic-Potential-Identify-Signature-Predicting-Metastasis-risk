########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/figure/"

dir_cmap="G:/My Drive/DB/CMap/" # Dataset's directory The dataset directory was obtained from Connectivity Map (CMap) Genomics.
dir_omni="G:/My Drive/DB/OmniPath/homo_sapiens/" # Dataset's directory The dataset directory was obtained from OmniPath Genomics.
dir_dep="G:/My Drive/DB/DepMap/" # Dataset's directory The dataset directory was obtained from DepMap Genomics.
dir_otp="G:/My Drive/DB/Open Target Platform/" # Dataset's directory The dataset directory was obtained from <Open Traget Platform> Genomics.
dir_ensembl="G:/My Drive/DB/Ensembl/" # Dataset's directory The dataset directory was obtained from Ensemb Genomics.



########################################
## Library
######################################### 
# library(mcr)
library(dplyr)
# library(data.table)
library(metafor)
library(stringr)
library(rtracklayer)
library(igraph)
library(ggplot2)
library(ggrepel)



########################################
## Load Data
########################################
load(paste0(dir_data,"2. DEG analysis res.Rdata")) # variable name: ddf
load(paste0(dir_data,"2. Metastasis biomarker.Rdata")) # variable name: fdegs
load(paste0(dir_data,"2. Gene interaction Network analysis res.Rdata")) # variable name: graph

## cmap
load(paste0(dir_cmap,"level5(gold)_compound_perturbation(2021-11-23).Rdata")) 
# pm: Cmap level5 dataset
# ginfo: gene information
# sinfo: sample information
# tinfo: treatment information
# cinfo: compound information
cinfo1=cinfo

## omnipath
load(paste0(dir_omni,"ppi.Rdata")) # variable name: ppi; Protein protein interaction data

## depmap
load(paste0(dir_dep,"cinfo_24Q4.Rdata")) # variable name: cinfo; Compound information data
cinfo2=cinfo
load(paste0(dir_dep,"cinfo_CTD2_230723.Rdata")) # variable name: cinfo; Compound information data
cinfo3=cinfo

## open target
load(paste0(dir_otp,"compound_target.Rdata")) # variable name: ct; Compound-Target information

## 10x
gtf=import(paste0(dir_ensembl,"Homo_sapiens.GRCh38.110.gtf")) # Gene information data



########################################
## CMap dataset processing
########################################
sinfo=sinfo[sinfo$Pooling=="tumor",]
tinfo=tinfo[tinfo$CC_q75>=0.8 & tinfo$Replicate_self_rank_q25<=0.05 & !is.na(tinfo$CC_q75) & !is.na(tinfo$Replicate_self_rank_q25) & tinfo$Replicate_no>=3 & tinfo$QC_pass>=1,] # remove low quality; Replicate_self_rank_q25=pct_self_rank_q25 
tinfo=tinfo[tinfo$Perturbation_type=="trt_cp",] # only perturbation
tinfo=tinfo[tinfo$Cell_Line %in% sinfo$Cell_Line,]
length(unique(tinfo$Cell_Line)) # 43

pm=pm[,colnames(pm) %in% tinfo$Signature_id]

pm=pm[gsub(".*\\(|\\)","",rownames(pm)) %in% fdegs,]
dim(pm) # 9 2221
rownames(pm)=gsub(".*\\(|\\)","",rownames(pm))



########################################
## 1st Drug repurposing: based drug signature
########################################

dsig=ddf$log2FoldChange[match(rownames(pm),ddf$Gene)] # disease score

cdf=data.frame(treat=colnames(pm),rho=NA,pvalue=NA)

cdf$rho=apply(pm,2,function(psig) cor.test(psig,dsig,method="pearson")$estimate)
cdf$pvalue=apply(pm,2,function(psig) cor.test(psig,dsig,method="pearson")$p.value)
cdf$drug=tinfo$Perturbagen_cmap[match(cdf$treat,tinfo$Signature_id)]

cdf$n=length(dsig)
cdf$rho=pmin(pmax(cdf$rho,  -0.999999), 0.999999)
cdfl=split(cdf, cdf$drug)

cresl=lapply(names(cdfl), function(drug) {
  df=cdfl[[drug]]
  r=df$rho
  n=df$n

  keep=is.finite(r) & is.finite(n) & (n>3) & (abs(r)<1)
  if(!any(keep)) {
    return(data.frame(drug=drug, rho=NA_real_, pvalue=NA_real_, k_used=0L, stringsAsFactors=F))
  }

  r=r[keep]
  n=n[keep]
  z=atanh(r)
  vi=1/(n-3)

  dat = metafor::escalc(measure = "ZCOR", ri = r, ni = n)
  fit=metafor::rma(yi = dat$yi, vi = dat$vi, method = "FE") 
  pres=predict(fit, transf = tanh)

  data.frame(drug=drug, rho=as.numeric(pres$pred), pvalue=as.numeric(fit$pval), k_used=length(r), stringsAsFactors=F)
})
cdf=do.call(rbind, cresl)



########################################
## 2nd Drug repurposing: drug-target network analysis
########################################

ctl=setNames(vector("list",nrow(cdf)) , cdf$drug)


## cinfo2
ctl=lapply(names(ctl), function(drug) {
  targets=cinfo2$GeneSymbolOfTargets[grepl(drug,cinfo2$CompoundName) | grepl(drug,cinfo2$Synonyms) | grepl(str_to_title(drug),cinfo2$CompoundName) | grepl(str_to_title(drug),cinfo2$Synonyms) | grepl(toupper(drug),cinfo2$CompoundName) | grepl(toupper(drug),cinfo2$Synonyms)]
  unique(unlist(strsplit(targets,";",fixed=T)))
} )
names(ctl)=cdf$drug


## cinfo3
ctl=lapply(names(ctl), function(drug) {
  targets=cinfo3$gene_symbol_of_protein_target[grepl(drug,cinfo3$cpd_name) | grepl(str_to_title(drug),cinfo3$cpd_name) | grepl(toupper(drug),cinfo3$cpd_name)]
  unique(c(ctl[[drug]] , unlist(strsplit(targets,";",fixed=T))))
})
names(ctl)=cdf$drug


## ct
ct$drug=sapply(ct$drugs, function(c) c[[1]][[1]])
ct$target=gtf$gene_name[match(ct$targetFromSourceId,gtf$gene_id)]
ctl=lapply(names(ctl), function(drug) {
  targets=ct$target[grepl(drug,ct$drug) | grepl(str_to_title(drug),ct$drug) | grepl(toupper(drug),ct$drug) | grepl(tolower(drug),ct$drug)]
  unique(c(ctl[[drug]] , targets[!is.na(targets)]))
})
names(ctl)=cdf$drug


## drug-target network
degs=ddf$Gene[ddf$log2FoldChange>0 & ddf$deg==T]
ppi=rbind(ppi[ppi$source_genesymbol %in% unlist(ctl) & ppi$target_genesymbol %in% degs , c("source_genesymbol","target_genesymbol")],
          ppi[ppi$source_genesymbol %in% degs & ppi$target_genesymbol %in% unlist(ctl) , c("source_genesymbol","target_genesymbol")])
ctl=ctl[sapply(ctl,length)>0]
tdf=data.frame(target=rep(names(ctl), sapply(ctl,length)),source=unlist(ctl))
colnames(tdf)=colnames(ppi)

edf=as_edgelist(graph)
colnames(edf)=colnames(ppi)

ppi=do.call(rbind, list(tdf, edf, ppi))
ppi=ppi[apply(ppi,1, function(c) sum(is.na(c)))==0,]

graph=graph_from_data_frame(ppi, directed=F)
V(graph)$type[V(graph)$name %in% names(ctl)]="drug"
V(graph)$type[V(graph)$name %in% fdegs]="hub_signature"
V(graph)$type[is.na(V(graph)$type)]="gene"


## network analysis
nmode = ifelse(is_directed(graph), "all", "out") # directed일 경우 out
dst=distances(graph, v=V(graph)$name[V(graph)$type=="drug"], to=V(graph)$name[V(graph)$type=="hub_signature"], mode=nmode, weights=NA) # distance matrix
reach_bool = is.finite(dst) # 각 13개 gene reach 여부
coverage = rowSums(reach_bool) # coverage 총 몇 개의 유전자에 도달??  => range: 0~13
adst = apply(dst, 1, function(x) { # average distance
  x=x[is.finite(x)]
  if (length(x) == 0) Inf else mean(x)
})
ndf=data.frame(drug=V(graph)$name[V(graph)$type=="drug"], coverage=coverage, average_distance=adst, stringsAsFactors=F)



########################################
## Final Durg repurposing score compiling
########################################
cdrugs=unique(intersect(cdf$drug,ndf$drug))

cdf=cdf[cdf$drug %in% cdrugs,]
ndf=ndf[ndf$drug %in% cdrugs,]

cdf$padj=p.adjust(cdf$pvalue, method="BH")
cdf$score=log10(cdf$pvalue)*sign(cdf$rho)*(-1)*(-1)
cdf=cdf[order(cdf$pvalue),]
cdf=cdf[order(cdf$score, decreasing=T),]
# save(cdf, file=paste0(dir_data,"6. 1st Drug repurposing res (based drug signature).Rdata"))

ndf=ndf[order(ndf$average_distance),]
ndf=ndf[order(ndf$coverage, decreasing = T),]
ndf$re_average_distance=1/(ndf$average_distance)
# save(ndf, file=paste0(dir_data,"6. 2nd Drug repurposing res (based DTI network).Rdata"))


## Save Drug list used in analysis
drugs=unique(intersect(ndf$drug,cdf$drug))
sdf=cbind(ndf[match(drugs,ndf$drug),c("coverage","average_distance")] , cdf[match(drugs,cdf$drug),c("rho","pvalue","padj")])
sdf=cbind(data.frame(`Drug`=drugs) , sdf)
rownames(sdf)=1:nrow(sdf)
colnames(sdf)=c("Drug","Drug-target gene coverage","Average distance on Network","Disease signature pearson rho","pearson rho p","pearson rho adjusted p")
# write.csv(sdf, file="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/Supplementary Data/Supplementary Data 5.csv", row.names=F, quote=F)



########################################
## DTI network viz
########################################
ggraph=graph
degs=ddf$Gene[ddf$log2FoldChange<log2(2) & !is.na(ddf$log2FoldChange) & !is.na(ddf$padj) & ddf$padj<0.05]

load(paste0(dir_data,"2. Gene interaction Network analysis res.Rdata"))
V(ggraph)$cluster=ndf$cluster[match(V(ggraph)$name,ndf$Gene)]
V(ggraph)$type[V(ggraph)$name %in% setdiff(degs,fdegs)]="DEG"

createNetworkFromIgraph(ggraph, title = "Main Figure 6(A)", collection = "Main Figure 6(A)")



########################################
## Final Durg repurposing rank Viz
########################################

gdf=left_join(cdf,ndf)

pcols=setNames(c("#ff0000","#ff2424ff","#ff3434ff","#ff5555ff","#ff6f6fff","#ff8383ff","#000dff","#222dffff","#414bffff","#5d65ffff","#787fffff","#8d92ffff"), c(rep("Up",6),rep("Down",6)))

xcut=log10(0.05)*(-1)
ycut=0.3823529
gp=ggplot(gdf, aes(x=score, y = re_average_distance)) +
      theme_bw() +
      theme(plot.title = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(color="black",  size=20),
            axis.text.y = element_text(color="black",  size=20),
            axis.ticks.y = element_line(linewidth = 0.5,  color="black"), 
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5)) +
      geom_point(data=gdf, shape = 21, col = "#6f6f6f", fill = "darkgrey", alpha = 0.5, size = 3, stroke = 1) + 
      geom_point(data=gdf[gdf$score>0,], shape = 21, col = "#ffbfd2ff", fill = "#fff0f4ff", size = 3, stroke = 1) + 
      geom_point(data=gdf[gdf$score>=xcut,], shape = 21, col = "#ff417d", fill = "#ffd4e1", size = 4, stroke = 1) + 
      geom_point(data=gdf[gdf$re_average_distance>=ycut,], shape = 21, col = "#119901", fill = "#9eff5d", size = 4, stroke = 1) + 
      geom_point(data=gdf[gdf$score>xcut & gdf$re_average_distance>=ycut,], shape = 21, col = "#dda415", fill = "#ffe100", size = 5, stroke = 1) + 
      geom_vline(xintercept=c(xcut,0), linetype='dashed', color="black", linewidth = 0.7) +
      geom_hline(yintercept=ycut, linetype='dashed', color="black", linewidth = 0.7) +
      geom_text_repel(data=gdf[setdiff(which(gdf$score>xcut | gdf$re_average_distance>=ycut) , which(gdf$score>xcut & gdf$re_average_distance>=ycut)),], aes(label = drug), size=4.5, color = 'black', segment.size=0.4, force=10, nudge_x=rep(c(1,0.5,-0.5,-1),length(setdiff(which(gdf$score>xcut | gdf$re_average_distance>=ycut) , which(gdf$score>xcut & gdf$re_average_distance>=ycut))))) +
      geom_text_repel(data=gdf[gdf$score>xcut & gdf$re_average_distance>=ycut,], aes(label = drug), size=6, color = 'red', segment.size=0.4, force=10, nudge_x=2, nudge_y=0.02) 
# ggsave(filename=paste0(dir_fig, "Main Figure 6(B).tiff"), plot = gp, width=15, height=14, units = "cm", dpi = 300, limitsize = TRUE)


## legend
tiff(filename=paste0(dir_fig, "Main Figure 6(B)_legend.tif"), width=10, height = 7, units = 'cm',res=300)
par(mar=c(0,0,0,0),mgp=c(2,0.1,0), cex.lab=0.8, cex.axis=0.5, cex.main=0.8, tck=-0.02, las=1,bty="l")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", c("Both Significant","Significant Positive Correlation","Positive Correlation","Hub Signature reachability"), nc=1, pch=21, pt.bg=c("#ffe100","#ffd4e1","#fff0f4ff","#9eff5d"), col=c("#dda415","#ff417d","#ffbfd2ff","#119901"), bty="n", pt.lwd=2.2, pt.cex=3, cex=1.5, xpd=T)
dev.off()


################### network cut-off #####################
tiff(filename=paste0(dir_fig, "6. Main Figure 6(B)_Network-cut-off.tif"), width=12, height = 12, units = 'cm',res=300)
par(mar=c(1.5,1.5,0,0),mgp=c(2,0.1,0), cex.lab=0.8, cex.axis=1, cex.main=0.8, tck=-0.005, las=1,bty="l")
plot(sort(ndf$re_average_distance, decreasing=T), pch=21, cex=1, lwd=1, col='black', bg='#616161')
abline(h=ycut,lty=1,col="#119901")
dev.off()
#########################################################