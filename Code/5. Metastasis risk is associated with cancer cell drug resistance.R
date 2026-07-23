########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/figure/"

dir_dep="G:/My Drive/DB/DepMap/" # The location where the DepMap data is stored



########################################
## Library
########################################
library(GSVA)
library(dplyr)
library(ggplot2)



########################################
## Load Data
########################################
load(paste0(dir_dep,"tpm_24Q4.Rdata")) # DepMap cancer cell line gene expression (TPM) dataset; variable name: tpm

load(paste0(dir_dep,"auc_CTD2_230723.Rdata")) # DepMap Drug sensitivity dataset; variable name: auc 
auc[1:5,1:5] # row: Drug; column: cell-line
#               ACH-000002 ACH-000004 ACH-000005 ACH-000006 ACH-000007
# CIL56                 NA     11.733         NA         NA     14.562
# FQI-1                 NA     10.825         NA         NA     11.469
# BRD-K92856060     13.349     13.911     13.646         NA     14.418
# B02               12.987     13.445     12.313         NA     12.920
# BRD6340           14.904     14.193     13.739     14.016     14.742

length(intersect(colnames(tpm),colnames(auc))) # No. of Cancer Cell Lines Used in the Analysis: 800

load(paste0(dir_data,"2. Metastasis biomarker.Rdata")) # variable name: fdegs



########################################
## Drug sensitivity association analysis across 545 compounds
########################################

## Metastasis Risk calculation based <Orignial Metastasis Signature>
param=ssgseaParam(tpm , list(`Metastasis Signature`=fdegs))
mdf=t(as.data.frame(gsva(param,verbose=T)))


## Metastasis Risk calculation using <Random Gene set>
tpm = tpm[ , colnames(tpm) %in% colnames(auc)]
param=ssgseaParam(tpm, list('Metastasis Risk'=fdegs))
mdf=cbind(data.frame(Sample=colnames(tpm)),t(as.data.frame(gsva(param,verbose=T))))

gsl=c()
for(i in 1:1999) {
    print(i)

    gis=sample(1:nrow(tpm), size=length(fdegs), replace=F) # Random Gene set
    param=ssgseaParam(tpm , setNames(list(rownames(tpm)[gis]),paste0("Random ",i)))
    mdf=cbind(mdf , t(as.data.frame(gsva(param,verbose=T))))
    gsl[[paste0("Random ",i)]]=rownames(tpm)[gis]
}
# save(mdf, gsl, file=paste0(dir_data,"5. Predicted Metastasis risk score with random gene set.Rdata"))


## Wilcoxon Rank-Sum Tests for Individual Drug Sensitivities
wdfl=c()
for(gset in setdiff(colnames(mdf),"Sample")) {
    print(gset)

    wresl=lapply(rownames(auc), function(drug) {
        aucs=auc[drug,]
        cuts=quantile(aucs,c(1/4,3/4),na.rm=T)

        rid=colnames(auc)[aucs>cuts[2] & !is.na(aucs)] # Drug resistant ID
        sid=colnames(auc)[aucs<cuts[1] & !is.na(aucs)] # Drug sensitive ID

        rpots=mdf[mdf$Sample %in% rid,gset]
        spots=mdf[mdf$Sample %in% sid,gset]

        if(length(rpots)>2 & length(spots)>2) {
            wilcox.test(rpots,spots,alternative="greater")
        } else {
            NA
        }
    })

    pvs=sapply(wresl, function(wres) ifelse(class(wres)=="htest",wres$p.value,NA))
    wdf=data.frame(Drug=rownames(auc), pvalue=pvs, padj=p.adjust(pvs, method="BH"))
    wdf=wdf[order(wdf$pvalue),]
    wdfl[[gset]]=wdf[order(wdf$padj),]
}
# save(wdfl, file=paste0(dir_data,"5. Metastasis risk Difference for Individual Drug Sensitivities.Rdata"))


## Metastasis Signature Rank
wdf=wdfl[["Metastasis Signature"]]

sum(is.na(wdf$padj)) # Drugs that could not be tested: 0
sum(wdf$padj>=0.05,na.rm=T) # 158
sum(wdf$padj<0.05,na.rm=T) # 387
sum(wdf$padj<0.01,na.rm=T) # 363
sum(wdf$padj<0.001,na.rm=T) # 334


## Viz: Drug Ranking by Adjusted P-value
cutoff = -log10(0.05) # Significant cut-off
pcols = c("Significant" = "#ffe100", "Non-Significant" = "#e0e0e0") # plot color

pwdf=wdf
pwdf$logpadj = -log10(pwdf$padj)
pwdf$Drug=1:nrow(pwdf)
pwdf=mutate(pwdf, group=ifelse(pwdf$logpadj > (-log10(0.05)), "Significant", "Non-Significant"))

gp = ggplot(pwdf, aes(x = Drug, y = logpadj)) +
  geom_area(fill = "#d3d3d3", color = NA) +  # 기본 배경
  geom_area(aes(fill = group), alpha = 1, color = "#808080", linewidth = 1.2) + 
  scale_fill_manual(values = pcols) +
  labs(x = "Drugs", y = "-log10(Adjusted p-value)", fill = NULL) +
  theme_classic() +
  theme(
    plot.margin = margin(l = 20, b = 7, t = 1, r=0.3),
    legend.position = "none",
    axis.title = element_text(size = 30, color='black'),
    axis.title.y = element_text(margin = margin(r = 30)),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 30, margin = margin(t = -25), color='black'),
    axis.text.y = element_text(size = 30, margin = margin(r = -48), color='black'),
    axis.line = element_blank()
  )
# ggsave(paste0(dir_fig,"Main Figure 5(A)_waterfallplot.tiff"), plot=gp, width=40, height=22, units="cm", dpi=300)

# legend
tiff(filename=paste0(dir_fig, "Main Figure 5(A)_waterfallplot(legend).tif"), width=12, height =1.5, units = 'cm',res=300)
par(mar=c(0,0,0,0),mgp=c(2,0.1,0), cex.lab=0.8, cex.axis=0.5, cex.main=0.8, tck=-0.02, las=1,bty="l")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", c("Significant","Non-Significant"), nc=2, pch=22, col="#808080", pt.bg=pcols, bty="n", pt.lwd=2.4, pt.cex=2.7, cex=1.4, xpd=T)
dev.off()


## Viz: Example of Top 3 Drug
pcuts=c(0.01,0.001)
pcols=c("#6532f9","#ff417d")

wdf=wdfl[["Metastasis Signature"]]

for(drug in wdf$Drug[order(wdf$padj)][1:3]) {
    wpv=unique(wdf$padj[wdf$Drug==drug])

    aucs=auc[drug,]
    cuts=quantile(aucs,c(1/4,3/4),na.rm=T)

    rid=colnames(auc)[aucs>cuts[2] & !is.na(aucs)] # resistant ID
    sid=colnames(auc)[aucs<cuts[1] & !is.na(aucs)] # sensitive ID

    rpots=mdf$`Metastasis Signature`[match(rid,mdf$Sample)]
    spots=mdf$`Metastasis Signature`[match(sid,mdf$Sample)]

    gdf = data.frame(`Sample`=c(rid,sid),
                    `Group`=c(rep("Resistant",length(rid)),rep("Sensitive",length(sid))),
                    `Metastasis Potential`=c(rpots,spots),
                    check.names=F, row.names=1:length(c(rid,sid)))
    gdf=gdf[apply(gdf,1,function(v) sum(is.na(v)))==0,]

    ldf=summarise(group_by(gdf, Group), mean=mean(`Metastasis Potential`))
    ldf=mutate(rowwise(ldf), density=list(density(gdf$`Metastasis Potential`[gdf$Group == Group])), dy=density$y[which.min(abs(density$x-mean))])

    ymax=max(unlist(lapply(ldf$density, function(den) den$y))) * 1.1

    gp=ggplot(gdf, aes(x=`Metastasis Potential`, color=Group, fill=Group)) +
        geom_density(position="identity", linewidth=1.2) +
        geom_segment(data=ldf,aes(x=mean, xend=mean, y=0, yend=dy, color=Group), linetype="dashed", linewidth=0.7) +
        geom_segment(aes(x=min(ldf$mean), xend=max(ldf$mean), y=ymax, yend=ymax), color="black", linewidth=0.5) +
        geom_text(aes(x = mean(c(min(ldf$mean),max(ldf$mean))), y = ymax*1.02, label = ifelse(wpv<pcuts[2],"**","*")), color="black", hjust=0.5, size=10)+
        labs(title=drug, y="Density")+
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5, face="bold", size=19, color="black", margin=margin(b=20)),
            axis.title = element_text(size=17), axis.title.x = element_text(margin = margin(t = 25)), axis.title.y = element_text(margin = margin(r = 25)),
            axis.text = element_text(size=15), legend.position="none",
            plot.margin=margin(20,15,20,20)) +
        scale_color_manual(values=pcols) +
        scale_fill_manual(values=paste0(pcols,"3a"))
    ggsave(paste0(dir_fig,"Main Figure 5(A)_Drug sensitivity-Metastasis Risk association(",drug,").tiff"), plot=gp, width=14, height=14, units="cm", dpi=300)
}

# legend
tiff(filename=paste0(dir_fig, "Main Figure 5(A)_Drug sensitivity-Metastasis Risk association_legend.tif"), width=6, height =0.8, units = 'cm',res=300)
par(mar=c(0,0,0,0),mgp=c(2,0.1,0), cex.lab=0.8, cex.axis=0.5, cex.main=0.8, tck=-0.02, las=1,bty="l")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", c("Resistant","Sensitive"), nc=2, pch=22, col=pcols, pt.bg=paste0(pcols,"3a"), bty="n", pt.lwd=2.4, pt.cex=2.5, cex=1.2, xpd=T)
dev.off()



########################################
## Performance 'Metastasis Signature' vs. 'Random gene sets'
########################################

## performance
sum(sapply(wdfl, function(wdf) sum(wdf$padj<0.05)>=387)) # 34-1=33
sum(sapply(wdfl, function(wdf) sum(wdf$padj<0.05)<387)) # 1966
# 34/2000*100=1.7 => Top 1.7% (0.017)


## Viz
alpha=0.05 # Significant cut-off

pdf=data.frame(Gene_Set=names(wdfl), No._of_Significant_Drug=sapply(wdfl, function(wdf) sum(wdf$padj<0.05)))
pdf$Group=ifelse(pdf$No._of_Significant_Drug>=pdf$No._of_Significant_Drug[pdf$Gene_Set=="Metastasis Signature"] & pdf$Gene_Set!="Metastasis Signature","Outperforming Random Signature",ifelse(pdf$No._of_Significant_Drug<pdf$No._of_Significant_Drug[pdf$Gene_Set=="Metastasis Signature"],"Underperforming Random Signature","Metastasis Signature"))
groups=c("Underperforming Random Signature","Outperforming Random Signature")

pdf=pdf[order(pdf$No._of_Significant_Drug),]
width=0.3
ymin=0
ymax=545
yseq=seq(ymin,ymax,100)
ylim=range(ymin,ymax)
breaks=seq(0,2000,500)
xaxis=breaks
sxseq=max(which(pdf$No._of_Significant_Drug==387))
counts=pdf$No._of_Significant_Drug

lpos=pdf$No._of_Significant_Drug[pdf$Gene_Set=="Metastasis Signature"]
den=density(pdf$No._of_Significant_Drug)
dy=den$y[ which.min(abs(den$x - lpos)) ]
apos=sort(pdf$No._of_Significant_Drug,decreasing=T)[nrow(pdf)*alpha]
ady=den$y[ which.min(abs(den$x - apos)) ]

gp=ggplot(pdf, aes(x=No._of_Significant_Drug)) +
    geom_density(position="identity", linewidth=1.2, color="#a8cbff", fill="#a8cbff") +
    geom_segment(aes(x=lpos, xend=lpos, y=0, yend=dy), inherit.aes = FALSE, color="red", linewidth=0.5, linetype="dashed") +
    geom_segment(aes(x=apos, xend=apos, y=0, yend=ady), inherit.aes = FALSE, color="#0011ff", linewidth=0.5, linetype="solid") +
    theme_classic() +
    theme(axis.title = element_blank(),
        axis.text.x = element_text(size=15, color='black', margin=margin(t=1)),
        axis.text.y = element_text(size=15, color='black', margin=margin(r=1)),
        legend.position="none",
        plot.margin = margin(l = 20, b = 7, t = 1, r=0.3),
        axis.line = element_blank())
# ggsave(paste0(dir_fig,"5. Main Figure 5(B).tiff"), plot=gp, width=14, height=14, units="cm", dpi=300)