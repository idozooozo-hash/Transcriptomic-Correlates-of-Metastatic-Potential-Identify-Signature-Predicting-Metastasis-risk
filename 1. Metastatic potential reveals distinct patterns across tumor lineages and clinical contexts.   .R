########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/github/figure/"



########################################
## Library
########################################
library(ComplexHeatmap)
library(colorRamp2)
library(stringr)
library(car)



########################################
## Load Data
########################################
load(file=paste0(dir_data,"0. filtered_sample_information(DepMap).Rdata")) # variable name: sinfo
fsinfo=sinfo
load(file=paste0(dir_data,"0. all_sample_information(DepMap).Rdata")) # variable name: sinfo



########################################
## Metastasis Potential distribution by target organ
########################################
ttypes=colnames(fsinfo)[grep("meta to ",colnames(fsinfo))]

for(ttype in ttypes) {
    breaks=seq(floor(min(fsinfo[,ttype])) , ceiling(max(fsinfo[,ttype])), by=0.5)
    freqs=table(cut(fsinfo[,ttype], breaks=breaks, right=F))
    if(min(freqs)==0) {
        freqs=freqs[-1]
        breaks=breaks[-1]
    }
    if(max(freqs)==0) {
        freqs=freqs[-length(freqs)]
        breaks=breaks[-length(breaks)]
    }
    counts=as.numeric(freqs)

    width=1.7
    ymin=0
    ymax=250
    yseq=seq(ymin,ymax,50)
    ylim=range(ymin,ymax)
    tcols=setNames(c("#3d3aff","#ff3a3a") , c(-4,-2))

    tiff(filename=paste0(dir_fig,"Supplementary Figure 1. ",gsub("meta to ","",ttype),".tif"), width=23, height=20, units="cm", res=300)
    par(mfrow=c(1,1), mar=c(10,5,5,0.3), las=1, bty="l") 
    bp=barplot(freqs, main=gsub("meta to ","",ttype), xaxt="n", yaxt='n', xlab="", col="#ffe100", border="#000000", width=width, space=0, lwd=4, ylim=ylim, cex.main=2.5)
    xaxis=breaks
    sxseq=c(bp[breaks==-4]-width/2 , bp[breaks==-2]-width/2)
    text("Metastasis Potential", x=mean(bp), y=ymin-(ymax-ymin)/5, cex=2, xpd=T)
    text(xaxis, x=c(bp-width/2,max(bp)+width/2), y=ymin-(ymax-ymin)/23, xpd=T, srt=0, cex=ifelse(breaks %in% c(-4,-2) , 1.8 , 1.5), col=ifelse(breaks %in% c(-4,-2) , tcols[match(breaks,names(tcols))] , "black"))
    axis(side=2, at=yseq, label=NA, xpd=T, adj=1, line=-1.1, lwd=1.2, lwd.ticks=1.2, tcl=-0.3)
    mtext(yseq, at=yseq, side=2, line=-0.6, cex=1.7)
    text("Frequency", x=min(bp)-4.5, y=mean(yseq), cex=2, srt=90, xpd=T)
    segments(x0=sxseq , x1=sxseq , y0=ymin, y1=ymax, lty=2, lwd=3.5, col=c('#000dff','#ff0000'), xpd=T)
    text(counts, x=bp, y=freqs+(ymax-ymin)/25, cex=1.8, xpd=T)
    dev.off()
}




########################################
## Metastasis Potential Summary per target organ
########################################

## filter cancer type
ectypes=names(table(sinfo$`Cancer Type`)[table(sinfo$`Cancer Type`)>=10])
sinfo=sinfo[sinfo$`Cancer Type` %in% ectypes,]
sinfo$`Cancer Type`=factor(sinfo$`Cancer Type`, levels=names(sort(table(sinfo$`Cancer Type`)[table(sinfo$`Cancer Type`)>=10],decreasing = T)))


## correlation test
ttypes=colnames(sinfo)[grepl("meta to",colnames(sinfo))]
ctypes=names(sort(table(sinfo$`Cancer Type`),decreasing = T))

hmtx=t(as.data.frame(lapply(split(sinfo[,ttypes] , sinfo$`Cancer Type`), function(df) apply(df,2,mean))))
rownames(hmtx)=gsub("\\."," ",rownames(hmtx))
colnames(hmtx)=str_to_title(gsub("meta to ","",colnames(hmtx)))

tmtx=round(hmtx,3)

counts=table(sinfo$`Cancer Type`)
counts=as.numeric(counts[match(rownames(hmtx),names(counts))])


## viz
color.hp=colorRamp2(c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100])), c('#403dff','#d4d4d4','#ff3d3d'))
ran=rowAnnotation(`Count of\nSample` = anno_barplot(counts, width = unit(5, "cm"), gp=gpar(fill="#fff239", col="#ffee00", lwd=2.5), border=F, axis=F, add_numbers=T, numbers_gp=gpar(fontsize=30, col='#000000'), numbers_rot=0, numbers_offset=unit(0.3,"cm")),
                            annotation_name_offset=c(`Count of\nSample`="0.5cm"), annotation_name_gp=gpar(fontsize=30, col="#000000"), annotation_name_side='top', annotation_name_rot=0)
ran.gap = rowAnnotation(gap_dummy = anno_barplot(rep(0, length(counts)), height = unit(2.5, "cm"), width = unit(1, "cm"), gp = gpar(fill = "#ffffff", col = "#ffffff"), border = FALSE, axis=F),show_annotation_name = FALSE)
hp=Heatmap(hmtx, show_heatmap_legend=F, col=color.hp, border=F, cluster_columns=T, cluster_rows=T, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_row_dend=T, column_names_rot=90, column_gap=unit(0.7,"cm"), row_labels=rownames(hmtx), column_labels=gsub("meta to ","",colnames(hmtx)),
        width=ncol(hmtx)*unit(3.5,"cm"), height=nrow(hmtx)*unit(2,"cm"), rect_gp = gpar(col = "#ffffff", lwd = 2.2, lty=1), 
        row_names_gp = gpar(fontsize=42), column_title_gp = gpar(fontsize=40, fontface="bold"), column_names_gp = gpar(fontsize=47),
        column_dend_height=unit(3,"cm"), row_dend_width=unit(3,"cm"), row_dend_gp=gpar(lwd=3, col="#808080"), column_dend_gp=gpar(lwd=3, col="#808080"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(tmtx[i, j], x, y, gp = gpar(fontsize = 28, col = "#ffffff")
        )
        },
        right_annotation = c(ran,ran.gap))
tiff(filename=paste0(dir_fig,"Main Figure 1(a).tif"), width=50, height=50, units = 'cm',res=300)
draw(hp, ht_gap=unit(0.7,"cm"), padding=unit(c(0, 0, 0, 7), "cm")) 
dev.off()


# legend
lg.hp=Legend(at=c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100])), labels=c("Minimum",-3,"Maximum"), col_fun = color.hp, title = "Average\n metastasis potential", border='#808080', title_position="topcenter", size=unit(4,"cm"))
tiff(filename=paste0(dir_fig,"Main Figure 1(a)_legend.tiff"), width=5, height=8, units = 'cm',res=300)
draw(packLegend(list=list(lg.hp))) 
dev.off() 



########################################
## Difference in metastasis potential by target organ
########################################
psinfo=as.data.frame(sinfo)

ctypes=names(table(psinfo$`Cancer Type`))[table(psinfo$`Cancer Type`)>=10]
psinfo=psinfo[psinfo$`Cancer Type` %in% ctypes,]

cols=colnames(psinfo)[grepl("^meta to ",colnames(psinfo))]
psinfo=psinfo[,cols,drop=F]

pdf=do.call(rbind,lapply(colnames(psinfo), function(ttype) data.frame(target=str_to_title(gsub("meta to ","",ttype)),mp=psinfo[,ttype])))
pdf$target=ifelse(pdf$target=="All5","all5",pdf$target)
pdf=pdf[!is.na(pdf$mp),]
pdf$mtype=ifelse(pdf$mp<=(-4),"Non Metastatic",ifelse(pdf$mp>=(-2),"Metastatic","Weakly Metastatic"))
mtypes=unique(pdf$mtype)

# sorting
mres=aggregate(mp ~ target, data = pdf, median)
targets=mres$target[order(mres$mp)]
pdf$target=factor(pdf$target, levels=targets)
pdf=pdf[order(pdf$target),]
targets=levels(pdf$target)
xseq=seq(1,length(targets),1)


## test
lm=lm(mp~target, data=pdf)
ares=Anova(lm, type=2)
pv=ares[1,"Pr(>F)"]

tres=TukeyHSD(aov(mp~target, data=pdf))


## viz
pcols=setNames(c("#000dff","#a4a4a4","#ff0000"),c("Non Metastatic","Weakly Metastatic","Metastatic"))
pfils=setNames(c("#403dff73","#d4d4d473","#ff3d3d73"),c("Non Metastatic","Weakly Metastatic","Metastatic"))

ymin=min(pdf$mp)
ymax=max(pdf$mp)
yrange=ymax-ymin
ylim=c(ymin,ymax+yrange*0.05)

tiff(filename=paste0(dir_fig,"Main Figure 1(b).tif"), width=7.5, height=9.5, units="cm", res=300)
par(mfrow=c(1,1), mar=c(5,4.5,5.5,1), mgp=c(2.8,0.5,0), tck=-0.03, las=1, bty="l") 
bp=boxplot(pdf$mp~pdf$target, cex.lab=1.1, cex.main=1.3, cex.axis=0.9, lwd=1.2, boxwex=0.6, frame=T, xaxt='n', yaxt='n', border="#5a5a5aff", col=NA, outline=F, main="", xlab="", ylab="Metastasis potential", ylim=ylim)
title(main="Target organ", cex.main=0.85, line=4.3)
axis(2, at=c(-4,-2,0,2), labels=c(-4,-2,0,2), las=1, cex.axis=0.9)
# for(target in targets) {
#     for(mtype in mtypes) {
#         stripchart(pdf$mp[pdf$target==target & pdf$mtype==mtype]~pdf$target[pdf$target==target & pdf$mtype==mtype], at=which(targets==target), cex=1, lwd=1.4, pch=21, col=pcols[mtype], bg=pfils[mtype], method="jitter", jitter=0.24, vertical=TRUE, add = TRUE)
#     }
# }
pvlab=if(pv<2.2e-16 | pv==0) expression(italic(p) < 2.2 %*% 10^-16) else if(pv<0.00001) bquote(italic(p) == .(round(pv / 10^floor(log10(pv)), 4)) %*% 10^.(floor(log10(pv)))) else bquote(italic(p) == .(round(pv,4)))
ytop=par("usr")[4]
segments(y0=ytop+yrange*0.3,y1=ytop+yrange*0.3,x0=min(xseq),x1=max(xseq),col="#000000ff",lty=1,lwd=1.2,xpd=NA)
segments(y0=rep(ytop+yrange*0.27,2),y1=rep(ytop+yrange*0.3,2),x0=c(min(xseq),max(xseq)),x1=c(min(xseq),max(xseq)),col="#000000ff",lty=1,lwd=1.2,xpd=NA)
text(labels=pvlab, x=mean(xseq), adj=0.5, cex=0.7, y=ytop+yrange*0.37, xpd=NA, srt=0)
segments(y0=ytop+yrange*0.17,y1=ytop+yrange*0.17,x0=min(xseq),x1=max(xseq),col="#000000ff",lty=1,lwd=1.2,xpd=NA)
text(labels="**", x=mean(xseq), adj=0.5, cex=0.85, y=ytop+yrange*0.2, xpd=NA, srt=0)
segments(y0=ytop+yrange*0.05,y1=ytop+yrange*0.05,x0=min(xseq),x1=max(xseq)-1,col="#000000ff",lty=1,lwd=1.2,xpd=NA)
text(labels="**", x=mean(c(min(xseq),max(xseq)-1)), adj=0.5, cex=0.85, y=ytop+yrange*0.08, xpd=NA, srt=0)
text(labels=targets, x=xseq, adj=1, cex=0.7, y=ymin-yrange*0.08, xpd=NA, srt=60)
dev.off()



########################################
## Metastasis potential category by organ
########################################
cats = colnames(sinfo)[grep("meta type to ", colnames(sinfo))]
cats=c(setdiff(cats,"meta type to all5"),"meta type to all5")

meta_levels = c("Non Metastatic", "Weakly Metastatic (Low Confidence)", "Metastatic")
pcols = setNames(c("#403dff", "#c7c7c7", "#ff3d3d"), meta_levels)

pm = matrix(0, nrow=length(meta_levels), ncol=length(cats))
rownames(pm) = meta_levels
colnames(pm) = gsub("meta type to ", "", cats)

for(i in seq_along(cats)) {
    cat = cats[i]
    tmp = factor(sinfo[,cat], levels=meta_levels)
    gdf = table(tmp)
    pm[,i] = as.numeric(gdf / sum(gdf) * 100)
}


## Viz
organs_ordered=names(sort(sapply(gsub("type ","",cats) , function(type) median(sinfo[,type]))))
pm=pm[,gsub("meta to ","",organs_ordered)]
tiff(filename=paste0(dir_fig, "Main Figure 1(c).tif"), width=6, height=9, units="cm", res=300)
par(mar=c(5,1,2,1), las=1, bty="l")
bp = barplot(pm, main="Potential category", col=pcols[rownames(pm)], border="black", width=1, space=0, ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="", lwd=0.5)
text(x=bp, y=-3, labels=colnames(pm), srt=60, adj=1, xpd=TRUE, cex=1)
dev.off()


## fraction
pm=round(pm,1)[c("Metastatic","Weakly Metastatic (Low Confidence)","Non Metastatic"),]
#                                    brain bone liver kidney lung all5
# Metastatic                          23.6 21.7  24.6   31.0 27.9 53.3
# Weakly Metastatic (Low Confidence)  64.5 69.2  65.5   62.4 65.7 44.2
# Non Metastatic                      12.0  9.1   9.9    6.6  6.4  2.5



########################################
## Number of metastatic target organs
########################################

ctypes=names(table(sinfo$`Cancer Type`)[table(sinfo$`Cancer Type`)>=10])
psinfo=sinfo[sinfo$`Cancer Type` %in% ctypes , ] 

ctypes=names(sort(table(psinfo$`Cancer Type`),decreasing = T))
ttypes=c(0:5)

hmtx=matrix(NA, nrow=length(ctypes),ncol=length(ttypes), dimnames=list(ctypes,ttypes))

for(ctype in ctypes) {
        for(ttype in ttypes) {
                esinfo=psinfo[psinfo$`Cancer Type`==ctype , ]
                counts=apply(esinfo[,grep("meta type to ",colnames(psinfo))], 1, function(t) sum(t=="Metastatic"))  
                values=esinfo$`meta to all5`[counts==as.character(ttype)]                     
                hmtx[ctype,as.character(ttype)]=ifelse(sum(!is.na(values))==0,100,mean(values,na.rm=T))
        }
}

tmtx=round(hmtx,3)
tmtx[tmtx==100]="None"
tmtx_b=tmtx
tmtx_w=tmtx
tmtx_b[tmtx_b!="None"]=""
tmtx_w[tmtx_w=="None"]=""

color.hp=colorRamp2(c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100]),100), c('#403dff','#d4d4d4','#ff3d3d','#ffffffff'))
counts=table(psinfo$`Cancer Type`)
counts=as.numeric(counts[match(rownames(hmtx),names(counts))])
ran=rowAnnotation(`Count of\nSample` = anno_barplot(counts, width = unit(5, "cm"), gp=gpar(fill="#fff239", col="#ffee00", lwd=2.5), border=F, axis=F, add_numbers=T, numbers_gp=gpar(fontsize=30, col='#000000'), numbers_rot=0, numbers_offset=unit(0.3,"cm")),
                            annotation_name_offset=c(`Couunt of\nSample`="0.5cm"), annotation_name_gp=gpar(fontsize=26, col="#000000"), annotation_name_side='top', annotation_name_rot=0)
ran.gap = rowAnnotation(gap_dummy = anno_barplot(rep(0, length(counts)), height = unit(2.5, "cm"), width = unit(1, "cm"),, gp = gpar(fill = "#ffffff", col = "#ffffff"), border = FALSE, axis=F),show_annotation_name = FALSE)
hp=Heatmap(hmtx, show_heatmap_legend=F, col=color.hp, border=F, cluster_columns=F, cluster_rows=F, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_column_dend=F, row_title=NULL, column_title="Metastasis to",
        width=ncol(hmtx)*unit(3.5,"cm"), height=nrow(hmtx)*unit(2,"cm"), rect_gp = gpar(col = "#808080", lwd = 2.2, lty=1), column_names_rot=0,
        row_names_gp = gpar(fontsize=42), column_names_gp = gpar(fontsize=47), column_title_gp = gpar(fontsize=40),
        row_dend_gp=gpar(lwd=3, col="#808080"), row_dend_width=unit(4,"cm"), column_dend_gp=gpar(lwd=3, col="#808080"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(tmtx_b[i, j], x, y, gp = gpar(fontsize = 33, col = "#000000ff"))
                grid.text(tmtx_w[i, j], x, y, gp = gpar(fontsize = 33, col = "#ffffffff"))
                },
        right_annotation = c(ran,ran.gap)
        )

tiff(filename=paste0(dir_fig,"Main Figure 1(d).tif"), width=50, height=45, units='cm',res=300)
draw(hp, ht_gap=unit(0.7,"cm"), padding=unit(c(0, 0, 0, 6), "cm")) 
dev.off()


## legend
lg.hp=Legend(at=c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100])), labels=c(round(min(hmtx[hmtx!=100]),1),"-3",round(max(hmtx[hmtx!=100]),1)), col_fun = color.hp, title = "Average\n metastasis potential", border='#808080', title_position="topcenter", size=unit(4,"cm"))
tiff(filename=paste0(dir_fig,"Main Figure 1(d)_legend.tif"), width=5, height=8, units = 'cm',res=300)
draw(packLegend(list=list(lg.hp))) 
dev.off() 



########################################
## Confounding factor Identification 
########################################
psinfo=sinfo

ctypes=names(table(psinfo$`Cancer Type`))[table(psinfo$`Cancer Type`)>=9]
psinfo=psinfo[psinfo$`Cancer Type` %in% ctypes,]

psinfo$`Age Group`=floor(psinfo$Age*0.1)*10
psinfo$`Age Group`=ifelse(is.na(psinfo$`Age Group`),psinfo$`Age Group`,paste0(psinfo$`Age Group`,"s"))
cfs=c("Cell Line Origin","Age Group","Sex","Race")


## viz
mdf=summarise(group_by(psinfo,`Cancer Type`), median=median(`meta to all5`))
psinfo$`Cancer Type`=factor(psinfo$`Cancer Type`, levels=mdf$`Cancer Type`[order(mdf$median, decreasing=T)])
psinfo$`meta type to all5`=factor(psinfo$`meta type to all5`, levels=c("Non Metastatic","Metastatic","Weakly Metastatic (Low Confidence)"))
psinfo$`Cell Line Origin`=factor(psinfo$`Cell Line Origin`, levels=c("Primary","Metastatic"))
mtypes=unique(psinfo$`meta type to all5`)
pcols=setNames(c("#000dff","#a4a4a4","#ff0000"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))
pfils=setNames(c("#403dff73","#d4d4d473","#ff3d3d73"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))

for(cf in cfs) {
    for(ctype in ctypes) {
        gdf=psinfo[psinfo$`Cancer Type`==ctype,]
        gdf=gdf[,c("meta to all5","meta type to all5",cf)]
        gdf=gdf[!is.na(gdf[,cf]) & gdf[,cf]!="Unknown" & gdf[,cf]!="",]
        gdf[,cf]=factor(gdf[,cf], levels=sort(unique(gdf[,cf])))    
        gdf=gdf[order(gdf[,cf]),]
        ftypes=unique(gdf[,cf])

        if(length(unique(gdf[,cf]))<2) next

        ymin=min(gdf$`meta to all5`)
        ymax=max(gdf$`meta to all5`)*1.6
        ylim=c(ymin,ymax)
        xseq=seq(1,length(ftypes),1)

        tiff(filename=paste0(dir_fig,"1.14. EDA(pan)/1.14.2.1. confounding_factor_identification_per_cancer type/",cf,"_",ctype,".tif"), width=5+length(ftypes)*2.2, height=27, units="cm", res=300)
        par(mfrow=c(1,1), mar=c(25,10,10,2), mgp=c(5,1,0), las=1, bty="l") 
        bp=boxplot(gdf$`meta to all5`~gdf[,cf],  cex.lab=2.5, cex.main=2.9,  cex.axis=2.2, lwd=3, boxwex=0.47, frame=T, xaxt='n', border="#5a5a5aff", col=NA, outline=F, main=cf, xlab="", ylab="Metastasis Potential", ylim=ylim)
        for(ftype in ftypes) {
            for(mtype in mtypes) {
                stripchart(gdf$`meta to all5`[gdf[,cf]==ftype & gdf$`meta type to all5`==mtype]~gdf[,cf][gdf[,cf]==ftype & gdf$`meta type to all5`==mtype], cex=1.5, lwd=1.4, pch=21, col=pcols[mtype], bg=pfils[mtype], method="jitter", jitter=0.19, vertical=TRUE, add = TRUE)
            }
        }
        lm=lm(gdf[,"meta to all5"]~gdf[,cf])
        ares=Anova(lm,type=2)    
        pv=ares[1,"Pr(>F)"]
        if (pv<0.0001) pv=bquote(.(round( pv / 10^floor(log10(pv)), 4)) %*% 10^.(floor(log10(pv)))) else pv=round(pv,4)
        segments(y0=ymax+(ymax-ymin)*0.08, y1=ymax+(ymax-ymin)*0.08 , x0=min(xseq) , x1=max(xseq),col="#000000ff",lty=1,lwd=3, xpd=T)
        segments(y0=rep(c(ymax+(ymax-ymin)*0.08),2),y1=rep(ymax+(ymax-ymin)*0.06,2),x0=c(min(xseq),max(xseq)),x1=c(min(xseq),max(xseq)),col="#000000ff",lty=1,lwd=3, xpd=T)
        text(labels=pv, x=mean(xseq), adj=0.5, cex=2.5, y=ymax+(ymax-ymin)*0.17, xpd=T, srt=0)
        text(labels=gsub("_"," ",ftypes), x=xseq, adj=1, cex=2, y=ymin-0.5, xpd=T, srt=60)
        dev.off()
    }
}