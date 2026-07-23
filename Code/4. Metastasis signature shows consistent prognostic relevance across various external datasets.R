########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/figure/"



########################################
## Library
########################################
library(survival)
library(survminer)
library(dplyr)
library(clinfun)



########################################
## Load Data
########################################
load(paste0(dir_data,"0.3. external dataset(GEO).Rdata")) # variable name: sinfol, eml
load(paste0(dir_data,"0.2. sample information(TCGA-GDC).Rdata")) # variable name: sinfo
gsinfo=sinfo
load(paste0(dir_data,"3. Predicted Metastasis risk.Rdata")) # variable name: mdfl



########################################
## Continuous analysis:  external validation
########################################
## Analysis settings
atypes=c("mfs","rfs","os")
occurences=setNames(c("Metastasis Occurrence","Recurrence Occurrence","Vital Status"),atypes)
free_survivals=setNames(c("Metastasis Free Survival","Relapse Free Survival","Overall Survival"),atypes)
events=setNames(c("Metastasis","Recurrence","Dead"),atypes)
nons=setNames(c("No Metastasis","No Recurrence","Alive"),atypes)
endpoint_labels=setNames(c("Metastasis-free survival","Relapse-free survival","Overall survival"),atypes)
endpoint_short=setNames(c("MFS","RFS","OS"),atypes)

horizons=c(1,5,10)
dsets=names(sinfol)

time_divisor=setNames(rep(1,length(dsets)),dsets) # Survival-time units


## Output directories
dir_external=paste0(dir_fig,"Continuous external validation/")
dir_km_table=paste0(dir_external,"KM with risk table/")
dir_km_notable=paste0(dir_external,"KM without risk table/")

dir.create(dir_external,recursive=T,showWarnings=F)
dir.create(dir_km_table,recursive=T,showWarnings=F)
dir.create(dir_km_notable,recursive=T,showWarnings=F)


## Utility functions
safe_filename=function(x) {
    x=gsub("[\\\\/:*?\"<>|]","_",x)
    x=gsub("[[:space:]]+"," ",x)
    trimws(x)
}

format_p=function(x) {
    if(is.na(x)) return("NA")
    if(x<0.001) return(formatC(x,format="e",digits=2))
    formatC(x,format="f",digits=3)
}

find_first_column=function(x,candidates) {
    detected=candidates[candidates %in% colnames(x)]
    if(length(detected)==0) return(NA_character_)
    detected[1]
}

detect_covariate_columns=function(sinfo) {
    c(
        Age=find_first_column(sinfo,c("Age","age","AGE","age:ch1","Age:ch1")),
        Sex=find_first_column(sinfo,c("Sex","sex","SEX","Gender","gender","gender:ch1","Sex:ch1")),
        Stage=find_first_column(sinfo,c("Stage","stage","STAGE","Pathologic Stage","Pathological Stage","AJCC Stage","stage:ch1","Stage:ch1","pathologic stage:ch1")),
        Histology=find_first_column(sinfo,c("histology:ch1","Histology type:ch1","Histology Type:ch1","Histology:ch1","Histology","histology"))
    )
}

clean_age=function(x) {
    if(is.numeric(x)) return(as.numeric(x))
    x=trimws(as.character(x))
    x[x %in% c("","NA","N/A","Unknown","unknown","Not Reported","not reported")]=NA
    suppressWarnings(as.numeric(gsub("[^0-9.+-]","",x)))
}

clean_sex=function(x) {
    x=tolower(trimws(as.character(x)))
    x[x %in% c("m","male","man")]="Male"
    x[x %in% c("f","female","woman")]="Female"
    x[!x %in% c("Male","Female")]=NA
    factor(x,levels=c("Female","Male"))
}

clean_stage=function(x) {
    x=toupper(trimws(as.character(x)))
    x=gsub("PATHOLOGIC|PATHOLOGICAL|AJCC|STAGE","",x)
    x=gsub("[[:space:]_-]","",x)
    x[x %in% c("","NA","N/A","X","IS","UNKNOWN","NOTREPORTED","NOTAVAILABLE")]=NA
    output=rep(NA_character_,length(x))
    output[grepl("^4$|^IV",x)]="IV"
    output[grepl("^3$|^III",x)]="III"
    output[grepl("^2$|^II",x) & is.na(output)]="II"
    output[grepl("^1$|^I",x) & is.na(output)]="I"
    factor(output,levels=c("I","II","III","IV"))
}

clean_histology=function(x) {
    x=trimws(as.character(x))
    x[tolower(x) %in% c("","na","n/a","unknown","not reported","not available")]=NA
    factor(x)
}


## Construct endpoint data
make_external_endpoint_data=function(sinfo,dset,atype,horizon,time_divisor=1) {
    event_col=occurences[atype]
    time_col=free_survivals[atype]
    required_columns=c("Metastasis Risk",event_col,time_col)
    missing_endpoint=setdiff(required_columns,colnames(sinfo))

    if(length(missing_endpoint)>0) return(list(ok=F,reason=paste0("Missing endpoint column(s): ",paste(missing_endpoint,collapse=", "))))

    covariate_columns=detect_covariate_columns(sinfo)
    if(!is.finite(time_divisor) || time_divisor<=0) time_divisor=1
    sample_id=if("Sample" %in% colnames(sinfo)) as.character(sinfo$Sample) else rownames(sinfo)

    dat=data.frame(
        Sample=sample_id,
        RiskScore=suppressWarnings(as.numeric(sinfo[["Metastasis Risk"]])),
        Time_raw=suppressWarnings(as.numeric(sinfo[[time_col]])),
        Event_text=as.character(sinfo[[event_col]]),
        stringsAsFactors=F,check.names=F
    )

    if(!is.na(covariate_columns["Age"])) dat$Age=clean_age(sinfo[[covariate_columns["Age"]]])
    if(!is.na(covariate_columns["Sex"])) dat$Sex=clean_sex(sinfo[[covariate_columns["Sex"]]])
    if(!is.na(covariate_columns["Stage"])) dat$Stage=clean_stage(sinfo[[covariate_columns["Stage"]]])
    if(!is.na(covariate_columns["Histology"])) dat$Histology=clean_histology(sinfo[[covariate_columns["Histology"]]])

    event_text=tolower(trimws(dat$Event_text))
    dat$event_original=ifelse(event_text==tolower(events[atype]),1,ifelse(event_text==tolower(nons[atype]),0,NA))
    dat$time_year=dat$Time_raw/time_divisor
    dat=dat[!is.na(dat$RiskScore) & is.finite(dat$RiskScore) & !is.na(dat$time_year) & is.finite(dat$time_year) & dat$time_year>=0 & !is.na(dat$event_original),,drop=F]

    if(nrow(dat)==0) return(list(ok=F,reason="No complete score, time, and event observations",covariate_columns=covariate_columns))

    available_covariates=intersect(c("Age","Sex","Stage","Histology"),colnames(dat))
    usable_covariates=available_covariates

    repeat {
        vars=c("RiskScore","time_year","event_original",usable_covariates)
        x=dat[complete.cases(dat[,vars,drop=F]),,drop=F]

        if(nrow(x)==0) return(list(ok=F,reason="No complete observations after clinical-covariate filtering",covariate_columns=covariate_columns))

        invalid=usable_covariates[vapply(usable_covariates,function(v) {
            values=x[[v]]
            if(v=="Age") {
                values=values[is.finite(values)]
                length(values)<2 || !is.finite(sd(values)) || sd(values)==0
            } else {
                length(unique(values[!is.na(values)]))<2
            }
        },logical(1))]

        if(length(invalid)==0) break
        usable_covariates=setdiff(usable_covariates,invalid)
    }

    dat=x
    dropped_covariates=setdiff(available_covariates,usable_covariates)

    if("Sex" %in% colnames(dat)) dat$Sex=droplevels(factor(dat$Sex))
    if("Stage" %in% colnames(dat)) dat$Stage=droplevels(factor(dat$Stage,levels=c("I","II","III","IV")))
    if("Histology" %in% colnames(dat)) dat$Histology=droplevels(factor(dat$Histology))

    dat$time=pmin(dat$time_year,horizon)
    dat$event=as.integer(dat$time_year<=horizon & dat$event_original==1)

    if(length(unique(dat$event))<2) return(list(ok=F,reason="Only one event status after administrative censoring",covariate_columns=covariate_columns))
    if(sum(dat$event)<3) return(list(ok=F,reason="Fewer than three events",covariate_columns=covariate_columns))
    if(sum(dat$event==0)<3) return(list(ok=F,reason="Fewer than three non-events",covariate_columns=covariate_columns))

    list(ok=T,data=dat,covariates=usable_covariates,dropped_covariates=dropped_covariates,covariate_columns=covariate_columns)
}


## Adjusted HR per SD and score-only C-index
analyze_external_adjusted_score=function(dat,covariates) {
    score_sd=sd(dat$RiskScore,na.rm=T)
    if(!is.finite(score_sd) || score_sd==0) return(list(ok=F,reason="Risk score has zero or undefined standard deviation"))

    dat$ScoreSD=as.numeric(scale(dat$RiskScore))
    rhs=c("ScoreSD",covariates)
    cox_formula=as.formula(paste("Surv(time,event) ~",paste(rhs,collapse=" + ")))

    cox_fit=tryCatch(
        coxph(cox_formula,data=dat,ties="efron",x=T,y=T,model=T),
        error=function(e) NULL
    )

    if(is.null(cox_fit)) return(list(ok=F,reason="Adjusted Cox model failed"))

    cox_summary=summary(cox_fit)
    if(!"ScoreSD" %in% rownames(cox_summary$coefficients)) return(list(ok=F,reason="Score coefficient was not estimated"))

    hr=cox_summary$coefficients["ScoreSD","exp(coef)"]
    hr_p=cox_summary$coefficients["ScoreSD","Pr(>|z|)"]
    hr_low=cox_summary$conf.int["ScoreSD","lower .95"]
    hr_high=cox_summary$conf.int["ScoreSD","upper .95"]

    # Score-only C-index, but calculated in the same complete-case population
    cfit=tryCatch(
        concordance(Surv(time,event)~ScoreSD,data=dat,reverse=T,timewt="n"),
        error=function(e) NULL
    )

    if(is.null(cfit)) return(list(ok=F,reason="Score-only C-index estimation failed"))

    c_index=as.numeric(cfit$concordance)[1]
    c_variance=as.numeric(cfit$var)[1]

    if(is.finite(c_variance) && c_variance>=0) {
        c_se=sqrt(c_variance)
        c_low=max(0,c_index-1.96*c_se)
        c_high=min(1,c_index+1.96*c_se)
    } else {
        c_se=NA_real_
        c_low=NA_real_
        c_high=NA_real_
    }

    # KM curve uses all complete-case patients and median split
    median_cut=median(dat$RiskScore,na.rm=T)
    dat$Level=ifelse(dat$RiskScore<=median_cut,"Low","High")
    dat$Level=factor(dat$Level,levels=c("Low","High"))

    level_count=table(dat$Level)
    km_available=length(level_count)==2 && all(level_count>=3)
    sres=NULL
    logrank_p=NA_real_

    if(km_available) {
        sres=tryCatch(survfit(Surv(time,event)~Level,data=dat),error=function(e) NULL)
        logrank_test=tryCatch(survdiff(Surv(time,event)~Level,data=dat),error=function(e) NULL)
        if(!is.null(logrank_test)) logrank_p=pchisq(logrank_test$chisq,df=length(logrank_test$n)-1,lower.tail=F)
        if(is.null(sres)) km_available=F
    }

    result=data.frame( N=nrow(dat), Events=sum(dat$event), Non_events=sum(dat$event==0), Adjustment_variables=if(length(covariates)>0) paste(covariates,collapse=" + ") else "None", Score_mean=mean(dat$RiskScore), Score_SD=score_sd, Adjusted_HR_per_SD=hr, Adjusted_HR_low=hr_low, Adjusted_HR_high=hr_high, Adjusted_HR_P=hr_p, C_score=c_index, C_score_low=c_low, C_score_high=c_high, C_score_SE=c_se, KM_cut=median_cut, KM_low_N=as.integer(level_count["Low"]), KM_high_N=as.integer(level_count["High"]), KM_logrank_P=logrank_p, KM_available=km_available, stringsAsFactors=F, check.names=F )

    list(ok=T,data=dat,cox_fit=cox_fit,concordance=cfit,survfit=sres,summary=result)
}


## KM plot function
plot_external_adjusted_km=function(result,dset,atype,horizon,save_plot=T) {
    if(!result$summary$KM_available || is.null(result$survfit)) return(NULL)

    dat=result$data
    sres=result$survfit
    rdf=result$summary
    minimum_survival=min(sres$surv,na.rm=T)
    ymin=max(0,min(0.80,floor(minimum_survival*10)/10-0.05))
    if(!is.finite(ymin)) ymin=0

    annotation_x=0.03*horizon
    annotation_y=ymin+0.24*(1-ymin)
    break_value=ifelse(horizon==1,0.25,1)

    annotation_label=paste0(
        "Adjusted HR per SD = ",sprintf("%.3f",rdf$Adjusted_HR_per_SD),
        " (95% CI ",sprintf("%.3f",rdf$Adjusted_HR_low),"\u2013",sprintf("%.3f",rdf$Adjusted_HR_high),")\n",
        "Adjusted Cox P = ",format_p(rdf$Adjusted_HR_P),"\n",
        "Score-only C-index = ",sprintf("%.3f",rdf$C_score),
        " (95% CI ",sprintf("%.3f",rdf$C_score_low),"\u2013",sprintf("%.3f",rdf$C_score_high),")\n",
        "KM log-rank P = ",format_p(rdf$KM_logrank_P),"\n"
    )

    km_cols=c(Low="#403dff",High="#ff3d3d")
    plot_title=paste0(dset," | ",horizon,"-year ",endpoint_short[atype])

    gplot=ggsurvplot(
        sres,data=dat,conf.int=T,risk.table=T,
        legend.labs=c("Low","High"),
        legend.title="Median-split risk level",
        palette=unname(km_cols[c("Low","High")]),
        title=plot_title,xlab="Years",break.time.by=break_value,
        risk.table.height=0.30,surv.median.line="none",
        size=0.70,fontsize=4,censor=T,censor.shape="+",censor.size=4,
        ylim=c(ymin,1),xlim=c(0,horizon),pval=F
    )

    gplot$plot=gplot$plot +
        labs(y=paste0(endpoint_labels[atype]," probability")) +
        annotate("text",x=annotation_x,y=annotation_y,label=annotation_label,hjust=0,vjust=1,color="black",size=4.1,fontface="plain") +
        theme(
            plot.title=element_text(size=13,face="bold",hjust=0.5),
            axis.title=element_text(size=13),
            axis.text=element_text(size=12),
            legend.title=element_text(size=11),
            legend.text=element_text(size=11),
            panel.grid=element_blank()
        )

    gplot$table=gplot$table +
        labs(title="Number at risk",y="Risk level") +
        theme(plot.title=element_text(size=11,face="italic"),axis.title=element_text(size=11),axis.text=element_text(size=10))

    if(save_plot) {
        file_base=safe_filename(paste(dset,toupper(atype),paste0(horizon,"year"),sep="_"))

        tiff(
            filename=paste0(dir_km_table,file_base,".tiff"),
            width=16,height=17,units="cm",res=300,compression="lzw"
        )
        print(gplot)
        dev.off()

        ggsave(
            filename=paste0(dir_km_notable,file_base,".tiff"),
            plot=gplot$plot,width=15,height=13,units="cm",
            dpi=300,limitsize=F,compression="lzw"
        )
    }

    gplot
}


## Run function
external_result_list=list()
external_skip_list=list()
result_index=1
skip_index=1

for(dset in dsets) {
    cat("\n========================================\n")
    cat("Dataset:",dset,"\n")
    cat("Detected columns:\n")
    print(detect_covariate_columns(sinfol[[dset]]))
    cat("========================================\n")

    sinfo=sinfol[[dset]]
    divisor=time_divisor[dset]
    if(is.na(divisor)) divisor=1

    for(atype in atypes) {
        for(horizon in horizons) {
            cat(dset,"|",toupper(atype),"|",horizon,"year\n")

            endpoint_object=make_external_endpoint_data(
    sinfo=sinfo,dset=dset,atype=atype,horizon=horizon,time_divisor=divisor
)

            if(!endpoint_object$ok) {
                detected=if(!is.null(endpoint_object$covariate_columns)) paste(names(endpoint_object$covariate_columns),endpoint_object$covariate_columns,sep="=",collapse="; ") else NA_character_
                external_skip_list[[skip_index]]=data.frame(
                    Dataset=dset,Endpoint=atype,Horizon_year=horizon,
                    Reason=endpoint_object$reason,Detected_columns=detected,
                    stringsAsFactors=F,check.names=F
                )
                skip_index=skip_index+1
                next
            }

            analysis_object=analyze_external_adjusted_score(
                dat=endpoint_object$data,
                covariates=endpoint_object$covariates
            )

            if(!analysis_object$ok) {
                external_skip_list[[skip_index]]=data.frame(
                    Dataset=dset,Endpoint=atype,Horizon_year=horizon,
                    Reason=analysis_object$reason,
                    Detected_columns=paste(names(endpoint_object$covariate_columns),endpoint_object$covariate_columns,sep="=",collapse="; "),
                    stringsAsFactors=F,check.names=F
                )
                skip_index=skip_index+1
                next
            }

            temp=analysis_object$summary
            temp$Dataset=dset
            temp$Endpoint=atype
            temp$Horizon_year=horizon
            temp$Dropped_invariant_covariates=if(length(endpoint_object$dropped_covariates)>0) paste(endpoint_object$dropped_covariates,collapse=" + ") else "None"
            temp$Age_column=endpoint_object$covariate_columns["Age"]
            temp$Sex_column=endpoint_object$covariate_columns["Sex"]
            temp$Stage_column=endpoint_object$covariate_columns["Stage"]
            temp$Histology_column=endpoint_object$covariate_columns["Histology"]

            temp=temp[,c(
                "Dataset","Endpoint","Horizon_year",
                setdiff(colnames(temp),c("Dataset","Endpoint","Horizon_year"))
            )]

            external_result_list[[result_index]]=temp
            result_index=result_index+1

            if(analysis_object$summary$KM_available) {
                gplot=plot_external_adjusted_km(
                    result=analysis_object,dset=dset,
                    atype=atype,horizon=horizon,save_plot=T
                )
                if(!is.null(gplot)) print(gplot)
            }
        }
    }
}


## Combine results

if(length(external_result_list)>0) {
    external_adjusted_results=do.call(rbind,external_result_list)
    rownames(external_adjusted_results)=NULL
} else {
    external_adjusted_results=NULL
}

if(length(external_skip_list)>0) {
    external_adjusted_skip_log=do.call(rbind,external_skip_list)
    rownames(external_adjusted_skip_log)=NULL
} else {
    external_adjusted_skip_log=NULL
}


## Figure 4A and 4B candidates
if(!is.null(external_adjusted_results)) {
    figure4A_candidates=external_adjusted_results[
        external_adjusted_results$Endpoint=="rfs" &
        external_adjusted_results$Horizon_year==5,,drop=F
    ]

    figure4B_candidates=external_adjusted_results[
        external_adjusted_results$Endpoint=="os" &
        external_adjusted_results$Horizon_year==10,,drop=F
    ]
} else {
    figure4A_candidates=NULL
    figure4B_candidates=NULL
}


## Save
# save( external_adjusted_results, external_adjusted_skip_log, figure4A_candidates, figure4B_candidates, file=paste0(dir_data,"4. External adjusted HR score-only C-index and KM results.Rdata") )

if(!is.null(external_adjusted_results)) {
    write.csv(
        external_adjusted_results,
        paste0(dir_data,"4. External adjusted HR and score-only C-index results.csv"),
        row.names=F
    )
}

if(!is.null(figure4A_candidates)) {
    write.csv(
        figure4A_candidates,
        paste0(dir_data,"4. Figure 4A candidate adjusted 5-year RFS results.csv"),
        row.names=F
    )
}

if(!is.null(figure4B_candidates)) {
    write.csv(
        figure4B_candidates,
        paste0(dir_data,"4. Figure 4B candidate adjusted 10-year OS results.csv"),
        row.names=F
    )
}

if(!is.null(external_adjusted_skip_log)) {
    write.csv(
        external_adjusted_skip_log,
        paste0(dir_data,"4. External adjusted analysis skip log.csv"),
        row.names=F
    )
}



########################################
## Metastasis Risk score across clinical stages in LUAD
########################################
dsets=c("LUAD (GSE31210)","NSCLC (GSE50081)","NSCLC (GSE14814)","NSCLC (GSE42127)","Pan cancer")
column_stages=setNames(c('Pathologic Stage', 'Stage', 'Stage', 'Pathologic Stage', 'AJCC Stage') , dsets)

## Jonckheere-Terpstra Trend Test
jresl=c()
for(set in dsets) {
    # data processing
    if(set=='Pan cancer') {
        sinfo=gsinfo
        sinfo=left_join(sinfo,mdfl[['GDC']])
        colnames(sinfo)[colnames(sinfo)=='Metastasis Signature']='Metastasis Risk'
    } else
        {sinfo=sinfol[[set]]}
    column_stage=column_stages[set]
    sinfo=sinfo[!is.na(sinfo[,column_stage]) & !(sinfo[,column_stage] %in% c("IS","X","unknown",0,"","Not Reported")),]
    # test
    sinfo[,column_stage]=factor(sinfo[,column_stage] , levels=sort(unique(sinfo[,column_stage])), ordered=T)
    jresl[[set]]=jonckheere.test(sinfo$`Metastasis Risk`, sinfo[,column_stage])
}


## Viz
pcols=c("#6bd652","#ff84c1","#d67deb","#621a86")
xseq_multi=1.5

for(set in names(jresl)) {
    # data processing
    if(set=='Pan cancer') {
        data=gsinfo
        data=left_join(data,mdfl[['GDC']])
        colnames(data)[colnames(data)=='Metastasis Signature']='Metastasis Risk'
    } else
        {data=sinfol[[set]]}
    column_stage=column_stages[set]
    data=data[!is.na(data[,column_stage]) & !(data[,column_stage] %in% c("IS","X","unknown",0,"","Not Reported")),]

    pv=jresl[[set]]$p.value
    xseq=setNames(c(1:length(unique(data[,column_stage])))*xseq_multi , sort(unique(data[,column_stage])))
    data$xseq=xseq[match(data[,column_stage],names(xseq))]

    data$stage=data[,column_stage]
    data_halv=data_halve=filter(group_by(data,stage), n()>1)
    gp=ggplot(data, aes(x = xseq, y = `Metastasis Risk`, color=stage, fill = stage)) +
        geom_violin(width = .7) +
        geom_boxplot(width = .2, outlier.shape = NA, fill='white') +
        stat_summary(aes(group = 1), fun = median, geom = "line", color = "#ff00005d", linewidth = 0.45) +
        stat_summary(aes(group = 1), fun = median, geom = "point", color = "#ff0000", size=1.5) +
        geom_segment(aes(x=min(xseq), xend=max(xseq), y=1.7, yend=1.7), color="black", linewidth=0.5) +
        geom_segment(aes(x=min(xseq), xend=min(xseq), y=1.665, yend=1.7), color="black", linewidth=0.5) +
        geom_segment(aes(x=max(xseq), xend=max(xseq), y=1.665, yend=1.7), color="black", linewidth=0.5) +
        geom_text(x = mean(xseq), y = 1.82, label = ifelse(pv<(2.2*10^-16),"< 2.2e-16",ifelse(pv<0.00001,sprintf("%.3e", pv),round(pv,3))), color="black", hjust=0.5, size=5) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5, face="bold", size=17, color='black', margin=margin(b=20)), axis.title = element_text(size=17), axis.title.x = element_text(margin = margin(t = 25)), axis.title.y = element_text(margin = margin(r = 25)), axis.text = element_text(size=17, color='black'), legend.position = "none") +
        labs(title=set, y = "Metastasis Risk", x = "Stage") +
        scale_color_manual(values=pcols) +
        scale_fill_manual(values=paste0(pcols,"5d")) +
        scale_x_continuous(breaks = sort(unique(data$xseq)), labels = sort(unique(data$stage))) + coord_cartesian(clip = "off") +
        scale_y_continuous(limits = c(0, 1.7)) 
    ggsave(paste0(dir_fig,"Main Figure 4(C)/",set,".tiff"), plot=gp, width=ifelse(length(unique(data$stage))==4,10,5.8), height=12, units="cm", dpi=300)
}

# legend
tiff(filename=paste0(dir_fig, "Main Figure 4(C)/Legend.tiff"), width=15, height =4, units = 'cm',res=300)
par(mar=c(0,0,0,0),mgp=c(2,0.1,0), cex.lab=0.8, cex.axis=0.5, cex.main=0.8, tck=-0.02, las=1,bty="l")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", c("I","II","III","IV"), nc=4, pch=22, col=pcols, pt.bg=paste0(pcols,"5d"), bty="n", pt.lwd=2.7, pt.cex=4, cex=2.5, xpd=T)
dev.off()