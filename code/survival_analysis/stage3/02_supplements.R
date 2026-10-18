# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Stage3_finalization");options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
des<-readRDS(.release_path("ODKD_WORK_ROOT", "Stage1/data/survey_design.rds"));v<-des$variables
des<-update(des,log2alb=log2(ualb_mgl),after24=time-24)
fb<-Surv(time,death)~group+ns(age10,df=3)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
out<-list();fits<-list()
extract<-function(f,name,ds){b<-coef(f);V<-vcov(f);cc<-list(G2_G0=c(groupG2=1),G1_G0=c(groupG1=1),G3_G0=c(groupG3=1),G2_G1=c(groupG2=1,groupG1=-1));do.call(rbind,lapply(names(cc),function(k){if(!all(names(cc[[k]])%in%names(b)))return(NULL);a<-setNames(rep(0,length(b)),names(b));a[names(cc[[k]])]<-cc[[k]];e<-sum(a*b);se<-sqrt(drop(t(a)%*%V%*%a));data.frame(analysis=name,contrast=k,n=nrow(ds$variables),events=sum(ds$variables$death),HR=exp(e),low=exp(e-1.96*se),high=exp(e+1.96*se),p=2*pnorm(-abs(e/se)))}))}
base<-v$domain==1&v$ccB==1
for(k in c('negative_domain_B','negative_domain_B_log2alb','landmark24','sex1','sex2')){
 keep<-base
 if(grepl('negative',k))keep<-keep&!is.na(v$ualb_mgl)&v$ualb_mgl<30
 if(k=='landmark24')keep<-keep&v$time>24
 if(k=='sex1')keep<-keep&v$sex=='1'
 if(k=='sex2')keep<-keep&v$sex=='2'
 ds<-des[which(keep),];ds$variables<-droplevels(ds$variables);f<-fb
 if(k=='negative_domain_B_log2alb')f<-update(f,.~.+log2alb)
 if(k=='landmark24')f<-update(f,Surv(after24,death)~.)
 if(grepl('sex',k))f<-update(f,.~.-sex)
 fit<-svycoxph(f,design=ds,ties='breslow',x=TRUE,model=TRUE);fits[[k]]<-fit;out[[k]]<-extract(fit,k,ds);cat(k,'complete\n')
}
ds<-des[which(base),];ds$variables<-droplevels(ds$variables)
fi<-update(fb,.~.+group:sex); interaction<-svycoxph(fi,design=ds,ties='breslow',model=TRUE,x=TRUE);jt<-regTermTest(interaction,~group:sex)
capture.output(jt,file=file.path(S,'audit/sex_interaction_test.txt'));write.csv(data.frame(test='group_by_sex_joint',df=jt$df,denominator_df=jt$ddf,F=unname(jt$Ftest),p=jt$p),file.path(S,'audit/sex_interaction.csv'),row.names=FALSE)
write.csv(do.call(rbind,out),file.path(S,'audit/supplemental_Cox.csv'),row.names=FALSE);saveRDS(c(fits,list(interaction=interaction)),file.path(S,'models/supplemental_Cox.rds'))
dd<-v[v$domain==1,];counts<-aggregate(cbind(n=rep(1,nrow(dd)),events=dd$death,completeB=dd$ccB,events_B=dd$death*dd$ccB),by=list(sex=dd$sex,group=dd$group),FUN=sum);write.csv(counts,file.path(S,'audit/sex_counts.csv'),row.names=FALSE)
lk<-dd$time>24;write.csv(data.frame(n_total=nrow(dd),eligible=sum(lk),excluded_death_by24=sum(dd$death==1&dd$time<=24),excluded_censored_by24=sum(dd$death==0&dd$time<=24),B_eligible=sum(lk&dd$ccB==1),B_events=sum(dd$death[lk&dd$ccB==1])),file.path(S,'audit/landmark_flow.csv'),row.names=FALSE)
res<-list()
for(t in c(30,20))for(dn in c('concentration_negative','ACR_positive')){
 keep<-v$domain==1&!is.na(v$ualb_mgl)&if(dn=='concentration_negative')v$ualb_mgl<t else v$acr_adjusted>=30
 q<-des[which(keep),];q<-update(q,g2=as.numeric(ualb_mgl<t&acr_adjusted>=30));est<-svyciprop(~g2,q,method='logit',na.rm=TRUE);ci<-confint(est)
 res[[length(res)+1]]<-data.frame(source='NHANES5905',threshold=t,denominator=dn,n=nrow(q$variables),numerator=sum(q$variables$g2),unweighted=mean(q$variables$g2),weighted=coef(est)[1],low=ci[1],high=ci[2],design_df=degf(q))
}
write.csv(do.call(rbind,res),file.path(S,'audit/weighted_classification.csv'),row.names=FALSE)
cat('Stage3 supplemental analyses complete\n')

capture.output(warnings(),file=file.path(S,'audit/R_warnings.txt'))
