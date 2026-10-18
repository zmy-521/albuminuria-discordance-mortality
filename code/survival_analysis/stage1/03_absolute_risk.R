# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Stage1");options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
des<-readRDS(file.path(S,'data/survey_design.rds'))
rd<-as.svrepdesign(des,type='JKn',mse=TRUE)
rd<-subset(rd,domain==1 & year<=2013)
# No outcome-dependent choice. spline basis and factor levels fixed by this baseline subcohort.
dat<-rd$variables;dat<-droplevels(dat)
kn<-ns(dat$age10[dat$ccB==1],df=3);kb<-attr(kn,'knots');bb<-attr(kn,'Boundary.knots')
f<-Surv(time,death)~group+ns(age10,knots=kb,Boundary.knots=bb)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
rd$variables<-dat
counter<-0
theta<-function(w,data){
 counter<<-counter+1
 if(counter%%30==0){cat('replicate',counter,'\n');flush.console()}
 data$.w<-w
 km<-sapply(paste0('G',0:3),function(g){d<-data[data$group==g&w>0,]; sf<-survfit(Surv(time,death)~1,data=d,weights=.w);1-summary(sf,times=60,extend=TRUE)$surv})
 names(km)<-paste0('KM_',paste0('G',0:3))
 dd<-data[data$ccB==1&w>0,]
 cf<-coxph(f,data=dd,weights=.w,ties='breslow',x=TRUE,model=TRUE,control=coxph.control(iter.max=50))
 stopifnot(all(is.finite(coef(cf))))
 bh<-basehaz(cf,centered=FALSE);h<-max(c(0,bh$hazard[bh$time<=60]))
 risks<-sapply(paste0('G',0:3),function(g){nd<-dd;nd$group<-factor(g,levels=levels(data$group));lp<-predict(cf,newdata=nd,type='lp',reference='zero');sum(dd$.w*(-expm1(-h*exp(lp))))/sum(dd$.w)})
 names(risks)<-paste0('AdjustedB_',paste0('G',0:3))
 c(km,risks,KM_RD_G2_G0=km[3]-km[1],AdjustedB_RD_G2_G0=risks[3]-risks[1],AdjustedB_RD_G2_G1=risks[3]-risks[2])
}
ans<-withReplicates(rd,theta,return.replicates=TRUE)
v<-coef(ans);se<-SE(ans);lo<-v-1.96*se;hi<-v+1.96*se
pr<-seq_along(v)<=8;lo[pr]<-plogis(qlogis(v[pr])-1.96*se[pr]/(v[pr]*(1-v[pr])));hi[pr]<-plogis(qlogis(v[pr])+1.96*se[pr]/(v[pr]*(1-v[pr])))
write.csv(data.frame(estimand=names(v),estimate=as.numeric(v),se=as.numeric(se),low=as.numeric(lo),high=as.numeric(hi),variance_method='stratified_delete_one_PSU_JKn'),file.path(S,'audit/risk60.csv'),row.names=FALSE)
saveRDS(ans,file.path(S,'data/risk60_replicates.rds'))
write.csv(data.frame(group=paste0('G',0:3),n=sapply(paste0('G',0:3),function(g)sum(dat$group==g)),deaths60=sapply(paste0('G',0:3),function(g)sum(dat$group==g&dat$death==1&dat$time<=60)),censored_before60=sapply(paste0('G',0:3),function(g)sum(dat$group==g&dat$death==0&dat$time<60)),completeB=sapply(paste0('G',0:3),function(g)sum(dat$group==g&dat$ccB==1))),file.path(S,'audit/risk60_counts.csv'),row.names=FALSE)
capture.output(str(ans,max.level=1),file=file.path(S,'audit/risk60_replication_log.txt'))
cat('absolute risks complete: ',counter,' fits\n')
