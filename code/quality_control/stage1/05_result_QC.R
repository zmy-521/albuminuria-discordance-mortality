# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines);S<-.release_path("ODKD_WORK_ROOT", "Stage1")
des<-readRDS(file.path(S,'data/survey_design.rds'));z<-subset(des,domain==1)
z<-update(z,female=as.numeric(sex=='2'))
out<-list()
for(g in paste0('G',0:3)){
 dz<-subset(z,group==g)
 for(v in c('age','female','bmi','hba1c','egfr','ucr_adjusted','ualb_mgl','pir','htn','cvd')){
  m<-svymean(as.formula(paste0('~',v)),dz,na.rm=TRUE)
  out[[length(out)+1]]<-data.frame(group=g,variable=v,weighted_mean=as.numeric(coef(m)),SE=as.numeric(SE(m)))
 }
}
write.csv(do.call(rbind,out),file.path(S,'audit/weighted_description.csv'),row.names=FALSE)
ans<-readRDS(file.path(S,'data/risk60_replicates.rds'));V<-vcov(ans);v<-coef(ans)
stopifnot(all(is.finite(ans$replicates)),max(abs(ans$replicates[,9]-(ans$replicates[,3]-ans$replicates[,1])))<1e-12,max(abs(ans$replicates[,10]-(ans$replicates[,7]-ans$replicates[,5])))<1e-12)
stopifnot(abs(V[10,10]-(V[7,7]+V[5,5]-2*V[7,5]))<1e-12)
# Independent cumulative weighted Breslow baseline and explicit model matrix standardization.
d<-droplevels(subset(z,year<=2013&ccB==1)$variables)
fa<-Surv(time,death)~group+ns(age10,df=3)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
f<-coxph(fa,data=d,weights=weight,ties='breslow',x=TRUE,model=TRUE)
eta<-drop(f$x%*%coef(f));times<-sort(unique(d$time[d$death==1&d$time<=60]));hh<-sum(sapply(times,function(t)sum(d$weight[d$death==1&d$time==t])/sum(d$weight[d$time>=t]*exp(eta[d$time>=t]))))
rr<-sapply(paste0('G',0:3),function(g){nd<-d;nd$group<-factor(g,levels=levels(d$group));X<-model.matrix(f,data=nd);sum(d$weight*(-expm1(-hh*exp(drop(X%*%coef(f))))))/sum(d$weight)})
stopifnot(max(abs(rr-v[5:8]))<1e-9)
write.csv(data.frame(check=c('JKn_replicates_all_finite','paired_RD_identity_every_replicate','RD_variance_covariance_identity','independent_Breslow_standardization'),result='PASS',n=c(nrow(ans$replicates),nrow(ans$replicates),1,length(d$SEQN))),file.path(S,'audit/risk_independent_QC.csv'),row.names=FALSE)
# Five-year risk estimate models are descriptive, not causal; retain fitted model and PH coefficients.
saveRDS(f,file.path(S,'data/risk60_standardization_model.rds'))
pf<-readRDS(file.path(S,'data/PH_model.rds'));b<-coef(pf);vv<-vcov(pf);rows<-list()
for(g in 1:3)for(period in c('0_60','after60')){
 aa<-rep(0,length(b));names(aa)<-names(b);aa[paste0('groupG',g)]<-1
 if(period=='after60')aa[paste0('g',g,'late')]<-1
 val<-sum(aa*b);se<-sqrt(drop(t(aa)%*%vv%*%aa));rows[[length(rows)+1]]<-data.frame(group=paste0('G',g),period,HR=exp(val),low=exp(val-1.96*se),high=exp(val+1.96*se))
}
write.csv(do.call(rbind,rows),file.path(S,'audit/PH_period_HR.csv'),row.names=FALSE)
cat('Risk independent QC PASS\n')
