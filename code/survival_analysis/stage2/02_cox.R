# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Stage2");P<-.release_path("ODKD_WORK_ROOT", "Stage1");options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
des<-readRDS(file.path(P,'data/survey_design.rds'));d<-des$variables
d$group20<-factor(ifelse(d$ualb_mgl>=20,ifelse(d$acr_adjusted>=30,'G3','G1'),ifelse(d$acr_adjusted>=30,'G2','G0')),levels=paste0('G',0:3))
d$log2ucr<-log2(d$ucr_adjusted);d$log2alb<-log2(d$ualb_mgl)
d$ccNoPIR<-complete.cases(d[,c('age10','sex','race','cycle','bmi5','hba1c','egfr10','smoking','htn','cvd')]);des$variables<-d
fa<-Surv(time,death)~group+ns(age10,df=3)+sex+race+cycle
fb<-update(fa,.~.+bmi5+hba1c+egfr10+smoking+htn+cvd+pir)
extract<-function(f,threshold,model,dat){
 b<-coef(f);V<-vcov(f);cs<-list(G1_G0=c(groupG1=1),G2_G0=c(groupG2=1),G3_G0=c(groupG3=1),G2_G1=c(groupG2=1,groupG1=-1))
 if(model=='components')cs<-list(UAlb_per_doubling=c(log2alb=1),UCr_per_doubling=c(log2ucr=1))
 do.call(rbind,lapply(names(cs),function(nm){a<-setNames(rep(0,length(b)),names(b));a[names(cs[[nm]])]<-cs[[nm]];q<-sum(a*b);se<-sqrt(drop(t(a)%*%V%*%a));data.frame(threshold,model,contrast=nm,n=nrow(dat),events=sum(dat$death),logHR=q,se,HR=exp(q),low=exp(q-1.96*se),high=exp(q+1.96*se),p=2*pnorm(-abs(q/se)))}))
}
rows<-list();fits<-list()
for(t in c(30,20))for(mod in c('A','A_Bsample','B')){
 ds<-des;ds$variables$group<-if(t==30)d$group else d$group20
 ds<-ds[d$domain==1 & (if(mod=='A')d$ccA==1 else d$ccB==1),];ds$variables<-droplevels(ds$variables)
 f<-svycoxph(if(mod=='B')fb else fa,ds,ties='breslow',x=TRUE,model=TRUE)
 rows[[length(rows)+1]]<-extract(f,t,mod,ds$variables);fits[[paste(t,mod,sep='_')]]<-f
}
for(mod in c('B_plus_UCr','components','NoPIR_sameB','NoPIR_expanded')){
 ds<-des[d$domain==1 & (if(mod=='NoPIR_expanded')d$ccNoPIR else d$ccB==1),];ds$variables<-droplevels(ds$variables)
 form<-switch(mod,B_plus_UCr=update(fb,.~.+log2ucr),components=update(fb,.~.-group+log2alb+log2ucr),update(fb,.~.-pir))
 f<-svycoxph(form,ds,ties='breslow',x=TRUE,model=TRUE);rows[[length(rows)+1]]<-extract(f,30,mod,ds$variables);fits[[mod]]<-f
}
write.csv(do.call(rbind,rows),file.path(S,'audit/cox_results.csv'),row.names=FALSE)
saveRDS(fits,file.path(S,'data/cox_models.rds'));saveRDS(des,file.path(S,'data/survey_design.rds'))
capture.output(sessionInfo(),file=file.path(S,'audit/R_session.txt'))
cat('Cox thresholds/composition/noPIR complete\n')
