# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Stage1");options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
sink(file.path(S,'audit/R_session.txt'));sessionInfo();sink()
d<-read.csv(file.path(S,'data/survey_population.csv'),stringsAsFactors=FALSE,na.strings=c('','NA'))
for(v in c('group','group_raw','sex','race','cycle','smoking')) d[[v]]<-factor(d[[v]])
des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~weight,data=d,nest=TRUE)
fa<-Surv(time,death)~group+ns(age10,df=3)+sex+race+cycle
fb<-update(fa,.~.+bmi5+hba1c+egfr10+smoking+htn+cvd+pir)
rows<-list();phrows<-list();checks<-list();fits<-list()
extract<-function(f,scenario,model,design,n,events){
 b<-coef(f);V<-vcov(f);cons<-list(G1_G0=c(groupG1=1),G2_G0=c(groupG2=1),G3_G0=c(groupG3=1),G2_G1=c(groupG2=1,groupG1=-1))
 do.call(rbind,lapply(names(cons),function(nm){a<-setNames(rep(0,length(b)),names(b));a[names(cons[[nm]])]<-cons[[nm]];v<-sum(a*b);se<-sqrt(drop(t(a)%*%V%*%a));data.frame(scenario,model,design,contrast=nm,n,events,logHR=v,se,HR=exp(v),low=exp(v-1.96*se),high=exp(v+1.96*se),p=2*pnorm(-abs(v/se)))}))
}
for(sc in c('primary','post2007','baseline1999_2014','urine_unadjusted')){
 dd<-des
 if(sc=='urine_unadjusted')dd<-update(dd,group=group_raw)
 idx<-dd$variables$domain==1
 if(sc=='post2007')idx<-idx&dd$variables$year>=2007
 if(sc=='baseline1999_2014')idx<-idx&dd$variables$year<=2013
 for(mod in c('A','A_Bsample','B')){
  keep<-idx & (if(mod=='A')dd$variables$ccA==1 else dd$variables$ccB==1)
  ds<-dd[keep,];dat<-droplevels(ds$variables);ds$variables<-dat
  f<-if(mod=='B')fb else fa
  fit<-svycoxph(f,design=ds,ties='breslow',x=TRUE,model=TRUE)
  key<-paste(sc,mod,sep='_');fits[[key]]<-fit
  rows[[length(rows)+1]]<-extract(fit,sc,mod,'survey_Taylor',nrow(dat),sum(dat$death))
  # ordinary weighted fit only verifies the point estimates, not survey inference
  wfit<-coxph(f,data=dat,weights=weight,ties='breslow')
  checks[[key]]<-data.frame(check=key,max_coef_difference=max(abs(coef(fit)-coef(wfit))))
  uf<-coxph(f,data=dat,ties='breslow');rows[[length(rows)+1]]<-extract(uf,sc,mod,'unweighted',nrow(dat),sum(dat$death))
  if(sc=='primary'&&mod=='B'){
   zz<-cox.zph(fit,transform='km');write.csv(cbind(term=rownames(zz$table),as.data.frame(zz$table)),file.path(S,'audit/PH_schoenfeld_diagnostic_not_survey.csv'),row.names=FALSE)
  }
  cat(key,nrow(dat),sum(dat$death),'done\n');flush.console()
 }
}
write.csv(do.call(rbind,rows),file.path(S,'audit/cox_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(S,'audit/point_estimate_checks.csv'),row.names=FALSE)
saveRDS(fits,file.path(S,'data/cox_models.rds'))
# Predetermined time-varying group diagnostic, survey design covariance.
dd<-d[d$domain==1&d$ccB==1,];sp<-survSplit(Surv(time,death)~.,data=dd,cut=60,episode='period',start='start')
for(g in 1:3)sp[[paste0('g',g,'late')]]<-as.integer(sp$group==paste0('G',g)&sp$start>=60)
pdes<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~weight,data=sp,nest=TRUE)
pform<-update(fb,Surv(start,time,death)~.+g1late+g2late+g3late)
pfit<-svycoxph(pform,design=pdes,ties='breslow',x=TRUE,model=TRUE)
ptest<-regTermTest(pfit,~g1late+g2late+g3late)
capture.output(ptest,file=file.path(S,'audit/PH_design_test.txt'))
write.csv(data.frame(term=c('G1_after60','G2_after60','G3_after60'),logHR_change=coef(pfit)[c('g1late','g2late','g3late')],SE=sqrt(diag(vcov(pfit)))[c('g1late','g2late','g3late')]),file.path(S,'audit/PH_design_coefficients.csv'),row.names=FALSE)
saveRDS(des,file.path(S,'data/survey_design.rds'));saveRDS(pfit,file.path(S,'data/PH_model.rds'))
cat('Cox analyses complete\n')
