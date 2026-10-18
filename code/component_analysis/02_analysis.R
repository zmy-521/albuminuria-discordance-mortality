# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Review_supplement_urine_creatinine");options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE,warn=1)
des<-readRDS(.release_path("ODKD_WORK_ROOT", "Stage2/data/survey_design.rds"));d<-des$variables
old<-readRDS(.release_path("ODKD_WORK_ROOT", "Stage2/data/cox_models.rds"));stopifnot(!any(grepl('group',names(coef(old$components)))))
full<-subset(des,domain==1);ds<-subset(des,domain==1&ccB==1);stopifnot(nrow(full$variables)==5905,nrow(ds$variables)==5155,sum(ds$variables$death)==1115)
fb<-Surv(time,death)~group+ns(age10,df=3)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
fc<-update(fb,.~.-group+log2alb+log2ucr)
fit<-svycoxph(fc,design=ds,ties='breslow',model=TRUE,x=TRUE)
stopifnot(max(abs(coef(fit)-coef(old$components)))<1e-10,max(abs(vcov(fit)-vcov(old$components)))<1e-10)
write.csv(data.frame(term=names(coef(fit)),coefficient=coef(fit)),file.path(S,'audit/component_coefficients.csv'),row.names=FALSE)
write.csv(vcov(fit),file.path(S,'audit/component_design_covariance.csv'))
rows<-list();fits<-list(components_linear=fit);checks<-list()
contrast<-function(f,a,name,scenario,dd){b<-coef(f);V<-vcov(f);z<-setNames(rep(0,length(b)),names(b));z[names(a)]<-a;e<-sum(z*b);var<-drop(t(z)%*%V%*%z);se<-sqrt(var);df<-degf(dd)-length(b)+1;w<-e^2/var;data.frame(scenario,contrast=name,n=nrow(dd$variables),events=sum(dd$variables$death),logHR=e,variance=var,SE=se,HR=exp(e),low=exp(e-1.96*se),high=exp(e+1.96*se),Wald_chisq=w,p=2*pnorm(-abs(e/se)),design_df=df,F_p=pf(w,1,df,lower.tail=FALSE),t_low=exp(e-qt(.975,df)*se),t_high=exp(e+qt(.975,df)*se))}
for(sc in c('original_linear_eGFR','fixed_RCS3_eGFR')){
 dd<-des
 if(sc=='fixed_RCS3_eGFR'){
  k<-quantile(ds$variables$egfr10,c(.05,.35,.65,.95));xx<-dd$variables$egfr10;keep<-is.finite(xx);mm<-matrix(NA_real_,length(xx),3);mm[keep,]<-ns(xx[keep],knots=k[2:3],Boundary.knots=k[c(1,4)])
  for(j in 1:3)dd$variables[[paste0('egfrS',j)]]<-mm[,j]
  write.csv(data.frame(percentile=c(.05,.35,.65,.95),eGFR=k*10),file.path(S,'audit/fixed_eGFR_knots.csv'),row.names=FALSE)
 }
 dfull<-subset(dd,domain==1&ccB==1);form<-fc
 if(sc=='fixed_RCS3_eGFR')form<-update(form,.~.-egfr10+egfrS1+egfrS2+egfrS3)
 f<-if(sc=='original_linear_eGFR')fit else svycoxph(form,dfull,ties='breslow',model=TRUE,x=TRUE);fits[[paste0('component_',sc)]]<-f
 rows[[length(rows)+1]]<-contrast(f,c(log2alb=1),'Alb_per_doubling',sc,dfull)
 rows[[length(rows)+1]]<-contrast(f,c(log2ucr=1),'UCr_per_doubling_at_fixed_Alb',sc,dfull)
 rows[[length(rows)+1]]<-contrast(f,c(log2alb=1,log2ucr=1),'UCr_per_doubling_at_fixed_ACR',sc,dfull)
 dfull<-update(dfull,log2acr=log2(acr_adjusted));reform<-update(form,.~.-log2alb+log2acr)
 rf<-svycoxph(reform,dfull,ties='breslow',model=TRUE,x=TRUE);sumb<-coef(f)['log2alb']+coef(f)['log2ucr'];sumv<-sum(vcov(f)[c('log2alb','log2ucr'),c('log2alb','log2ucr')]);stopifnot(abs(sumb-coef(rf)['log2ucr'])<1e-9,abs(sumv-vcov(rf)['log2ucr','log2ucr'])<1e-9)
 jt<-regTermTest(rf,~log2ucr);capture.output(jt,file=file.path(S,paste0('audit/reparameterized_test_',sc,'.txt')))
 checks[[sc]]<-data.frame(check=sc,beta_sum=sumb,reparameterized=coef(rf)['log2ucr'],variance_sum=sumv,reparameterized_variance=vcov(rf)['log2ucr','log2ucr'])
 aa<-subset(dd,domain==1&ccB==1&acr_adjusted>=30&acr_adjusted<300);aa<-update(aa,g2=as.numeric(group=='G2'),log2acr=log2(acr_adjusted));af<-update(form,.~.-log2alb-log2ucr+g2+log2acr)
 afit<-svycoxph(af,aa,ties='breslow',x=TRUE,model=TRUE);fits[[paste0('A2_',sc)]]<-afit;rows[[length(rows)+1]]<-contrast(afit,c(g2=1),'G2_vs_G3_A2_adjust_log2ACR',sc,aa)
 if(sc=='fixed_RCS3_eGFR'){
  gf<-svycoxph(update(fb,.~.-egfr10+egfrS1+egfrS2+egfrS3),dfull,ties='breslow',model=TRUE,x=TRUE);fits$groupB_RCS<-gf;rows[[length(rows)+1]]<-contrast(gf,c(groupG2=1),'original_G2_G0_RCS_check',sc,dfull)
 }
}
write.csv(do.call(rbind,rows),file.path(S,'audit/effect_estimates.csv'),row.names=FALSE);write.csv(do.call(rbind,checks),file.path(S,'audit/reparameterization_QC.csv'),row.names=FALSE);saveRDS(fits,file.path(S,'models/supplemental_models.rds'))
basecheck<-data.frame(n=nrow(ds$variables),events=sum(ds$variables$death),strata=length(unique(ds$variables$SDMVSTRA)),PSU=nrow(unique(ds$variables[,c('SDMVSTRA','SDMVPSU')])),design_df=degf(ds),components_contains_group=FALSE,max_old_coef_difference=max(abs(coef(fit)-coef(old$components))),max_old_cov_difference=max(abs(vcov(fit)-vcov(old$components))))
write.csv(basecheck,file.path(S,'audit/frozen_model_audit.csv'),row.names=FALSE)
capture.output(list(formula=fc,variables=names(ds$variables),design=des$call,weights_formula=des$call$weights),file=file.path(S,'audit/model_definition.txt'))
# Descriptive survey functions: no new outcome models in flow subset.
describe<-function(dd,vars,label){out<-list();for(v in vars){xx<-dd$variables[[v]];ok<-is.finite(xx);if(!sum(ok))next;ss<-dd[which(ok),];form<-as.formula(paste0('~',v));m<-svymean(form,ss);ci<-confint(m);qq<-as.numeric(coef(svyquantile(form,ss,quantiles=c(.25,.5,.75),ci=FALSE)));out[[v]]<-data.frame(population=label,variable=v,n=sum(ok),missing=sum(!ok),mean=coef(m)[1],low=ci[1],high=ci[2],q1=qq[1],median=qq[2],q3=qq[3],min=min(xx[ok]),max=max(xx[ok]));};do.call(rbind,out)}
adr<-list();cnt<-list();bins<-list();support<-list()
for(pop in c('A2_all','A2_B')){
 aa<-subset(des,domain==1&acr_adjusted>=30&acr_adjusted<300);if(pop=='A2_B')aa<-subset(aa,ccB==1)
 for(g in c('G2','G3')){
  q<-subset(aa,group==g);q<-update(q,female=as.numeric(sex=='2'));cnt[[paste(pop,g)]]<-data.frame(population=pop,group=g,n=nrow(q$variables),events=sum(q$variables$death),PSU=nrow(unique(q$variables[,c('SDMVSTRA','SDMVPSU')])),strata=length(unique(q$variables$SDMVSTRA)),design_df=degf(q))
  adr[[paste(pop,g)]]<-describe(q,c('age','female','bmi','hba1c','egfr','ucr_adjusted','ualb_mgl','acr_adjusted','scr_mgdl','pir'),paste(pop,g))
  for(i in 1:4){lo<-c(30,60,100,200)[i];hi<-c(60,100,200,300)[i];bins[[length(bins)+1]]<-data.frame(population=pop,group=g,interval=paste0(lo,'-<',hi),n=sum(q$variables$acr_adjusted>=lo&q$variables$acr_adjusted<hi))}
 }
 ar<-aa$variables$acr_adjusted;gr<-aa$variables$group;lo<-max(min(ar[gr=='G2']),min(ar[gr=='G3']));hi<-min(max(ar[gr=='G2']),max(ar[gr=='G3']))
 for(g in c('G2','G3')){q<-subset(aa,group==g);q<-update(q,in_shared_range=as.numeric(acr_adjusted>=lo&acr_adjusted<=hi));support[[length(support)+1]]<-data.frame(population=pop,group=g,shared_min=lo,shared_max=hi,n=nrow(q$variables),outside=sum(q$variables$in_shared_range==0),weighted_inside=coef(svymean(~in_shared_range,q))[1])}
}
write.csv(do.call(rbind,cnt),file.path(S,'audit/A2_counts.csv'),row.names=FALSE);write.csv(do.call(rbind,adr),file.path(S,'audit/A2_weighted_descriptives.csv'),row.names=FALSE);write.csv(do.call(rbind,bins),file.path(S,'audit/A2_support_bins.csv'),row.names=FALSE);write.csv(do.call(rbind,support),file.path(S,'audit/A2_common_support.csv'),row.names=FALSE)
flow<-read.csv(file.path(S,'data/flow_to_R.csv'),na.strings=c('','NA'));ix<-match(des$variables$SEQN,flow$SEQN)
for(v in setdiff(names(flow),'SEQN'))des$variables[[v]]<-flow[[v]][ix]
fr<-subset(des,domain==1&year>=2009&flow_valid==1);dr<-list();cn<-list();gp<-list()
for(g in paste0('G',0:3)){
 q<-subset(fr,group==g);cn[[g]]<-data.frame(group=g,n=nrow(q$variables),events=sum(q$variables$death),heart_failure_missing=sum(is.na(q$variables$heart_failure)),diuretic_missing=sum(is.na(q$variables$diuretic)),BMI_missing=sum(is.na(q$variables$bmi)))
 dr[[g]]<-describe(q,c('flow','AER','CER','ucr_adjusted','ualb_mgl','acr_adjusted','age','bmi','female','scr_mgdl','heart_failure','diuretic','total_volume','total_time'),g)
}
q<-subset(fr,group=='G2');r<-svyciprop(~AER20,q,method='logit');ci<-confint(r)
write.csv(data.frame(group='G2',n=nrow(q$variables),count=sum(q$variables$AER20),unweighted=mean(q$variables$AER20),weighted=coef(r)[1],low=ci[1],high=ci[2],design_df=degf(q)),file.path(S,'audit/G2_AER20.csv'),row.names=FALSE)
write.csv(do.call(rbind,cn),file.path(S,'audit/flow_counts.csv'),row.names=FALSE);write.csv(do.call(rbind,dr),file.path(S,'audit/flow_weighted_descriptives.csv'),row.names=FALSE)
# Selection description, same original2009–2018 domain; not selection-corrected inference.
sel<-list();for(z in c(0,1)){q<-subset(des,domain==1&year>=2009&flow_valid==z);sel[[as.character(z)]]<-describe(q,c('age','female','bmi','ucr_adjusted','ualb_mgl','egfr','death'),paste0('flow_valid=',z))}
write.csv(do.call(rbind,sel),file.path(S,'audit/flow_selection_description.csv'),row.names=FALSE)
saveRDS(des,file.path(S,'data/review_survey_design.rds'));capture.output(sessionInfo(),file=file.path(S,'audit/R_session.txt'))
cat('Finite review analyses complete\n')
