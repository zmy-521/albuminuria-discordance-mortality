# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=FALSE,digits=17)
R<-.release_path("ODKD_WORK_ROOT", "");O<-file.path(R,'Table1_flow_completion');dir.create(file.path(O,'internal'),recursive=TRUE,showWarnings=FALSE)
des<-readRDS(file.path(R,'Stage1/data/survey_design.rds'));z<-subset(des,domain==1);stopifnot(nrow(z$variables)==5905)
continuous<-c(age='years',bmi='kg/m2',hba1c='%',scr_mgdl='mg/dL',egfr='mL/min/1.73m2',pir='ratio',ualb_mgl='mg/L',ucr_adjusted='mg/dL',acr_adjusted='mg/g')
categorical<-list(sex=c('1','2'),race=as.character(1:5),smoking=c('never','former','current'),htn=c('0','1'),cvd=c('0','1'))
res<-list()
for(g in c('Overall',paste0('G',0:3))){
 dd<-if(g=='Overall')z else subset(z,group==g);N<-nrow(dd$variables)
 for(v in names(continuous)){
  ok<-is.finite(dd$variables[[v]]);ss<-dd[ok,];fm<-as.formula(paste0('~',v));m<-svymean(fm,ss);qq<-rep(NA_real_,3)
  if(v%in%c('ualb_mgl','ucr_adjusted','acr_adjusted'))qq<-as.numeric(coef(svyquantile(fm,ss,quantiles=c(.25,.5,.75),qrule='math',ci=FALSE)))
  res[[length(res)+1]]<-data.frame(group=g,variable=v,category='',unit=continuous[v],group_N=N,valid_n=sum(ok),missing_n=sum(!ok),category_n=NA_integer_,mean=as.numeric(coef(m)),mean_SE=as.numeric(SE(m)),percent=NA_real_,percent_SE=NA_real_,Q1=qq[1],median=qq[2],Q3=qq[3])
 }
 for(v in names(categorical)){
  raw<-as.character(dd$variables[[v]]);ok<-!is.na(raw);stopifnot(all(raw[ok]%in%categorical[[v]]));ss<-dd[ok,]
  for(level in categorical[[v]]){
   ss$variables$.indicator<-as.numeric(as.character(ss$variables[[v]])==level);m<-svymean(~.indicator,ss)
   res[[length(res)+1]]<-data.frame(group=g,variable=v,category=level,unit='%',group_N=N,valid_n=sum(ok),missing_n=sum(!ok),category_n=sum(raw[ok]==level),mean=NA_real_,mean_SE=NA_real_,percent=100*as.numeric(coef(m)),percent_SE=100*as.numeric(SE(m)),Q1=NA_real_,median=NA_real_,Q3=NA_real_)
  }
 }
}
write.csv(do.call(rbind,res),file.path(O,'internal/baseline_full_precision.csv'),row.names=FALSE)
# Actual saved model members; no fitting, prediction or complete-case substitution.
co<-readRDS(file.path(R,'Stage1/data/cox_models.rds'));mb<-co$primary_B$survey.design$variables;stopifnot(identical(rownames(co$primary_B$model),rownames(mb)))
su<-readRDS(file.path(R,'Review_supplement_urine_creatinine/models/supplemental_models.rds'));aa<-su$A2_original_linear_eGFR$survey.design$variables
lo<-readRDS(file.path(R,'Stage2/data/logistic_survey_models.rds'))[['30']];lv<-lo$survey.design$variables;stopifnot(identical(rownames(lo$model),rownames(lv)))
rv<-readRDS(file.path(R,'Review_supplement_urine_creatinine/data/review_survey_design.rds'))$variables
sets<-list(main=z$variables$SEQN,ModelB=mb$SEQN,A2_actual=aa$SEQN,five_year_source=z$variables$SEQN[z$variables$year<=2013],five_year_target=co$baseline1999_2014_B$survey.design$variables$SEQN,five_year_logistic_fit=lv$SEQN,flow_source=rv$SEQN[rv$domain==1&rv$year>=2009],flow_valid=rv$SEQN[rv$domain==1&rv$year>=2009&rv$flow_valid==1])
sets<-lapply(sets,function(x)x[!is.na(x)])
stopifnot(length(sets$ModelB)==5155,sum(mb$death)==1115,length(sets$A2_actual)==1120,length(sets$five_year_source)==4315,length(sets$five_year_target)==3789,length(sets$five_year_logistic_fit)==3785,length(sets$flow_source)==3539,length(sets$flow_valid)==3261)
std<-readRDS(file.path(R,'Stage1/data/risk60_standardization_model.rds'));target<-co$baseline1999_2014_B$survey.design$variables
stopifnot(nrow(std$model)==3789,setequal(rownames(std$model),rownames(target)))
stopifnot(all.equal(unname(as.matrix(std$y)),unname(as.matrix(target[match(rownames(std$model),rownames(target)),c('time','death')])),check.attributes=FALSE)==TRUE)
write.csv(data.frame(SEQN=z$variables$SEQN,lapply(sets,function(s)as.integer(z$variables$SEQN%in%s))),file.path(O,'internal/member_flags_LOCAL_ONLY.csv'),row.names=FALSE)
relations<-list();for(a in names(sets))for(b in names(sets))if(a<b)relations[[length(relations)+1]]<-data.frame(set_A=a,set_B=b,n_A=length(sets[[a]]),n_B=length(sets[[b]]),intersection=length(intersect(sets[[a]],sets[[b]])),A_not_B=length(setdiff(sets[[a]],sets[[b]])),B_not_A=length(setdiff(sets[[b]],sets[[a]])))
write.csv(do.call(rbind,relations),file.path(O,'internal/member_relations.csv'),row.names=FALSE)
capture.output(list(design_class=class(des),design_call=des$call,full_design_n=nrow(des$variables),main_n=nrow(z$variables),strata=length(unique(z$variables$SDMVSTRA)),design_df=degf(z),quantile_rule='survey::svyquantile qrule=math; inverse weighted empirical CDF; CI not requested',Cox_standardization_members_verified=TRUE),file=file.path(O,'internal/design_QC.txt'))
capture.output(sessionInfo(),file=file.path(O,'internal/R_session.txt'))
capture.output(warnings(),file=file.path(O,'internal/R_warnings.txt'))
cat('Baseline descriptives and actual membership checks completed; no model fit.\n')
