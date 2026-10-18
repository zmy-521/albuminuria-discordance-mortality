# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
root<-.release_path("ODKD_WORK_ROOT", "");out<-file.path(root,'Figure_v3_PlotData_export/internal');dir.create(out,recursive=TRUE,showWarnings=FALSE)
models<-readRDS(file.path(root,'Stage1/data/cox_models.rds'));f<-models$primary_B
capture.output(list(names=names(f),model_n=nrow(f$model),model_columns=names(f$model),survey_columns=names(f$survey.design$variables),model_nobs=f$n,events=f$nevent),file=file.path(out,'model_structure.txt'))
v<-f$survey.design$variables
stopifnot(!is.null(v$SEQN),!anyDuplicated(v$SEQN),nrow(v)==5155,sum(v$death)==1115)
stopifnot(identical(rownames(f$model),rownames(v)))
write.csv(v[,c('SEQN','cycle','group','weight','death','SDMVPSU','SDMVSTRA')],file.path(out,'actual_ModelB_members.csv'),row.names=FALSE)
review<-readRDS(file.path(root,'Review_supplement_urine_creatinine/models/supplemental_models.rds'))
a<-review$A2_original_linear_eGFR;av<-a$survey.design$variables
stopifnot(identical(rownames(a$model),rownames(av)),nrow(av)==1120)
write.csv(av[,c('SEQN','group','acr_adjusted','death')],file.path(out,'actual_A2_members.csv'),row.names=FALSE)
cf<-review$components_linear
cc<-read.csv(file.path(root,'Review_supplement_urine_creatinine/audit/component_coefficients.csv'),check.names=FALSE)
cv<-as.matrix(read.csv(file.path(root,'Review_supplement_urine_creatinine/audit/component_design_covariance.csv'),row.names=1,check.names=FALSE))
stopifnot(identical(cc$term,names(coef(cf))),!any(grepl('group',cc$term)),max(abs(cc$coefficient-coef(cf)))<1e-12,max(abs(cv-vcov(cf)))<1e-12)
capture.output(list(component_formula=formula(cf),n=cf$n,events=cf$nevent,max_coefficient_export_difference=max(abs(cc$coefficient-coef(cf))),max_covariance_export_difference=max(abs(cv-vcov(cf))),contains_group=any(grepl('group',cc$term))),file=file.path(out,'component_QC.txt'))
files<-unlist(lapply(c('Stage1','Stage2','Stage3_finalization','Review_supplement_urine_creatinine'),function(stage)unlist(lapply(c('data','models'),function(folder)list.files(file.path(root,stage,folder),pattern='\\.rds$',recursive=FALSE,full.names=TRUE)))))
inventory<-list();curves<-list()
visit<-function(x,p,path,depth=0){
 if(inherits(x,c('survfit','svykm','svykmlist')))curves[[length(curves)+1]]<<-data.frame(file=p,object=path,class=paste(class(x),collapse=';'))
 if(depth<2 && is.list(x) && !inherits(x,c('coxph','survey.design','survey.design2','svyglm','svrep.design'))){for(n in names(x))visit(x[[n]],p,paste0(path,'$',n),depth+1)}
}
for(p in files){x<-readRDS(p);inventory[[p]]<-data.frame(file=p,class=paste(class(x),collapse=';'),members=paste(head(names(x),60),collapse=';'));visit(x,p,'root')}
write.csv(do.call(rbind,inventory),file.path(out,'RDS_inventory.csv'),row.names=FALSE)
writeLines(if(length(curves))capture.output(do.call(rbind,curves)) else 'No saved survfit/svykm/svykmlist objects found in scoped frozen RDS files.',file.path(out,'curve_object_search.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
cat('Saved-model membership verified: Model B 5155/1115; A2',nrow(av),'; curves',length(curves),'\n')
