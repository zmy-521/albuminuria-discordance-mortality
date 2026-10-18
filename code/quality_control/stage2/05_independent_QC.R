# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines);S<-.release_path("ODKD_WORK_ROOT", "Stage2")
a<-readRDS(file.path(S,'data/fixed_risk_replicates.rds'));v<-coef(a);V<-vcov(a);stopifnot(all(is.finite(a$replicates)))
for(t in c(30,20))for(mode in c('known','all0','all1')){
 prefix<-paste(t,mode,sep='_');nm<-paste0(prefix,c('_risk_G0','_risk_G1','_risk_G2','_RD_G2_G0','_RD_G2_G1'));z<-a$replicates[,nm]
 stopifnot(max(abs(z[,4]-(z[,3]-z[,1])))<1e-12,max(abs(z[,5]-(z[,3]-z[,2])))<1e-12)
 stopifnot(abs(V[nm[4],nm[4]]-(V[nm[3],nm[3]]+V[nm[1],nm[1]]-2*V[nm[1],nm[3]]))<1e-12)
}
des<-readRDS(file.path(S,'data/survey_design.rds'));d<-droplevels(subset(des,domain==1&year<=2013&ccB==1)$variables);kn<-ns(d$age10,df=3)
for(j in 1:3)d[[paste0('ageS',j)]]<-kn[,j]
fits<-readRDS(file.path(S,'data/logistic_survey_models.rds'));rows<-list()
for(t in c(30,20))for(g in paste0('G',0:3)){
 nd<-d;nd$group<-factor(g,levels=paste0('G',0:3));pp<-as.numeric(predict(fits[[as.character(t)]],newdata=nd,type='response'));risk<-weighted.mean(pp,d$weight)
 ref<-v[paste(t,'known','risk',g,sep='_')];stopifnot(abs(risk-ref)<1e-10);rows[[length(rows)+1]]<-data.frame(check=paste('independent svyglm prediction',t,g),difference=abs(risk-ref),result='PASS')
}
for(file in c('cox20_risk_replicates.rds','selection_replicates.rds')){
 b<-readRDS(file.path(S,'data',file));stopifnot(all(is.finite(b$replicates)));rows[[length(rows)+1]]<-data.frame(check=file,difference=0,result='PASS')
}
rows[[length(rows)+1]]<-data.frame(check='All301 paired RD identities and covariance',difference=0,result='PASS')
write.csv(do.call(rbind,rows),file.path(S,'audit/independent_R_QC.csv'),row.names=FALSE);cat('R independent checks PASS\n')
