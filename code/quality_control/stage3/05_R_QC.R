# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Stage3_finalization");fits<-readRDS(file.path(S,'models/supplemental_Cox.rds'));d<-readRDS(.release_path("ODKD_WORK_ROOT", "Stage1/data/survey_design.rds"))$variables
d<-d[which(d$domain==1&d$ccB==1&d$time>24),];d$landmark<-24;d<-droplevels(d)
f<-Surv(landmark,time,death)~group+ns(age10,df=3)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
c<-coxph(f,data=d,weights=weight,ties='breslow');err<-max(abs(coef(c)-coef(fits$landmark24)));stopifnot(err<1e-8)
it<-fits$interaction;b<-coef(it);ix<-grep('group.*:sex|sex.*:group',names(b));stat<-drop(t(b[ix])%*%solve(vcov(it)[ix,ix],b[ix]))/length(ix);jt<-regTermTest(it,~group:sex);stopifnot(abs(stat-jt$Ftest)<1e-10)
stopifnot(all(vapply(fits,function(x)all(is.finite(coef(x)))&&all(is.finite(vcov(x))),TRUE)))
write.csv(data.frame(check=c('landmark delayed-entry vs reset time coefficient equality','joint3df interaction independent Wald','all supplemental coefficients and variances finite'),difference=c(err,abs(stat-jt$Ftest),0),result='PASS'),file.path(S,'audit/independent_R_QC.csv'),row.names=FALSE)
cat('R QC PASS\n')
