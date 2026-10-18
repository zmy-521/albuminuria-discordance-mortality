# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines);S<-.release_path("ODKD_WORK_ROOT", "Stage2")
des<-readRDS(file.path(S,'data/survey_design.rds'));dom<-subset(des,domain==1);d<-droplevels(dom$variables);dom$variables<-d
kn<-ns(d$age10,df=3);ks<-attr(kn,'knots');bs<-attr(kn,'Boundary.knots')
for(j in 1:3)d[[paste0('ageS',j)]]<-kn[,j]
dom$variables<-d
fs<-ccB~group+ageS1+ageS2+ageS3+sex+race+cycle
sf<-svyglm(fs,dom,family=quasibinomial());prob<-as.numeric(predict(sf,type='response'))
stopifnot(min(prob)>.05)
cc<-d$ccB==1;kk<-ns(d$age10[cc],df=3);kb<-attr(kk,'knots');bb<-attr(kk,'Boundary.knots')
fc<-Surv(time,death)~group+ns(age10,knots=kb,Boundary.knots=bb)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
diag<-data.frame(n=nrow(d),completeB=sum(cc),min_prob=min(prob),max_prob=max(prob),max_inverse_weight=max(1/prob),ESS_before=sum(d$weight[cc])^2/sum(d$weight[cc]^2),ESS_after=sum(d$weight[cc]/prob[cc])^2/sum((d$weight[cc]/prob[cc])^2))
write.csv(diag,file.path(S,'audit/selection_weights.csv'),row.names=FALSE)
rd<-as.svrepdesign(des,type='JKn',mse=TRUE);rd<-subset(rd,domain==1);rd$variables<-d
X<-model.matrix(sf);counter<-0
theta<-function(w,data){
 counter<<-counter+1;if(counter%%50==0){cat('selection replicate',counter,'\n');flush.console()}
 use<-w>0;g<-glm.fit(X[use,],data$ccB[use],weights=w[use]/mean(w[use]),family=quasibinomial(),control=glm.control(maxit=50,epsilon=1e-9));stopifnot(g$converged)
 pp<-plogis(drop(X%*%g$coefficients));stopifnot(min(pp)>.05)
 dd<-data[cc&use,];dd$.w<-w[cc&use]/pp[cc&use]
 cf<-coxph(fc,data=dd,weights=.w,ties='breslow',x=TRUE,model=TRUE)
 b<-coef(cf);c(G2_G0=unname(b['groupG2']),G1_G0=unname(b['groupG1']),G3_G0=unname(b['groupG3']),G2_G1=unname(b['groupG2']-b['groupG1']))
}
a<-withReplicates(rd,theta,return.replicates=TRUE);v<-coef(a);se<-SE(a)
write.csv(data.frame(contrast=names(v),n=sum(cc),events=sum(d$death[cc]),logHR=as.numeric(v),se=as.numeric(se),HR=exp(v),low=exp(v-1.96*se),high=exp(v+1.96*se),p=2*pnorm(-abs(v/se))),file.path(S,'audit/selection_cox.csv'),row.names=FALSE)
saveRDS(a,file.path(S,'data/selection_replicates.rds'));saveRDS(sf,file.path(S,'data/selection_probability_model.rds'));cat('Selection completed\n')
