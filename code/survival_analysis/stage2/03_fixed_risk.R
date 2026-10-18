# Release adaptation: use the installed R libraries; local runtime paths removed.
.release_path <- function(env, relative="") { root <- Sys.getenv(env); if (!nzchar(root)) stop(paste("Set",env)); file.path(root,relative) }
library(survey);library(survival);library(splines)
S<-.release_path("ODKD_WORK_ROOT", "Stage2");options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
des<-readRDS(file.path(S,'data/survey_design.rds'));rd<-as.svrepdesign(des,type='JKn',mse=TRUE)
rd<-subset(rd,domain==1 & year<=2013 & ccB==1);dat<-droplevels(rd$variables)
dat$y60<-ifelse(dat$death==0&dat$months<60,NA,as.integer(dat$death==1&dat$months<=60))
kn<-ns(dat$age10,df=3);kb<-attr(kn,'knots');bb<-attr(kn,'Boundary.knots')
for(j in 1:3)dat[[paste0('ageS',j)]]<-kn[,j]
rd$variables<-dat
f<-y60~group+ageS1+ageS2+ageS3+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
xs<-list();xgs<-list();checks<-list();svyfits<-list()
for(t in c(30,20)){
 dd<-dat;if(t==20)dd$group<-dd$group20
 mf<-model.frame(delete.response(terms(f)),data=dd,na.action=na.pass);X<-model.matrix(delete.response(terms(f)),mf)
 xs[[as.character(t)]]<-X
 xgs[[as.character(t)]]<-lapply(paste0('G',0:3),function(g){ng<-dd;ng$group<-factor(g,levels=paste0('G',0:3));model.matrix(delete.response(terms(f)),ng)})
 sd<-subset(des,domain==1&year<=2013&ccB==1);sd$variables<-dd
 sf<-svyglm(f,design=subset(sd,!is.na(y60)),family=quasibinomial());svyfits[[as.character(t)]]<-sf
 known<-!is.na(dd$y60);gf<-glm.fit(X[known,],dd$y60[known],weights=dd$weight[known]/mean(dd$weight[known]),family=quasibinomial(),control=glm.control(maxit=50,epsilon=1e-9))
 checks[[as.character(t)]]<-data.frame(threshold=t,max_coef_difference=max(abs(coef(sf)-gf$coefficients)))
}
saveRDS(svyfits,file.path(S,'data/logistic_survey_models.rds'));write.csv(do.call(rbind,checks),file.path(S,'audit/logistic_point_checks.csv'),row.names=FALSE)
counter<-0
calc<-function(w,X,Xg,y){
 use<-w>0&!is.na(y);fit<-glm.fit(X[use,],y[use],weights=w[use]/mean(w[use]),family=quasibinomial(),control=glm.control(maxit=50,epsilon=1e-9))
 stopifnot(fit$converged,all(is.finite(fit$coefficients)))
 r<-sapply(Xg,function(x)sum(w*plogis(drop(x%*%fit$coefficients)))/sum(w))
 setNames(c(r,r[3]-r[1],r[3]-r[2]),c(paste0('risk_G',0:3),'RD_G2_G0','RD_G2_G1'))
}
theta<-function(w,data){
 counter<<-counter+1;if(counter%%40==0){cat('risk replicate',counter,'\n');flush.console()}
 out<-c()
 for(t in c(30,20))for(mode in c('known','all0','all1')){
  y<-data$y60;if(mode!='known')y[is.na(y)]<-if(mode=='all0')0 else 1
  r<-calc(w,xs[[as.character(t)]],xgs[[as.character(t)]],y);names(r)<-paste(t,mode,names(r),sep='_');out<-c(out,r)
 }
 # Paired threshold differences, not a claim of causal threshold effect.
 c(out,threshold20_minus30_RD20=unname(out['20_known_RD_G2_G0']-out['30_known_RD_G2_G0']),threshold20_minus30_RD21=unname(out['20_known_RD_G2_G1']-out['30_known_RD_G2_G1']))
}
ans<-withReplicates(rd,theta,return.replicates=TRUE);saveRDS(ans,file.path(S,'data/fixed_risk_replicates.rds'))
v<-coef(ans);se<-SE(ans);lo<-v-1.96*se;hi<-v+1.96*se;pr<-grepl('_risk_',names(v));lo[pr]<-plogis(qlogis(v[pr])-1.96*se[pr]/(v[pr]*(1-v[pr])));hi[pr]<-plogis(qlogis(v[pr])+1.96*se[pr]/(v[pr]*(1-v[pr])))
write.csv(data.frame(estimand=names(v),estimate=as.numeric(v),se=as.numeric(se),low=as.numeric(lo),high=as.numeric(hi),inference='JKn301'),file.path(S,'audit/fixed_risk.csv'),row.names=FALSE)
# All unknown outcome assignments: four in B target; fifth outside B cannot change these fits.
unknown<-which(is.na(dat$y60));comb<-expand.grid(rep(list(0:1),length(unknown)));enum<-list()
for(t in c(30,20))for(i in 1:nrow(comb)){
 y<-dat$y60;y[unknown]<-as.numeric(comb[i,]);r<-calc(dat$weight,xs[[as.character(t)]],xgs[[as.character(t)]],y)
 enum[[length(enum)+1]]<-data.frame(threshold=t,assignment=i,unknown_events=sum(comb[i,]),estimand=names(r),estimate=as.numeric(r))
}
write.csv(do.call(rbind,enum),file.path(S,'audit/censoring_enumeration.csv'),row.names=FALSE)
# Stage1-compatible Cox standardized risk at20, explicitly distinct from PH-free logistic.
fc<-Surv(time,death)~group+ns(age10,knots=kb,Boundary.knots=bb)+sex+race+cycle+bmi5+hba1c+egfr10+smoking+htn+cvd+pir
ct<-function(w,data){dd<-data[w>0,];dd$.w<-w[w>0];dd$group<-dd$group20;fit<-coxph(fc,data=dd,weights=.w,ties='breslow',x=TRUE,model=TRUE);bh<-basehaz(fit,centered=FALSE);h<-max(bh$hazard[bh$time<=60]);r<-sapply(paste0('G',0:3),function(g){nd<-dd;nd$group<-factor(g,levels=paste0('G',0:3));sum(dd$.w*(-expm1(-h*exp(predict(fit,newdata=nd,type='lp',reference='zero')))))/sum(dd$.w)});setNames(c(r,r[3]-r[1],r[3]-r[2]),c(paste0('risk_G',0:3),'RD_G2_G0','RD_G2_G1'))}
cat('Logistic complete; Cox20 risk starting\n');flush.console()
ca<-withReplicates(rd,ct,return.replicates=TRUE);saveRDS(ca,file.path(S,'data/cox20_risk_replicates.rds'));v<-coef(ca);se<-SE(ca);lo<-v-1.96*se;hi<-v+1.96*se;pr<-1:4;lo[pr]<-plogis(qlogis(v[pr])-1.96*se[pr]/(v[pr]*(1-v[pr])));hi[pr]<-plogis(qlogis(v[pr])+1.96*se[pr]/(v[pr]*(1-v[pr])))
write.csv(data.frame(estimand=names(v),estimate=as.numeric(v),se=as.numeric(se),low=as.numeric(lo),high=as.numeric(hi)),file.path(S,'audit/cox20_risk.csv'),row.names=FALSE)
cat('Fixed risk complete\n')
