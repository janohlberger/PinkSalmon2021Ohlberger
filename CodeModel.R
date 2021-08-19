##=======================================================================##
##                                                                       ##
## Wild PWS pink salmon linear Ricker stock-recruit model and covariates ##
##                                                                       ##
##=======================================================================##
pkgs<-c("readxl","tidyverse","ggplot2","viridis","MuMIn","effects","Hmisc", "visreg","RColorBrewer","car","fmsb","mgcv","ggplot2","fabricatr", "performance","R.utils","ggiraphExtra","writexl","asbio","jtools","lattice","grid","gridExtra","nlme","pedometrics","gsl")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) { install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=T) }
invisible(lapply(pkgs,library,character.only=T)) ## install and load 
if(exists("mainDir")) { } else { mainDir<-setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) } ## set home/main directory
is.even<-function(x) { x %% 2 == 0 }
'%!in%'<-function(x,y)!('%in%'(x,y))

##=================================================================## data
data<-read.csv("DataModel.csv")
##-----------------------------------------------------## rename for model
data$esc<-data$Spawners ## brood year (two years prior)
data$runs<-data$prevYearRunSize ## previous year run size
data$totr<-data$prev_tot_PWS_return ## previous year total return to PWS
data$prel<-data$PWS_released_prev_Total ## previous year hatchery releases
data$srel<-data$PWS_released_Total ## same year hatchery releases
data$comp<-data$sameYrcomp_regional ## smae year regional competitors
data$compNP<-data$sameYrcompetitors ## same year North Pacific competitors
data$pdo<-data$PDOann ## annual PDO
data$npgo<-data$NPGOann ## annual NPGO
data$temp<-data$egoa.spr.sst ## eatern GoA spring SST
##---------------------------------------------## only selected covariates
data<-dplyr::select(data,Year,BL,lnRS,regime,esc,temp,totr,comp,srel)
data<-data[complete.cases(data),] ## drop years with NAs
nY<-dim(data)[1]
rownames(data)<-seq(1,dim(data)[1],1)
##------------------------------------------------------------## model data
# write.csv(data,"DataModel.csv")
##-----------------------------------------------------------## for scaling
data_unscaled<-dplyr::select(data,-Year,-BL,-lnRS,-regime) 
data_other<-dplyr::select(data,Year,BL,lnRS,regime)
data_means<-sapply(data_unscaled,mean)
data_sds<-sapply(data_unscaled,sd)
##-------------------------------------## center/standardize to mean=0/sd=1
data<-data.frame(cbind(data_other,scale(data_unscaled)))
data_unscaled_all<-data.frame(cbind(data_other,data_unscaled))
##----------------------------------------------------------------## notes
## pCO2 time series shorter >> data and model selection need to be re-run
## multiple temperature or competitor metrics not to be included together

##=======================================================================##
##======================================## multiple linear regression model 
##=======================================================================##
options(na.action="na.fail")	## avoid fits to different datasets if NAs
rank<-"AICc" ## model selection criterion (AICc or BIC)
##-------------------------------------------------------## model selection
# mod_form<-formula(lnRS~esc*BL+esc:BL:regime+temp+I(temp^2)+temp:regime +totr+I(totr^2)+srel+I(srel^2)+comp+I(comp^2)+pdo+npgo)
mod_form<-formula(lnRS~esc*BL+esc:BL:regime+temp+I(temp^2)+temp:regime+totr+I(totr^2)+srel+I(srel^2)+comp+I(comp^2)) 
fixed<-c("esc","BL","BL:esc") ## stock-recruit relationship each broodline
mod_full<-lm(mod_form,data=data) 
mod_select<-dredge(mod_full,fixed=fixed,trace=F,rank=rank) 
# mod_select<-dredge(mod_full,trace=F,rank=rank) ## for table
##--------------------------------------------------------## selected model
mod<-get.models(mod_select,subset=1)[[1]]
out<-summary(mod) ## adjusted r2: 0.53 (full dataset)
out$r.squared ## r2: ~62%
AICc(mod) ## ~82.4 (full dataset)
anova(mod)
car::vif(mod)
##----------------------------------------------## residual autocorrelation
## check for autocorrelation (threshold=0.25, lags 1+)
acf_out<-as.numeric(acf(residuals(mod),lag=9,plot=F)$acf)[-1] # no 0-lag
acf_out>=0.25
pacf_out<-as.numeric(pacf(residuals(mod),lag=9,plot=F)$acf)
pacf_out>=0.25
##------------------------------------------------------## effects per unit
mycoeffs<-data.frame(out$coefficients[-1,1])
spawners_per_million_even_since1989<-mycoeffs[6,]/data_sds[names(data_sds)=="esc"]
spawners_per_million_even_upto1988<-(mycoeffs[6,]+mycoeffs[9,])/data_sds[names(data_sds)=="esc"]
spawners_per_million_odd_since1989<-(mycoeffs[6,]+mycoeffs[8,])/ data_sds[names(data_sds)=="esc"]
spawners_per_million_odd_upto1989<-(mycoeffs[6,]+mycoeffs[8,]+mycoeffs[10,])/ data_sds[names(data_sds)=="esc"]
temperature_per_degC_since1989<-mycoeffs[3,]/data_sds[names(data_sds)=="temp"]
temperature_per_degC_upto1988<-(mycoeffs[3,]+mycoeffs[7,])/data_sds[names(data_sds)=="temp"]
releases_per_million<-mycoeffs[2,]/data_sds[names(data_sds)=="srel"]
total_run_per_million<-mycoeffs[4,]/data_sds[names(data_sds)=="totr"]
competitors_per_million<-mycoeffs[1,]/data_sds[names(data_sds)=="comp"]

##=========================================================## save results 
data$resid<-residuals(mod) 
# save(mod,file="Model.rda")
##----------------------------------------------------------## model terms
terms<-attr(mod$terms,"term.labels")
addterms<-paste(terms,collapse="+")       
##------------------------------------------------## selected model formula
mod_final<-sel_mod<-formula(paste("lnRS~",addterms,sep=""))
##---------------------------------------------------## save anova results
anova_results<-data.frame(anova(mod))
anova_results$Variable<-rownames(anova_results)
write_xlsx(anova_results,"model_ANOVA_table.xlsx")

##================================================## model selection table
num_of_mods<-20 ## top X models
mod_select_topX<-mod_select[1:num_of_mods,]
mod_select_topX$cum_weight<-cumsum(mod_select_topX$weight)
##---------------------------------------------## model forms
mod_forms<-get.models(mod_select,subset=delta<5)
mod_forms<-mod_forms[1:num_of_mods]
mod_forms_topX<-NA
for(i in 1:num_of_mods) {
  use_mod<-mod_forms[[i]]	
  form<-as.character(use_mod$call)[2]
  mod_forms_topX[i]<-form
}
##---------------------------------------------## edit vars in model forms
mod_forms_topX<-unlist(mod_forms_topX)
mod_forms_topX<-gsub("lnRS","ln(R/S)",mod_forms_topX)
mod_forms_topX<-gsub("esc","S",mod_forms_topX)
mod_forms_topX<-gsub("prel","H",mod_forms_topX)
mod_forms_topX<-gsub("srel","H",mod_forms_topX)
mod_forms_topX<-gsub("comp","C",mod_forms_topX)
mod_forms_topX<-gsub("temp","T",mod_forms_topX)
mod_forms_topX<-gsub("runs","R",mod_forms_topX)
mod_forms_topX<-gsub("totr","R",mod_forms_topX)
mod_forms_topX<-gsub("regime","D",mod_forms_topX)
mod_forms_topX<-gsub("BL","B",mod_forms_topX)
mod_forms_topX<-gsub(" ","",mod_forms_topX)
mod_covars<-mod_forms_topX
# mod_covars<-lapply(mod_covars,function(x) substr(x,1,nchar(x)-2))
mod_covars<-lapply(mod_covars,function(x) substr(x,1,nchar(x)-10))
mod_covars<-lapply(mod_covars,function(x) substr(x,9,nchar(x)))
mod_covars<-unlist(mod_covars)
##-----------------------------------------------------## add back to table
mod_select_topX$mod_covars<-mod_covars
mod_select_topX$mod_forms<-mod_forms_topX
fn_round<-function(x) { x<-round(x,digits=3) } ## round numeric values
mod_select_topX<-mod_select_topX %>% mutate_if(is.numeric,fn_round)
mod_select_topX$delta<-round(mod_select_topX$delta,1)
mod_selection_table<-mod_select_topX
write_xlsx(mod_selection_table,"model_selection_table.xlsx")

##=======================================================================##
##==================================================================## plot
##=======================================================================##
alphap<-0.05 ## 90% CIs (alpha: 1-coverage)
pdf("PWS_wild_pinks_lnRS_vs_covars_selected.pdf",width=9,height=6.2)
layout(matrix(c(1:6),nrow=2,byrow=T))
par(mar=c(4,4,1,1),oma=c(0,0,0.5,0),mgp=c(2.2,0.5,0),cex.axis=1,cex.lab=1.2,tcl=-0.3)
col<-"gray25";col_pt<-"gray50" ## no interaction plots
cols_R<-c("orangered","dodgerblue3") ## regimes
cols_BL<-c("goldenrod","forestgreen") ## brood lines
ltys<-c(1,3)
cexll<-1.5 ## letter size
ll<-c(-0.3,-0.1) ## inset for letters
leg<-c(0.1,0) ## inset for legend
cex.pt<-1.2 ## point size
by<-"regime"
ylab<-"Partial effect on ln(recruits/spawner)"
min_y<--1.1;max_y<-2.6;ylim<-c(min_y,max_y)
##====================================================## spawner abundance
xtrans<-function(x) { x*data_sds[names(data_sds)=="esc"]+data_means[names(data_means)=="esc"] }
##----------------------------------------------------## even
obj_esc_even<-visreg(mod,xvar="esc",by=by,cond=list(BL="even"),partial=T,alpha=alphap,scale="response",plot=F)
obj_esc_even_unscaled<-visreg(mod,xvar="esc",by=by,cond=list(BL="even"), partial=T,alpha=alphap,scale="response",plot=F,xtrans=xtrans)
obj_esc_even$res$BL<-NA ## create new BL column with actual BL by year
for(i in 1:dim(obj_esc_even$res)[1]) { obj_esc_even$res$BL[i]<-data$BL[data$esc==obj_esc_even$res$esc[i]] }
min_esc<-min(obj_esc_even$res$esc) ## minimum escapement for this BL
max_esc<-max(obj_esc_even$res$esc) ## maximum escapement for this BL
index<-which(obj_esc_even$res$BL=="even") ## get row index 
obj_esc_even$res<-obj_esc_even$res[index,] ## drop residuals of other BL
obj_esc_even$res$esc<-obj_esc_even_unscaled$res$esc[index] ## replace with unscaled escapement
obj_esc_even$fit<-obj_esc_even$fit[obj_esc_even$fit$esc>=min_esc & obj_esc_even$fit$esc<=max_esc,] ## use only range predicted for this BL
obj_esc_even$fit$esc<-obj_esc_even$fit$esc*data_sds[names(data_sds)=="esc"]+data_means[names(data_means)=="esc"]
obj_esc_even$fit<-obj_esc_even$fit[obj_esc_even$fit$esc>=0.95*min(data_unscaled_all$esc[data_unscaled_all$BL=="even"]) & obj_esc_even$fit$esc<=1.01*max(data_unscaled_all$esc[data_unscaled_all$BL=="even"]),]
min_x<-min(obj_esc_even$res$esc);max_x<-max(obj_esc_even$res$esc) 
plot(obj_esc_even,overlay=T,legend=F,partial=T,xlab="Spawner abundance (millions)",ylab=ylab,xlim=c(min_x,max_x),ylim=ylim,line=list(col=cols_R,lty=ltys),fill=list(col=alpha(cols_R,0.25)),points=list(pch=21,bg=cols_R,col=1,lwd=0.2,cex=cex.pt)) 
box()
##-----------------------------------------------------------## legend
legend("bottomleft",c("",""), col=alpha(rev(cols_R),0.25),cex=1,lwd=7,bty="n",seg.len=2,inset=leg)
labs<-c("up to 1988","since 1989")
legend("bottomleft",labs,col=rev(cols_R),lty=rev(ltys),cex=1,lwd=2,bty="n",inset=leg)
mtext("even-year broodline",side=3,cex=0.8,font=2,line=-1.4)
legend("topleft","a",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##----------------------------------------------------## even
obj_esc_odd<-visreg(mod,xvar="esc",by=by,cond=list(BL="odd"),partial=T,alpha=alphap,scale="response",plot=F)
obj_esc_odd_unscaled<-visreg(mod,xvar="esc",by=by,cond=list(BL="odd"), partial=T,alpha=alphap,scale="response",plot=F,xtrans=xtrans)
obj_esc_odd$res$BL<-NA ## create new BL column with actual BL by year
for(i in 1:dim(obj_esc_odd$res)[1]) { obj_esc_odd$res$BL[i]<-data$BL[data$esc==obj_esc_odd$res$esc[i]] }
index<-which(obj_esc_odd$res$BL=="odd") ## get row index 
obj_esc_odd$res<-obj_esc_odd$res[index,] ## drop residuals of other BL
obj_esc_odd$res$esc<-obj_esc_odd_unscaled$res$esc[index] ## replace with unscaled escapement
obj_esc_odd$fit$esc<-obj_esc_odd$fit$esc*data_sds[names(data_sds)=="esc"]+data_means[names(data_means)=="esc"]
plot(obj_esc_odd,overlay=T,legend=F,partial=T,xlab="Spawner abundance (millions)",ylab=ylab,ylim=ylim,line=list(col=cols_R,lty=ltys), fill=list(col=alpha(cols_R,0.25)),points=list(pch=21,bg=cols_R,col=1,lwd=0.2,cex=cex.pt)) 
box()
##-----------------------------------------------------------## legend
legend("bottomleft",c("",""), col=alpha(rev(cols_R),0.25),cex=1,lwd=7,bty="n",seg.len=2,inset=leg)
labs<-c("up to 1988","since 1989")
legend("bottomleft",labs,col=rev(cols_R),lty=rev(ltys),cex=1,lwd=2,bty="n",inset=leg)
mtext("odd-year broodline",side=3,cex=0.8,font=2,line=-1.4)
legend("topleft","b",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##============================================================## temperture
#xtrans<-function(x) { x*data_sds[names(data_sds)=="temp"]+data_means[names(data_means)=="temp"] }
obj_temp_odd<-visreg(mod,xvar="temp",by="regime",xlab="SST anomaly (ºC)", overlay=T, ylab=ylab,ylim=ylim,partial=T,alpha=alphap, legend=F,scale="response",main="", line=list(col=cols_R,lty=ltys), fill=list(col=alpha(cols_R,0.25)), points=list(pch=21,bg=cols_R,col=1,lwd=0.2,cex=cex.pt)) ## ,xtrans=xtrans)
##-----------------------------------------------------------## legend
legend("bottomleft",c("",""), col=alpha(rev(cols_R),0.25),cex=1,lwd=7,bty="n",seg.len=2,inset=leg)
legend("bottomleft",c("up to 1988","since 1989"),col=rev(cols_R),lty=rev(ltys),cex=1,lwd=2,bty="n",inset=leg)
legend("topleft","c",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##==========================================================## total return
xtrans<-function(x) { x*data_sds[names(data_sds)=="totr"]+data_means[names(data_means)=="totr"] }
obj_totr<-visreg(mod,xvar="totr",overlay=T,legend=F,partial=T,alpha=alphap, scale="response",main="",xlab="Total return (millions)",ylab=ylab,ylim=ylim, line=list(col=col), fill=list(col=alpha(col,0.25)), points=list(pch=21,bg=col_pt,col=1,lwd=0.2,cex=cex.pt),xtrans=xtrans)
legend("topleft","d",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##====================================================## hatchery releases
xtrans<-function(x) { x*data_sds[names(data_sds)=="srel"]+data_means[names(data_means)=="srel"] }
obj_rel<-visreg(mod,xvar="srel",overlay=T,legend=F,partial=T,alpha=alphap, scale="response",main="",xlab="Hatchery releases (millions)",ylab=ylab,ylim=ylim, line=list(col=col), fill=list(col=alpha(col,0.25)), points=list(pch=21,bg=col_pt,col=1,lwd=0.2,cex=cex.pt),xtrans=xtrans)
legend("topleft","e",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##====================================================## competitor biomass
xtrans<-function(x) { x*data_sds[names(data_sds)=="comp"]+data_means[names(data_means)=="comp"] }
obj_comp<-visreg(mod,xvar="comp",overlay=T,legend=F,partial=T, alpha=alphap,scale="response",main="", xlab="Competitor abundance (millions)",ylim=ylim,ylab=ylab,line=list(col=col), fill=list(col=alpha(col,0.25)), points=list(pch=21,bg=col_pt,col=1,lwd=0.2,cex=cex.pt), xtrans=xtrans)
legend("topleft","f",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##==================================================================## save
dev.off()

##================================================## broodline main effect
pdf("PWS_wild_pinks_lnRS_vs_covars_selected_BLeffect.pdf",width=4.2,height=4.2)
par(mar=c(4,4,1,1),oma=c(0,0,0.5,0),mgp=c(2.2,0.5,0), cex.axis=1.2,cex.lab=1.3,tcl=-0.3,pch=21)
visreg(mod,xvar="BL",xlab="Broodline",ylab="ln(recruits/spawner)",overlay=T,legend=T,partial=T,alpha=alphap,scale="response", main="", line=list(col=col), fill=list(col=alpha(col,0.25)), points=list(pch=21,bg=col_pt,col=1,lwd=0.2,cex=1.1))
dev.off()

##=======================================================================##
##======================================================## cross-validation
##=======================================================================##
## out-of-sample predictions (using train and test data)
##========================================================## models to test
mod_select<-dredge(mod,trace=F,rank=rank)
mod_list<-get.models(mod_select,subset=delta<100) ## sub-models of selected
##------------------------------------------------------## get model terms
nM<-length(mod_list)
test_forms<-list()
for(i in 1:nM) {
mymod<-mod_list[[i]]
terms<-attr(mymod$terms,"term.labels")
nterms<-length(terms)
if(nterms!=0) my_mod<-formula("lnRS~1") ## intercept only model
if(nterms!=0) my_mod<-formula(paste("lnRS~",paste0(terms,collapse="+")))
test_forms[[i]]<-my_mod
}
##-----------------------------------## randomly sample train and test data
nS<-1e3 ## number of runs
RMSE<-array(dim=c(nS,nM))
formulas<-list()
start<-Sys.time()
for(i in 1:nS) {
set.seed(i)
train<-sort(sample(seq(nY),round(0.75*nY),replace=F)) ## 75%
test<-seq(nY)[seq(nY) %!in% train]
traindata<-data[train,] 
testdata<-data[test,]
##-------------------------------------------------## loop model structrues
for(j in 1:nM) { 
test_mod<-test_forms[[j]]
trainmod<-lm(test_mod,data=traindata)
##-----------------------------------------------------------## predictions
predicted<-predict(trainmod,newdata=list(esc=testdata$esc, BL=testdata$BL,regime=testdata$regime,temp=testdata$temp, srel=testdata$srel,comp=testdata$comp,totr=testdata$totr), se.fit=T,se=T,type="response")
##-----------------------------------------------## root mean squared error
pred<-as.numeric(predicted$fit)
true<-as.numeric(testdata$lnRS)
RMSE[i,j]<-sqrt(sum((pred-true)^2)/mean(true)) 
##---------------------------------------------------------## save formulas
prev_mod<-trainmod
save_terms<-attr(prev_mod$terms,"term.labels")
allterms<-paste(save_terms,collapse="+") 
new_mod<-formula(paste("lnRS~",allterms,sep=""))
formulas[[j]]<-new_mod
} ## end loop over models
} ## end stochastic loop
end<-Sys.time()
print(end-start)

##=============================================================## plot RMSE
pdf("PWS_wild_pinks_lnRS_vs_covars_RMSE.pdf",height=4,width=6.5)
par(mar=c(4,4,1,1),mgp=c(2.25,0.5,0),xaxs="i",yaxs="i", cex.axis=1.1,cex.lab=1.2,tcl=-0.3,las=1)
plot_RMSEs<-apply(RMSE,2,function(x) quantile(x,prob=c(0.5,0.05,0.25,0.75,0.95)))
nx<-dim(RMSE)[2]
xvec<-seq(nx)
plot(NA,NA,xlim=c(-1,nx+2),ylim=c(1,4),xlab="Model #",ylab="RMSE")
segments(xvec,plot_RMSEs[2,],xvec,plot_RMSEs[5,],lwd=0.5,col="gray10")
segments(xvec,plot_RMSEs[3,],xvec,plot_RMSEs[4,],lwd=1.5,col="gray10")
points(xvec,plot_RMSEs[1,],pch=21,cex=0.7,lwd=0.5,bg="white")
points(1,plot_RMSEs[1,1],pch=21,cex=0.7,lwd=0.5,bg="firebrick")
dev.off()

##=======================================================================##
##=======================================================================##
##=======================================================================##