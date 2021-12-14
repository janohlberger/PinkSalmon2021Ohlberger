##=======================================================================##
##                                                                       ##
## Wild PWS pink salmon linear Ricker stock-recruit model and covariates ##
##                                                                       ##
##=======================================================================##
pkgs<-c("readxl","tidyverse","ggplot2","viridis","MuMIn","effects","Hmisc", "visreg","RColorBrewer","car","fmsb","mgcv","ggplot2","fabricatr","asbio", "performance","R.utils","ggiraphExtra","writexl","jtools","lattice","grid","gridExtra","nlme","pedometrics","gsl","modelsummary","report","randomForest")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) { install.packages(setdiff(pkgs,rownames(installed.packages())),dependencies=T) }
invisible(lapply(pkgs,library,character.only=T)) ## install and load 
##=============================================================## directory
if(exists("mainDir")) { } else { mainDir<-setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) } ## set home/main directory
##=============================================================## functions
is.even<-function(x) { x %% 2 == 0 }
'%!in%'<-function(x,y)!('%in%'(x,y))

##=======================================================================##
##==================================================================## data
##=======================================================================##
data<-read.csv("Data.csv") ## full dataset
data<-data[data$Year>1960,]
##===========================================## explanation of column names
## Year	- return year
## BL	- brood line (odd/even)
## regime - ocean climate regime (upto1988/since1989)
## Escapement	- wild PWS pink salmon spawner escapement (millions)
## Run.Size	- wild PWS pink salmon run size, i.e. total return (millions)
## Spawners	- wild PWS pink salmon spawners two years prior (millions)
## RS	- wild PWS pink salmon recruits/spawner
## lnRS	- wild PWS pink salmon ln(recruits/spawner)
## prevYearRunSize - wild PWS pink salmon previous year run size (millions)
## prev_tot_PWS_return - previous year return wild plus hatchery (millions)
## PWS_released_Total	- total releases of hatchery pink salmon (millions)
## PWS_released_prev_Total - hatchery releases previous year (millions)
## sameYrcomp_regional - regional sockeye+chum salmon abundance (millions)
## sameYrcompetitors - NE Pacific sockeye+chum salmon abundance (millions)
## egoa.win.sst	- eastern Gulf of Alaska winter temperature anomaly 
## egoa.spr.sst	- eastern Gulf of Alaska spring temperature anomaly 
## wgoa.win.sst	- western Gulf of Alaska winter temperature anomaly 
## wgoa.spr.sst	- western Gulf of Alaska spring temperature anomaly 
## PDOann	- annual Pacific Decadal Oscillation index
## NPGOann - annual North Pacific Gyre Oscillation index
## pCO2_PWS	- pCO2 Prince William Sound
## pCO2_nGoA - pCO2 northern Gulf of Alaska
## pH_PWS	- pH Prince William Sound
## pH_nGoA - pH northern Gulf of Alaska
## mean_pCO2 - pCO2 subpolar gyre 
## mean_pH - pCO2 subpolar gyre 
## mean_NGAO - NGAO index subpolar gyre 
##========================================================## rename columns
##-------------------------------------------------------## spawners/return
data$esc<-data$Spawners 
data$runs<-data$prevYearRunSize 
data$totr<-data$prev_tot_PWS_return 
##-----------------------------------------------------## hatchery releases
data$srel<-data$PWS_released_Total 
##------------------------------------------------------## competitor index
data$comp<-data$sameYrcomp_regional 
# data$comp<-data$sameYrcompetitors 
##------------------------------------------------------------## SST metric
data$temp<-data$egoa.spr.sst 
# data$temp<-data$wgoa.spr.sst
# data$temp<-data$egoa.win.sst
# data$temp<-data$wgoa.win.sst 
##--------------------------------------------------------------## OA proxy
data$OA<-data$pCO2_nGoA 
# data$OA<-data$mean_pCO2 
# data$OA<-data$pH_nGoA  
# data$OA<-data$mean_pH
# data$OA<-data$mean_NGAO
##--------------------------------------------------------------## PDO/NPGO
# data$npgo<-data$NPGOann ## annual NPGO
# data$pdo<-data$PDOann ## annual PDO
  
##===========## restrict to selected covariates when using full time series
data<-dplyr::select(data,Year,BL,lnRS,regime,esc,temp,totr,comp,srel)
data<-data[complete.cases(data),] ## drop years with NAs
nY<-dim(data)[1]
rownames(data)<-seq(1,dim(data)[1],1)
##----------------------------------------------------------------## notes
## OA time series shorter >> data and model selection need to be re-run
## multiple temperature or competitor metrics not to be included together

##==========================================================## data scaling
data_unscaled<-dplyr::select(data,-Year,-BL,-lnRS,-regime) 
data_other<-dplyr::select(data,Year,BL,lnRS,regime)
data_means<-sapply(data_unscaled,mean)
data_sds<-sapply(data_unscaled,sd)
##-------------------------------------## center/standardize to mean=0/sd=1
data<-data.frame(cbind(data_other,scale(data_unscaled)))
data_unscaled_all<-data.frame(cbind(data_other,data_unscaled))

##=======================================================================##
##======================================## multiple linear regression model 
##=======================================================================##
options(na.action="na.fail")	## avoid fits to different datasets if NAs
rank<-"AICc" ## model selection criterion (AICc or BIC)
##------------------------------------------------------------## model form
mod_form<-formula(lnRS~esc*BL+esc:BL:regime+temp+I(temp^2)+temp:regime+totr+I(totr^2)+srel+I(srel^2)+comp+I(comp^2)) 
##-------------------------------------------------------## model selection
fixed<-c("esc","BL","BL:esc") ## stock-recruit relationship each broodline
mod_full<-lm(mod_form,data=data) 
mod_select<-dredge(mod_full,fixed=fixed,trace=F,rank=rank) 
##--------------------------------------------------------## selected model
mod<-get.models(mod_select,subset=1)[[1]]
out<-summary(mod)
out$r.squared ## r2: ~62% (selected variables without OA and PDO/NPGO)
AICc(mod) 
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
anova_results<-data.frame(anova(mod))

##===============================================## reporting model results
report_model(mod)
report_performance(mod)
report_statistics(mod)
data.frame(tidy(mod))

##================================================## model selection table
num_of_mods<-10 ## top X models
mod_select_topX<-mod_select[1:num_of_mods,]
mod_select_topX$cum_weight<-cumsum(mod_select_topX$weight)
##-----------------------------------------------------------## model forms
mod_forms<-get.models(mod_select,subset=delta<100)
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
mod_covars<-lapply(mod_covars,function(x) substr(x,1,nchar(x)-2))
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
pdf("Plot_covariate_effects.pdf",width=9,height=6.2)
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
pdf("Plot_broodline_effect.pdf",width=4.2,height=4.2)
par(mar=c(4,4,1,1),oma=c(0,0,0.5,0),mgp=c(2.2,0.5,0), cex.axis=1.2,cex.lab=1.3,tcl=-0.3,pch=21)
visreg(mod,xvar="BL",xlab="Broodline",ylab="ln(recruits/spawner)",overlay=T,legend=T,partial=T,alpha=alphap,scale="response", main="", line=list(col=col), fill=list(col=alpha(col,0.25)), points=list(pch=21,bg=col_pt,col=1,lwd=0.2,cex=1.1))
dev.off()

##=========================================## parameter estimates with CIs
pdf("Plot_model_estimates_with_CIs.pdf",width=5,height=6)
par(mar=c(4,8,1,1),mgp=c(2.2,0.5,0),cex.lab=1,tcl=-0.3)
renamed<-c("(Intercept)"="Intercept","I(comp^2)"="Competitor abundance\n(quardatic)","srel"="Hatchery releases","temp"="SST anomaly","I(totr^2)"="Total return\n(quardatic)","BLodd"="Broodline odd","esc"="Spawner abundance","temp:regimeupto1988"="SST anomaly\n:regime up to 1988","BLodd:esc"="Spawner abundance\n:broodline odd","regimeupto1988:BLeven:esc"="Spawner abundance\n:regime up to 1988\n:broodline even","regimeupto1988:BLodd:esc"="Spawner abundance\n:regime up to 1988\n:broodline odd")
p<-modelplot(mod,coef_rename=renamed,coef_omit='Interc')
pdat<-data.frame(p$data)
pdat<-pdat[rev(c(5,3,1,2,8,4,7,9,10)),]
names<-pdat$term
nnames<-length(names)
y<-seq(nnames)
xlab="Coefficient estimates with 95% CIs"
plot(NA,NA,xlab=xlab,ylab="",xlim=c(-2.2,2.2),ylim=c(0,nnames),yaxt="n")
abline(v=0,lwd=0.5)
segments(pdat$conf.low,y,pdat$conf.high,y,lwd=1.2)
points(pdat$estimate,y,pch=21,cex=1.2,bg="darkgray",lwd=0.5)
axis(side=2,at=y,labels=pdat$term,las=2,cex.axis=0.8)
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
pdf("Plot_cross_validation_RMSE.pdf",height=4,width=6.5)
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
##================================## variable importance via random forests
##=======================================================================##
data$comp2 = data$comp^2
data$totr2 = data$totr^2
##--------------------------------------------------## type of measure used
type<-1
## 1=%IncMSE: % increase in mean squared error when variable is excluded
## 2=IncNodePurity: increase in node purity calculated based on reduction in sum of squared errors whenever a variable is chosen to split
##-----------------------------------------------------## fit random forest
ntrees<-1e4
rf<-randomForest(lnRS~comp2+srel+temp+totr2+regime:temp+BL:esc:regime+1+BL+esc+BL:esc,data=data,importance=TRUE,ntree=ntrees,mtry=3,replace=FALSE)
##---------------------------------------------------## variable importance
variable_importance<-importance(rf,type=type)
# varImpPlot(rf)
##----------------------------------------------------## get variable names
index<-sort(as.numeric(variable_importance),index.return=T)
ind<-index$ix
varnames<-rownames(variable_importance)[ind]
varnames<-gsub("esc","Spawners",varnames)
varnames<-gsub("srel","Hatchery releases",varnames)
varnames<-gsub("comp2","Competitors",varnames)
varnames<-gsub("temp","SST anomaly",varnames)
varnames<-gsub("BL","Broodline",varnames)
varnames<-gsub("totr2","Total return",varnames)
varnames<-gsub("regime","Regime",varnames)
nvars<-length(varnames)
##----------------------------------------------------------## plot results
pdf("Plot_variable_importance_random_forests.pdf",width=5,height=4)
par(mar=c(4,8,1,1),mgp=c(2,0.5,0),cex.axis=1,cex.lab=1.2,tcl=-0.3)
x<-index$x
y<-seq(nvars)
cols<-rep("black",nvars)
cols[varnames %in% c(c("Spawners","Broodline"))]<-"darkgray"
if(type==1) xlab<-"Percent increase in MSE"
if(type==2) xlab<-"Increase in node purity"
plot(x,y,yaxt="n",xlab=xlab,pch=21,cex=1.5,ylab="",xlim=c(0,max(x)*1.05),ylim=c(0.9,0.1+nvars),bg=cols,lwd=0.5)
axis(2,at=y,labels=varnames,las=2)
dev.off()

##=======================================================================##
##=====================================================## model diagnostics
##=======================================================================#$
summary(mod)$call
##==============================================## extract model components
response<-data$lnRS
time<-data$Year
fitted<-fitted(mod)
covariate1<-data$esc
covariate2<-data$prel
covariate3<-data$comp
factor1<-data$BL
factor2<-data$regime
residuals<-residuals(mod) ## standardized
r2(mod) 
check_singularity(mod)
model_performance(mod,rank=TRUE)
## residual standard error 
k<-length(mod$coefficients)-1 # -1 to ignore intercept
SSE<-sum(mod$residuals^2)
n<-length(mod$residuals)
res_std_err<-sqrt(SSE/(n-(1+k)))
round(res_std_err,2)
## adjusted r-squared
SSyy<-sum((response-mean(response))^2)
adjR2<-1-(SSE/SSyy)*(n-1)/(n-(k+1))
round(adjR2,2) ## variance explained

##==============================================## diagnostics plots for SI
pdf("Plot_model_diagnostics.pdf",width=6.5,height=6.5)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0), cex.lab=1.3,cex.axis=1,tcl=-0.3,pch=16)
layout(matrix(c(1:4),nrow=2,byrow=T))
##----------------------------------------## check for temporal correlation
plot(time,residuals,xlab="Year",ylab="Residuals", xlim=c(1962,2015),ylim=c(-1,1))
abline(h=0,lty=3)
#acf(residuals,lag.max=50,xlab="Time lag",ylab="Autocorrelation residuals")
##----------------------------------------------## Check for heterogeneiety
plot(fitted,residuals,xlab="Fitted values",ylab="Residuals")
abline(0,0,lty=3)
##-------------------------------------------------## observed vs predicted
plot(fitted,response,ylab="Observed",xlab="Predicted", xlim=c(-0.75,2.25),ylim=c(-0.75,2.25));abline(0,1,lty=3)
##------------------------------------------------## normality of residuals
qqnorm(residuals,main="");qqline(residuals)
dev.off()


##=======================================================================##
##=====================================================## model predictions
##=======================================================================##

################################################################### by year
predicted<-predict(mod,newdata=list(esc=data$esc,BL=data$BL,temp=data$temp, totr=data$totr,srel=data$srel,comp=data$comp,regime=data$regime),se.fit=T,se=T,type="response")
pred<-as.vector(predicted$fit)
upper<-as.vector(predicted$fit+2*predicted$se.fit)
lower<-as.vector(predicted$fit-2*predicted$se.fit)

##=============================================## plot ln(recruits/spawner)
pdf("Plot_lnRS_predictions.pdf",width=4.5,height=3.5)
par(mar=c(3.5,3.5,0.5,1),mgp=c(2,0.5,0),cex.axis=1,cex.lab=1.25,tcl=-0.3)
plot(NA,NA,xlab="Year",ylab="ln(recruits/spawner)", xlim=c(min(data$Year),max(data$Year)),ylim=c(-1.2,2.7))
abline(h=0,lty=1,col="grey50")
##-----------------------------------------------------## plot shading area 
X.Vec<-c(data$Year,tail(data$Year,1),rev(data$Year),data$Year[1])
Y.Vec<-c(lower,tail(upper,1),rev(upper),lower[1])
polygon(X.Vec,Y.Vec,col="grey90",border=NA) 
##------------------------------------------------------## full predictions
lines(data$Year,pred,lwd=1.5,lty=1,col="grey50")
lines(data$Year,upper,lwd=0.5,lty=1,col="grey50")
lines(data$Year,lower,lwd=0.5,lty=1,col="grey50")
##--------------------------------------------------## observed time-series
lines(data$Year,data$lnRS,type="o",pch=16,cex=0.8,lwd=1,col=1)
dev.off()

##=======================================================================##
##=====================================## wild ln(R/S) vs hatchery releases
##=======================================================================##
BL_scen<-unique(data$BL) # both broodlines
prob<-c(0.1,0.25,0.5,0.75,0.9)
##====================================================## release scenarios
nn<-1e6 ## number of samples when drawing normal dist
min_release<-min(data$srel)
max_release<-max(data$srel)
##----------------------------------------------------------## extrapolate
## calculate max release value equal two times the actual
true_max<-max(data$srel)*data_sds[names(data_sds)=="srel"]+data_means[names(data_means)=="srel"]
desired_max<-2*true_max
use_max<-(desired_max-data_means[names(data_means)=="srel"])/data_sds[names(data_sds)=="srel"]
##---------------------------------------------------------## set new range
release_scen<-seq(min_release,use_max,length=100)
nreleases<-length(release_scen)
##---------------------------------------------------------## which regime
use_regime<-"since1989" ## 'upto1988' or 'since1989' 

##==================================================================## even
newdat<-expand.grid(srel=release_scen,esc=median(data$esc[data$BL=="even"]), BL="even",regime=use_regime,temp=median(data$temp), comp=median(data$comp),totr=median(data$totr))
predicted<-predict(mod,newdata=newdat,se.fit=T,se=T,type="response")
predicted<-predict(mod,newdata=newdat,se.fit=T,se=T, interval="confidence",level=0.95,type="response")
nscen<-dim(newdat)[1]
##-----------------------------------------------------------## predictions
newdat$pred<-data.frame(predicted$fit)[,1]
newdat$se<-data.frame(predicted$se.fit)[,1]
newdat$lower<-data.frame(predicted$fit)[,2]
newdat$upper<-data.frame(predicted$fit)[,3]
##-------------------------------------------------## wild recruits/spawner
newdat$RperS_fit<-exp(newdat$pred) 	
newdat$RperS_low<-exp(newdat$lower) 	
newdat$RperS_upp<-exp(newdat$upper) 	
##-----------------------------------------------------------## rename data
newdat_even<-newdat

##===================================================================## odd
newdat<-expand.grid(srel=release_scen,esc=median(data$esc[data$BL=="odd"]), BL="odd",regime=use_regime,temp=median(data$temp), comp=median(data$comp),totr=median(data$totr))
predicted<-predict(mod,newdata=newdat,se.fit=T,se=T,interval="confidence",level=0.95,type="response")
nscen<-dim(newdat)[1]
##-----------------------------------------------------------## predictions
newdat$pred<-data.frame(predicted$fit)[,1]
newdat$se<-data.frame(predicted$se.fit)[,1]
newdat$lower<-data.frame(predicted$fit)[,2]
newdat$upper<-data.frame(predicted$fit)[,3]
##-------------------------------------------------## wild recruits/spawner
newdat$RperS_fit<-exp(newdat$pred) 	
newdat$RperS_low<-exp(newdat$lower) 	
newdat$RperS_upp<-exp(newdat$upper) 	
##-----------------------------------------------------------## rename data
newdat_odd<-newdat

##==================================================================## plot
newdat<-data.frame(rbind(newdat_even,newdat_odd))
pdf("Plot_predicted_RperS_vs_releases.pdf", width=4.2,height=4.2)
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.5,0),cex.lab=1.1,cex.axis=0.8, tcl=-0.3,xaxs="i",yaxs="i")
cols<-c("darkorchid4","goldenrod1") 
lnames<-c("even","odd")
ltys<-c(1,1)
cexll<-1.5 ## letter size
ll<-c(-0.28,-0.14) ## inset for letters
ymax<-7
maxval<-max(data_unscaled$srel)
xlim<-c(min(data_unscaled$srel), max(data_unscaled$srel))
xlim<-c(xlim[1],xlim[2]*2)
plot(NA,NA,xlim=xlim,ylim=c(0,ymax),xlab="Hatchery releases (millions)", ylab="Wild recruits/spawner")
abline(h=1,lty=1,lwd=1,col="gray")
for(i in 1:2) {
if(i==1) plotD<-newdat[newdat$BL=="even",] 
if(i==2) plotD<-newdat[newdat$BL=="odd",]
x<-plotD$srel
x<-x*data_sds[names(data_sds)=="srel"]+data_means[names(data_means)=="srel"]
y<-plotD$RperS_fit
lines(x,y,lty=ltys[1],lwd=1,col=cols[i]) 
y1<-plotD$RperS_low
y2<-plotD$RperS_upp
X.Vec<-c(x,tail(x,1),rev(x),x[1])
Y.Vec<-c(y1,tail(y2,1),rev(y2),y1[1])
polygon(X.Vec,Y.Vec,col=alpha(cols[i],0.1),border=NA) 
lines(x,y1,lty=ltys[1],lwd=0.5,col=cols[i])
lines(x,y2,lty=ltys[1],lwd=0.5,col=cols[i])
## extrapolated part of the relationship
x<-x[1:length(x)/2];y<-y[1:length(y)/2]
lines(x,y,lty=ltys[1],lwd=3,col=cols[i]) 
y1<-y1[1:length(y1)/2];y2<-y2[1:length(y2)/2]
X.Vec<-c(x,tail(x,1),rev(x),x[1])
Y.Vec<-c(y1,tail(y2,1),rev(y2),y1[1])
polygon(X.Vec,Y.Vec,col=alpha(cols[i],0.1),border=NA) 
lines(x,y1,lty=ltys[1],lwd=1,col=cols[i])
lines(x,y2,lty=ltys[1],lwd=1,col=cols[i])
}
abline(v=maxval,lty=1,lwd=3,col="gray")
box()
legend("topright",c("        ","        "), col=alpha(cols,0.2), cex=0.8,lwd=7,bty="n",seg.len=1.8) 
legend("topright",lnames, col=cols,lty=ltys,cex=0.8,lwd=2,bty="n")
dev.off()

##=======================================================================##
##=======================================================================##
##=======================================================================##