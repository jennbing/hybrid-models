# this script uses hybrid mixture models for parameter rating

# install packages
install.packages(stringr)
install.packages(matrixStats)

# library
library(stringr)
library(matrixStats)

# filename
filename1="pathway to rmm_predict.R"
filename2="pathway to kddcup.data_10_percent_corrected"
filename3="pathway to image output"

# function definition
source(filename1)

# data extraction
delim = ","  # or is it "\t" ?
dec = "."    # or is it "," ?
para=read.csv(filename2,sep=",",header=FALSE)
para=para[rowSums(para=="?")==0,]
level42=levels(para[,42])
para[]=lapply(para, as.numeric)
paraA=para[para[,ncol(para)]==12,]
paraB=para[para[,ncol(para)]!=12,]
paraA=data.matrix(paraA,rownames.force = NA)
paraB=data.matrix(paraB,rownames.force = NA)
para=data.matrix(para,rownames.force = NA)
para=para[sample(nrow(para)),]
paraC=paraA

# parameter initialization
y1=matrix(0,nrow=0,ncol=0)
v=matrix(0,nrow=0,ncol=0)
f1=c(1,1,1,1)
N1=0
N2=0
N1o=0
N1a=0
N1b=0
index=0.1
index0=0
index1=0
index2=0
GMM_gamma_adapto=0
GMM_gamma_adapt=0
SMM_gamma_adapt=0

# whitened principal component analysis
y1=wpca(paraC[,1:ncol(paraC)-1],0)
y1$x=y1$x[,1:20]
y1$sdev=y1$sdev[1:20]
y1$rotation=y1$rotation[,1:20]

# model training in batches
for (m in seq(1,nrow(y1$x),10000)){

	if (m+9999 <= nrow(y1$x)) y1a=y1$x[seq(m,m+9999,1),]
	else y1a=y1$x[seq(m,nrow(y1$x),1),]	
	
	gamma_o=matrix(0,nrow=0,ncol=0)
	GMM_gamma=matrix(0,nrow=0,ncol=0)
	SMM_gamma=matrix(0,nrow=0,ncol=0)

	# model training
	gamma_o=Gmeans_clust(y1a,0.0001)
	GMM_gamma=GMM_ML(y1a,0.05,100)
	while ((GMM_gamma$tol==0 | is.infinite(GMM_gamma$tol)==1) & (index<=0.2)) {
		GMM_gamma=GMM_ML(y1a,index,100)
		index=index+0.1}
	index=0.1

	SMM_gamma=SMM_ML(y1a,0.05,100) # Student t distribution Mixture Model (SMM) parameters estimation
	while ((SMM_gamma$tol==0 | is.infinite(SMM_gamma$tol)==1) & (index<=0.2)) {
		SMM_gamma=SMM_ML(y1a,index,100)
		index=index+0.1}
	index=0.1

	# model adaptation
	if (index0==0 & length(gamma_o$w)>0) {
		if (sum(gamma_o$sigma==0)==0) {
			GMM_gamma_adapto=gamma_o
			index0=1
			N1=nrow(y1a)}}
	if (GMM_gamma$tol>=0 & GMM_gamma$tol<=0.2 & index1==0 & length(GMM_gamma$w)>0) {
		if (sum(GMM_gamma$sigma==0)==0) { 
			GMM_gamma_adapt=GMM_gamma
			index1=1
			N1a=nrow(y1a)}}
	if (SMM_gamma$tol>=0 & SMM_gamma$tol<=0.2 & index2==0 & length(SMM_gamma$w)>0) {
		if (sum(SMM_gamma$sigma==0)==0) { 
			SMM_gamma_adapt=SMM_gamma
			index2=1
			N1b=nrow(y1a)}}
	N2=nrow(y1a)
	if (index0==1 & length(gamma_o$w)>0) {
		if (sum(gamma_o$sigma==0)==0) {
			GMM_gamma_adapto=model_adapt(f1,N1,N2,GMM_gamma_adapto$w,gamma_o$w,GMM_gamma_adapto$mu,gamma_o$mu,GMM_gamma_adapto$sigma,gamma_o$sigma,v,v)
			N1=N1+N2}}
	if (GMM_gamma$tol>=0 & GMM_gamma$tol<=0.2 & index1==1 & length(GMM_gamma$w)>0) {
		if (sum(GMM_gamma$sigma==0)==0) { 
			GMM_gamma_adapt=model_adapt(f1,N1a,N2,GMM_gamma_adapt$w,GMM_gamma$w,GMM_gamma_adapt$mu,GMM_gamma$mu,GMM_gamma_adapt$sigma,GMM_gamma$sigma,v,v)
			N1a=N1a+N2}}
	if (SMM_gamma$tol>=0 & SMM_gamma$tol<=0.2 & index2==1 & length(SMM_gamma$w)>0) {
		if (sum(SMM_gamma$sigma==0)==0) { 
			SMM_gamma_adapt=model_adapt(f1,N1b,N2,SMM_gamma_adapt$w,SMM_gamma$w,SMM_gamma_adapt$mu,SMM_gamma$mu,SMM_gamma_adapt$sigma,SMM_gamma$sigma,SMM_gamma_adapt$v,SMM_gamma$v)
			N1b=N1b+N2}}
}	
	# whitened principal component analysis
	para1=(para-repmat(colMeans(paraC),nrow(para),1))
	y3=para1[,1:41]%*%y1$rotation/((nrow(paraC)-1)*repmat(y1$sdev,nrow(para1),1))
	y3=y3[,1:20]

	# parameter rating
	para_rate=matrix(0,0,0)
	color=c("magenta","black","red","blue")
	png(filename3,width=6,height=8,units="in",res=300)
	par(mfrow=c(3,2),oma = c(0,0,2,0),mar=c(2,2,1,1),mgp=c(2,1,0))
for (n in 1:max(para[,ncol(para)])){
	xbar_KS=par_rate(paraC[,1:41],y1$x,para[para[,ncol(para)]==n,1:41],y3[para[,ncol(para)]==n,],GMM_gamma_adapto$w,GMM_gamma_adapto$mu,GMM_gamma_adapto$sigma,matrix(1,nrow=nrow(GMM_gamma_adapto$w),ncol=1),2)
	xbar_KS[is.infinite(xbar_KS)]=NaN
	xbar_KS=apply(xbar_KS,2,mean,na.rm=TRUE)
	if (n==1) para_rate=matrix(xbar_KS,ncol=1)
	if (n>1) para_rate=cbind(para_rate,matrix(xbar_KS,ncol=1))
	xbar_KS=xbar_KS/max(xbar_KS)
	if ((n%%4)==1) plot(seq(1,length(xbar_KS),1),xbar_KS,type="b",lty=1,pch=1,col=color[(n%%4)+1])
	else points(seq(1,length(xbar_KS),1),xbar_KS,type="b",lty=1,pch=1,col=color[(n%%4)+1])#,main=bquote("KS, f="~.(f1)))
	if ((n%%4)==1) legend("topright",legend=level42[n:(n+3)],lty=1,col=c("black","red","blue","magenta"))
}
	mtext("Network Intrusion Parameter Rating", outer = TRUE, cex = 1)
	dev.off()