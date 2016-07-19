# this script uses hybrid mixture models for predictive modeling

# install packages
install.packages(stringr)
install.packages(matrixStats)

# library
library(stringr)
library(matrixStats)

# data filename
filename1="pathway to rmm_predict.R"
filename2="pathway to kddcup.data_10_percent_corrected"

# function definition
source(filename1)

# data preprocessing
delim = ","  # or is it "\t" ?
dec = "."    # or is it "," ?
para=read.csv(filename2,sep=",",header=FALSE)
para=para[rowSums(para=="?")==0,]
para[]=lapply(para, as.numeric)
para=data.matrix(para,rownames.force = NA)
para=para[sample(nrow(para)),]
para=para[1:(10*floor(nrow(para)/10)),]
# para=cbind(para[,3:ncol(para)],para[,2])

# ------------------------------ 10-fold cross-validation --------------------

# parameter initialization
sample_size=nrow(para)/10
acc1=matrix(0,nrow=0,ncol=0)
acc2=matrix(0,nrow=0,ncol=0)
acc3=matrix(0,nrow=0,ncol=0)
acc4=matrix(0,nrow=0,ncol=0)
acc5=matrix(0,nrow=0,ncol=0)
acc6=matrix(0,nrow=0,ncol=0)
acc7a=matrix(0,nrow=0,ncol=0)
acc7b=matrix(0,nrow=0,ncol=0)
acc7=matrix(0,nrow=0,ncol=0)
acc8=matrix(0,nrow=0,ncol=0)
acc9=matrix(0,nrow=0,ncol=0)
acc10=matrix(0,nrow=0,ncol=0)
y1=matrix(0,nrow=0,ncol=0)
v=matrix(0,nrow=0,ncol=0)
f1=c(0.5,0.5,0.5,0.5)

for (j in seq(1,nrow(para),sample_size)){

	# data sampling
	para1=para[seq(-j,-(j+sample_size-1),-1),]
	para2=para[seq(j,j+sample_size-1,1),]
	para3=para2

	# whitened principal component analysis
	y1=wpca(para1,0)

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

# model training in batches
for (m in seq(1,nrow(y1$x),1e5)){

	if (m+99999 <= nrow(y1$x)) y1a=y1$x[seq(m,m+99999,1),]
	else y1a=y1$x[seq(m,nrow(y1$x),1),]	
	
	gamma_o=matrix(0,nrow=0,ncol=0)
	GMM_gamma=matrix(0,nrow=0,ncol=0)
	SMM_gamma=matrix(0,nrow=0,ncol=0)

	# model training
	gamma_o=Gmeans_clust(y1a,1e-4)
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
	if (SMM_gamma$tol>=0 & SMM_gamma$tol<=0.5 & index2==0 & length(SMM_gamma$w)>0) {
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
	if (SMM_gamma$tol>=0 & SMM_gamma$tol<=0.5 & index2==1 & length(SMM_gamma$w)>0) {
		if (sum(SMM_gamma$sigma==0)==0) { 
			SMM_gamma_adapt=model_adapt(f1,N1b,N2,SMM_gamma_adapt$w,SMM_gamma$w,SMM_gamma_adapt$mu,SMM_gamma$mu,SMM_gamma_adapt$sigma,SMM_gamma$sigma,SMM_gamma_adapt$v,SMM_gamma$v)
			N1b=N1b+N2}}
}
	prob1=matrix(0,nrow=nrow(para3),ncol=1)
	prob2=matrix(0,nrow=nrow(para3),ncol=1)
	prob3=matrix(0,nrow=nrow(para3),ncol=1)
	prob4=matrix(0,nrow=nrow(para3),ncol=1)
	prob5=matrix(0,nrow=nrow(para3),ncol=1)
	prob6=matrix(0,nrow=nrow(para3),ncol=1)

	class1=matrix(0,nrow=nrow(para3),ncol=1)
	class2=matrix(0,nrow=nrow(para3),ncol=1)
	class3=matrix(0,nrow=nrow(para3),ncol=1)
	class4=matrix(0,nrow=nrow(para3),ncol=1)
	class5=matrix(0,nrow=nrow(para3),ncol=1)
	class6=matrix(0,nrow=nrow(para3),ncol=1)

for (k in 1:max(para[,ncol(para)])){

	# whitened principal component analysis
	para3[,ncol(para3)]=k
	para3a=(para3-repmat(colMeans(para1),nrow(para3),1))
	y3=para3a%*%y1$rotation/(sqrt(nrow(para1-1))*repmat(y1$sdev,nrow(para3a),1))

	# predictive density
	# KG
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapto$w)){
		p1=p1+matrix(GMM_gamma_adapto$w[i]*gauss(y3,GMM_gamma_adapto$mu[i,],GMM_gamma_adapto$sigma[i,]),ncol=1)}
	p1[is.na(p1)]=0
	class1[p1>prob1]=k
	prob1[p1>prob1]=p1[p1>prob1]

	# KS
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapto$w)){ 
		p1=p1+matrix(GMM_gamma_adapto$w[i]*student(y3,GMM_gamma_adapto$mu[i,],GMM_gamma_adapto$sigma[i,],1),ncol=1)}
	p1[is.na(p1)]=0
	class2[p1>prob2]=k
	prob2[p1>prob2]=p1[p1>prob2]

	# KEG
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapt$w)){
		p1=p1+matrix(GMM_gamma_adapt$w[i]*gauss(y3,GMM_gamma_adapt$mu[i,],GMM_gamma_adapt$sigma[i,]),ncol=1)}
	p1[is.na(p1)]=0
	class3[p1>prob3]=k
	prob3[p1>prob3]=p1[p1>prob3]
	
	# KES
	#p1=matrix(0,nrow=nrow(y3),ncol=1)
	#for (i in 1:length(SMM_gamma_adapt$w)){
	#	p1=p1+matrix(SMM_gamma_adapt$w[i]*student(y3,SMM_gamma_adapt$mu[i,],SMM_gamma_adapt$sigma[i,],SMM_gamma_adapt$sigma[i]),ncol=1)}
	#p1[is.na(p1)]=0
	#class4[p1>prob4]=k
	#prob4[p1>prob4]=p1[p1>prob4]

	# KEGS
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapt$w)){ 
		p1=p1+matrix(GMM_gamma_adapt$w[i]*student(y3,GMM_gamma_adapt$mu[i,],GMM_gamma_adapt$sigma[i,],1),ncol=1)}
	p1[is.na(p1)]=0
	class5[p1>prob5]=k
	prob5[p1>prob5]=p1[p1>prob5]

	# KESG
	#p1=matrix(0,nrow=nrow(y3),ncol=1)
	#for (i in 1:length(SMM_gamma_adapt$w)){
	# 	p1=p1+matrix(SMM_gamma_adapt$w[i]*gauss(y3,SMM_gamma_adapt$mu[i,],SMM_gamma_adapt$sigma[i,]),ncol=1)}
	#p1[is.na(p1)]=0
	#class6[p1>prob6]=k
	#prob6[p1>prob6]=p1[p1>prob6]
}
	# prediction accuracies
	acc1=c(acc1,apply(matrix(class1,nrow=1),1,"==",matrix(para2[,ncol(para2)],nrow=1)))
	acc2=c(acc2,apply(matrix(class2,nrow=1),1,"==",matrix(para2[,ncol(para2)],nrow=1)))
	acc3=c(acc3,apply(matrix(class3,nrow=1),1,"==",matrix(para2[,ncol(para2)],nrow=1)))
	acc4=c(acc4,apply(matrix(class4,nrow=1),1,"==",matrix(para2[,ncol(para2)],nrow=1)))
	acc5=c(acc5,apply(matrix(class5,nrow=1),1,"==",matrix(para2[,ncol(para2)],nrow=1)))
	acc6=c(acc6,apply(matrix(class6,nrow=1),1,"==",matrix(para2[,ncol(para2)],nrow=1)))
}
acc=c(sum(acc1)/nrow(para),sum(acc2)/nrow(para),sum(acc3)/nrow(para),sum(acc4)/nrow(para),sum(acc5)/nrow(para),sum(acc6)/nrow(para))
print(acc)