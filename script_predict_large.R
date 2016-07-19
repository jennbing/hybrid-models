# this script uses hybrid mixture models for large-scale predictive modeling

# install packages
install.packages("stringr")
install.packages("matrixStats")
install.packages("far")

# library
library(stringr)
library(matrixStats)
library(far)

# filename
filename1="pathway to rmm_predict.R"
filename2="pathway to corrected"
filename3="pathway to kddcup.data.corrected"
filename4="pathway to picture output"

# function code
source(filename1)

# data extraction
delim = ","  # or is it "\t" ?
dec = "."    # or is it "," ?
para=read.csv(filename2,sep=",",header=FALSE)
level2=levels(para[,2])
level3=levels(para[,3])
level4=levels(para[,4])
level42=levels(para[,42])
rm(para)

# count number of lines in file
testcon <- file(filename3,open="r")
readsizeof <- 20000
nooflines <- 0
( while((linesread <- length(readLines(testcon,readsizeof))) > 0 ) 
nooflines <- nooflines+linesread )
close(testcon)

# parameter initialization
v=matrix(0,nrow=0,ncol=0)
f1=c(1,1,1,1)
acc=matrix(0,nrow=0,ncol=0)
accA=matrix(0,nrow=0,ncol=0)
acc1=matrix(0,nrow=0,ncol=0)
acc2=matrix(0,nrow=0,ncol=0)
acc3=matrix(0,nrow=0,ncol=0)
acc4=matrix(0,nrow=0,ncol=0)
acc5=matrix(0,nrow=0,ncol=0)
acc6=matrix(0,nrow=0,ncol=0)
N1=0
N2=0
N1o=0
N1a=0
N1b=0
N2=0
index=0.1
index0=0
index1=0
index2=0
GMM_gamma_adapto=0
GMM_gamma_adapt=0
SMM_gamma_adapt=0
sample_size=1e5
label_max=length(level42)
pc_num=matrix(0,0,0)

# read data from file
con = file(filename,open="r")
# for (j in seq(1,nooflines,sample_size)){
for (j in seq(1,20e5,sample_size)){
	print(j)
	# data preprocessing
	delim = ","  # or is it "\t" ?
	dec = "."    # or is it "," ?
	para=matrix(0,0,0)
	para=read.table(con,sep=",",header=FALSE,nrows=sample_size)
	para[,2]=match(para[,2],level2)
	para[,3]=match(para[,3],level3)
	para[,4]=match(para[,4],level4)
	para[,42]=match(para[,42],level42)
	para=para[rowSums(para=="?")==0,]
	para=para[is.na(rowSums(para))==0,]
	para[]=lapply(para, as.numeric)
	para=data.matrix(para,rownames.force = NA)
	para=para[sample(nrow(para)),]

# ---------------------------- model testing ----------------------------------

if (j>1) {

	prob1=matrix(0,nrow=nrow(para),ncol=1)
	prob2=matrix(0,nrow=nrow(para),ncol=1)
	prob3=matrix(0,nrow=nrow(para),ncol=1)
	prob4=matrix(0,nrow=nrow(para),ncol=1)
	prob5=matrix(0,nrow=nrow(para),ncol=1)
	prob6=matrix(0,nrow=nrow(para),ncol=1)

	class1=matrix(0,nrow=nrow(para),ncol=1)
	class2=matrix(0,nrow=nrow(para),ncol=1)
	class3=matrix(0,nrow=nrow(para),ncol=1)
	class4=matrix(0,nrow=nrow(para),ncol=1)
	class5=matrix(0,nrow=nrow(para),ncol=1)
	class6=matrix(0,nrow=nrow(para),ncol=1)

for (k in 1:label_max){

	# whitened principal component analysis
	para3=para
	para3[,ncol(para3)]=k
	para3a=para3-repmat(para_mean,nrow(para3),1)
	y3=para3a%*%y2$rotation/(sqrt(sample_size-1)*repmat(y2$sdev,nrow(para3a),1))
	if (ncol(GMM_gamma_adapto$mu)<ncol(y3)) y3a=y3[,1:ncol(GMM_gamma_adapto$mu)]
	else y3a=y3
	if (ncol(GMM_gamma_adapt$mu)<ncol(y3)) y3b=y3[,1:ncol(GMM_gamma_adapt$mu)]
	else y3b=y3
	if (ncol(SMM_gamma_adapt$mu)<ncol(y3)) y3c=y3[,1:ncol(SMM_gamma_adapt$mu)]
	else y3c=y3

	# predictive density
	# KG
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapto$w)){
		p1=p1+matrix(GMM_gamma_adapto$w[i]*gauss(y3a,GMM_gamma_adapto$mu[i,],GMM_gamma_adapto$sigma[i,]),ncol=1)}
	p1[is.na(p1)]=0
	class1[p1>prob1]=k
	prob1[p1>prob1]=p1[p1>prob1]

	# KS
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapto$w)){ 
		p1=p1+matrix(GMM_gamma_adapto$w[i]*student(y3a,GMM_gamma_adapto$mu[i,],GMM_gamma_adapto$sigma[i,],1),ncol=1)}
	p1[is.na(p1)]=0
	class2[p1>prob2]=k
	prob2[p1>prob2]=p1[p1>prob2]

	# KEG
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapt$w)){
		p1=p1+matrix(GMM_gamma_adapt$w[i]*gauss(y3b,GMM_gamma_adapt$mu[i,],GMM_gamma_adapt$sigma[i,]),ncol=1)}
	p1[is.na(p1)]=0
	class3[p1>prob3]=k
	prob3[p1>prob3]=p1[p1>prob3]
	
	# KES
	p1=matrix(0,nrow=nrow(y3a),ncol=1)
	for (i in 1:length(SMM_gamma_adapt$w)){ 
		p1=p1+matrix(SMM_gamma_adapt$w[i]*student(y3c,SMM_gamma_adapt$mu[i,],SMM_gamma_adapt$sigma[i,],SMM_gamma_adapt$sigma[i]),ncol=1)}
	p1[is.na(p1)]=0
	class4[p1>prob4]=k
	prob4[p1>prob4]=p1[p1>prob4]

	# KEGS
	p1=matrix(0,nrow=nrow(y3),ncol=1)
	for (i in 1:length(GMM_gamma_adapt$w)){ 
		p1=p1+matrix(GMM_gamma_adapt$w[i]*student(y3b,GMM_gamma_adapt$mu[i,],GMM_gamma_adapt$sigma[i,],1),ncol=1)}
	p1[is.na(p1)]=0
	class5[p1>prob5]=k
	prob5[p1>prob5]=p1[p1>prob5]

	# KESG
	p1=matrix(0,nrow=nrow(y3a),ncol=1)
	for (i in 1:length(SMM_gamma_adapt$w)){
		p1=p1+matrix(SMM_gamma_adapt$w[i]*gauss(y3c,SMM_gamma_adapt$mu[i,],SMM_gamma_adapt$sigma[i,]),ncol=1)}
	p1[is.na(p1)]=0
	class6[p1>prob6]=k
	prob6[p1>prob6]=p1[p1>prob6]
}
	# prediction accuracy
	acc1=c(apply(matrix(class1,nrow=1),1,"==",matrix(para[,ncol(para)],nrow=1)))
	acc2=c(apply(matrix(class2,nrow=1),1,"==",matrix(para[,ncol(para)],nrow=1)))
	acc3=c(apply(matrix(class3,nrow=1),1,"==",matrix(para[,ncol(para)],nrow=1)))
	acc4=c(apply(matrix(class4,nrow=1),1,"==",matrix(para[,ncol(para)],nrow=1)))
	acc5=c(apply(matrix(class5,nrow=1),1,"==",matrix(para[,ncol(para)],nrow=1)))
	acc6=c(apply(matrix(class6,nrow=1),1,"==",matrix(para[,ncol(para)],nrow=1)))
}
accA=matrix(c(sum(acc1)/length(acc1),sum(acc2)/length(acc2),sum(acc3)/length(acc3),sum(acc4)/length(acc4),sum(acc5)/length(acc5),sum(acc6)/length(acc6)),nrow=1)
if (length(acc)==0) acc=accA
else acc=rbind(acc,accA)

# ----------------------------- model training -----------------------------------

	# whitened principal component analysis
	y1=matrix(0,nrow=0,ncol=0)
	para1=para
	y1=wpca(para1,0)
	y1$x=y1$x[,1:15]
	y1$sdev=y1$sdev[1:15]
	y1$rotation=y1$rotation[,1:15]
if (j==1) {
	y2=y1
	para_mean=colMeans(para1)
	y1b=y2$x}
else {
	rot_vec=y1$rotation
rot_vec1=matrix(0,0,0)
for (p in seq(1,ncol(rot_vec),1)){
for (n in seq(1,ncol(y2$rotation),1)){
	tryerror=try(orthonormalization(cbind(y2$rotation[,n],rot_vec[,p]),basis=FALSE,norm=FALSE)[,-1],silent=FALSE)
	if (length(tryerror)>1) rot_vec[,p]=tryerror
	else rot_vec[,p]=0}
if (sum(rot_vec[,p]^2)>0.9) {
	if (length(rot_vec1)==0) rot_vec1=rot_vec[,p]
	else rot_vec1=cbind(rot_vec1,rot_vec[,p])}
}
rot_vec=rot_vec1
if (length(rot_vec)>0){
	rot_vec=orthonormalization(rot_vec,basis=FALSE,norm=TRUE)
	y2$rotation=cbind(y2$rotation,rot_vec)
	y2$sdev=cbind(matrix(y2$sdev,nrow=1),matrix(colSums(para1%*%rot_vec/sqrt(sample_size-1)),nrow=1))
}
	para1A=para1-repmat(para_mean,nrow(para1),1)
	y1b=para1A%*%y2$rotation/(sqrt(sample_size-1)*repmat(y2$sdev,nrow(para1A),1))
}
if (length(pc_num)==0) pc_num=ncol(y2$rotation)
else pc_num=c(pc_num,ncol(y2$rotation))
print(pc_num)

# model training in batches
for (m in seq(1,nrow(y1b),100000)){

	if (m+99999 <= nrow(y1b)) y1a=y1b[seq(m,m+99999,1),]
	else y1a=y1b[seq(m,nrow(y1b),1),]	
	
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
}}
close(con)

# plot prediction accuracy
png(filename4,width=4,height=4,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(2,2,2,1),mgp=c(2,1,0))
plot(seq(0,19e5,1e5),acc[,1],main="KG (shaded) and KS", xlab="Number of Samples Analyzed",pch=16,ylim=c(0,1))
points(seq(0,19e5,1e5),acc[,2])
plot(seq(0,19e5,1e5),acc[,3],main="KEG", xlab="Number of Samples Analyzed",pch=16,ylim=c(0,1))
points(seq(0,19e5,1e5),acc[,4])
plot(seq(0,19e5,1e5),acc[,5],main="KEGS", xlab="Number of Samples Analyzed",pch=16,ylim=c(0,1))
points(seq(0,19e5,1e5),acc[,6])
plot(seq(0,19e5,1e5),pc_num,main="No. of PCs", xlab="Number of Samples Analyzed",pch=16)
mtext("Large-scale Predictive Modeling", outer = TRUE, cex = 1)
dev.off()