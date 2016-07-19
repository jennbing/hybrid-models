# this script uses hybrid mixture models for anomaly detection

# install packages
install.packages(stringr)
install.packages(matrixStats)

# library
library(stringr)
library(matrixStats)

# filename
filename1="pathway to rmm_predict.R"
filename2="pathway to kddcup.data_10_percent_corrected"

# function code
source(filename1)

# data extraction
delim = ","  # or is it "\t" ?
dec = "."    # or is it "," ?
para=read.csv(filename2,sep=",",header=FALSE)
para=para[rowSums(para=="?")==0,]
para[]=lapply(para, as.numeric)
paraA=para[para[,ncol(para)]==12,]
paraB=para[para[,ncol(para)]!=12,]
paraA=data.matrix(paraA,rownames.force = NA)
paraB=data.matrix(paraB,rownames.force = NA)
rm(para)

# network attack types
# dos=c(1,7,10,15,19,21)
# u2r=c(2,8,13,17)
# r2l=c(3,4,5,9,14,20,22,23)
# probe=c(6,11,16,18)
# attack_types=dos
# paraB1=matrix(0,0,0)
# for (i in 1:length(attack_types)){
#	if (i==1) paraB1=paraB[paraB[,ncol(paraB)]==attack_types[i],]
#	else paraB1=rbind(paraB1,paraB[paraB[,ncol(paraB)]==attack_types[i],])}
# paraB=paraB1

# inject anomalies
prob=0.1
N=95000
paraA=paraA[sample(nrow(paraA)),]
paraB=paraB[sample(nrow(paraB)),]
paraC=rbind(paraA[sample(N*(1-prob)),],paraB[sample(N*prob),])
paraC=paraC[sample(nrow(paraC)),]

# parameter initialization
y1=matrix(0,nrow=0,ncol=0)
v=matrix(0,nrow=0,ncol=0)
f1=c(1,1,1,1)
N1=0
N2=0
N1o=0
N1oo=0
N1a=0
N1b=0
N1c=0
N1d=0
index=0.1
index0=0
index00=0
index1=0
index2=0
index3=0
index4=0
GMM_gamma_adapto1=0
GMM_gamma_adapt1=0
SMM_gamma_adapt1=0

# whitened principal component analysis
y1=wpca(paraC[,1:ncol(paraC)-1],0)
y1$x=y1$x[,1:7]
y1$sdev=y1$sdev[1:7]
y1$rotation=y1$rotation[,7]

# model training in batches
for (m in seq(1,nrow(y1$x),1000)){

	if (m+999 <= nrow(y1$x)) y1a=y1$x[seq(m,m+999,1),]
	else y1a=y1$x[seq(m,nrow(y1$x),1),]
	
	gamma_o1=matrix(0,nrow=0,ncol=0)
	GMM_gamma1=matrix(0,nrow=0,ncol=0)
	SMM_gamma1=matrix(0,nrow=0,ncol=0)

	# model training
	y1a_b=diffMaps(y1a) # noise removal
	if (length(y1a_b)==ncol(y1a)) next
	gamma_o1=Gmeans_clust(y1a_b,0.0001)
	GMM_gamma1=GMM_ML(y1a_b,0.05,100)
	while ((GMM_gamma1$tol==0 | is.infinite(GMM_gamma1$tol)==1) & (index<=0.2)) {
		GMM_gamma1=GMM_ML(y1a_b,index,100)
		index=index+0.1}
	index=0.1
	# while (GMM_gamma1$tol==0) GMM_gamma1=GMM_ML(y1a_b,0.1,100)
	SMM_gamma1=SMM_ML(y1a_b,0.05,100) # Student t distribution Mixture Model (SMM) parameters estimation
	while ((SMM_gamma1$tol==0 | is.infinite(SMM_gamma1$tol)==1) & (index<=0.2)) {
		SMM_gamma1=SMM_ML(y1a_b,index,100)
		index=index+0.1}
	index=0.1
	# while (SMM_gamma1$tol==0) SMM_gamma1=SMM_ML(y1a_b,0.1,100)
	
	# model adaptation
	if (index00==0 & length(gamma_o1$w)>0) {
		if (sum(gamma_o1$sigma==0)==0) {
			GMM_gamma_adapto1=gamma_o1
			index00=1
			N1o=nrow(y1a)}}
	if (GMM_gamma1$tol>=0 & GMM_gamma1$tol<=0.2 & index3==0 & length(GMM_gamma1$w)>0) {
		if (sum(GMM_gamma1$sigma==0)==0) { 
			GMM_gamma_adapt1=GMM_gamma1
			index3=1
			N1c=nrow(y1a)}}
	if (SMM_gamma1$tol>=0 & SMM_gamma1$tol<=0.5 & index4==0 & length(SMM_gamma1$w)>0) {
		if (sum(SMM_gamma1$sigma==0)==0) { 
			SMM_gamma_adapt1=SMM_gamma1
			index4=1
			N1d=nrow(y1a)}}

	N2=nrow(y1a)
	if (index00==1 & length(gamma_o1$w)>0) {
		if (sum(gamma_o1$sigma==0)==0) {
			GMM_gamma_adapto1=model_adapt(f1,N1o,N2,GMM_gamma_adapto1$w,gamma_o1$w,GMM_gamma_adapto1$mu,gamma_o1$mu,GMM_gamma_adapto1$sigma,gamma_o1$sigma,v,v)
			N1o=N1o+N2}}
	if (GMM_gamma1$tol>=0 & GMM_gamma1$tol<=0.2 & index3==1 & length(GMM_gamma1$w)>0) {
		if (sum(GMM_gamma1$sigma==0)==0) { 		
			GMM_gamma_adapt1=model_adapt(f1,N1c,N2,GMM_gamma_adapt1$w,GMM_gamma1$w,GMM_gamma_adapt1$mu,GMM_gamma1$mu,GMM_gamma_adapt1$sigma,GMM_gamma1$sigma,v,v)
			N1c=N1c+N2}}
	if (SMM_gamma1$tol>=0 & SMM_gamma1$tol<=0.5 & index4==1 & length(SMM_gamma1$w)>0) {
		if (sum(SMM_gamma1$sigma==0)==0) { 
			SMM_gamma_adapt1=model_adapt(f1,N1d,N2,SMM_gamma_adapt1$w,SMM_gamma1$w,SMM_gamma_adapt1$mu,SMM_gamma1$mu,SMM_gamma_adapt1$sigma,SMM_gamma1$sigma,SMM_gamma_adapt1$v,SMM_gamma1$v)
			N1d=N1d+N2}}
}	
	# predictive density
	# DKG
	prob7a=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(GMM_gamma_adapto1$w)){
		prob7a=prob7a+matrix(GMM_gamma_adapto1$w[i]*gauss(y1$x,GMM_gamma_adapto1$mu[i,],GMM_gamma_adapto1$sigma[i,]),ncol=1)}

	# DKS
	prob7b=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(GMM_gamma_adapto1$w)){
		prob7b=prob7b+matrix(GMM_gamma_adapto1$w[i]*student(y1$x,GMM_gamma_adapto1$mu[i,],GMM_gamma_adapto1$sigma[i,],1),ncol=1)}

	# DKEG
	prob7=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(GMM_gamma_adapt1$w)){
		prob7=prob7+matrix(GMM_gamma_adapt1$w[i]*gauss(y1$x,GMM_gamma_adapt1$mu[i,],GMM_gamma_adapt1$sigma[i,]),ncol=1)}
	
	# DKES
	prob8=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(SMM_gamma_adapt1$w)){ 
		prob8=prob8+matrix(SMM_gamma_adapt1$w[i]*student(y1$x,SMM_gamma_adapt1$mu[i,],SMM_gamma_adapt1$sigma[i,],SMM_gamma_adapt1$sigma[i]),ncol=1)}

	# DKEGS
	prob9=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(GMM_gamma_adapt1$w)){ 
		prob9=prob9+matrix(GMM_gamma_adapt1$w[i]*student(y1$x,GMM_gamma_adapt1$mu[i,],GMM_gamma_adapt1$sigma[i,],1),ncol=1)}

	# DKESG
	prob10=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(SMM_gamma_adapt1$w)){
		prob10=prob10+matrix(SMM_gamma_adapt1$w[i]*gauss(y1$x,SMM_gamma_adapt1$mu[i,],SMM_gamma_adapt1$sigma[i,]),ncol=1)}

	# true positive rate
	tp7=sum(paraC[log(prob7)<quantile(log(prob7),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)
	tp7a=sum(paraC[log(prob7a)<quantile(log(prob7a),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)
	tp7b=sum(paraC[log(prob7b)<quantile(log(prob7b),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)
	tp8=sum(paraC[log(prob8)<quantile(log(prob8),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)
	tp9=sum(paraC[log(prob9)<quantile(log(prob9),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)
	tp10=sum(paraC[log(prob10)<quantile(log(prob10),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)

	# false positive rate
	fp7=sum(paraC[log(prob7)<quantile(log(prob7),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))
	fp7a=sum(paraC[log(prob7a)<quantile(log(prob7a),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))
	fp7b=sum(paraC[log(prob7b)<quantile(log(prob7b),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))
	fp8=sum(paraC[log(prob8)<quantile(log(prob8),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))
	fp9=sum(paraC[log(prob9)<quantile(log(prob9),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))
	fp10=sum(paraC[log(prob10)<quantile(log(prob10),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))

	print(tp7a)
	print(tp7b)
	print(tp7)
	print(tp8)
	print(tp9)
	print(tp10)

	print(fp7a)
	print(fp7b)
	print(fp7)
	print(fp8)
	print(fp9)
	print(fp10)
