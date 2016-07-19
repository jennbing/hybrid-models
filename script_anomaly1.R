# this script uses hybrid mixture models for anomaly detection

# install packages
install.packages(stringr)
install.packages(matrixStats)

# library
library(stringr)
library(matrixStats)

# filename
filename1="pathway to rmm_anomaly.R"
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
GMM_gamma_adapto=0
GMM_gamma_adapto1=0
GMM_gamma_adapt=0
SMM_gamma_adapt=0
GMM_gamma_adapt1=0
SMM_gamma_adapt1=0

# whitened principal component analysis
y1=wpca(paraC[,1:ncol(paraC)-1],0)
y1$x=y1$x[,1:7]
y1$sdev=y1$sdev[1:7]
y1$rotation=y1$rotation[,7]

# model training in batches
for (m in seq(1,nrow(y1$x),100000)){

	if (m+99999 <= nrow(y1$x)) y1a=y1$x[seq(m,m+99999,1),]
	else y1a=y1$x[seq(m,nrow(y1$x),1),]
	
	# model training
	gamma_o=matrix(0,nrow=0,ncol=0)
	gamma_o=Gmeans_clust(y1a,1e-4)
	
	mu_min=min(abs(gamma_o$mu))
	gamma_o$w=gamma_o$w[rowSums(gamma_o$sigma<(mu_min/100))==0]
	gamma_o$mu=gamma_o$mu[rowSums(gamma_o$sigma<(mu_min/100))==0,]
	gamma_o$sigma=gamma_o$sigma[rowSums(gamma_o$sigma<(mu_min/100))==0,]
	gamma_o$w=gamma_o$w/sum(gamma_o$w)
	GMM_gamma_adapto=gamma_o
	
}	
	# predictive density
	# KS
	prob2=matrix(0,nrow=nrow(y1$x),ncol=1)
	for (i in 1:length(GMM_gamma_adapto$w)){ 
		prob2A=matrix(GMM_gamma_adapto$w[i]*student(y1$x,GMM_gamma_adapto$mu[i,],GMM_gamma_adapto$sigma[i,],1),ncol=1)
		if (sum(is.na(prob2A) | is.infinite(prob2A))==0) prob2=prob2+prob2A}

	# true positive rate
	tp2=sum(paraC[log(prob2)<quantile(log(prob2),probs=prob,na.rm=TRUE),ncol(paraC)]!=12)/((N)*prob)

	# false positive rate
	fp2=sum(paraC[log(prob2)<quantile(log(prob2),probs=prob,na.rm=TRUE),ncol(paraC)]==12)/((N)*(1-prob))

	print(tp2)
	print(fp2)
