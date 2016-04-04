# this script validates the written algorithms by Barkan & Averbuch (2015).
# for more information, please refer to Barkan & Averbuch (2015), 
# Robust Subspace Mixing Models for Anomaly Detection in High Dimensions.

# install packages
install.packages("abind")
install.packages("pracma")
install.packages("nortest")
install.packages("ggm")
install.packages("mvtnorm")
install.packages("fitdistrplus")

# define file path
filepath1="file path to robust mixture model (RMM) 1.R"
filepath2="file path to OnlineNewsPopularity plots/simulation1.png"
filepath3="file path to OnlineNewsPopularity plots/par_rate1.png"
filepath4="file path to OnlineNewsPopularity plots/model_adapt1.png"

# functions definition
source(filepath1)

# Gaussian Mixture Model (GMM) simulation with noise
sigma1=matrix(0,nrow=10,ncol=10)
diag(sigma1)=seq(1,5.5,0.5)
sigma2=matrix(0,nrow=10,ncol=10)
diag(sigma2)=seq(1,5.5,0.5)
xa=rmvnorm(100,mean=seq(1,10,1),sigma1)
xb=rmvnorm(100,mean=seq(11,20,1),sigma2)
z=matrix(runif(1000,0,20),ncol=10)
z1=matrix(5,nrow=1,ncol=10)
z1a=matrix(15,nrow=1,ncol=10)
z1a[1]=2
x1z=abind(xa,z,xb,along=1)
x1z1=abind(xa,z1,z1a,xb,along=1)
x1a=abind(xa,xb,along=1)

# coarse filtering using Diffusion Maps (DM)
x1z_b=diffMaps(x1z)

# dimensionality reduction using Whitened Principal Component Analysis (WPCA)
y1z=wpca(x1z,0.01)
y1z1=wpca(x1z1,0.01)
y1a=wpca(x1a,0.01)
y1z_b=wpca(x1z_b,0.01)

# Gaussian Mixture Model (GMM) parameters estimation
KG_gamma=Gmeans_clust(y1z$x,0.0001)
KG_gamma1=Gmeans_clust(y1z1$x,0.0001)
KEGS_gamma1=GMM_ML(y1z1$x,0.05,100)
while (KEGS_gamma1$tol==0) KEGS_gamma1=GMM_ML(y1z1$x,0.05,100)

# Student t distribution Mixture Model (SMM) parameters estimation
KESG_gamma=SMM_ML(y1z$x,0.05,100)
while (KESG_gamma$tol==0) KESG_gamma=SMM_ML(y1z$x,0.05,100)
DKESG_gamma=SMM_ML(y1z_b$x,0.05,100)
while (DKESG_gamma$tol==0) DKESG_gamma=SMM_ML(y1z_b$x,0.05,100)
KS_gamma1=Gmeans_clust(y1z1$x,0.0001)
KESG_gamma1=SMM_ML(y1z1$x,0.05,100)
while (KESG_gamma1$tol==0) KESG_gamma1=SMM_ML(y1z1$x,0.05,100)

# probabilities of detected anomalies
p1=matrix(0,nrow=nrow(y1z$x),ncol=1)
p2=matrix(0,nrow=nrow(y1z$x),ncol=1)
p3=matrix(0,nrow=nrow(y1z$x),ncol=1)
p4=matrix(0,nrow=nrow(y1z1$x),ncol=1)
p4a=matrix(0,nrow=nrow(y1z1$x),ncol=1)
p4b=matrix(0,nrow=nrow(y1z1$x),ncol=1)
p4d=matrix(0,nrow=nrow(y1z1$x),ncol=1)
p4e=matrix(0,nrow=nrow(y1z1$x),ncol=1)
for (i in 1:length(KG_gamma$w)){ # KG
	p1=p1+matrix(KG_gamma$w[i]*gauss(y1z$x,KG_gamma$mu[i,],KG_gamma$sigma[i,]),ncol=1)}
for (i in 1:length(KESG_gamma$w)){ # KESG
	p2=p2+matrix(KESG_gamma$w[i]*student(y1z$x,KESG_gamma$mu[i,],KESG_gamma$sigma[i,],KESG_gamma$v[i]),ncol=1)}
for (i in 1:length(DKESG_gamma$w)){ # DKESG
	p3=p3+matrix(DKESG_gamma$w[i]*gauss(y1z$x,DKESG_gamma$mu[i,],DKESG_gamma$sigma[i,]),ncol=1)}
logp1=log(p1)
logp2=log(p2)
logp3=log(p3)
logp1_coff=quantile(log(p1),probs=0.1,na.rm=TRUE) # KG
logp1_coffA=quantile(log(p1),probs=0.4,na.rm=TRUE) # KG
logp2_coff=quantile(log(p2),probs=0.1,na.rm=TRUE) # KESG
logp2_coffA=quantile(log(p2),probs=0.5,na.rm=TRUE) # KESG
logp3_coff=quantile(log(p3),probs=0.1,na.rm=TRUE) # DKESG
logp3_coffA=quantile(log(p3),probs=0.4,na.rm=TRUE) # DKESG

# plotting (noise removal)
png(filepath2,width=4,height=4,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,0,0),mar=c(2,2,2,2),mgp=c(2,1,0))
plot(x1a[,1],x1a[,2],main="Data (Dim. 2 vs 1)", xlab="", ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2])))  # data points
plot(x1z[,1],x1z[,2],main="Data (with noise)", xlab="", ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2])))  # data points with white noise added
#plot(x1z_b[,1],x1z_b[,2],main="Data (DM coarse-filtered)", xlab="", ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2])))  # data points (coarse-filtered)
#plot(x1z[logp3>logp3_coff,1],x1z[logp3>logp3_coff,2],main="DKESG (low threshold)",xlab="",ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2]))) # DKESG
#plot(x1z[logp1>logp1_coff,1],x1z[logp1>logp1_coff,2],main="KESG (low threshold)",xlab="",ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2]))) # KESG
#plot(x1z[,1],x1z[,2],main="KG (low threshold)",xlab="",ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2])))	# KG
plot(x1z[logp3>logp3_coffA,1],x1z[logp3>logp3_coffA,2],main="DKESG",xlab="",ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2]))) # DKESG
#plot(x1z[logp2>logp2_coffA,1],x1z[logp2>logp2_coffA,2],main="KESG (high threshold)",xlab="",ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2]))) # KESG
plot(x1z[logp1>logp1_coffA,1],x1z[logp1>logp1_coffA,2],main="KG",xlab="",ylab="",pch=16,xlim=c(min(x1z[,1]),max(x1z[,1])),ylim=c(min(x1z[,2]),max(x1z[,2])))	# KG
dev.off()

# plotting (parameter rating)
png(filepath3,width=6,height=4,units="in",res=300)
par(mfrow=c(2,3),oma = c(0,0,0,0),mar=c(2,2,2,2),mgp=c(2,1,0))
plot(x1a[,1],x1a[,2],main="Data (Dim. 2 vs 1)", xlab="", ylab="",pch=16,xlim=c(0,20),ylim=c(0,20))  # data points
points(z1[1],z1[2],pch=8,col="red")
points(z1a[1],z1a[2],pch=8,col="red")
v=matrix(0,nrow=0,ncol=0)
xbar_a1=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KESG_gamma1$w,KESG_gamma1$mu,KESG_gamma1$sigma,v,1)
plot(xbar_a1[1,],type="p",col="red",main="Par. Rating (DKESG)",ylim=c(0,15))
points(xbar_a1[2,],type="p",col="red")
xbar_a2=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KESG_gamma1$w,KESG_gamma1$mu,KESG_gamma1$sigma,v,2)
points(xbar_a2[1,],type="l",col="blue")
points(xbar_a2[2,],type="l",col="blue",lty=2)
xbar_a1=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KG_gamma1$w,KG_gamma1$mu,KG_gamma1$sigma,v,1)
plot(xbar_a1[1,],type="p",col="red",main="Par. Rating (KG)",ylim=c(0,15))
points(xbar_a1[2,],type="p",col="red")
xbar_a2=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KG_gamma1$w,KG_gamma1$mu,KG_gamma1$sigma,v,2)
points(xbar_a2[1,],type="l",col="blue")
points(xbar_a2[2,],type="l",col="blue",lty=2)
plot(x1a[,5],x1a[,6],main="Data (Dim. 6 vs 5)", xlab="", ylab="",pch=16,xlim=c(0,20),ylim=c(0,20))  # data points
points(z1[5],z1[6],pch=8,col="red")
points(z1a[5],z1a[6],pch=8,col="red")
xbar_a1=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KEGS_gamma1$w,KEGS_gamma1$mu,KEGS_gamma1$sigma,matrix(1,nrow=length(KEGS_gamma1$w),ncol=1),1)
plot(xbar_a1[1,],type="p",col="red",main="Par. Rating (DKEGS)",ylim=c(0,15))
points(xbar_a1[2,],type="p",col="red")
xbar_a2=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KEGS_gamma1$w,KEGS_gamma1$mu,KEGS_gamma1$sigma,matrix(1,nrow=length(KEGS_gamma1$w),ncol=1),2)
points(xbar_a2[1,],type="l",col="blue")
points(xbar_a2[2,],type="l",col="blue",lty=2)
xbar_a1=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KS_gamma1$w,KS_gamma1$mu,KS_gamma1$sigma,matrix(1,nrow=length(KS_gamma1$w),ncol=1),1)
plot(xbar_a1[1,],type="p",col="red",main="Par. Rating (KS)",ylim=c(0,15))
points(xbar_a1[2,],type="p",col="red")
xbar_a2=par_rate(x1z1,y1z1$x,x1z1[101:102,],y1z1$x[101:102,],KS_gamma1$w,KS_gamma1$mu,KS_gamma1$sigma,matrix(1,nrow=length(KEGS_gamma1$w),ncol=1),2)
points(xbar_a2[1,],type="l",col="blue")
points(xbar_a2[2,],type="l",col="blue",lty=2)
dev.off()

# model adaptation
sigma1=matrix(0,nrow=2,ncol=2)
diag(sigma1)=9
sigma2=matrix(0,nrow=2,ncol=2)
diag(sigma2)=9
xa=rmvnorm(100,mean=matrix(c(1,3),nrow=1),sigma1)
xb=rmvnorm(100,mean=matrix(c(11,20),nrow=1),sigma2)
x1=abind(xa,xb,along=1)
xa=rmvnorm(100,mean=matrix(c(1.5,10.5),nrow=1),sigma1)
xb=rmvnorm(100,mean=matrix(c(11.5,20.5),nrow=1),sigma2)
x2=abind(xa,xb,along=1)
N1=nrow(x1)
N2=nrow(x2)
KS_gamma1=Gmeans_clust(x1,0.0001)
KS_gamma2=Gmeans_clust(x2,0.0001)
v=matrix(1,nrow=length(KS_gamma1$w),ncol=1)
f1=c(1,1,1,1)
f2=c(0.5,0.5,0.5,0.5)
gamma_adapt1=model_adapt(f1,N1,N2,KS_gamma1$w,KS_gamma2$w,KS_gamma1$mu,KS_gamma2$mu,KS_gamma1$sigma,KS_gamma2$sigma,v,v)
gamma_adapt2=model_adapt(f2,N1,N2,KS_gamma1$w,KS_gamma2$w,KS_gamma1$mu,KS_gamma2$mu,KS_gamma1$sigma,KS_gamma2$sigma,v,v)

# plotting (model adaptation)
prob1=matrix(0,nrow=length(seq(-5,30,1)),ncol=length(seq(-5,30,1)))
for (i in seq(-5,30,1)){
for (j in seq(-5,30,1)){	
	for (k in 1:length(gamma_adapt1$w)){
		prob1[i+6,j+6]=prob1[i+6,j+6]+gamma_adapt1$w[k]*gauss(matrix(c(i,j),nrow=1),gamma_adapt1$mu[k,],gamma_adapt1$sigma[k,])
	}}}
prob2=matrix(0,nrow=length(seq(-5,30,1)),ncol=length(seq(-5,30,1)))
for (i in seq(-5,30,1)){
for (j in seq(-5,30,1)){
	for (k in 1:length(gamma_adapt2$w)){
		prob2[i+6,j+6]=prob2[i+6,j+6]+gamma_adapt2$w[k]*gauss(matrix(c(i,j),nrow=1),gamma_adapt2$mu[k,],gamma_adapt2$sigma[k,])
	}}}
png(filepath4,width=4,height=4,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,0,0),mar=c(2,2,2,2),mgp=c(2,1,0))
plot(x1[,1],x1[,2],main="Data 1 (Dim. 2 vs 1)", xlab="", ylab="",pch=16,xlim=c(-5,20),ylim=c(-5,30))
points(1,3,pch=3,col="red")
points(11,20,pch=3,col="red")
plot(x2[,1],x2[,2],main="Data 2 (Dim. 2 vs 1)", xlab="", ylab="",pch=16,xlim=c(-5,20),ylim=c(-5,30))
points(1.5,10.5,pch=3,col="red")
points(11.5,20.5,pch=3,col="red")
contour(seq(-5,30,1),seq(-5,30,1),prob1,nlevels=7,main=bquote("PDF (KS), f=" ~ .(f1)),xlim=c(-5,20),ylim=c(-5,30))
contour(seq(-5,30,1),seq(-5,30,1),prob2,nlevels=7,main=bquote("PDF (KS), f=" ~ .(f2)),xlim=c(-5,20),ylim=c(-5,30))
dev.off()