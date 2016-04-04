# this script computes the results in (Ong & Ng 2015), 
# Hybrid Subspace Mixing Models for Pattern Mining and Anomaly Detection 
# in High Dimensions.

# install packages
install.packages("abind")
install.packages("pracma")
install.packages("nortest")
install.packages("ggm")
install.packages("mvtnorm")
install.packages("fitdistrplus")

# define file path
filepath1="path to robust mixture model (RMM) 1.R"
filepath2="path to OnlineNewsPopularity/OnlineNewsPopularity.csv"
filepath3="path to OnlineNewsPopularity plots/news_par_rate1.png"

# functions definition
source(filepath1)

# read online news popularity .csv file
delim = ","  # or is it "\t" ?
dec = "."    # or is it "," ?
news_par=as.matrix(read.csv(filepath2,header=TRUE,sep=delim,dec=dec,stringsAsFactors=FALSE))

# dimensionality reduction using principal component analysis
news_para=matrix(as.numeric(news_par[,2:61]),ncol=60)
#news_para=matrix(as.numeric(news_par[matrix(as.numeric(news_par[,2]),ncol=1)<38,2:61]),ncol=60)
y1=wpca(news_para[,2:60],0.01)

# model adaptation
day1=min(as.numeric(news_par[,2]))
count=0
v=matrix(0,nrow=0,ncol=0)
while (1){
	if (count==0){
		data=matrix(y1$x[news_para[,1]<day1+10 & news_para[,1]>=day1,],ncol=ncol(y1$x))
		N1=nrow(data)

		# KS
		KS_gamma=Gmeans_clust(data,0.0001)
		w1a1=KS_gamma$w
		mu1a1=KS_gamma$mu
		sigma1a1=KS_gamma$sigma
		#print("KS")
		#print(KS_gamma)		
		
		# KES		
		KES_gamma=SMM_ML(data,0.05,200) # KES
		while (KES_gamma$tol==0 | is.infinite(KES_gamma$tol)==1) KES_gamma=SMM_ML(data,0.05,200)
		w1b1=KES_gamma$w
		mu1b1=KES_gamma$mu
		sigma1b1=KES_gamma$sigma
		v1b1=KES_gamma$v
		#print("KES")
		#print(KES_gamma)

		# KEGS
		KEGS_gamma=GMM_ML(data,0.05,200) # KEGS
		while (KEGS_gamma$tol==0 | is.infinite(KEGS_gamma$tol)==1) KEGS_gamma=GMM_ML(data,0.05,200)
		w1c1=KEGS_gamma$w
		mu1c1=KEGS_gamma$mu
		sigma1c1=KEGS_gamma$sigma
		#print("KEGS")
		#print(KEGS_gamma)

		# DKES
		data=diffMaps(data)
		DKES_gamma=SMM_ML(data,0.05,200) 
		while (DKES_gamma$tol==0 | is.infinite(DKES_gamma$tol)==1) DKES_gamma=SMM_ML(data,0.05,200)
		w1d1=DKES_gamma$w
		mu1d1=DKES_gamma$mu
		sigma1d1=DKES_gamma$sigma
		v1d1=DKES_gamma$v
		#print("DKES")
		#print(DKES_gamma)

		# DKEGS
		DKEGS_gamma=GMM_ML(data,0.05,200) 
		while (DKEGS_gamma$tol==0 | is.infinite(DKEGS_gamma$tol)==1) DKEGS_gamma=GMM_ML(data,0.05,200)
		w1e1=DKEGS_gamma$w
		mu1e1=DKEGS_gamma$mu
		sigma1e1=DKEGS_gamma$sigma
		v1e1=DKEGS_gamma$v
		#print("DKEGS")
		#print(DKEGS_gamma)

		day1=day1+10
		count=count+1
		next		
	}
	data=matrix(y1$x[news_para[,1]<day1+10 & news_para[,1]>=day1,],ncol=ncol(y1$x))
	N2=nrow(data)
	f1=c(1,1,1,1)

	# KS
	KS_gamma=Gmeans_clust(data,0.0001)
	gamma_adapt=model_adapt(f1,N1,N2,w1a1,KS_gamma$w,mu1a1,KS_gamma$mu,sigma1a1,KS_gamma$sigma,v,v)
	w1a1=gamma_adapt$w
	mu1a1=gamma_adapt$mu
	sigma1a1=gamma_adapt$sigma
	#print("KS")
	#print(KS_gamma)
	
	# KES
	KES_gamma=SMM_ML(data,0.05,200)
	while (KES_gamma$tol==0 | is.infinite(KES_gamma$tol)==1) KES_gamma=SMM_ML(data,0.05,200)
	gamma_adapt=model_adapt(f1,N1,N2,w1b1,KES_gamma$w,mu1b1,KES_gamma$mu,sigma1b1,KES_gamma$sigma,v1b1,KES_gamma$v)
	w1b1=gamma_adapt$w
	mu1b1=gamma_adapt$mu
	sigma1b1=gamma_adapt$sigma
	v1b1=gamma_adapt$v
	#print("KES")
	#print(KES_gamma)
	
	# KEGS
	KEGS_gamma=GMM_ML(data,0.05,200)
	while (KEGS_gamma$tol==0 | is.infinite(KEGS_gamma$tol)==1) KEGS_gamma=GMM_ML(data,0.05,200)
	gamma_adapt=model_adapt(f1,N1,N2,w1c1,KEGS_gamma$w,mu1c1,KEGS_gamma$mu,sigma1c1,KEGS_gamma$sigma,v,v)
	w1c1=gamma_adapt$w
	mu1c1=gamma_adapt$mu
	sigma1c1=gamma_adapt$sigma
	#print("KEGS")
	#print(KEGS_gamma)	

	# DKES
	data=diffMaps(data)
	DKES_gamma=SMM_ML(data,0.05,200)
	while (DKES_gamma$tol==0 | is.infinite(DKES_gamma$tol)==1) DKES_gamma=SMM_ML(data,0.05,200)
	gamma_adapt=model_adapt(f1,N1,N2,w1d1,DKES_gamma$w,mu1d1,DKES_gamma$mu,sigma1d1,DKES_gamma$sigma,v1d1,DKES_gamma$v)
	w1d1=gamma_adapt$w
	mu1d1=gamma_adapt$mu
	sigma1d1=gamma_adapt$sigma
	v1d1=gamma_adapt$v
	#print("DKES")
	#print(DKES_gamma)	

	# DKEGS
	DKEGS_gamma=GMM_ML(data,0.05,200)
	while (DKEGS_gamma$tol==0 | is.infinite(DKEGS_gamma$tol)==1) DKEGS_gamma=GMM_ML(data,0.05,200)
	gamma_adapt=model_adapt(f1,N1,N2,w1e1,DKEGS_gamma$w,mu1e1,DKEGS_gamma$mu,sigma1e1,DKEGS_gamma$sigma,v1e1,DKEGS_gamma$v)
	w1e1=gamma_adapt$w
	mu1e1=gamma_adapt$mu
	sigma1e1=gamma_adapt$sigma
	v1e1=gamma_adapt$v
	#print("DKEGS")
	#print(DKEGS_gamma)	

	N1=N1+N2
	day1=day1+10
	print(day1)
	if (day1>max(news_para[,1])) break
}

# KS, f=1
xbar_KS_f1a=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]<1400,2:60],y1$x[news_para[,60]<1400,],w1a1,mu1a1,sigma1a1,matrix(0.1,nrow=nrow(w1a1),ncol=1),2)
xbar_KS_f1b=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]>=1400,2:60],y1$x[news_para[,60]>=1400,],w1a1,mu1a1,sigma1a1,matrix(0.1,nrow=nrow(w1a1),ncol=1),2)
xbar_KS_f1a[is.infinite(xbar_KS_f1a)]=NaN
xbar_KS_f1b[is.infinite(xbar_KS_f1b)]=NaN
xbar_KS_f1a=apply(xbar_KS_f1a,2,mean,na.rm=TRUE)
xbar_KS_f1b=apply(xbar_KS_f1b,2,mean,na.rm=TRUE)

# KES, f=1
xbar_KES_f1a=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]<1400,2:60],y1$x[news_para[,60]<1400,],w1b1,mu1b1,sigma1b1,v1b1,2)
xbar_KES_f1b=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]>=1400,2:60],y1$x[news_para[,60]>=1400,],w1b1,mu1b1,sigma1b1,v1b1,2)
xbar_KES_f1a[is.infinite(xbar_KES_f1a)]=NaN
xbar_KES_f1b[is.infinite(xbar_KES_f1b)]=NaN
xbar_KES_f1a=apply(xbar_KES_f1a,2,mean,na.rm=TRUE)
xbar_KES_f1b=apply(xbar_KES_f1b,2,mean,na.rm=TRUE)

# KEGS, f=1
xbar_KEGS_f1a=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]<1400,2:60],y1$x[news_para[,60]<1400,],w1c1,mu1c1,sigma1c1,matrix(0.1,nrow=nrow(w1c1),ncol=1),2)
xbar_KEGS_f1b=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]>=1400,2:60],y1$x[news_para[,60]>=1400,],w1c1,mu1c1,sigma1c1,matrix(0.1,nrow=nrow(w1c1),ncol=1),2)
xbar_KEGS_f1a[is.infinite(xbar_KEGS_f1a)]=NaN
xbar_KEGS_f1b[is.infinite(xbar_KEGS_f1b)]=NaN
xbar_KEGS_f1a=apply(xbar_KEGS_f1a,2,mean,na.rm=TRUE)
xbar_KEGS_f1b=apply(xbar_KEGS_f1b,2,mean,na.rm=TRUE)

# DKES, f=1
xbar_DKES_f1a=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]<1400,2:60],y1$x[news_para[,60]<1400,],w1d1,mu1d1,sigma1d1,v1d1,2)
xbar_DKES_f1b=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]>=1400,2:60],y1$x[news_para[,60]>=1400,],w1d1,mu1d1,sigma1d1,v1d1,2)
xbar_DKES_f1a[is.infinite(xbar_DKES_f1a)]=NaN
xbar_DKES_f1b[is.infinite(xbar_DKES_f1b)]=NaN
xbar_DKES_f1a=apply(xbar_DKES_f1a,2,mean,na.rm=TRUE)
xbar_DKES_f1b=apply(xbar_DKES_f1b,2,mean,na.rm=TRUE)

# DKEGS, f=1
xbar_DKEGS_f1a=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]<1400,2:60],y1$x[news_para[,60]<1400,],w1e1,mu1e1,sigma1e1,matrix(0.1,nrow=nrow(w1e1),ncol=1),2)
xbar_DKEGS_f1b=par_rate(news_para[,2:60],y1$x,news_para[news_para[,60]>=1400,2:60],y1$x[news_para[,60]>=1400,],w1e1,mu1e1,sigma1e1,matrix(0.1,nrow=nrow(w1e1),ncol=1),2)
xbar_DKEGS_f1a[is.infinite(xbar_DKEGS_f1a)]=NaN
xbar_DKEGS_f1b[is.infinite(xbar_DKEGS_f1b)]=NaN
xbar_DKEGS_f1a=apply(xbar_DKEGS_f1a,2,mean,na.rm=TRUE)
xbar_DKEGS_f1b=apply(xbar_DKEGS_f1b,2,mean,na.rm=TRUE)

png(filepath3,width=4,height=6,units="in",res=300)
par(mfrow=c(3,2),oma = c(0,0,0,0),mar=c(2,2,2,2),mgp=c(2,1,0))
news_par1=as.numeric(news_par[,61])
news_par1=news_par1[news_par1<3000]
hist(news_par1,freq=FALSE,main="Article Popularity",ylab="probability density",xlab="number of shares")
logfit <- fitdist(news_par1, "lnorm")
prob=dlnorm(seq(0,3000,1),logfit$estimate[1],logfit$estimate[2],log=FALSE)
points(seq(0,3000,1),prob,type="l",col="red")
points(c(1400,1400),c(0,1),type="l",lty=2,col="blue")
ylim=c(log10(min(c(xbar_KS_f1a,xbar_KS_f1b),na.rm=TRUE)),log10(max(c(xbar_KS_f1a,xbar_KS_f1b),na.rm=TRUE)))
plot(seq(2,60,1),log10(xbar_KS_f1a),type="b",lty=1,pch=1,col="black",main=bquote("KS"),ylim=ylim)
points(seq(2,60,1),log10(xbar_KS_f1b),type="b",lty=2,pch=1,col="red")
ylim=c(log10(min(c(xbar_KES_f1a,xbar_KES_f1b),na.rm=TRUE)),log10(max(c(xbar_KES_f1a,xbar_KES_f1b),na.rm=TRUE)))
plot(seq(2,60,1),log10(xbar_KES_f1a),type="b",lty=1,pch=1,col="black",main=bquote("KES"),ylim=ylim)
points(seq(2,60,1),log10(xbar_KES_f1b),type="b",lty=2,pch=1,col="red")
ylim=c(log10(min(c(xbar_KEGS_f1a,xbar_KEGS_f1b),na.rm=TRUE)),log10(max(c(xbar_KEGS_f1a,xbar_KEGS_f1b),na.rm=TRUE)))
plot(seq(2,60,1),log10(xbar_KEGS_f1a),type="b",lty=1,pch=1,col="black",main=bquote("KEGS"),ylim=ylim)
points(seq(2,60,1),log10(xbar_KEGS_f1b),type="b",lty=2,pch=1,col="red")
ylim=c(log10(min(c(xbar_DKES_f1a,xbar_DKES_f1b),na.rm=TRUE)),log10(max(c(xbar_DKES_f1a,xbar_DKES_f1b),na.rm=TRUE)))
plot(seq(2,60,1),log10(xbar_DKES_f1a),type="b",lty=1,pch=1,col="black",main=bquote("DKES"),ylim=ylim)
points(seq(2,60,1),log10(xbar_DKES_f1b),type="b",lty=2,pch=1,col="red")
ylim=c(log10(min(c(xbar_DKEGS_f1a,xbar_DKEGS_f1b),na.rm=TRUE)),log10(max(c(xbar_DKEGS_f1a,xbar_DKEGS_f1b),na.rm=TRUE)))
plot(seq(2,60,1),log10(xbar_DKEGS_f1a),type="b",lty=1,pch=1,col="black",main=bquote("DKEGS"),ylim=ylim)
points(seq(2,60,1),log10(xbar_DKEGS_f1b),type="b",lty=2,pch=1,col="red")
dev.off()