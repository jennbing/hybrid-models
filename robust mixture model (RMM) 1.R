# functions for computational study in Barkan & Averbuch (2015) are included here.
# for more information, please refer to Barkan & Averbuch (2015), 
# Robust Subspace Mixing Models for Anomaly Detection in High Dimensions.

# load R packages
library(abind)
library(pracma)
library(nortest)
library(ggm)
library(mvtnorm)
library(fitdistrplus)

# Mixture Models - model adaptation
# input:
# f - rate of loss of relevance of the old data
# n1 - total probability of the old data points
# n2 - total probability of the new data points
# N1 - total number of old data points 
# N2 - total number of new data points
# w1 - set of weightages of the old model parameters
# w2 - set of weightages of the new model parameters
# mu1 - set of means ... old ...
# mu2 - set of means ... new ...
# sigma1 - set of variance vectors ... old ...
# sigma2 - set of variance vectors ... new ...
# v1 - set of degrees of freedom of old SMM parameters
# v2 - set of degrees of freedom of new SMM parameters
# output:
# w - adapted set of weightages of the mixture models
# mu - adapted set of means ...
# sigma - adapted set of variance vectors ...
# v - adapted set of degrees of freedom of the SMM model
model_adapt <- function(f,N1,N2,w1,w2,mu1,mu2,sigma1,sigma2,v1,v2){
# data-dependent adaptation coefficient
w_adapt=w1
mu_adapt=mu1
sigma_adapt=sigma1

# truncate the parameters with low weightage
mu1=matrix(mu1[w1>0.01,],ncol=ncol(mu1))
sigma1=matrix(sigma1[w1>0.01,],ncol=ncol(sigma1))
if (length(v1)!=0) v1=matrix(v1[w1>0.01],ncol=1)
w1=matrix(w1[w1>0.01],ncol=1)
w1=w1/sum(w1)
mu2=matrix(mu2[w2>0.01,],ncol=ncol(mu2))
sigma2=matrix(sigma2[w2>0.01,],ncol=ncol(sigma2))
if (length(v1)!=0) v2=matrix(v2[w2>0.01],ncol=1)
w2=matrix(w2[w2>0.01],ncol=1)
w2=w2/sum(w2)

alpha=matrix(0,nrow=3,ncol=1)
if (length(v1)!=0){
	v_adapt=v1
	alpha=matrix(0,nrow=4,ncol=1)
}
for (i in 1:length(w1)){
	index1=seq(1,nrow(mu2),1)
	mu_diff=sqrt(rowSums((mu2-repmat(mu1[i,],nrow(mu2),1))^2))
	while (1){
		if (length(mu_diff)==0) break
		index=which.min(mu_diff)
		# F test of the equality of two covariance matrices
		test1=F_test_covariance(sigma1[i,],sigma2[index1[index],],w1[i]*N1,w2[index1[index]]*N2)
		# Hotelling's T-square test of two sample means with equal covariance matrices
		test2=Hotel_T2_test_mean(mu1[i,],mu2[index1[index],],sigma1[i,],sigma2[index1[index],],w1[i]*N1,w2[index1[index]]*N2)
		if (test1 & test2){
			index1=index1[index]
			break
		}
	 	else {
			index1=index1[-index]
	 		mu_diff=mu_diff[-index]}}	
	if (length(mu_diff)==0) next
	# adapted set of GMM parameters
	alpha[1]=(N2*w2[index1])/((N2*w2[index1])+(N1*w1[i]*f[1]))
	alpha[2]=(N2*w2[index1])/((N2*w2[index1])+(N1*w1[i]*f[2]))
	alpha[3]=(N2*w2[index1])/((N2*w2[index1])+(N1*w1[i]*f[3]))
	if (length(v1)!=0) alpha[4]=(N2*w2[index1])/((N2*w2[index1])+(N1*w1[i]*f[4]))
	w_adapt[i]=alpha[1]*w2[index1]+(1-alpha[1])*w1[i]
	mu_adapt[i,]=alpha[2]*mu2[index1,]+(1-alpha[2])*mu1[i,]
	sigma_adapt[i,]=alpha[3]*(sigma2[index1,]+mu2[index1,]^2)+(1-alpha[3])*(sigma1[i,]+mu1[i,]^2)-mu_adapt[i,]^2
	if (length(v1)!=0) v_adapt[i]=alpha[4]*v2[index1]+(1-alpha[4])*v1[i]
	w2=matrix(w2[-index1],ncol=1)
	mu2=matrix(mu2[-index1,],ncol=ncol(mu2))
	sigma2=matrix(sigma2[-index1,],ncol=ncol(sigma2))
	if (length(v1)!=0) v2=matrix(v2[-index1,],ncol=ncol(v2))
	if (length(w2)==0) break}
if (length(w2)!=0){
	w_adapt=abind(w_adapt,w2,along=1)
	mu_adapt=abind(mu_adapt,mu2,along=1)
	sigma_adapt=abind(sigma_adapt,sigma2,along=1)
	if (length(v1)!=0) v_adapt=abind(v_adapt,v2,along=1)}
# normalization
w_adapt=w_adapt/sum(w_adapt)
mu_adapt=matrix(mu_adapt[w_adapt>0.01,],ncol=ncol(mu_adapt))
sigma_adapt=matrix(sigma_adapt[w_adapt>0.01,],ncol=ncol(sigma_adapt))
if (length(v1)!=0) v_adapt=matrix(v_adapt[w_adapt>0.01],ncol=1)
w_adapt=matrix(w_adapt[w_adapt>0.01],ncol=1)
w_adapt=w_adapt/sum(w_adapt)
if (length(v1)!=0) gamma_adapt=list("w"=w_adapt,"mu"=mu_adapt,"sigma"=sigma_adapt,"v"=v_adapt)
else gamma_adapt=list("w"=w_adapt,"mu"=mu_adapt,"sigma"=sigma_adapt)
return(gamma_adapt)}
 
# soft parameter rating
# input:
# x - data points in high dimensional space
# y - data points in embedded space
# x_a - detected anomaly in high dimensional space
# y_a - corresponding detected anomaly in embedded space
# w - set of weightages of GMM
# mu - set of means ...
# sigma - set of variance vectors ... 
# method - 1: Equation 5.3 in [Barkan & Averbuch 2015]
#	     2: soft parameter rating, Equation 5.4 ... 
# output:
# xbar_a - parametric vector for the detected anomaly
par_rate <- function(x,y,x_a,y_a,w,mu,sigma,v,method){
E_q=matrix(0,nrow=1,ncol=ncol(x))
sigma_S=matrix(0,nrow=1,ncol=ncol(x))
x_a=matrix(x_a,ncol=ncol(x))
y_a=matrix(y_a,ncol=ncol(y))
xbar_a=matrix(0,nrow=nrow(x_a),ncol=ncol(x_a))
sum_prob1=matrix(0,nrow=nrow(x_a),ncol=1)
prob=matrix(0,nrow=nrow(x_a),ncol=1)
k1=matrix(0,nrow=nrow(x_a),ncol=1)
sigma_S=matrix(0,nrow=nrow(x_a),ncol=ncol(x_a))
if (method==1){
	E_q=matrix(0,nrow=nrow(x_a),ncol=ncol(x_a))
	for (k in 1:length(w)){
		if (length(v)==0) prob1=w[k]*matrix(gauss(y_a,mu[k,],sigma[k,]),ncol=1)
		if (length(v)!=0)	prob1=w[k]*matrix(student(y_a,mu[k,],sigma[k,],v[k]),ncol=1)
		k1[prob1>prob]=k
		prob[prob1>prob]=prob1[prob1>prob]
	}
	for (k in 1:length(k1)){
		k2=k1[k]
		if (length(v)==0){
			sum_prob=colSums(w[k2]*matrix(gauss(y,mu[k2,],sigma[k2,]),ncol=1),na.rm=TRUE,1)
			E_q[k,]=colSums(repmat(w[k2]*matrix(gauss(y,mu[k2,],sigma[k2,]),ncol=1),1,ncol(x))*x/sum_prob,na.rm=TRUE,1)
 		      sigma_S[k,]=colSums(repmat(w[k2]*matrix(gauss(y,mu[k2,],sigma[k2,]),ncol=1),1,ncol(x))*(x-repmat(E_q[k,],nrow(x),1))^2/sum_prob,na.rm=TRUE,1)
		}
		if (length(v)!=0){
			sum_prob=colSums(w[k2]*matrix(student(y,mu[k2,],sigma[k2,],v[k2,]),ncol=1),na.rm=TRUE,1)
			E_q[k,]=colSums(repmat(w[k2]*matrix(student(y,mu[k2,],sigma[k2,],v[k2,]),ncol=1),1,ncol(x))*x/sum_prob,na.rm=TRUE,1)
 		      sigma_S[k,]=colSums(repmat(w[k2]*matrix(student(y,mu[k2,],sigma[k2,],v[k2,]),ncol=1),1,ncol(x))*(x-repmat(E_q[k,],nrow(x),1))^2/sum_prob,na.rm=TRUE,1)
		}
	}
	xbar_a=abs(x_a-E_q)/sqrt(sigma_S)
}
if (method==2){
	for (k in 1:length(w)){
		if (length(v)==0){
			sum_prob=colSums(w[k]*matrix(gauss(y,mu[k,],sigma[k,]),ncol=1),na.rm=TRUE,1)
		      sum_prob1a=w[k]*matrix(gauss(y_a,mu[k,],sigma[k,]),ncol=1)
			E_q=colSums(repmat(w[k]*matrix(gauss(y,mu[k,],sigma[k,]),ncol=1),1,ncol(x))*x/sum_prob,na.rm=TRUE,1)
	 	      sigma_S=colSums(repmat(w[k]*matrix(gauss(y,mu[k,],sigma[k,]),ncol=1),1,ncol(x))*(x-repmat(E_q,nrow(x),1))^2/sum_prob,na.rm=TRUE,1)
		}
		if (length(v)!=0){
			sum_prob=colSums(w[k]*matrix(student(y,mu[k,],sigma[k,],v[k,]),ncol=1),na.rm=TRUE,1)
		      sum_prob1a=w[k]*matrix(student(y_a,mu[k,],sigma[k,],v[k,]),ncol=1)
			E_q=colSums(repmat(w[k]*matrix(student(y,mu[k,],sigma[k,],v[k,]),ncol=1),1,ncol(x))*x/sum_prob,na.rm=TRUE,1)
	 	      sigma_S=colSums(repmat(w[k]*matrix(student(y,mu[k,],sigma[k,],v[k,]),ncol=1),1,ncol(x))*(x-repmat(E_q,nrow(x),1))^2/sum_prob,na.rm=TRUE,1)
		}
		sum_prob1=sum_prob1a+sum_prob1
		xbar_a=xbar_a+repmat(sum_prob1a,1,ncol(x_a))*abs(x_a-repmat(E_q,nrow(x_a),1))/sqrt(repmat(sigma_S,nrow(x_a),1))
	}
	xbar_a=xbar_a/repmat(sum_prob1,1,ncol(x_a))
}
return(xbar_a)}

# G-means clustering of observed data
# [Hamerly (2004), Learning the k in k-means]
# input:
# x - set of d-dimensional vectors in the row
# alpha - confidence level of Anderson-Darling test, e.g. alpha=0.0001
# output:
# w - set of estimated weightages of GMM using G-means clustering
# mu - set of estimated means ...
# sigma - set of estimated variance vectors ...
Gmeans_clust <- function(x,alpha){
k=2
split_lev=1
split_index=1
while (split_index==1){
	split_index=0
	if (split_lev==1){
		x_centers=matrix(apply(x,2,mean),nrow=1)
		x_clust=matrix(1,nrow=nrow(x),ncol=1)
		clust=list("cluster"=x_clust,"centers"=x_centers)
		split_lev=split_lev+1}
	else {
		clust=kmeans(x,centers,iter.max=10,algorithm=c("Hartigan-Wong"), trace=FALSE)}
	centers=clust$centers
	clust_index=0
	for (i in sort(unique(clust$cluster),decreasing=FALSE)){
		if (length(x[clust$cluster==i,1])<8) next
		clust1=kmeans(x[clust$cluster==i,],2,iter.max=10,algorithm=c("Hartigan-Wong"), trace=FALSE)
		v=clust1$centers[1,]-clust1$centers[2,]
		x_v=rowSums(x[clust$cluster==i,]*repmat(v,length(x[clust$cluster==i,1]),1))/sum(v^2)
		if (ad.test(x_v)$p.value<alpha){
			centers=abind(centers,clust1$centers,along=1)
			centers=centers[-(i-clust_index),]
			clust_index=clust_index+1
			split_index=1}}}
index=1
w=matrix(0,nrow=nrow(clust$centers))
mu=matrix(0,nrow=nrow(clust$centers),ncol=ncol(clust$centers))
sigma=matrix(0,nrow=nrow(clust$centers),ncol=ncol(clust$centers))
for (i in unique(clust$cluster)){
	w[index]=length(x[clust$cluster==i,1])/length(x[,1])
	mu[index,]=apply(matrix(x[clust$cluster==i,],ncol=ncol(mu)),2,mean)
	sigma[index,]=apply(matrix(x[clust$cluster==i,],ncol=ncol(mu)),2,var)
	index=index+1}	
gamma_o=list("w"=w,"mu"=mu,"sigma"=sigma)
return(gamma_o)}

# F test of two sample covariance matrices
# [Charles Zaiontz (dated 30 Jan 2016), Box?????????????¡§?????o???????????¡§???????????????¡§?????s M Test Basic Concepts] 
# (http://www.real-statistics.com/multivariate-statistics/boxs-test-equality-covariance-matrices/boxs-test-basic-concepts/)
# input:
# sigma1,sigma2 - sample variance vectors ...
# n1,n2 - sample sizes ...
# output:
# test - logical (TRUE if the null hypothesis is accepted, i.e. the two sample means are equal, FALSE if vice versa)
F_test_covariance <- function(sigma1,sigma2,n1,n2){
n=n1+n2
m=2
sigma=((n1-1)*sigma1+(n2-1)*sigma2)/(n-m) # pooled covariance matrix
k=length(sigma)
M=(n-m)*log(prod(sigma))-((n1-1)*log(prod(sigma1))+(n2-1)*log(prod(sigma2)))
c=(2*k^2+3*k-1)/(6*(k+1)*(m-1))*((1/(n1-1))+(1/(n2-1))-(1/(n-m)))
c2=((k-1)*(k+2))/(6*(m-1))*((1/(n1-1)^2)+(1/(n2-1)^2)-(1/(n-m)^2))
df=k*(k+1)*(m-1)/2
df2=(df+2)/abs(c2-c^2)
a_plus=df/(1-c-(df/df2))
F_plus=M/a_plus
a_minus=df2/(1-c+(2/df2))
F_minus=(df2*M)/(df*(a_minus-M))
if (c2>c^2) F=F_plus
if (c2<c^2) F=F_minus
F_crit=qf(0.0001,df,df2,lower.tail=FALSE,log.p=FALSE)
test=(F<F_crit)
return(test)}

# Hotelling's T-square test of two sample means with equal covariance matrices
# [Charles Zaiontz (dated 30 Jan 2016), Hotelling's T-square Test for Two Independent Samples]
# (http://www.real-statistics.com/multivariate-statistics/hotellings-t-square-statistic/hotellings-t-square-independent-samples/)
# input:
# mu1,mu2 - sample means of two datasets
# sigma1,sigma2 - sample variance vectors ...
# n1,n2 - sample sizes ...
# output:
# test - logical (TRUE if the null hypothesis is accepted, i.e. the two sample means are equal, FALSE if vice versa)
Hotel_T2_test_mean <- function(mu1,mu2,sigma1,sigma2,n1,n2){
# pooled covariance matrix
sigma=((n1-1)*sigma1+(n2-1)*sigma2)/(n1+n2-2)
# Hotelling's T-square
t2=sum((mu1-mu2)^2/(sigma*((1/n1)+(1/n2))))
# F test
n=n1+n2-1
k=length(mu1)
f=(n-k)/(k*(n-1))*t2
f_crit=qf(0.0001,k,n-k,lower.tail=FALSE,log.p=FALSE)
test=(f<f_crit)
return(test)}

# Gaussian Mixture Model (GMM) parameters estimation using Expectation-Maximization (EM) algorithm
# based on Maximum Likelihood (ML) estimation 
# input:
# x - set of d-dimensional vectors in the row
# tol - a value indicating the magnitude below which the convergence happens. 
# (Convergence happens if their values are less than or equal to tol times the value one time step before.)
# iter - maximum number of iterations
# output:
# gamma$w - set of estimated weightages of GMM
# gamma$mu - set of estimated means of GMM
# gamma$sigma - set of estimated variance vectors of GMM
# gamma$tol - similar to tol mentioned above
GMM_ML <- function(x,tol,iter){
tol1=tol
# initial estimates of w, mu, and sigma using G-means clustering 
# [Hamerly (2004), Learning the k in k-means]
gamma_o=Gmeans_clust(x,0.0001)
w=gamma_o$w
mu=gamma_o$mu
sigma=gamma_o$sigma
mu1=matrix(mu[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0,],ncol=ncol(mu))
sigma1=matrix(sigma[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0,],ncol=ncol(sigma))
w1=matrix(w[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0],ncol=1)
mu=mu1
sigma=sigma1
w=w1
tol_w=matrix(0,nrow(w),ncol(w))
tol_mu=matrix(0,nrow(mu),ncol(mu))
tol_sigma=matrix(0,nrow(sigma),ncol(sigma))
for (t in 1:iter){
	# expectation step
	r=matrix(0,nrow=nrow(x),ncol=length(w))
	for (i in 1:length(w)){
		r[,i]=w[i]*gauss(x,mu[i,],sigma[i,])}
	r=r/repmat(matrix(rowSums(r,na.rm=TRUE),ncol=1),1,ncol(r))
	r[is.nan(r)]=0
	# maximization step
	w1=matrix(0,nrow=nrow(w))
	mu1=matrix(0,nrow=nrow(mu),ncol=ncol(mu))
	sigma1=matrix(0,nrow=nrow(sigma),ncol=ncol(sigma))
	for (i in 1:length(w)){
		w1[i]=sum(r[,i])/sum(r)
		mu1[i,]=colSums(repmat(matrix(r[,i],ncol=1),1,ncol(x))*x)/sum(r[,i])
		sigma1[i,]=colSums(repmat(matrix(r[,i],ncol=1),1,ncol(x))*(x-repmat(mu1[i,],nrow(x),1))^2)/sum(r[,i])}
	tol_w=abs(2*(w1-w)/(w+w1))
	tol_mu=abs(2*(mu1-mu)/(mu+mu1))
	tol_sigma=abs(2*(sigma1-sigma)/(sigma+sigma1))
	w=w1
	mu=mu1
	sigma=sigma1
	tol=max(matrix(c(max(tol_w),max(tol_mu),max(tol_sigma)),ncol=1),na.rm=TRUE)
	if (tol<tol1){
		break}}
print(tol)
mu1=matrix(mu[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0,],ncol=ncol(mu))
sigma1=matrix(sigma[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0,],ncol=ncol(sigma))
w1=matrix(w[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0],ncol=ncol(w))
mu=mu1
sigma=sigma1
w=w1
gamma=list("w"=w,"mu"=mu,"sigma"=sigma,"tol"=tol)
return(gamma)}

# Student t distribution Mixture Model (SMM) parameters estimation using Expectation-Maximization (EM) algorithm
# based on Maximum Likelihood (ML) estimation 
# input:
# x - set of d-dimensional vectors in the row
# tol - a value indicating the magnitude below which the convergence happens. 
# (Convergence happens if their values are less than or equal to tol times the value one time step before.)
# iter - maximum number of iterations
# output:
# gamma$w - set of estimated weightages of GMM
# gamma$mu - set of estimated means ...
# gamma$sigma - set of estimated variance vectors ...
# gamma$tol - similar to tol as mentioned above
SMM_ML <- function(x,tol,iter){
tol1=tol
# initial estimates of w, mu, and sigma using G-means clustering 
# [Hamerly (2004), Learning the k in k-means]
gamma_o=Gmeans_clust(x,0.0001)
w=gamma_o$w
mu=gamma_o$mu
sigma=gamma_o$sigma
mu1=matrix(mu[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0,],ncol=ncol(mu))
sigma1=matrix(sigma[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0,],ncol=ncol(sigma))
w1=matrix(w[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0],ncol=1)
mu=mu1
sigma=sigma1
w=w1
v=matrix(40,nrow=nrow(mu),ncol=1)
for (t in 1:iter){
	# expectation step
	r=matrix(0,nrow=nrow(x),ncol=length(w))
	y=matrix(0,nrow=nrow(x),ncol=length(w))
	for (i in 1:length(w)){
		r[,i]=w[i]*student(x,mu[i,],sigma[i,],v[i])
		y[,i]=(v[i]+ncol(x))/(v[i]+ rowSums((x-repmat(matrix(mu[i,],nrow=1),nrow(x),1))^2/repmat(matrix(sigma[i,],nrow=1),nrow(x),1)))}
	r=r/repmat(matrix(rowSums(r,na.rm=TRUE),ncol=1),1,ncol(r))
	r[is.nan(r)]=0
	# maximization step
	w1=matrix(0,nrow=nrow(w))
	mu1=matrix(0,nrow=nrow(mu),ncol=ncol(mu))
	sigma1=matrix(0,nrow=nrow(sigma),ncol=ncol(sigma))
	v1=matrix(0,nrow=nrow(v),ncol=1)
	for (i in 1:length(w)){
		w1[i]=sum(r[,i])/sum(r)
		mu1[i,]=colSums(repmat(matrix(r[,i],ncol=1),1,ncol(x))*repmat(matrix(y[,i],ncol=1),1,ncol(x))*x)/sum(r[,i]*y[,i])
		sigma1[i,]=colSums(repmat(matrix(r[,i],ncol=1),1,ncol(x))*repmat(matrix(y[,i],ncol=1),1,ncol(x))*(x-repmat(matrix(mu1[i,],nrow=1),nrow(x),1))^2)/sum(r[,i])	
		v1[i]=tryCatch(uniroot(function(z) log(z/2)-digamma(z/2)+1-log((v[i]+ncol(x))/2)+digamma((v[i]+ncol(x))/2)+sum(r[,i]*(log(y[,i])-y[,i]))/sum(r[,i]),lower=1,upper=300,tol=0.001,maxiter=100)$root, error=function(e) v[i])}
	tol_w=abs(2*(w1-w)/(w+w1))
	tol_mu=abs(2*(mu1-mu)/(mu+mu1))
	tol_sigma=abs(2*(sigma1-sigma)/(sigma+sigma1))
	tol_v=abs(2*(v1-v)/(v1+v))
	w=w1
	mu=mu1
	sigma=sigma1
	v=v1
	tol=max(matrix(c(max(tol_w),max(tol_mu),max(tol_sigma),max(tol_v)),ncol=1),na.rm=TRUE)
	if (tol<tol1) break}
print(tol)
mu1=matrix(mu[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0 & is.na(v)==0,],ncol=ncol(mu))
sigma1=matrix(sigma[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0 & is.na(v)==0,],ncol=ncol(sigma))
w1=matrix(w[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0 & is.na(v)==0],ncol=ncol(w))
v1=matrix(w[w>0.01 & is.na(w)==0 & is.na(apply(mu,1,mean))==0 & is.na(apply(sigma,1,mean))==0 & is.na(v)==0],ncol=ncol(v))
mu=mu1
sigma=sigma1
w=w1
v=v1
gamma=list("w"=w,"mu"=mu,"sigma"=sigma,"v"=v,"tol"=tol)
return(gamma)}

# dimensionality reduction using whitened principal component analysis
# input:
# x - set of d-dimensional vectors in the row
# tol - a value indicating the magnitude below which components should be omitted. 
# (Components are omitted if their standard deviations are less than or equal to tol times the standard deviation of the first component.)
# output:
# y$sdev - standard deviations of the principal components  (i.e., the square roots of the eigenvalues of the covariance/correlation matrix)
# y$u - set of eigenvectors in the columns
# y$x - projected data (whitened)
wpca <- function(x,tol){
y=prcomp(x,retx=TRUE,center=TRUE,tol=tol)
y$x=y$x/repmat(y$sdev,nrow(y$x),1)
return(y)}

# dimensionality reduction using diffusion maps
# input:
# x - set of d-dimensional vectors in the row
# output:
# x_b - biggest connected components in the embedding space 
diffMaps <- function(x){
# Euclidean distance matrix
chi=as.matrix(dist(x,"euclidean",diag=TRUE,upper=TRUE))^2
chi=matrix(chi,nrow=nrow(x))
chi1=chi
chi1[chi1==0]=NaN
# kernel scaling parameter
eps=mean(apply(chi1,2,min,na.rm=TRUE))
chi=exp(-chi/eps)
rm(chi1)
# transition probability matrix
P=chi/repmat(colSums(chi,na.rm=TRUE),nrow(chi),1)
# singular value decomposition
D=svd(P,nu=0,nv=ncol(P),LINPACK=FALSE)
lambda=D$d[2:length(D$d)]
v=D$v[,2:length(D$d)]
v=v[,lambda>(0.1*lambda[1])]
lambda=lambda[lambda>(0.1*lambda[1])]
# embedding space
y=repmat(matrix(lambda,nrow=1),nrow(v),1)*v
# adjacency matrix
y_dist=matrix(as.matrix(dist(y,"euclidean",diag=TRUE,upper=TRUE)),nrow=nrow(y))
zeta=4
del_zeta=apply(y_dist,2,sort,decreasing=FALSE)
del_zeta=colSums(del_zeta[2:zeta+1,],na.rm=TRUE)
mu_zeta=sum(del_zeta)/length(del_zeta)
sigma_zeta=sqrt(sum((del_zeta-mu_zeta)^2)/length(del_zeta))
a_ij=abs(y_dist-mu_zeta)/sigma_zeta
tao=quantile(a_ij,probs=0.6,na.rm=TRUE)
a_ij[a_ij<=tao]=0
a_ij[a_ij>tao]=1
diag(a_ij)=0
# biggest connected component
con=conComp(a_ij,method=2)
index=which.max(rowSums(repmat(con,length(unique(con)),1)==repmat(matrix(unique(con),ncol=1),1,length(con))))
x_b=x[con==unique(con)[index],]
return(x_b)}

# multivariate Gaussian distribution
# input:
# x - d-dimensional vectors in the row
# mu - mean vector
# sigma - variance vector
# output:
# p - probability density
gauss <- function(x,mu,sigma){
p=((2*pi)^(ncol(x)/2)*prod(sigma)^(0.5))^(-1)*exp(-(0.5)*rowSums((x-repmat(mu,nrow(x),1))^2/repmat(sigma,nrow(x),1)))
return(p)}

# multivariate Student t distribution
# input: 
# x - d-dimensional vectors in the row
# mu - mean vector
# sigma - variance vector
# v - degrees of freedom
# output:
# p - probability density
student <- function(x,mu,sigma,v){
p=gamma((v+ncol(x))/2)/((v*pi)^(ncol(x)/2)*prod(sigma)^(0.5)*gamma(v/2))*(rowSums((x-repmat(mu,nrow(x),1))^2/repmat(sigma,nrow(x),1))/v+1)^(-(v+ncol(x))/2)
return(p)}