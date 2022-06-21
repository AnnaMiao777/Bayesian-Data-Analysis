library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)


############Data Analysis#############
#------------------------------------#

data.re<- read_excel("/Users/annabellgao/Library/Mobile Documents/com~apple~CloudDocs/Spring 2022/8310_Bayesian Data Analysis /Final Project/Real estate valuation data set.xlsx")

colnames(data.re)[colnames(data.re) == "X1 transaction date"] <- "x1"
colnames(data.re)[colnames(data.re) == "X2 house age"] <- "x2"
colnames(data.re)[colnames(data.re) == "X3 distance to the nearest MRT station"] <- "x3"
colnames(data.re)[colnames(data.re) == "X4 number of convenience stores"] <- "x4"
colnames(data.re)[colnames(data.re) == "X5 latitude"] <- "x5"
colnames(data.re)[colnames(data.re) == "X6 longitude"] <- "x6"
colnames(data.re)[colnames(data.re) == "Y house price of unit area"] <- "y"

head(data.re) 
colSums(is.na(data.re)) #There is no missing data

#Descriptive statistics
summary(data.re) 

##Relation between 2 variables
ggplot(aes(x=x2,y=y),data=data.re)+
  geom_point()+
  ggtitle("relation between house age and house price") 

ggplot(aes(x=x3,y=y),data=data.re)+
  geom_point()+
  ggtitle("relation between distance to the nearest MRT station and house price") 

ggplot(aes(x=x4,y=y),data=data.re)+
  geom_point()+
  ggtitle("relation between number of convenience stores and house price") 

ggplot(aes(x=x5,y=y),data=data.re)+
  geom_point()+
  ggtitle("relation between latitude and house price") 

ggplot(aes(x=x6,y=y),data=data.re)+
  geom_point()+
  ggtitle("relation between longitude and house price") 

pairs(data.re[,1:4], col= "blue", pch=18, main= "Relationship between x1, x2, x3, x4")


###############Frequentist Analysis##############
##----------------------------------------------#

# Regress on Y
fre.re <- lm(y ~ x5 + x2 + x3 + x4, data=data.re)
summary(fre.re)
confint(fre.re)
# Checking fit 
par(mfrow = c(2,2)) 
plot(fre.re)  
# Create regression table
library(moderndive)
get_regression_table(fre.re,digits = 3)


# Check multicollinearity  

vif(fre.re)
#Desire VIFs < 10.

###############Bayesian Analysis##############
##--------Prior 1----------------------------#
# Inputs:
# Y = response vector
# X = design matrix
# a = prior mean for beta
# R = prior covariance matrix for beta (i.e., phi^{-1})   
# a0 = prior parameter for phi 
# b0 = prior parameter for phi, where phi~Gamma(a0,b0)
# G  = Number of Gibbs itreates
# beta = initial value of regression coefficients
# phi  = initial value of precision parameter 
Gibbs.MLR<-function(Y,X,a,R,a0=0,b0=0,G,beta,phi,verbose=TRUE){
  library(mvtnorm)
  p<-dim(X)[2]
  n<-dim(X)[1]
  beta.MCMC<-matrix(-99,nrow=G,ncol=p)
  phi.MCMC<-rep(-99,G)
  # this is for the non-informative Jeffreys prior 
  Ri<-solve(R)
  Ria<-Ri%*%a
  XTXRI<-solve(t(X)%*%X + Ri)
  XTY<-t(X)%*%Y
  for(g in 1:G){
    Cv.beta<- XTXRI/phi
    mu.beta<- XTXRI%*%(Ria+XTY)
    beta<-as.vector(rmvnorm(1,mu.beta,Cv.beta))
    a0star <- (n+p)/2 + a0
    b0star <- (sum( (Y-X%*%beta)^2 )+(beta-a)%*%Ri%*%(beta-a))/2 + b0
    phi<-rgamma(1,a0star,b0star)
    # Saving the results
    beta.MCMC[g,]<-beta
    phi.MCMC[g]<-phi
    if(verbose==TRUE){print(g)}
  }
  return(list("beta"=beta.MCMC,"phi"=phi.MCMC))
}
###################################################################################
# Generate data: Multiple Linear Regression
n<-nrow(data.re)
X<-cbind(rep(1, n), data.re$x5, data.re$x2,data.re$x3,data.re$x4)
Y<-data.re$y
a<-c(0,0,0,0,0)
R<-diag(rep(1500,5))
a0<-0
b0<-0

beta0 <- solve(t(X)%*%X)%*%t(X)%*%Y    
phi0 <- n/sum((Y-X%*%beta0)^2) 
set.seed(1)
res1<-Gibbs.MLR(Y=Y,X=X,a=a,R=R,a0=a0,b0=b0,G=5000,beta=beta0,phi=phi0,verbose=FALSE)
###################################
# Summarizing the results
(bayescoef <- apply(res1$beta[2500:5000,],2,mean))


######################
# Summarize results
plot(res1$beta[,1], col="blue")
plot(res1$beta[,2], col="blue")
plot(res1$beta[,3], col="blue")
plot(res1$beta[,4], col="blue")
plot(res1$beta[,5], col="blue")
plot(res1$phi, col="blue")


# Calculating sigma from phi
sqrt(1/mean(res1$phi[2500:5000]))
#Calculate Root MSE
yhat <- X %*% bayescoef
sqrt(sum((yhat-Y)^2)/(n-4-1))

# HPD intervals 
library(coda)
beta.mcmc = as.mcmc(res1$beta[2500:5000,])           # Coerce the vector into a MCMC object
HPDinterval(beta.mcmc, prob = 0.95)                  # Find 95% HPD interval for theta using the CODA function
phi.mcmc = as.mcmc(res1$phi[2500:5000])             # Coerce the vector into a MCMC object
HPDinterval(sqrt(1/phi.mcmc), prob = 0.95)                  # Find 95% HPD interval for sigma2 using the CODA function

###### second prior
######-------------------------------------------####
n<-nrow(data.re)
X<-cbind(rep(1, n), data.re$x5, data.re$x2,data.re$x3,data.re$x4)
Y<-data.re$y
a<-c(-5916, 239, -0.269, -0.004, 1.16)     
R<-diag(rep(750,5))
a0<-2
b0<-4

beta0 <- solve(t(X)%*%X)%*%t(X)%*%Y    
phi0  <- n/sum((Y-X%*%beta0)^2) 
res2<-Gibbs.MLR(Y=Y,X=X,a=a,R=R,a0=a0,b0=b0,G=5000,beta=beta0,phi=phi0,verbose=FALSE)
######################
# Summarize results
par(mfrow=c(3,1))
plot(res2$beta[,1], col="blue")
plot(res2$beta[,2], col="blue")
plot(res2$beta[,3], col="blue")
plot(res2$beta[,4], col="blue")
plot(res2$beta[,5], col="blue")
plot(res2$phi, col="blue")
###################################
# Summarizing the results
(bayescoef2 <- apply(res2$beta[2500:5000,],2,mean))
# Calculating sigma from phi
sqrt(1/mean(res2$phi[2500:5000]))
#Calculate Root MSE
yhat2 <- X %*% bayescoef2
sqrt(sum((yhat2-Y)^2)/(n-4-1))
#HPD intervals 
beta2.mcmc = as.mcmc(res2$beta[2500:5000,])           
HPDinterval(beta2.mcmc, prob = 0.95)                  
phi2.mcmc = as.mcmc(res2$phi[2500:5000])             
HPDinterval(sqrt(1/phi2.mcmc), prob = 0.95)               

