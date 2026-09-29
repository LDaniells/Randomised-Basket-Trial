#BMA Application From Liu Paper

#library(RBesT)
library(reshape2)
library(gridExtra)
library(foreach)
#library(bhmbasket)
## approxima
#library(cubature)
#library(arm)
require(LearnBayes)
#library(rjags)
library(mvtnorm)
library(xtable)
#Calculate logit transformation
#library(car)

library(rjags)
library(textmineR)
library(matrixStats)
library(parallel)
library(MASS)


#source("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/Functions.R")
source("~/Randomized Basket Trial/TryAgain/Final Code/Functions.R")

#number of simulations
nsim <- 10000
r <- 0
# Noninferiority margin
ni_mar <- 0.05

R <- 2

p0 <- c(0.60, 0.88, 0.70, 0.51, 0.67)
# day42 death rates for the treatment:
p1 <- p0 - 0.30

# B: the number of indications 
B <- length(p0)
# ind_size: the total sample size of each indication
ind_size <- c(30, 22, 22, 45, 30)
n_soc <- c(10, 7, 7, 15, 10)
n_trt <- c(20, 15, 15, 30, 20)
## ----------------------------------------------------------------
## historical control:
## ----------------------------------------------------------------
## sample size of relevant historical controls:
nh <- c(112, 89, 20, 37, 21)
## number of day42 deaths in historical controls
rh <- c(66, 78, 14, 19, 14)
## prior ess for borrowing from historical trial
prior_n <- n_soc

ph <- rh/nh
## unit information of the log odds for SOC
u <- ph*(1-ph)
## sd of the prior for log odds
sd0 <- sqrt((u*n_soc)^(-1))
m0 <- log(ph/(1-ph))

pmat0 <- rbind(gen.scenario(p0 = p0, p1 = p0))
pmat1 <- rbind(gen.scenario(p0 = p0, p1 = p1))
##non-inferiority test
pmat1[1, ] <- pmat1[1, ] + 0.05
for(i in 2:5){
  pmat1[i, i:B] <- pmat1[i, i:B] + 0.05
}
for(i in 2:5){
  pmat1[i+5, i:B] <- pmat1[i+5, i:B] + 0.05
}
for(i in 2:5){
  pmat1[i+10, i:B] <- pmat1[i+10, i:B] + 0.05
}


#S: the number of scenarios
pmat0 <- rbind(pmat0[1:6, ], c(0.6, 0.88, 0.7, 0.51, 0.67), 
               pmat0[7:11, ], c(0.6, 0.88, 0.7, 0.51, 0.67),
               c(0.6, 0.88, 0.7, 0.51, 0.67),c(0.6, 0.88, 0.7, 0.51, 0.67),c(0.6, 0.88, 0.7, 0.51, 0.67))
pmat1 <- rbind(pmat1[1:6, ], pmat1[7:11, ], 
               c(0.6, 0.88, 0.7, 0.51, 0.67),
               c(0.53, 0.76, 0.53, 0.29, 0.40),
               c(0.43, 0.66, 0.43, 0.21, 0.52),c(0.75,0.93,0.85,0.76,0.82),c(0.55,0.83,0.65,0.46,0.62))
S <- nrow(pmat0)


# array to store results from S scenarios
ref_pp <- bma_pp <- bma2D_pp <- nex_pp <- nex2D_pp <- array(NA, c(S, nsim, B))

pool_zval <- test_pool_zval <- array(NA, c(S, nsim, B))

# Matrix for saving the power of each basket in each scenario
ref_pow <- bma_pow <- bma2D_pow <- nex_pow <- nex2D_pow <- matrix(NA, nrow = S, ncol = B) 

# Informative prior for the control arm using exnex
infor_prior_ab <- cbind(prior_n * rh/nh, prior_n * (nh-rh)/nh)

# array to store bias from S scenarios
ref_bias <- pool_bias <- pool_test_bias <- bma_bias <- bma2D_bias <- nex_bias <- nex2D_bias <- array(NA, c(S, nsim, B))


Simulation <- function(s,pmat1,pmat0,B,ind_size,Tk,Basket,delta_NI,pi,nu.zeta,m.eta,nu.eta,m.new,nu.nex,a,b,yh,nh,power,pi.alpha,nsim){
  exnex_pp <- matrix(NA,ncol=B,nrow=nsim)
  for(j in 1:nsim){
    set.seed(j)
    y <- c()
    for(i in 1:sum(ind_size)){ #Generate Data for each individual i 
      if(Tk[i]==0){
        y[i] <- rbinom(1,1,pmat0[s,Basket[i]]+delta_NI)
      }else{
        y[i] <- rbinom(1,1,pmat1[s,Basket[i]])
      }
    }
    
    Nobs <- sum(ind_size)
    jags.data <- list('K'=K,'Nobs'=Nobs,'y'=y,'Basket'=Basket,'T'=Tk,'pi'=pi,'m.eta'=m.eta,'nu.eta'=nu.eta,
                      'nu.zeta'=nu.zeta,'m.nex'=m.nex,'nu.nex'=nu.nex,'a'=a,'b'=b,'yh'=yh,'nh'=nh,'power'=power,'pi.alpha'=pi.alpha)
    jags.fit <- jags.model(file='EXppNEX_Robustified.txt',data=jags.data,n.adapt=1000,n.chains=4,quiet=T)
    
    
    samples <- coda.samples(jags.fit, variable.names = c('alpha', 'theta'), n.iter = 10000, thin = 5, silent = TRUE)
    samples_mat <- as.matrix(samples)
    
    
    #Inference on Probability Scale
    for (k in 1:B) {
      a_samps <- samples_mat[, paste0('alpha[', k, ']')]
      t_samps <- samples_mat[, paste0('theta[', k, ']')]
      
      p0_mcmc <- 1 / (1 + exp(-a_samps))
      p1_mcmc <- 1 / (1 + exp(-(a_samps + t_samps)))
      
      
      NI_indicator <- (p1_mcmc < (p0_mcmc + delta_NI))
      
      exnex_pp[j, k] <- mean(NI_indicator)
    }
    
    print(j)
  }
  return(exnex_pp)
}

s <- list()
for(i in 1:16){
  s[[i]] <- i
}
SS <- ind_size
Basket <- c(rep(1,SS[1]),rep(2,SS[2]),rep(3,SS[3]),rep(4,SS[4]),rep(5,SS[5]))
Tk <- c(rep(0,10),rep(1,20),rep(0,7),rep(1,15),rep(0,7),rep(1,15),rep(0,15),rep(1,30),rep(0,10),rep(1,20))
delta_NI <- ni_mar 
K <- B

pi <- matrix(0.5,nrow=K,ncol=2)
pi.alpha <- matrix(0.5,nrow=K,ncol=2)
Nobs <- 149

m.alpha <- rep(1,K) #was 0
nu.alpha <- rep(2,K) #was 100
m.eta <- log(4) #was log(4)
nu.eta <- 1 
m.nex <- rep(0,K)
nu.nex <- rep(1,K) #was 2
nu.zeta <- 0.5

power <- 0.5
a <- b <- 1

yh <- c(66,78,14,19,14)
nh <- c(112,89,20,37,21)

nsim <- 10000



Test <- mclapply(s,Simulation,pmat1,pmat0,B,ind_size,Tk,Basket,delta_NI,pi,nu.zeta,m.eta,nu.eta,m.new,nu.nex,a,b,yh,nh,power,pi.alpha,nsim,mc.cores=17)
save(Test,file='EXppNEX_Rubustified_10000Sims.RData')


#
#Single Criteria-------
fwer1 <- function(post_prob, fwer_target = 0.10, epsilon = 0.0001, gam = NULL){
  if(is.null(gam)){
    FWER <- 0
    gam <- 1
    while(FWER<=fwer_target){
      gam <- gam - epsilon
      FWER <- mean(apply(post_prob>gam, 1, any))
    }
    gam <- gam + epsilon
    FWER <- mean(apply(post_prob>gam, 1, any))
  }else{
    FWER <- mean(apply(post_prob>gam, 1, any))
  }
  return(list("gam"=gam, "FWER"=FWER))
}
post_prob <- Test[[1]]
fwer1(post_prob,fwer_target=0.1,epsilon=0.0001,gam=NULL)


#Basket-Wise Criteria-------
quant.try <- seq(0.95,0.99,0.0001)
fwer <- c()
Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
for(j in 1:length(quant.try)){
  PostProb <- Test[[1]]
  Q <- colQuantiles(PostProb,probs=quant.try[j])
  #Q[2] <- 0.99
  reject <- c()
  hypo <- matrix(NA,nrow=10000,ncol=5)
  for(i in 1:5){
    hypo[,i] <- PostProb[,i]>Q[i]
    reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
  }
  fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/10000 #FWER Under the Null
  Q.mat[j,] <- Q
}
quant.try[which(fwer==10)]
Q.mat[which(fwer==10),]
quant.try[which.min(abs(fwer-10))]
fwer[which.min(abs(fwer-10))]
Q.mat[which.min(abs(fwer-10)),]

#
#Q <- rep(0.9125,5) #If using a single efficacy criteria
Q <- c(0.9167500, 0.8538832, 0.8703805, 0.9553777, 0.9103777)
PostProb <- Test[[16]]
reject <- c()
hypo <- matrix(NA,nrow=10000,ncol=5)
for(i in 1:5){
  hypo[,i] <- PostProb[,i]>Q[i]
  reject[i] <- 100*(length(which(PostProb[,i]>Q[i]))/10000)
}
100*(length(which(apply(hypo[,1:5],1,sum)!=0)))/10000 #FWER Under the Null
reject #Type I error/Power
mean(reject[1:5])


EXNEXR <- Test
EXNEXNR <- Test


hist(EXNEXR[[1]][,1],ylim=c(0,2000),main='Basket 1',xlab='Post Prob')
hist(EXNEXNR[[1]][,1],add=T,col=rgb(1,0,0,0.2))
legend('topleft',legend=c('EXppNEX Robust','EXppNEX Non-Robust'),fill=c('grey',rgb(1,0,0,0.5)))
hist(EXNEXR[[1]][,2],ylim=c(0,2000),main='Basket 2',xlab='Post Prob')
hist(EXNEXNR[[1]][,2],add=T,col=rgb(1,0,0,0.2))
legend('topleft',legend=c('EXppNEX Robust','EXppNEX Non-Robust'),fill=c('grey',rgb(1,0,0,0.5)))
hist(EXNEXR[[1]][,3],ylim=c(0,2000),main='Basket 3',xlab='Post Prob')
hist(EXNEXNR[[1]][,3],add=T,col=rgb(1,0,0,0.2))
legend('topleft',legend=c('EXppNEX Robust','EXppNEX Non-Robust'),fill=c('grey',rgb(1,0,0,0.5)))
