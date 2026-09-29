#BMA Application From Liu Paper

library(RBesT)
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

library(randomizr)
source("~/Randomized Basket Trial/TryAgain/Final Code/Functions.R")
#source("Functions.R")
#number of simulations
nsim <- 2000
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

pess_prior <- enth_prior <- matrix(NA, nrow = 2, ncol = length(p1))

for(b in 1:B){
  n11 <- ind_size[b]*(R/(R+1))*p1[b]
  n10 <- ind_size[b]*(R/(R+1))*(1-p1[b])
  n01 <- ind_size[b]*(1/(R+1))*(p0[b])
  n00 <- ind_size[b]*(1/(R+1))*(1-p0[b])
  v <- sum(c(n11, n10, n01, n00)^(-1))
  #pess_prior[, b] <- c(0, sqrt(v))
  #assume it is worse
  pess_prior[, b] <- c(logit(p0[b]+ni_mar)-logit(p0[b]), sqrt(v))
  enth_prior[, b] <- c(logit(p1[b])-logit(p0[b]), sqrt(v))
}

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

bma_mods <- init.mod(B, prior_indicator = c(0, 1))

bma_prior <- gen.prior(mods = bma_mods, r = r)$prior

# array to store results from S scenarios
ref_pp <- bma_pp <- bma2D_pp <- nex_pp <- nex2D_pp <- array(NA, c(S, nsim, B))

pool_zval <- test_pool_zval <- array(NA, c(S, nsim, B))

# Matrix for saving the power of each basket in each scenario
ref_pow <- bma_pow <- bma2D_pow <- nex_pow <- nex2D_pow <- matrix(NA, nrow = S, ncol = B) 

# Informative prior for the control arm using exnex
infor_prior_ab <- cbind(prior_n * rh/nh, prior_n * (nh-rh)/nh)

# array to store bias from S scenarios
ref_bias <- pool_bias <- pool_test_bias <- bma_bias <- bma2D_bias <- nex_bias <- nex2D_bias <- array(NA, c(S, nsim, B))


for(s in 1:16){
  for(i in 1:nsim){
    set.seed(i)
    sim_out <- sim.dat(n_soc = n_soc, n_trt = n_trt, 
                       p_soc = pmat0[s, ] + ni_mar, p_trt = pmat1[s, ])
    
    trial_dat <- sim_out$raw_data
    
    ## ----------------------------------------------------
    ## BMA: learning across indications
    ## ----------------------------------------------------
    b0_prior_mat <- cbind(rep(0, B), rep(5, B), 1)
    
    post_mod_wt <- update.mod.wt(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                                 enth_prior = enth_prior, mod_prior = bma_prior)$mods
    
    post_enth_wt <- colSums(post_mod_wt[, 1:B]*post_mod_wt$p)
    
    beta_ni <- log((p0+ni_mar)/(1-(p0+ni_mar)))-log(p0/(1-p0))
    
    fit <- fit.bma(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                   enth_prior = enth_prior, mod_prior = bma_prior,
                   post_enth_wt = post_enth_wt, beta_ni = beta_ni)
    bma_pp[s, i, ] <- fit$pp
    
    
    ## ----------------------------------------------------
    ## two dimensional learning 
    ## ----------------------------------------------------
    b0_prior_mat <- cbind(m0, sd0, 0.5)
    
    post_mod_wt <- update.mod.wt(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                                 enth_prior = enth_prior, mod_prior = bma_prior)$mods
    
    post_enth_wt <- colSums(post_mod_wt[, 1:B]*post_mod_wt$p)
    
    fit <- fit.bma(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                   enth_prior = enth_prior, mod_prior = bma_prior,
                   post_enth_wt = post_enth_wt, beta_ni = beta_ni)
    
    bma2D_pp[s, i, ] <- fit$pp
    
    
    print(c(i, s))
  }
}

save(bma_pp,file='BMA_AllScenarios_2000.RData')
save(bma2D_pp,file='BMA-2D_AllScenarios_2000.RData')




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
post_prob <- bma2D_pp[1,,]
fwer1(post_prob,fwer_target=0.1,epsilon=0.00001,gam=NULL)


Q <- rep(0.9648,5) #If using a single efficacy criteria
PostProb <- bma2D_pp[6,,]
reject <- c()
hypo <- matrix(NA,nrow=2000,ncol=5)
for(i in 1:5){
  hypo[,i] <- PostProb[,i]>Q[i]
  reject[i] <- 100*(length(which(PostProb[,i]>Q[i]))/2000)
}
100*(length(which(apply(hypo[,4:5],1,sum)!=0)))/2000 #FWER Under the Null
reject #Type I error/Power
mean(reject[1:5])


quant.try <- seq(0.95,0.99,0.0001)
fwer <- c()
Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
for(j in 1:length(quant.try)){
  PostProb <- bma_pp[1,,]
  Q <- colQuantiles(PostProb,probs=quant.try[j])
  reject <- c()
  hypo <- matrix(NA,nrow=2000,ncol=5)
  for(i in 1:5){
    hypo[,i] <- PostProb[,i]>Q[i]
    reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/2000
  }
  fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/2000 #FWER Under the Null
  Q.mat[j,] <- Q
}
fwer[which.min(abs(fwer-10))]
Q.mat[which.min(abs(fwer-10)),]

Q <- c(0.9848046, 0.9088138 ,0.9706000, 0.9912138 ,0.9840000) #If using a single efficacy criteria
PostProb <- bma_pp[6,,]
reject <- c()
hypo <- matrix(NA,nrow=2000,ncol=5)
for(i in 1:5){
  hypo[,i] <- PostProb[,i]>Q[i]
  reject[i] <- 100*(length(which(PostProb[,i]>Q[i]))/2000)
}
100*(length(which(apply(hypo[,4:5],1,sum)!=0)))/2000 #FWER Under the Null
reject #Type I error/Power
mean(reject[1:5])



par(mfrow=c(1,3))
hist(bma2D_pp[1,,][,1],ylim=c(0,2500),xlim=c(0,1),main='Basket 1',xlab='Post Prob')
hist(bma_pp[1,,][,1],add=T,col=rgb(1,0,0,0.2))
legend('topleft',legend=c('BMA','BMA-2D'),fill=c(rgb(1,0,0,0.5),'grey'))

hist(bma2D_pp[1,,][,2],ylim=c(0,2500),xlim=c(0,1),main='Basket 2',xlab='Post Prob')
hist(bma_pp[1,,][,2],add=T,col=rgb(1,0,0,0.2))
legend('topleft',legend=c('BMA','BMA-2D'),fill=c(rgb(1,0,0,0.5),'grey'))
hist(bma2D_pp[1,,][,3],xlim=c(0,1),main='Basket 3',xlab='Post Prob',ylim=c(0,2500))
hist(bma_pp[1,,][,3],add=T,col=rgb(1,0,0,0.2))
legend('topleft',legend=c('BMA','BMA-2D'),fill=c(rgb(1,0,0,0.5),'grey'))


