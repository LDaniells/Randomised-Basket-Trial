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

#library(randomizr)

#source("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/Functions.R")
source("Functions.R")

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
               c(0.6, 0.88, 0.7, 0.51, 0.67))
pmat1 <- rbind(pmat1[1:6, ], pmat1[7:11, ], 
               c(0.6, 0.88, 0.7, 0.51, 0.67),
               c(0.53, 0.76, 0.53, 0.29, 0.40),
               c(0.43, 0.66, 0.43, 0.21, 0.52))
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


Simulation <- function(s,pmat1,pmat0,B,ind_size,Tk,Basket,delta_NI,pi,nexmu,nexsigma,zeta,yh,nh,nsim,n_soc,n_trt){
  exnex_pp <- matrix(NA,ncol=B,nrow=nsim)
  K <- B
  for(j in 1:nsim){
    #set.seed(j)
    # sim_out <- sim.dat(n_soc = n_soc, n_trt = n_trt, 
    #                    p_soc = pmat0[s, ] + ni_mar, p_trt = pmat1[s, ])
    # 
    # trial_dat <- sim_out$raw_data
    # y.t <- sim_out$summary_data$y_trt
    # n.t <- sim_out$summary_data$n_trt
    # y.c <- sim_out$summary_data$y_soc
    # n.c <- sim_out$summary_data$n_soc
    
    y.t <- y.c <- n.t <- n.c <- c()
    for(b in 1:B){
      y.t[b] <- rbinom(1,n_trt[b],pmat1[s,b])
      y.c[b] <- rbinom(1,n_soc[b],pmat0[s,b]+delta_NI)
      n.t[b] <- n_trt[b]
      n.c[b] <- n_soc[b]
    }
    
    Nobs <- sum(ind_size)
    
    exnex2.data <- list('K'=K,'n'=n.t,'y'=y.t,'MU'=log(4),'nexmu'=rep(nexmu,K),'nexsigma'=rep(nexsigma,K),'pi'=pi,'zeta'=zeta)
    jags.exnex2 <- jags.model(file='LiuImplementation.txt',data=exnex2.data,n.adapt=25000,n.chains=4,quiet=T)
    samples <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=20000,thin=10,silent=T)
    pt_mcmc <- as.matrix(samples)
    
    pc_mcmc <- NULL
    for(b in 1:K){
      ## beta mixture priors 
      bm_prior <- mixbeta(rob = c(0.5, 1, 1), inf = c(0.5, cbind(n.c[b] * yh[b]/nh[b], n.c[b]* (nh[b]-yh[b])/nh[b])))
      bm_post <- postmix(priormix = bm_prior, n = n.c[b], 
                         r = y.c[b])
      
      pc_mcmc <- cbind(pc_mcmc, rmix(bm_post, 8000))
    }
    
    exnex_pp[j,] <- colMeans(pt_mcmc < pc_mcmc + delta_NI)
    
    print(j)
  }
  return(exnex_pp)
}

s <- list()
for(i in 1:6){
  s[[i]] <- i
}
SS <- ind_size
Basket <- c(rep(1,SS[1]),rep(2,SS[2]),rep(3,SS[3]),rep(4,SS[4]),rep(5,SS[5]))
Tk <- c(rep(0,10),rep(1,20),rep(0,7),rep(1,15),rep(0,7),rep(1,15),rep(0,15),rep(1,30),rep(0,10),rep(1,20))
delta_NI <- ni_mar 
K <- B
pi <- matrix(0.5,nrow=K,ncol=2)
nexmu <- 0
nexsigma <- 2
zeta <- 0.5

yh <- rh
nsim <- 1000


task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=1))

#x <- Simulation(s,pmat1,pmat0,B,ind_size,Tk,Basket,delta_NI,pi,nexmu,nexsigma,zeta,nsim)
Test <- mclapply(s,Simulation,pmat1,pmat0,B,ind_size,Tk,Basket,delta_NI,pi,nexmu,nexsigma,zeta,yh,nh,nsim,n_soc,n_trt,mc.cores=7)
saveRDS(Test,file=paste0('EXNEX2D_Liu_MyWay_seeds_100Sims',task_id,'.rds'))
#save(Test,file='EXNEX_Liu.RData')

