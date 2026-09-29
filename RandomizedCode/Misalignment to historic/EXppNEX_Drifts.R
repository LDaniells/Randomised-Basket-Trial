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

#library(randomizr)

#source("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/Functions.R")
source("Functions.R")

#number of simulations
nsim <- 10000
r <- 0
# Noninferiority margin
ni_mar <- 0.05

R <- 2

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=1))
drift.vec <- seq(-0.25,0.05,0.01)
drift <- drift.vec[task_id]

p0 <- c(0.60+drift, 0.88+drift, 0.70+drift, 0.51+drift, 0.67+drift)
#p0 <- c(0.60, 0.88, 0.70, 0.51, 0.67)
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
for(i in 1:6){
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

Test <- mclapply(s,Simulation,pmat1,pmat0,B,ind_size,Tk,Basket,delta_NI,pi,nu.zeta,m.eta,nu.eta,m.new,nu.nex,a,b,yh,nh,power,pi.alpha,nsim,mc.cores=7)
#save(Test,file='EXppNEX_Rubustified.RData')
saveRDS(Test,file=paste0('EXppNEX_Robustified_Drifts_10000Sims_',task_id,'.rds'))





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

#NEW------
library(ggplot2)
library(patchwork)

drift.vec <- seq(-0.25,0.05,0.01)
files <- list.files(pattern="EXppNEX_*.", full.names=TRUE)

n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- rep(0.908,5)

for(k in seq_along(drift.vec)){
  
  file <- files[grep(paste0("Drifts_10000Sims_",k,".rds"), files)]
  Results <- readRDS(file)
  
  for(j in 1:n_scen){
    
    PostProb <- Results[[j]]
    
    reject <- colMeans(PostProb > Q)*100
    hypo <- PostProb > Q
    
    if(j==6){
      FWER <- NA
    } else if(j==5){
      FWER <- reject[5]
    } else{
      FWER <- mean(rowSums(hypo[,j:5])!=0)*100
    }
    
    Power_array[k,j,] <- reject
    FWER_mat[k,j] <- FWER
  }
}

idx <- drift.vec >= -0.2 & drift.vec <= 0.05

Power_array2 <- Power_array[idx,,]
drift_sub <- drift.vec[idx]
Power_array <- Power_array[idx,,]
FWER_sub <- FWER_mat[idx,]


Scenarios <- paste0("S",2:6)

Power_array[,2,2:5] <- NA
Power_array[,3,3:5] <- NA
Power_array[,4,4:5] <- NA
Power_array[,5,5]   <- NA

df_power <- expand.grid(
  Drift = drift_sub,
  Scenario = Scenarios,
  Basket = paste("Basket",1:5)
)

df_power$Power <- as.vector(Power_array[,2:6,])
Power_plot <- ggplot(df_power,
                     aes(Drift, Power, color=Scenario)) +
  
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  
  facet_wrap(~Basket, ncol=5) +
  ylim(25,100) +
  
  labs(
    x="Drift",
    y="Power",
    title="EXppNEX-Robust Power Across Drift by Basket and Scenario"
  ) +
  
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  )

df_fwer <- data.frame(
  Drift = rep(drift_sub,6),
  Scenario = rep(paste0("S",1:6),each=length(drift_sub)),
  FWER = as.vector(FWER_sub)
)

df_fwer <- subset(df_fwer, Scenario!="S6")

FWER_plot <- ggplot(df_fwer,
                    aes(Drift,FWER)) +
  
  geom_line(size=1,color="steelblue") +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=10,linetype="dotted") +
  
  facet_wrap(~Scenario,ncol=5) +
  
  labs(
    x="Drift",
    y="FWER (%)",
    title="FWER Across Drift by Scenario"
  ) +
  
  theme_bw() +
  theme(
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  )


Power_S1 <- Power_array2[,1,]
Power_S2 <- Power_array2[,2,]
Power_S3 <- Power_array2[,3,]
Power_S4 <- Power_array2[,4,]
Power_S5 <- Power_array2[,5,]

Scenarios <- c("S1","S2","S3","S4","S5")

# Apply null masking
Power_S2[,1]   <- NA
Power_S3[,1:2] <- NA
Power_S4[,1:3] <- NA
Power_S5[,1:4] <- NA

Power_list <- list(
  S1 = Power_S1,
  S2 = Power_S2,
  S3 = Power_S3,
  S4 = Power_S4,
  S5 = Power_S5
)

n_drift  <- length(drift_sub)
n_basket <- 5
n_scen   <- length(Scenarios)

Drift    <- rep(drift_sub, times = n_basket * n_scen)
Basket   <- rep(rep(1:5, each = n_drift), times = n_scen)
Scenario <- rep(Scenarios, each = n_drift * n_basket)

Power <- numeric(length(Drift))

idx_fill <- 1

for(s in Scenarios){
  
  mat <- Power_list[[s]]
  
  for(b in 1:5){
    
    Power[idx_fill:(idx_fill + n_drift - 1)] <- mat[,b]
    idx_fill <- idx_fill + n_drift
    
  }
}

df_t1e <- data.frame(
  Drift = Drift,
  Basket = factor(Basket, labels=paste("Basket",1:5)),
  Scenario = Scenario,
  T1E = Power
)
T1E_plot <- ggplot(df_t1e,
                   aes(Drift, T1E, color=Scenario)) +
  
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  
  facet_wrap(~Basket, ncol=5) +
  ylim(0,30) +
  
  labs(
    x="Drift",
    y="Type I Error (%)",
    title="EXppNEX-Robust Type I Error Across Drift by Basket and Scenario"
  ) +
  
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  )

# Define a common color palette for the 5 scenarios
scenario_levels1 <- c("S1","S2","S3","S4","S5")
scenario_colors1 <- c("S1"="#1b9e77",
                      "S2"="#d95f02",
                      "S3"="#7570b3",
                      "S4"="#e7298a",
                      "S5"="#66a61e"
)
# Define a common color palette for the 5 scenarios
scenario_levels2 <- c("S2","S3","S4","S5","S6")
scenario_colors2 <- c(
  "S2"="#d95f02",
  "S3"="#7570b3",
  "S4"="#e7298a",
  "S5"="#66a61e",
  "S6"="grey"
)

# Make Scenario a factor with the same levels
df_t1e$Scenario   <- factor(df_t1e$Scenario, levels = scenario_levels1)
df_power$Scenario <- factor(df_power$Scenario, levels = scenario_levels2)

# Power plot
Power_plot <- ggplot(df_power,
                     aes(Drift, Power, color=Scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  facet_wrap(~Basket, ncol=5) +
  ylim(25,100) +
  labs(
    x="Drift",
    y="Power",
    title="EXppNEX-Robust Power Across Drift by Basket and Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  ) +
  scale_color_manual(values=scenario_colors2)

# Type I Error plot
T1E_plot <- ggplot(df_t1e,
                   aes(Drift, T1E, color=Scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  facet_wrap(~Basket, ncol=5) +
  ylim(0,60) +
  labs(
    x="Drift",
    y="Type I Error (%)",
    title="EXppNEX-Robust Type I Error Across Drift by Basket and Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  ) +
  scale_color_manual(values=scenario_colors1)
(Power_plot / T1E_plot) / FWER_plot




#Basket Specific-----
library(ggplot2)
library(patchwork)

drift.vec <- seq(-0.25,0.05,0.01)
files <- list.files(pattern="EXppNEX_*.", full.names=TRUE)

n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- c(0.917,0.845,0.857,0.950,0.892)

for(k in seq_along(drift.vec)){
  
  file <- files[grep(paste0("Drifts_10000Sims_",k,".rds"), files)]
  Results <- readRDS(file)
  
  for(j in 1:n_scen){
    
    PostProb <- Results[[j]]
    hypo <- matrix(NA,ncol=5,nrow=10000)
    reject<- c()
    for(z in 1:5){
      reject[z] <- mean(PostProb[,z] > Q[z])*100
      hypo[,z] <- PostProb[,z] > Q[z]
    }
    if(j==6){
      FWER <- NA
    } else if(j==5){
      FWER <- reject[5]
    } else{
      FWER <- mean(rowSums(hypo[,j:5])!=0)*100
    }
    
    Power_array[k,j,] <- reject
    FWER_mat[k,j] <- FWER
  }
}

idx <- drift.vec >= -0.2 & drift.vec <= 0.05

Power_array2 <- Power_array[idx,,]
drift_sub <- drift.vec[idx]
Power_array <- Power_array[idx,,]
FWER_sub <- FWER_mat[idx,]


Scenarios <- paste0("S",2:6)

Power_array[,2,2:5] <- NA
Power_array[,3,3:5] <- NA
Power_array[,4,4:5] <- NA
Power_array[,5,5]   <- NA

df_power <- expand.grid(
  Drift = drift_sub,
  Scenario = Scenarios,
  Basket = paste("Basket",1:5)
)

df_power$Power <- as.vector(Power_array[,2:6,])
Power_plot <- ggplot(df_power,
                     aes(Drift, Power, color=Scenario)) +
  
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  
  facet_wrap(~Basket, ncol=5) +
  ylim(25,100) +
  
  labs(
    x="Drift",
    y="Power",
    title="EXppNEX-Robust Power Across Drift by Basket and Scenario"
  ) +
  
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  )

df_fwer <- data.frame(
  Drift = rep(drift_sub,6),
  Scenario = rep(paste0("S",1:6),each=length(drift_sub)),
  FWER = as.vector(FWER_sub)
)

df_fwer <- subset(df_fwer, Scenario!="S6")

FWER_plot <- ggplot(df_fwer,
                    aes(Drift,FWER)) +
  
  geom_line(size=1,color="steelblue") +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=10,linetype="dotted") +
  
  facet_wrap(~Scenario,ncol=5) +
  
  labs(
    x="Drift",
    y="FWER (%)",
    title="FWER Across Drift by Scenario"
  ) +
  
  theme_bw() +
  theme(
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  )


Power_S1 <- Power_array2[,1,]
Power_S2 <- Power_array2[,2,]
Power_S3 <- Power_array2[,3,]
Power_S4 <- Power_array2[,4,]
Power_S5 <- Power_array2[,5,]

Scenarios <- c("S1","S2","S3","S4","S5")

# Apply null masking
Power_S2[,1]   <- NA
Power_S3[,1:2] <- NA
Power_S4[,1:3] <- NA
Power_S5[,1:4] <- NA

Power_list <- list(
  S1 = Power_S1,
  S2 = Power_S2,
  S3 = Power_S3,
  S4 = Power_S4,
  S5 = Power_S5
)

n_drift  <- length(drift_sub)
n_basket <- 5
n_scen   <- length(Scenarios)

Drift    <- rep(drift_sub, times = n_basket * n_scen)
Basket   <- rep(rep(1:5, each = n_drift), times = n_scen)
Scenario <- rep(Scenarios, each = n_drift * n_basket)

Power <- numeric(length(Drift))

idx_fill <- 1

for(s in Scenarios){
  
  mat <- Power_list[[s]]
  
  for(b in 1:5){
    
    Power[idx_fill:(idx_fill + n_drift - 1)] <- mat[,b]
    idx_fill <- idx_fill + n_drift
    
  }
}

df_t1e <- data.frame(
  Drift = Drift,
  Basket = factor(Basket, labels=paste("Basket",1:5)),
  Scenario = Scenario,
  T1E = Power
)
T1E_plot <- ggplot(df_t1e,
                   aes(Drift, T1E, color=Scenario)) +
  
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  
  facet_wrap(~Basket, ncol=5) +
  ylim(0,50) +
  
  labs(
    x="Drift",
    y="Type I Error (%)",
    title="EXppNEX-Robust Type I Error Across Drift by Basket and Scenario"
  ) +
  
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  )

# Define a common color palette for the 5 scenarios
scenario_levels1 <- c("S1","S2","S3","S4","S5")
scenario_colors1 <- c("S1"="#1b9e77",
                      "S2"="#d95f02",
                      "S3"="#7570b3",
                      "S4"="#e7298a",
                      "S5"="#66a61e"
)
# Define a common color palette for the 5 scenarios
scenario_levels2 <- c("S2","S3","S4","S5","S6")
scenario_colors2 <- c(
  "S2"="#d95f02",
  "S3"="#7570b3",
  "S4"="#e7298a",
  "S5"="#66a61e",
  "S6"="grey"
)

# Make Scenario a factor with the same levels
df_t1e$Scenario   <- factor(df_t1e$Scenario, levels = scenario_levels1)
df_power$Scenario <- factor(df_power$Scenario, levels = scenario_levels2)

# Power plot
Power_plot <- ggplot(df_power,
                     aes(Drift, Power, color=Scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  facet_wrap(~Basket, ncol=5) +
  ylim(25,100) +
  labs(
    x="Drift",
    y="Power",
    title="EXppNEX-Robust Power Across Drift by Basket and Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  ) +
  scale_color_manual(values=scenario_colors2)

# Type I Error plot
T1E_plot <- ggplot(df_t1e,
                   aes(Drift, T1E, color=Scenario)) +
  geom_line(size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  facet_wrap(~Basket, ncol=5) +
  ylim(0,60) +
  labs(
    x="Drift",
    y="Type I Error (%)",
    title="EXppNEX-Robust Type I Error Across Drift by Basket and Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  ) +
  scale_color_manual(values=scenario_colors1)
(Power_plot / T1E_plot) / FWER_plot

