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

source("~/Randomized Basket Trial/TryAgain/Functions.R")
#source("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/Functions.R")
#number of simulations
nsim <- 10000
r <- 0
# Noninferiority margin
ni_mar <- 0.05

R <- 2
Drift <- function(drift){
  p0 <- c(0.60+drift, 0.88+drift, 0.70+drift, 0.51+drift, 0.67+drift)
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
  
  # pess_prior <- enth_prior <- matrix(NA, nrow = 2, ncol = length(p1))
  # 
  # for(b in 1:B){
  #   n11 <- ind_size[b]*(R/(R+1))*p1[b]
  #   n10 <- ind_size[b]*(R/(R+1))*(1-p1[b])
  #   n01 <- ind_size[b]*(1/(R+1))*(p0[b])
  #   n00 <- ind_size[b]*(1/(R+1))*(1-p0[b])
  #   v <- sum(c(n11, n10, n01, n00)^(-1))
  #   #pess_prior[, b] <- c(0, sqrt(v))
  #   #assume it is worse
  #   pess_prior[, b] <- c(logit(p0[b]+ni_mar)-logit(p0[b]), sqrt(v))
  #   enth_prior[, b] <- c(logit(p1[b])-logit(p0[b]), sqrt(v))
  # }
  
  pess_prior <- enth_prior <- matrix(NA, nrow = 2, ncol = length(p1))
  
  for(b in 1:B){
    # Replace p0 and p1 with historical control and promising treatment
    theta0_star <- rh[b]/nh[b]   # observed D42 mortality in historical control
    theta1_star <- theta0_star - 0.30
    
    # Keep same variance calculation structure
    n11 <- ind_size[b]*(R/(R+1))*theta1_star
    n10 <- ind_size[b]*(R/(R+1))*(1-theta1_star)
    n01 <- ind_size[b]*(1/(R+1))*(theta0_star)
    n00 <- ind_size[b]*(1/(R+1))*(1-theta0_star)
    v <- sum(c(n11, n10, n01, n00)^(-1))
    
    # Pessimistic prior
    pess_prior[, b] <- c(logit(theta0_star + ni_mar)-logit(theta0_star), sqrt(v))
    
    # Enthusiastic prior
    enth_prior[, b] <- c(logit(theta1_star)-logit(theta0_star), sqrt(v))
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
                 c(0.6, 0.88, 0.7, 0.51, 0.67))
  pmat1 <- rbind(pmat1[1:6, ], pmat1[7:11, ], 
                 c(0.6, 0.88, 0.7, 0.51, 0.67),
                 c(0.53, 0.76, 0.53, 0.29, 0.40),
                 c(0.43, 0.66, 0.43, 0.21, 0.52))
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
  
  
  for(s in 1:6){
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
      
      #beta_ni <- log((p0+ni_mar)/(1-(p0+ni_mar)))-log(p0/(1-p0))
      beta_ni <- logit(theta0_star + ni_mar) - logit(theta0_star)
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
  return(list(bma_pp,bma2D_pp))
  
  # save(bma_pp,file='BMA_V3.RData')
  # save(bma2D_pp,file='BMA-2D_V3.RData')
}


drift.vec <- seq(-0.2,0.05,0.01)
drift <- list()
for(i in 1:length(drift.vec)){
  drift[[i]] <- drift.vec[i]
}

BMA_Approaches_Drifts <- mclapply(drift,Drift,mc.cores = 14)
save(BMA_Approaches_Drifts,file='BMA_BMA2D_Drifts_10000.RData')







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

drift.vec <- seq(-0.2,0.05,0.01)

n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- rep(0.963,5)

for(k in seq_along(drift.vec)){
  
  Results <- BMA_Approaches_Drifts[[k]][[2]]
  
  for(j in 1:n_scen){
    
    PostProb <- Results[j,,]
    
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
    title="BMA-2D Power Across Drift by Basket and Scenario"
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
    title="BMA-2D Type I Error Across Drift by Basket and Scenario"
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

drift.vec <- seq(-0.2,0.05,0.01)

n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- c(0.967,0.88,0.954,0.974,0.965)

for(k in seq_along(drift.vec)){
  
 
  Results <- BMA_Approaches_Drifts[[k]][[2]]
  
  for(j in 1:n_scen){
    
    PostProb <- Results[j,,]
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
    title="BMA-2D Power Across Drift by Basket and Scenario"
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
    title="BMA-2D Type I Error Across Drift by Basket and Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  ) +
  scale_color_manual(values=scenario_colors1)
(Power_plot / T1E_plot) / FWER_plot

