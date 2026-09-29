library(rjags)
library(textmineR)
library(matrixStats)
library(parallel)
library(MASS)
library(coda)
library(R2WinBUGS)
library(rjags)
library(R2jags)
library(dplyr)

# p0: hypothesized response rates for SOCs
# p1: hypothesized response rates for experimental treatments under h1
gen.scenario <- function(p0, p1){
  # assume the treatment has half the hypothesized effect
  p1_low <- p0+(p1-p0)/2
  # assume the treatment is better than the hypothesized effect
  p1_high <- p1+(p1-p0)/2
  # the number of indications
  B <- length(p0)
  
  pmat1 <- pmat2 <- pmat3 <- matrix(NA, B+1, B)
  
  for(b in 1:(B+1)){
    pmat1[b, ] <- pmat2[b, ] <- pmat3[b, ] <- p0
  }
  
  for(b in 1:B){
    # create response rates same as p1
    pmat1[-1, ][b, 1:b] <- p1[1:b]
    # create response rates lower than p1
    pmat2[-1, ][b, 1:b] <- p1_low[1:b]
    # create response rates higher than p1
    pmat3[-1, ][b, 1:b] <- p1_high[1:b]
  }
  
  rbind(pmat1, pmat2[-1, ], pmat3[-1, ])
}


init.mod <- function(B, prior_indicator){
  zero_one <- list()
  #1 is promising, 0 is not promising
  for(b in 1:B) zero_one[[b]] <- prior_indicator
  mods <- expand.grid(zero_one)
  names(mods) <- LETTERS[1:B]
  mods
}

gen.prior <- function(mods, r){
  B <- ncol(mods)
  tab0 <- rbind(table(rowSums(mods==1)), (1+0:B)^r/sum((1+0:B)^r))
  prior <- tab0[2, 1+rowSums(mods==1)]/tab0[1, 1+rowSums(mods==1)]
  
  sim_mat <- diag(1, B, B)
  for(i in 1:B){
    for(j in 1:B){
      if(i==j){sim_mat[i, j] <- 1}
      if(i!=j){
        sim_mat[i, j] <- sim_mat[j, i] <- sum((mods[, i]*mods[, j]==1)*prior)
      }
    }
  }
  rownames(sim_mat) <- colnames(sim_mat) <- paste("IND", 1:B, sep="")
  
  list("tab" = tab0, "prior"=prior, "sim" = sim_mat)
}

post.den <- function(betas, b0_prior, b1_prior, x, y, log = T){
  b0 <- betas[1]
  b1 <- betas[2]
  lik <- prod(ifelse( y==1, exp(b0+b1*x)/(1+exp(b0+b1*x)), 
                      1/(1+exp(b0+b1*x)) ))
  ##prior density of intercept b0
  den0 <- b0_prior[3]*dnorm(b0, b0_prior[1], b0_prior[2])+
    (1-b0_prior[3])*dnorm(b0, 0, 100)
  ##prior density of slope b1
  den1 <- dnorm(b1, b1_prior[1], b1_prior[2])
  post <- lik*den0*den1
  if(log){
    log(post)
  }else{
    post
  }
}

## y: response to treatment
## x: treatment assignment

mlik <- function(y, x, b0_prior, b1_prior, approxi = TRUE){
  if(approxi){
    inits <- c(as.numeric(glm(y ~ 1, family = binomial)$coef), 0)
    logmp <- laplace(post.den, inits, b0_prior = b0_prior, b1_prior = b1_prior,
                     y = y, x = x)$int
    exp(logmp)
  }else{
    adaptIntegrate(post.den, b0_prior = b0_prior, b1_prior = b1_prior,
                   y = y, x = x, log = FALSE, 
                   lowerLimit=c(-100,-100),
                   upperLimit=c(100, 100))$integral
  }
}

update.mod.wt <- function(dat, b0_prior_mat, pess_prior, enth_prior, mod_prior, 
                          approxi = TRUE){
  B <- ncol(pess_prior)
  ## list all the models to average
  ## rows are models, columns are indications
  mods <- init.mod(B, prior_indicator = c(0, 1))
  ## marginal likelihood of each indication under 
  ## pessimisitc and enthusiastic priors
  mar_lik <- matrix(NA, nrow = 2, ncol = B)
  rownames(mar_lik) <- c("pess", "enth")
  
  for(b in 1:B){
    y <- dat[[b]]$outcome
    x <- dat[[b]]$trt
    
    if(approxi==FALSE){
      mar_lik[1, b] <- mlik(y = y, x = x, b0_prior = b0_prior_mat[b, ], 
                            b1_prior = pess_prior[, b], approxi = FALSE)
      
      mar_lik[2, b] <- mlik(y = y, x = x, b0_prior = b0_prior_mat[b, ], 
                            b1_prior = enth_prior[, b], approxi = FALSE)
    }
    if(approxi){
      mar_lik[1, b] <- mlik(y = y, x = x, b0_prior = b0_prior_mat[b, ], 
                            b1_prior = pess_prior[, b])
      
      mar_lik[2, b] <- mlik(y = y, x = x, b0_prior = b0_prior_mat[b, ], 
                            b1_prior = enth_prior[, b])
    }
  }
  
  p <- NULL
  for(k in 1:nrow(mods)){
    idx0 <- (mods[k, ]==0)
    idx1 <- (mods[k, ]==1)
    p[k] <- prod(ifelse(mods[k, ]==0, mar_lik[1, ], mar_lik[2, ]))*mod_prior[k]
  }
  ## calculate posterior prob of each grouping structure
  mods$p <- p/sum(p)
  
  list("mlik" = mar_lik, "mods" = mods)
}

fit.bma <- function(dat, b0_prior_mat, pess_prior, enth_prior, mod_prior, 
                    post_enth_wt, approxi = TRUE, beta_ni = NULL){
  result <- NULL
  for(b in 1:ncol(pess_prior)){
    y <- dat[[b]]$outcome
    x <- dat[[b]]$trt
    
    inits <- c(as.numeric(glm(y ~ 1, family = binomial)$coef), 0)
    
    post_fit0 <- laplace(post.den, inits, b0_prior = b0_prior_mat[b, ], 
                         b1_prior = pess_prior[, b], y = y, x = x)
    post_sim0 <- rmvnorm(5000, mean = post_fit0$mode, sigma = post_fit0$var)
    
    post_fit1 <- laplace(post.den, inits, b0_prior = b0_prior_mat[b, ], 
                         b1_prior = enth_prior[, b], y = y, x = x)
    post_sim1 <- rmvnorm(5000, mean = post_fit1$mode, sigma = post_fit1$var)
    ## the mixture of two models
    mod_ind <- rbinom(n = 5000, size = 1, prob = post_enth_wt[b])
    post_sim <- rbind(post_sim1[mod_ind==1, ], post_sim0[mod_ind==0, ])
    est <- quantile(post_sim[, 2], c(0.025, 0.5, 0.975))
    
    #Calculate estimated mortality rate
    logit_trt <- as.vector(post_sim[, 1]+post_sim[, 2])
    Diff      <- exp(logit_trt)/(1+exp(logit_trt)) - exp(post_sim[, 1])/(1+exp(post_sim[, 1]))
    est_diff  <- mean(Diff)
    
    if(is.null(beta_ni)){
      pp <- mean(post_sim[, 2]<0)
    }else{
      pp <- mean(post_sim[, 2]<beta_ni)
    }
    
    result <- rbind(result, c(est, pp, est_diff))
  }
  result <- data.frame(result)
  names(result) <- c("q2.5", "q50", "q97.5", "pp", "bias")
  result
}

ix_sel = function(sel,vnames,matrix2list=FALSE) {
  ns <- nchar(sel)
  # check for exact name matching
  match1 <- sapply(vnames, function(e)
    (substring(paste(e, sel, sep = ""), 1, ns) == sel)  &
      (nchar(e) == ns)
  )
  # check for matching matrix of vector, i.e., with "[" following sel
  match2 <- sapply(vnames, function(e) 
    (substring(paste(e, sel, sep = ""), 1, ns) == sel)  &
      (substring(e, ns + 1, ns + 1) == "[") 
    # &
    #   (!is.element(",",strsplit(e,split="")[[1]]))
  )
  if (matrix2list)
    match3 <- sapply(vnames, function(e) 
      (substring(paste(e, sel, sep = ""), 1, ns) == sel)  &
        (substring(e, ns + 1, ns + 1) == "[") &
        (is.element(",",strsplit(e,split="")[[1]]))
    )
  else
    match3 <- sapply(vnames,function(e) FALSE)
  
  out = NULL
  
  if (any(match1)) {
    out = (1:length(vnames))[match1]
    out = list(out)
    names(out) = sel
  }
  
  if (any(match2)) {
    out = (1:length(vnames))[match2]
    out = list(out)
    names(out) = sel
  }
  
  if (any(match3)) {
    sss = sapply(vnames,function(e) strsplit(e,split=",")[[1]][1])
    sss=sss[match3]
    sss = unique(sss)
    ssslab=lapply(sss,function(e) strsplit(e,split="")[[1]])
    ssslab=lapply(ssslab,function(e) paste(e[e != "["],collapse=""))
    
    out = lapply(sss,function(ee) {
      match3 <- sapply(vnames,
                       function(e)
                         substring(e,1,nchar(ee)) == ee
      )
      return(  (1:length(vnames))[match3])
    }
    )
    names(out) = ssslab
  }
  
  if (is.list(out)) {
    
    outvec = as.vector(sapply(out,function(e) vnames[e]))
    g1 = sapply(outvec,function(e) match(",",strsplit(e,"")[[1]]))
    g2 = sapply(1:length(outvec),function(e) substring(outvec[e],g1[e]+1,g1[e]+1))
    
    if (length(unique(g2))==1)
      out = as.vector(sapply(out,function(e) e))
  }
  
  return(out)
}

sim_jagswb <-  function(pars,
                        out,
                        matrix2list = TRUE) {
  allout = list()
  for (j in 1:length(pars)) {
    
    parsj = pars[j]
    out = cbind(out)
    vnames <- colnames(out)
    out.ix = ix_sel(parsj, vnames, matrix2list = matrix2list)
    
    if (is.list(out.ix))
      out.simsum = lapply(out.ix, function(e)
        out[, e])
    else
      out.simsum = out[, out.ix]
    
    allout[[j]] = out.simsum
  }
  
  names(allout) = pars
  
  return(allout)
}


sum_jagswb <-  function(pars,
                        out,
                        matrix2list = TRUE,
                        sum.cols =  c("mean", "sd", "2.5%", "50%", "97.5%")) {
  allout = list()
  
  for (j in 1:length(pars)) {
    parsj = pars[j]
    out = rbind(out)
    vnames <- rownames(out)
    out.ix = ix_sel(parsj, vnames, matrix2list = matrix2list)
    if (is.list(out.ix))
      out.simsum = lapply(out.ix, function(e)
        out[e, sum.cols])
    else
      out.simsum = rbind(out)[out.ix, sum.cols]
    
    allout[[j]] = out.simsum
    
  }
  names(allout) = pars
  
  allout = lapply(allout, function(e) {
    if (is.list(e) & length(e) == 1)
      return(e[[1]])
    else
      return(e)
  })
  
  return(allout)
}

###############################
#Binary Endpoint
###############################
Bin_rand_EXEXWB <- function(){
  
  ## EXNEX structure for treatment effect
  
  # priors for mu and tau(s)
  prec.mu.t <- pow(Prior.mu.t[2],-2)
  mu.t ~ dnorm(Prior.mu.t[1],prec.mu.t)
  
  for (j in 1:Ntau) {
    tau.prec.t[j] <- pow(Prior.tau.t[j,2],-2)
    tau.t[j] ~ dnorm(Prior.tau.t[j,1],tau.prec.t[j]) %_% I(0, )
  }
  
  # parameter model
  for (j in 1:Ngroups){
    Z.t[j] ~ dbin(p.exch.t[j],1)
    
    # theta from exch
    re.t[j] ~ dnorm(0,1)
    theta.ex.t[j] <- mu.t+re.t[j]*tau.t[tau.index[j]]
    
    # theta from non-exch
    prec.nex.t[j] <- pow(Prior.nex.t[j,2],-2)
    theta.nex.t[j]~ dnorm(Prior.nex.t[j,1], prec.nex.t[j])
    
    # select theta from latent variable
    theta.t[j] <- Z.t[j]*theta.ex.t[j]+(1-Z.t[j])*theta.nex.t[j]
  }
  
  # likelihood
  for (j in 1:Ngroups){
    y.t[j] ~ dbin(p.t[j], n.t[j])
    logit(p.t[j]) <- theta.t[j]
  }
}
#####################
#R-wrapper function 
#####################

Model_EXEX_rand <- function(y.t = NULL, n.t = NULL,
                            group       = NULL, 
                            group.nm    = NULL,
                            Prior.mu.t  = NULL,
                            Prior.tau.t = NULL,
                            p.exch.t    = NULL,
                            tau.index   = NULL,
                            Prior.nex.t = NULL,
                            pars        = c("p.t", "tau.t", "Z.t"),
                            MCMC        = c(3000, 7000, 4, 1),
                            DIC         = TRUE,
                            seed        = 19978
)
{
  #Data
  y.t <- y.t
  n.t <- n.t
  
  Ngroups <- length(unique(group))
  Ntau <- max(tau.index)
  
  #Priors
  ##EX parts
  
  Prior.mu.t = Prior.mu.t
  
  tau.index <- if(length(tau.index)==1) rep(1, Ngroups) else tau.index
  
  if(length(tau.index) != Ngroups){stop("tau.index must be same length as number of groups")}
  
  Prior.tau.t <- rbind(NULL,  Prior.tau.t)
  
  p.exch.t <- if(length(p.exch.t)==1) rep(p.exch.t, Ngroups) else p.exch.t
  
  if(length(p.exch.t) != Ngroups){stop("p.exch.t must be same length as number of groups")}
  
  
  #NEX part
  
  Prior.nex.t <- rbind(NULL, Prior.nex.t)
  
  
  Prior.nex.t <- if(nrow(Prior.nex.t) == 1) t(replicate(Ngroups, Prior.nex.t[1,])) else Prior.nex.t
  
  #JAGS call 
  model <- Bin_rand_EXEXWB 
  
  dat <- list("y.t" = y.t, "n.t" = n.t, 
              "Prior.mu.t" = Prior.mu.t, 
              "Prior.tau.t" = Prior.tau.t,
              "tau.index" = tau.index, "p.exch.t" = p.exch.t, 
              "Prior.nex.t" = Prior.nex.t, 
              "Ngroups" = Ngroups, "Ntau" = Ntau
  )
  
  initsfun = function(i){
    list(tau = rexp(Ntau, 1))
  }
  
  inits <- lapply(rep(1, MCMC[3]),  initsfun )
  
  fit.1 <- jags(
    data = dat,
    model = model,
    parameters.to.save = pars,
    n.burnin = MCMC[1],
    n.iter = MCMC[2],
    n.chains = MCMC[3],
    n.thin = MCMC[4],
    DIC = DIC,
    jags.seed = seed, 
    quiet = TRUE
  )
  
  #Collecting samples and summary
  
  sim.1 = sim_jagswb(pars,fit.1$BUGSoutput$sims.matrix)
  sum.1 = sum_jagswb(pars,fit.1$BUGSoutput$summary)
  
  # Convergence diagnostics
  #Rhat = fit.1$BUGSoutput$summary[,"Rhat"]
  #Rhat =Rhat[grep("theta", names(Rhat))]
  
  #if(max(Rhat) > 1.1)warning("At least one Rhat Statistics in > 1.1")
  
  return(list(fit=fit.1, sim=sim.1, sum=sum.1))
  
}

design.prior.bin <- function(theta, pt, lci, uci){
  fit <- get.beta.par(p = pt, q = c(lci, theta, uci),
                      show.output = FALSE, plot = FALSE)
  a <- fit[[1]]
  b <- fit[[2]]
  c(a, b)
}

fwer1 <- function(post_prob, fwer_target = 0.10, epsilon = 0.001, gam = NULL){
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


##x~beta(a,b), y~beta(c,d)
##P(X>y+delta)
beta.ineq <- function(a, b, c, d, delta)
{ 
  if (a <=0 | b <= 0 | c <=0 | d <= 0) 
    stop("paramters has to be positive")
  if (a <0.01 | b < 0.01 | c <0.01 | d < 0.01) 
    stop("paramters are to close to 0")
  if (delta>1) 
    stop("delta>1!")
  if (delta<0) 
    stop("delta<0!")
  
  integrand <- function(x) { dbeta(x, a, b)*pbeta(x-delta, c, d) }
  integrate(integrand, delta, 1, rel.tol=1e-4)$value
}

## simulate data for a randomized basket trial
sim.dat <- function(n_soc, n_trt, p_soc, p_trt){
  ## sample size of each indication
  ind_size <- n_soc + n_trt
  ## number of indications
  B <- length(p_soc)
  ## response rate (trt*indication)
  pmat <- rbind(p_soc, p_trt)
  ## sample size (trt*indication)
  rownames(pmat) <- c("soc", "trt")
  colnames(pmat) <- paste("IND-", LETTERS[1:B], sep = "")
  
  rawdat <- list()
  ## aggregated data
  y_soc <- y_trt <- NULL
  for(b in 1:B){
    ## randomize within each indication
    trt <- complete_ra(N = ind_size[b],  m = n_trt[b])
    ## response probability of each patient in indication b
    pvec <- pmat[1+trt, b]
    outcome <- rbinom(n = ind_size[b], size = 1, prob = pvec)
    rawdat[[b]] <- data.frame("trt"=trt, "outcome"=outcome)
    y_trt[b] <- sum(outcome[trt==1])
    y_soc[b] <- sum(outcome[trt==0])
  }
  
  summary_data <- data.frame(n_soc = n_soc, n_trt = n_trt, 
                             y_soc = y_soc, y_trt = y_trt)
  
  list("raw_data" = rawdat, "summary_data" = summary_data)
}

#eva.trt <- function(fit, dat, a_ctr, b_ctr, B, prior_n, 
#                    method = bma_method){
#  pp <- NULL
#  for(b in 1:B){
#    beta_par <- cbind(fit$w[, b], fit$a[, b], fit$b[, b])

#    if(method == "rob"){
#      trt_fit <- mixbeta(beta_par[1, ], beta_par[2, ], beta_par[3, ])
#    }else{
#     trt_fit <- mixbeta(beta_par[1, ], beta_par[2, ])
#   }

#    w <- c(prior_n/(a_ctr+b_ctr))[b]
#    mem_prior <- c(w, 1-w)
#    mem_fit0 <- set.mem.prior(num_study = 2, delta = 0)

#    mem_fit <- update.part.bin(x = c(a_ctr[b], dat$x[1, b]), 
#                           n = c(a_ctr[b]+b_ctr[b], dat$n[1, b]), 
#                           prior_part = mem_prior, 
#                           part = mem_fit0$part)

#ctr_prior <- mixbeta(c(0.5, 1, 1), c(0.5, a_ctr[b], b_ctr[b]))
#ctr_fit <- postmix(ctr_prior, n = dat$n[1, b], r = dat$x[1, b])
#    ctr_fit <- mixbeta(c(1, 1+mem_fit$post_sim*a_ctr[b]+dat$x[1, b], 
#                   1+mem_fit$post_sim*b_ctr[b]+dat$y[1, b]))
#    pp[b] <- pmixdiff(trt_fit, ctr_fit, 0, FALSE)
#  }
# pp
#}

no.borrow <- function(dat, ni_mar){
  pp <- rep(NA, nrow(dat))
  for(b in 1:nrow(dat)){
    yt <- dat$y_trt[b]
    yc <- dat$y_soc[b]
    nt <- dat$n_trt[b]
    nc <- dat$n_soc[b]
    
    pp[b] <- beta.ineq(a = yt + 1, b = nt-yt + 1, 
                       c = yc + 1, d = nc-yc + 1, 
                       delta = ni_mar )
  }
  pp
}

## post_prob: a vector of posterior probs
solve.gam <- function(post_prob, err1_target, epsilon = 0.001){
  err1 <- 0
  gam <- 1
  while(err1 <= err1_target){
    gam <- gam - epsilon
    err1 <- mean(post_prob>gam)
  }
  gam <- gam + epsilon
  err1 <- mean(post_prob>gam)
  return(list("gam" = gam, "err1" = err1))
}

cal.pow <- function(pp_mat, gam){
  pow <- NULL
  for(i in 1:ncol(pp_mat)){
    pow[i] <- mean(pp_mat[, i]>gam[i])
  }
  pow
}

## pp_mat: matrix of posterior probs (#simulations * #indication)
## gam: posterior thereshold of different indicatons
## true_neg: indicator of true negative (TRUE vs. False) of different indications
cal.fwer <- function(pp_mat, gam, true_neg){
  new_mat <- pp_mat[, true_neg]
  new_gam <- gam[true_neg]
  if(sum(true_neg)==0){
    FWER <- 0
  }
  if(sum(true_neg)==1){
    FWER <- mean(new_mat > new_gam)
  }
  if(sum(true_neg)>1){
    FWER <- mean(apply(new_mat, 1, function(x){any(x>new_gam)}))
  }
  FWER
}



library(RBesT)
library(randomizr)
library(ggplot2)
library(reshape2)
library(GGally)
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


#number of simulations
nsim <- 10000
r <- 0
# Noninferiority margin
ni_mar <- 0.05
## ----------------------------------------------------------------
## current trial:
## Participants will be randomized in a 2:1 ratio to an 
## experimental treatment or standard of care (SOC) based on 
## the confirmed or suspected fungal species. Below are the 
## details of comparator and sample size for each pathogen. 
##Fungal Pathogen	SOC	N (Trt/Control)
## Fusarium	L-AmB +/- VORI	20/10
## Scedosporium	VORI +/- other AF	15/7
## Lomentospora	VORI + TERB  (+/- other AF)	15/7
## Aspergillus	L-AmB +/-CAN	30/15
## Mucorales	L-AmB +/- other AF or ISA monoRx	20/10
## ----------------------------------------------------------------
# randomization ratio of trt vs. ctr
R <- 2
ind_name <- c("Fusarium", "Scedosporium", "Lomentospora", "Aspergillus", 
              "Mucorales")
# below are the assumed true death rates we simulate from:
# day42 death rates for SOC
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

xtable(data.frame(ind_name, rh, nh, rh/nh, n_soc, n_trt), 
       digits = c(rep(0,4), 2, rep(0, 2)))

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
               c(0.6, 0.88, 0.7, 0.51, 0.67))
pmat1 <- rbind(pmat1[1:6, ], pmat1[7:11, ], 
               c(0.6, 0.88, 0.7, 0.51, 0.67),
               c(0.53, 0.76, 0.53, 0.29, 0.40),
               c(0.43, 0.66, 0.43, 0.21, 0.52))
S <- nrow(pmat0)

## ----------------------------------------------------------------
## select the method for BMA
## ----------------------------------------------------------------
bma_mods <- init.mod(B, prior_indicator = c(0, 1))

bma_prior <- gen.prior(mods = bma_mods, r = r)$prior

aggregate(bma_prior, by = list(bma_mods[, 1]), sum)
aggregate(bma_prior, by = list(rowSums(bma_mods)), sum)

mat <- cbind(gen.prior(mods = bma_mods, r = 0)$prior, 
             gen.prior(mods = bma_mods, r = 1)$prior, 
             gen.prior(mods = bma_mods, r = 2)$prior)

sim0 <- gen.prior(mods = bma_mods, r = 0)$sim
sim1 <- gen.prior(mods = bma_mods, r = 1)$sim
sim2 <- gen.prior(mods = bma_mods, r = 2)$sim



# array to store results from S scenarios
ref_pp <- bma_pp <- bma2D_pp <- nex_pp <- nex2D_pp <- array(NA, c(nsim, B))

pool_zval <- test_pool_zval <- array(NA, c(S, nsim, B))

# Matrix for saving the power of each basket in each scenario
ref_pow <- bma_pow <- bma2D_pow <- nex_pow <- nex2D_pow <- matrix(NA, nrow = S, ncol = B) 

# Informative prior for the control arm using exnex
infor_prior_ab <- cbind(prior_n * rh/nh, prior_n * (nh-rh)/nh)

# array to store bias from S scenarios
ref_bias <- pool_bias <- pool_test_bias <- bma_bias <- bma2D_bias <- nex_bias <- nex2D_bias <- array(NA, c(S, nsim, B))

## Hypothesis of NI test:
## H0: p1 >= p0 + 0.05 (the D42 mortality rate is higher than p0+0.05)
## H1: p1 < p0 + 0.05
Simulation <- function(drift,nsim){
  p0 <- c(0.60+drift, 0.88+drift, 0.70+drift, 0.51+drift, 0.67+drift)
  # day42 death rates for the treatment:
  p1 <- p0 - 0.30
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
  my.list <- list()
  for(s in 1:6){
  for(i in 1:nsim){
    set.seed(i)
    sim_out <- sim.dat(n_soc = n_soc, n_trt = n_trt, 
                       p_soc = pmat0[s, ] + ni_mar, p_trt = pmat1[s, ])
    
    trial_dat <- sim_out$raw_data
    ## ----------------------------------------------------
    ## Reference method without any borrowing
    ## ----------------------------------------------------
    trial_summary <- sim_out$summary_data
    trial_summary$Disease <- c("Fusarium", "Scedosporium", "Lomentospora", "Aspergillus", "Mucorales")
    #ref_pp[s, i, ] <- 1-no.borrow(trial_summary, ni_mar = ni_mar)
    # #Calculate mortality rate est bias
    # phat0 <- (sim_out$summary_data$y_soc+1)/(sim_out$summary_data$n_soc+2)
    # phat1 <- (sim_out$summary_data$y_trt+1)/(sim_out$summary_data$n_trt+2)
    # bias.ref <- (phat1-phat0) - (pmat1[s, ]-pmat0[s, ])
    # ref_bias[s, i, ] <- bias.ref
    ## ----------------------------------------------------
    ## naive pooling
    ## ----------------------------------------------------
    #pooled <- data.frame(matrix(colSums(sim_out$summary_data), nrow = 1))
    #names(pooled) <- names(sim_out$summary_data)
    # pooled <- colSums(sim_out$summary_data)
    # phat0 <- pooled[3]/pooled[1]
    # phat1 <- pooled[4]/pooled[2]
    # se <- sqrt(phat0*(1-phat0)/pooled[1] + phat1*(1-phat1)/pooled[2])
    # pool_zval[s, i, ] <- (phat1 - phat0 - ni_mar )/se
    # #pool_pp[i, s] <- 1-no.borrow(pooled, ni_mar = ni_mar)
    # #Calculate mortality rate est bias
    # bias.pool <- (phat1-phat0) - (pmat1[s, ]-pmat0[s, ])
    # pool_bias[s, i, ] <- bias.pool
    
    ## ----------------------------------------------------
    ## if pooled result is significant, then proceed with individual tests:
    ## ----------------------------------------------------
    # x <- as.vector(sim_out$summary_data$n_trt - sim_out$summary_data$y_trt)
    # y <- as.vector(sim_out$summary_data$y_trt)
    # tab <- cbind(x, y)
    # #pval <- fisher.test(tab)$p.val
    # ##zvalue to pvalue
    # pval <- pnorm(unique(pool_zval[s, i, ]))
    # if(pval>0.1){
    #   test_pool_zval[s, i, ] <- pool_zval[s, i, ]
    #   #Calculate mortality rate est bias
    #   pool_test_bias[s, i, ] <- pool_bias[s, i, ]
    # }
    # if(pval<0.1){
    #   phat0 <- sim_out$summary_data$y_soc/sim_out$summary_data$n_soc
    #   phat1 <- sim_out$summary_data$y_trt/sim_out$summary_data$n_trt
    #   tot_n <- sim_out$summary_data$n_soc + sim_out$summary_data$n_trt
    #   tot_y <- sim_out$summary_data$y_soc + sim_out$summary_data$y_trt
    #   phat <- tot_y/tot_n
    #   se <- sqrt(phat0*(1-phat0)/sim_out$summary_data$n_soc + phat1*(1-phat1)/sim_out$summary_data$n_trt)
    #   
    #   test_pool_zval[s, i, ] <- (phat1 - phat0 - ni_mar)/se
    #   
    #   #Calculate mortality rate est bias
    #   phat0 <- sim_out$summary_data$y_soc/sim_out$summary_data$n_soc
    #   phat1 <- sim_out$summary_data$y_trt/sim_out$summary_data$n_trt
    #   bias.pool.test <- (phat1-phat0) - (pmat1[s, ]-pmat0[s, ])
    #   pool_test_bias[s, i, ] <- bias.pool.test
    # }
    
    
    # ## ----------------------------------------------------
    # ## BMA: learning across indications
    # ## ----------------------------------------------------
    # b0_prior_mat <- cbind(rep(0, B), rep(5, B), 1)
    # 
    # post_mod_wt <- update.mod.wt(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
    #                              enth_prior = enth_prior, mod_prior = bma_prior)$mods
    # 
    # post_enth_wt <- colSums(post_mod_wt[, 1:B]*post_mod_wt$p)
    # 
    # beta_ni <- log((p0+ni_mar)/(1-(p0+ni_mar)))-log(p0/(1-p0))
    # 
    # fit <- fit.bma(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
    #                enth_prior = enth_prior, mod_prior = bma_prior,
    #                post_enth_wt = post_enth_wt, beta_ni = beta_ni)
    # bma_pp[s, i, ] <- fit$pp
    # 
    # #Calculate mortality rate est bias
    # bias.bma <- fit$bias - (pmat1[s, ]-pmat0[s, ])
    # bma_bias[s, i, ] <- bias.bma 
    # #summary(glm(tmp[[3]]$outcome ~ tmp[[3]]$trt, family = binomial))$coef
    # ## ----------------------------------------------------
    # ## two dimensional learning 
    # ## ----------------------------------------------------
    # b0_prior_mat <- cbind(m0, sd0, 0.5)
    # 
    # post_mod_wt <- update.mod.wt(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
    #                              enth_prior = enth_prior, mod_prior = bma_prior)$mods
    # 
    # post_enth_wt <- colSums(post_mod_wt[, 1:B]*post_mod_wt$p)
    # 
    # fit <- fit.bma(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
    #                enth_prior = enth_prior, mod_prior = bma_prior,
    #                post_enth_wt = post_enth_wt, beta_ni = beta_ni)
    # 
    # bma2D_pp[s, i, ] <- fit$pp
    # 
    # #Calculate mortality rate est bias
    # bias.bma2D <- fit$bias - (pmat1[s, ]-pmat0[s, ])
    # bma2D_bias[s, i, ] <- bias.bma2D 
    ## ----------------------------------------------------
    ## EXNEX
    ## ----------------------------------------------------
    trial_summary$x_soc <- trial_summary$n_soc-trial_summary$y_soc
    
    nex_anal <- Model_EXEX_rand (y.t = trial_summary$y_trt,
                                 n.t = trial_summary$n_trt,
                                 group       = 1:B,
                                 group.nm    = trial_summary$Disease,
                                 Prior.mu.t  = c(log(4), 1),
                                 Prior.tau.t = c(0, 0.5),
                                 p.exch.t    = 0.5,
                                 tau.index   = 1,
                                 Prior.nex.t = cbind(rep(0, B), rep(2, B)),
                                 pars = c("p.t", "tau.t", "Z.t"),
                                 MCMC = c(5000, 10000, 1, 1),
                                 DIC =TRUE, seed = 19978)
    
    pt_mcmc <- nex_anal$sim$p.t
    
    pc_mcmc <- NULL
    for(basket in 1:B){
      pc_mcmc <- cbind(pc_mcmc, rbeta(5000, 1+trial_summary$y_soc[basket], 1+trial_summary$x_soc[basket]))
    }
    
    nex_pp[i, ] <- colMeans(pt_mcmc < pc_mcmc + ni_mar)
    
    
    ## ----------------------------------------------------
    ## EXNEX-2D
    ## ----------------------------------------------------
    pc_mcmc <- NULL
    for(basket in 1:B){
      ## beta mixture priors
      bm_prior <- mixbeta(rob = c(0.5, 1, 1), inf = c(0.5, infor_prior_ab[basket, ]))
      bm_post <- postmix(priormix = bm_prior, n = trial_summary$n_soc[basket],
                         r = trial_summary$y_soc[basket])
      
      pc_mcmc <- cbind(pc_mcmc, rmix(bm_post, 5000))
    }
    
    nex2D_pp[i, ] <- colMeans(pt_mcmc < pc_mcmc + ni_mar)
    
    print(c(i, s))
    
  }
  my.list[[s]] <- list('EXNEX'=nex_pp,'EXNEX2D'=nex2D_pp)
  }
  return(my.list)
}

drift.vec <- seq(-0.2,0.05,0.01)
drift <- list()
for(i in 1:length(drift.vec)){
  drift[[i]] <- drift.vec[i]
}

EXNEX2D_Drifts <-mclapply(drift,Simulation,nsim,mc.cores = 14)
save(EXNEX2D_Drifts,file='EXNEX2D_Drifts_10000.RData')





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

Q <- rep(0.97,5)

for(k in seq_along(drift.vec)){
  
  Results <- EXNEX2D_Drifts[[k]]
  
  for(j in 1:n_scen){
    
    PostProb <- Results[[j]]$EXNEX2D
    
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

Q <- c(0.968,0.969,0.967,0.975,0.974)

for(k in seq_along(drift.vec)){
  
  Results <- EXNEX2D_Drifts[[k]]
  
  for(j in 1:n_scen){
    
    PostProb <- Results[[j]]$EXNEX2D
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
    title="EXNEX-2D Power Across Drift by Basket and Scenario"
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
    title="EXNEX-2D Type I Error Across Drift by Basket and Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    strip.text=element_text(face="bold"),
    plot.title=element_text(hjust=0.5)
  ) +
  scale_color_manual(values=scenario_colors1)
(Power_plot / T1E_plot) / FWER_plot


