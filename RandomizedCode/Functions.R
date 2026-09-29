
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
    #trt <- sample(c(rep(1, n_trt[b]), rep(0, n_soc[b])))
    
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


## b0_prior: (mean, sd, weight of this component)
## b1_prior: (mean, sd)
## pess_prior: pess prior for each indication
## enth_prior: enth prior for each indication
## mod_prior: prior weight for each model of consideration

## set up all the models: 
## B: number of indications
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


