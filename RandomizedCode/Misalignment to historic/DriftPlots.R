
#EXppNEX Robust---------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/27thMarchResults/EXppNEX_Robust_HistoricDriftsMarch26")

library(ggplot2)
library(patchwork)

drift.vec <- seq(-0.25,0.05,0.01)
files <- list.files(pattern="EXppNEX_*.", full.names=TRUE)


n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
Reject.List <- list()
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- c(0.917,0.854,0.870,0.955,0.910)

for(k in seq_along(drift.vec)){
  
  file <- files[grep(paste0("Drifts_10000Sims_",k,".rds"), files)]
  Results <- readRDS(file)
  
  Reject <- matrix(NA,nrow=6,ncol=5)
  fwer <- c()
  for(j in 1:6){
    PostProb <- Results[[j]]
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    if(j==6){
      FWER <- NA
    } else if(j==5){
      FWER <- reject[5]
    } else{
      FWER <- mean(rowSums(hypo[,j:5])!=0)*100
    }
    fwer[j] <- FWER
    Reject[j,] <- reject
  }
  Reject.List[[k]] <- Reject
  FWER_mat[k,] <- fwer
  }
  
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts")

EXppNEXR <- list(
  FWER = FWER_mat,
  Reject = Reject.List
)
save(EXppNEXR,file='EXppNEX_Robust_Drifts.RData')


#EXppNEX Non-Robust---------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/27thMarchResults/EXppNEX_NonRobust_HistoricDriftsMarch26")
library(ggplot2)
library(patchwork)

drift.vec <- seq(-0.25,0.05,0.01)
files <- list.files(pattern="EXppNEX_*.", full.names=TRUE)

n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
Reject.List <- list()
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- c(0.922,0.883,0.875,0.946,0.906)

for(k in seq_along(drift.vec)){
  
  file <- files[grep(paste0("Drifts_10000Sims_",k,".rds"), files)]
  Results <- readRDS(file)
  
  Reject <- matrix(NA,nrow=6,ncol=5)
  for(j in 1:6){
    PostProb <- Results[[j]]
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    if(j==6){
      FWER <- NA
    } else if(j==5){
      FWER <- reject[5]
    } else{
      FWER <- mean(rowSums(hypo[,j:5])!=0)*100
    }
    fwer[j] <- FWER
    Reject[j,] <- reject
  }
  Reject.List[[k]] <- Reject
  FWER_mat[k,] <- fwer
}

setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts")

EXppNEX <- list(
  FWER = FWER_mat,
  Reject = Reject.List
)
save(EXppNEX,file='EXppNEX_NonRobust_Drifts.RData')



#BMA2D---------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts")
library(ggplot2)
library(patchwork)

drift.vec <- seq(-0.2,0.05,0.01)
load("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts/BMA_BMA2D_Drifts_10000.RData")

n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
Reject.List <- list()
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- c(0.967,0.88,0.954,0.974,0.965)

for(k in seq_along(drift.vec)){
  
  Results <- BMA_Approaches_Drifts[[k]][[2]]
  
  Reject <- matrix(NA,nrow=6,ncol=5)
  fwer <- c()
  for(j in 1:6){
    PostProb <- Results[j,,]
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    if(j==6){
      FWER <- NA
    } else if(j==5){
      FWER <- reject[5]
    } else{
      FWER <- mean(rowSums(hypo[,j:5])!=0)*100
    }
    fwer[j] <- FWER
    Reject[j,] <- reject
  }
  Reject.List[[k]] <- Reject
  FWER_mat[k,] <- fwer
}

setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts")

BMA2D <- list(
  FWER = FWER_mat,
  Reject = Reject.List
)
save(BMA2D,file='BMA2D_Drifts.RData')



#EXNEX2D---------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts")
library(ggplot2)
library(patchwork)

drift.vec <- seq(-0.2,0.05,0.01)
load("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts/EXNEX2D_Drifts_10000.RData")
n_drift <- length(drift.vec)
n_scen  <- 6
n_basket <- 5

Power_array <- array(NA, dim=c(n_drift,n_scen,n_basket))
Reject.List <- list()
FWER_mat <- matrix(NA,n_drift,n_scen)

Q <- c(0.968,0.969,0.967,0.975,0.974)

for(k in seq_along(drift.vec)){
  
  Results <- EXNEX2D_Drifts[[k]]
  
  Reject <- matrix(NA,nrow=6,ncol=5)
  fwer <- c()
  for(j in 1:6){
    PostProb <- Results[[j]]$EXNEX2D
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    if(j==6){
      FWER <- NA
    } else if(j==5){
      FWER <- reject[5]
    } else{
      FWER <- mean(rowSums(hypo[,j:5])!=0)*100
    }
    fwer[j] <- FWER
    Reject[j,] <- reject
  }
  Reject.List[[k]] <- Reject
  FWER_mat[k,] <- fwer
}

setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Paper_Plots/Drifts")

EXNEX2D <- list(
  FWER = FWER_mat,
  Reject = Reject.List
)
save(EXNEX2D,file='EXNEX2D_Drifts.RData')



MeanOC <- function(Data){
  MP <- c()
  ME <- c()
  for(i in 1:6){
    if(i==1|i==15){
      ME[i] <- mean(Data[i,])
      MP[i] <- 0
    }else if(i==2|i==3|i==4|i==5){
      ME[i] <- mean(Data[i,i:5])
      MP[i] <- mean(Data[i,1:(i-1)])
    }else if(i==6){
      ME[i] <- 0
      MP[i] <- mean(Data[i,])
    }
  }
  return(list(
    TypeIError = ME,
    Power = MP
  ))
}

computeOC <- function(RejectMat){
  # RejectMat: 6 (scenario) x 5 (basket) matrix of % rejection
  MP <- ME <- rep(NA, 6)
  for(i in 1:6){
    if(i == 1){
      ME[i] <- mean(RejectMat[i, ])
      MP[i] <- NA
    } else if(i %in% 2:5){
      ME[i] <- mean(RejectMat[i, i:5])
      MP[i] <- mean(RejectMat[i, 1:(i-1)])
    } else if(i == 6){
      ME[i] <- NA
      MP[i] <- mean(RejectMat[i, ])
    }
  }
  data.frame(Scenario = 1:6, Power = MP, TypeIError = ME)
}


library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# --- Load each method's saved results ---
load("EXppNEX_Robust_Drifts.RData")   # EXppNEXR: list(FWER, Reject)
load("EXppNEX_NonRobust_Drifts.RData")# EXppNEX
load("BMA2D_Drifts.RData")            # BMA2D
load("EXNEX2D_Drifts.RData")          # EXNEX2D

drift.vec.long  <- seq(-0.25, 0.05, 0.01)  # EXppNEX / EXppNEXR
drift.vec.short <- seq(-0.20, 0.05, 0.01)  # BMA2D / EXNEX2D

build_long <- function(method_name, obj, drift.vec){
  map2_dfr(seq_along(drift.vec), drift.vec, function(k, d){
    oc <- computeOC(obj$Reject[[k]])
    oc$FWER  <- obj$FWER[k, ]
    oc$Drift <- d
    oc$Method <- method_name
    oc
  })
}

plot_data <- bind_rows(
  build_long("EXNEX-2D",         EXNEX2D, drift.vec.short),
  build_long("EXppNEX",         EXppNEX, drift.vec.long),
  build_long("EXppNEX Robust",  EXppNEXR, drift.vec.long),
  build_long("BMA-2D",           BMA2D,   drift.vec.short)
) %>%
  filter(Drift >= -0.2)   # <-- restrict to the common range, -0.2 to 0.05


plot_data$Method <- factor(
  plot_data$Method,
  levels = c("EXNEX-2D", "EXppNEX", "EXppNEX Robust", "BMA-2D")
)

scenario_labels <- c(
  "1" = "Global Null", "2" = "Scenario 1", "3" = "Scenario 2",
  "4" = "Scenario 3",  "5" = "Scenario 4", "6" = "Global Alternative"
)
plot_data$Scenario <- factor(plot_data$Scenario, labels = scenario_labels[as.character(1:6)])

# Reshape to long format across the three metrics for faceting by row
plot_long <- plot_data %>%
  pivot_longer(cols = c(Power, TypeIError, FWER), names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(Value))

plot_long$Metric <- factor(
  plot_long$Metric,
  levels = c("Power", "TypeIError", "FWER"),
  labels = c("Mean Power (%)", "Mean Type I Error (%)", "FWER (%)")
)


hline_fwer <- data.frame(
  Metric = factor("FWER (%)", levels = levels(plot_long$Metric)),
  yint = 10
)

ggplot(plot_long, aes(x = Drift, y = Value, colour = Method)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(
    data = hline_fwer,
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "black"
  ) +
  facet_grid(Metric ~ Scenario) +
  labs(x = expression(Delta~"(drift)"), y = NULL, colour = "Method") +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  scale_colour_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 100))


#No Averaging
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# --- Build fully raw (unaveraged) rejection-rate data, now tagged with Metric ---
NullMat <- matrix(FALSE, nrow = 6, ncol = 5)
NullMat[1, ]    <- TRUE
NullMat[2, 2:5] <- TRUE
NullMat[3, 3:5] <- TRUE
NullMat[4, 4:5] <- TRUE
NullMat[5, 5]   <- TRUE
# Scenario 6 (Global Alternative): all baskets non-null (remain FALSE)

build_raw <- function(method_name, obj, drift.vec){
  map2_dfr(seq_along(drift.vec), drift.vec, function(k, d){
    RejectMat <- obj$Reject[[k]]  # 6 (scenario) x 5 (basket)
    df <- as.data.frame(RejectMat)
    names(df) <- paste0("Basket", 1:5)
    df$Scenario <- 1:6
    df$Drift <- d
    df$Method <- method_name
    df
  })
}



raw_data <- bind_rows(
  build_raw("EXNEX-2D",         EXNEX2D, drift.vec.short),
  build_raw("EXppNEX",         EXppNEX, drift.vec.long),
  build_raw("EXppNEX Robust",  EXppNEXR, drift.vec.long),
  build_raw("BMA-2D",           BMA2D,   drift.vec.short)
) %>%
  filter(Drift >= -0.2) %>%
  pivot_longer(cols = starts_with("Basket"), names_to = "Basket", values_to = "Reject") %>%
  mutate(BasketNum = as.integer(gsub("Basket", "", Basket)))

# Tag each row as Power or Type I Error, based on NullMat
raw_data$Metric <- mapply(function(s, b) if(NullMat[s, b]) "Type I Error" else "Power",
                          raw_data$Scenario, raw_data$BasketNum)
raw_data$Metric <- factor(raw_data$Metric, levels = c("Power", "Type I Error"))

raw_data$Method <- factor(raw_data$Method, levels = c("EXNEX-2D", "EXppNEX", "EXppNEX Robust", "BMA-2D"))
raw_data$Basket <- factor(raw_data$Basket, levels = paste0("Basket", 1:5), labels = paste("Basket", 1:5))

scenario_labels <- c(
  "1" = "Global Null", "2" = "Scenario 1", "3" = "Scenario 2",
  "4" = "Scenario 3",  "5" = "Scenario 4", "6" = "Scenario 5"
)
raw_data$Scenario <- factor(raw_data$Scenario, labels = scenario_labels[as.character(1:6)])

# --- Full grid: Scenario (rows) x Basket (columns), linetype = Metric, colour = Method ---
ggplot(raw_data, aes(x = Drift, y = Reject, colour = Method, linetype = Metric)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_grid(Scenario ~ Basket) +
  labs(x = expression(Delta~"(drift)"), y = "Rejection rate (%)", colour = "Method", linetype = "Metric") +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "top"
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("Power" = "solid", "Type I Error" = "dashed")) +
  coord_cartesian(ylim = c(0, 100))+scale_x_continuous(breaks = seq(-0.2, 0.05, by = 0.1)) 


scenario_labels <- c(
  "1" = "Global Null", "2" = "Scenario 1", "3" = "Scenario 2",
  "4" = "Scenario 3",  "5" = "Scenario 4", "6" = "Scenario 5"
)

fwer_scenario_long <- plot_long %>%
  filter(Metric == "FWER (%)")

ggplot(fwer_scenario_long, aes(x = Drift, y = Value, colour = Method)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = 10, linetype = "dashed", colour = "black") +
  facet_wrap(~ Scenario, nrow = 1) +
  labs(x = expression(Delta~"(drift)"), y = "FWER (%)", colour = "Method") +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(-0.2, 0.05, by = 0.05)) +
  coord_cartesian(ylim = c(0, 100))


#Just BMA2D
# --- Prepare Power / Type I Error data for BMA-2D, per basket ---
bma2d_power_error <- raw_data %>%
  filter(Method == "BMA-2D") %>%
  select(Scenario, Basket, Drift, Reject, Metric) %>%
  rename(Value = Reject)

# --- Prepare FWER data for BMA-2D (no basket breakdown) ---
bma2d_fwer <- fwer_scenario_long %>%
  filter(Method == "BMA-2D") %>%
  select(Scenario, Drift, Value) %>%
  mutate(Metric = "FWER")

bma2d_power_error$Metric <- factor(
  bma2d_power_error$Metric,
  levels = c("Power", "Type I Error", "FWER")
)
bma2d_fwer$Metric <- factor(bma2d_fwer$Metric, levels = levels(bma2d_power_error$Metric))

# --- Plot ---
ggplot() +
  geom_line(
    data = bma2d_power_error,
    aes(x = Drift, y = Value, colour = Basket, group = Basket),
    linewidth = 0.7
  ) +
  geom_line(
    data = bma2d_fwer,
    aes(x = Drift, y = Value, group = 1),
    colour = "black",
    linewidth = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_hline(
    data = data.frame(Metric = factor("FWER", levels = levels(bma2d_power_error$Metric)), yint = 10),
    aes(yintercept = yint),
    linetype = "dashed", colour = "black"
  ) +
  facet_grid(Metric ~ Scenario, scales = "free_y") +
  labs(
    x = expression(Delta~"(drift)"),
    y = "Rate (%)",
    colour = "Basket"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(-0.2, 0.05, by = 0.1)) +
  coord_cartesian(ylim = c(0, 100))
