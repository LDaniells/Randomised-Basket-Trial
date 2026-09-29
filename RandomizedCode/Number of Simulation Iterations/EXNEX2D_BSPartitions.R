library(textmineR)
library(matrixStats)
library(parallel)
library(MASS)


filenamesdelta_all <- list.files(pattern="EXNEX2D_Liu_MyWay_seeds_*.", 
                                 full.names=TRUE)
NSims <- function(partitions){
  NsimsResults_EXNEXLiu <- list()
  nsims <- partitions
  EXNEXLiuResults  <- list()
  
  for(i in 1:100){
    filenamesdelta <- filenamesdelta_all[grep(paste0('100Sims',i,".rds"), filenamesdelta_all)]
    EXNEX.Data <- readRDS(filenamesdelta)
    EXNEX <- EXNEX.Data[[1]]
    
    
    quant.try <- seq(0.95,0.99,0.001)
    fwer <- c()
    Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
    for(j in 1:length(quant.try)){
      PostProb <- EXNEX[1:nsims,]
      Q <- colQuantiles(PostProb,probs=quant.try[j])
      #Q[2] <- 0.99
      reject <- c()
      hypo <- matrix(NA,nrow=nsims,ncol=5)
      for(t in 1:5){
        hypo[,t] <- PostProb[,t]>Q[t]
        reject[t] <- 100*length(which(PostProb[,t]>Q[t]))/nsims
      }
      fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/nsims #FWER Under the Null
      Q.mat[j,] <- Q
    }
    Q <- Q.mat[which.min(abs(fwer-10)),] 

    Reject <- matrix(NA,nrow=6,ncol=5)
    MeanPower <- FWER <- c()
    for(j in 1:6){
      PostProb <- EXNEX.Data[[j]]
      PostProb <- PostProb[1:nsims,]
      reject <- c()
      hypo <- matrix(NA,nrow=nsims,ncol=5)
      for(k in 1:5){
        hypo[,k] <- PostProb[,k]>Q[k]
        reject[k] <- 100*(length(which(PostProb[,k]>Q[k]))/nsims)
      }
      Reject[j,] <- reject
      if(j==1){
        MeanPower[j] <- NA
      }else{
        MeanPower[j] <- mean(reject[1:(j-1)])
      }
      if(j==6){
        FWER[j] <- NA
      }else if(j==5){
        FWER[j] <- reject[j]
      }else{
        FWER[j] <- 100*(length(which(apply(hypo[,j:5],1,sum)!=0)))/nsims #FWER Under the Null
      }
    }
    
    xx <- list('Reject'=Reject,'FWER'=FWER,'Mean Power'=MeanPower)
    EXNEXLiuResults[[i]] <- xx
    
    
    print(i)
  }
  NsimsResults_EXNEXLiu <- EXNEXLiuResults
  
  saveRDS(EXNEXLiuResults,file=paste0('EXNEX2D_LiuMyWay_BSPartitions_',partitions,'.rds'))
  
  # save(NsimsResults_BMA,file='BMA_Nsims')
  # save(NsimsResults_BMA2D,file='BMA2D_Nsims')
}

Partitions <- seq(100,10000,10)
partitions <- list()
for(i in 1:length(Partitions)){
  #for(i in 1:2){
  partitions[[i]] <- Partitions[i]
}
Results <- mclapply(partitions,NSims,mc.cores =20)

# task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=1))
# XX <- NSims(partitions[[task_id]])




#For each basket how does power change
filenamesdelta_all <- list.files(pattern="EXNEX2D_LiuMyWay_BSPartitions_*.",
                                 full.names=TRUE)
Partitions <- seq(100,10000,10)
Mean.Power <- Median.Power <- sd.Power <- se.Power <- lower.ci.Power <- upper.ci.Power <- matrix(NA,nrow=length(Partitions),ncol=6)
for(i in 1:length(Partitions)){
  filenamesdelta <- filenamesdelta_all[grep(paste0("BSPartitions_",Partitions[i],".rds"), filenamesdelta_all)]
  Results <- readRDS(filenamesdelta)
  Power <- matrix(NA,nrow=20,ncol=6)
  for(j in 1:20){
    Power[j,] <- Results[[j]]$`Mean Power`
  }
  Mean.Power[i,] <- apply(Power,2,mean)
  Median.Power[i,] <- apply(Power,2,median)
  sd.Power[i,] <- apply(Power,2,sd)
  se.Power[i,] <- sd.Power[i,]/sqrt(20)
  alpha <- 0.05
  lower.ci.Power[i,] <- Mean.Power[i,]-qt(1-alpha/2,df=20-1)*se.Power[i,]
  upper.ci.Power[i,] <- Mean.Power[i,]+qt(1-alpha/2,df=20-1)*se.Power[i,]
}


library(tidyr)
library(dplyr)
library(ggplot2)

df <- data.frame(
  Partition = Partitions,
  Mean.Power,
  lower.ci.Power,
  upper.ci.Power
)

# Name columns
colnames(df) <- c(
  "Partition",
  paste0("Mean_", 1:6),
  paste0("Lower_", 1:6),
  paste0("Upper_", 1:6)
)

# Pivot to long format
df_long <- df %>%
  pivot_longer(
    cols = -Partition,
    names_to = c(".value", "Scenario"),
    names_sep = "_"
  )
# df_plot <- df_long %>%
#   dplyr::filter(Scenario != "1")
df_plot <- df_long %>%
  filter(Scenario != "1")

df_plot$Scenario <- factor(
  df_plot$Scenario,
  levels = c("2","3","4","5","6"),
  labels = c(
    "Scenario 2",
    "Scenario 3",
    "Scenario 4",
    "Scenario 5",
    "Scenario 6"
  )
)




# Initialize a long data frame to store all individual simulation-study means
df_raw <- data.frame()

for(i in seq_along(Partitions)){

  # Get the RDS file for this partition
  filenamesdelta <- filenamesdelta_all[grep(paste0("BSPartitions_",Partitions[i],".rds"), filenamesdelta_all)]
  Results <- readRDS(filenamesdelta)

  nsim <- length(Results)  # should be 20

  # Build Power matrix
  Power <- matrix(NA, nrow = nsim, ncol = 6)
  for(j in 1:nsim){
    Power[j,] <- Results[[j]]$`Mean Power`
  }

  # Convert to long format for this partition
  df_tmp <- as.data.frame(Power)
  df_tmp$SimulationStudy <- 1:nsim
  df_tmp$Partition <- Partitions[i]

  df_tmp_long <- df_tmp %>%
    pivot_longer(
      cols = starts_with("V"),   # columns 1:6 will be V1..V6
      names_to = "Scenario",
      names_prefix = "V",
      values_to = "Power"
    )

  df_raw <- bind_rows(df_raw, df_tmp_long)
}

# Remove scenario 1 if desired
df_raw <- df_raw %>% filter(Scenario != "1")



# Relabel scenarios consistently with df_plot
df_raw$Scenario <- factor(
  df_raw$Scenario,
  levels = c("2","3","4","5","6"),
  labels = c("Scenario 2","Scenario 3","Scenario 4","Scenario 5","Scenario 6")
)
#2000 sims
hline_df <- data.frame(
  Scenario = factor(
    c("Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6"),
    levels = c("Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6")
  ),
  y = c(59.8, 61.2, 59.1, 64.8, 65.9)  # your horizontal line values per scenario
)
#10000 sims
hline2_df <- data.frame(
  Scenario = factor(
    c("Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6"),
    levels = c("Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6")
  ),
  y = c(65.1, 63.3, 59.6, 66.4, 67.8)  # your horizontal line values per scenario
)
# hline2_df <- data.frame(
#   Scenario = factor(
#     c("Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6"),
#     levels = c("Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6")
#   ),
#   y = c(77.3, 75.6, 73.3, 84.6,89 )  # your horizontal line values per scenario
# )


ggplot() +
  # CI ribbon in dark pink


  # Individual simulation-study lines (grey)
  geom_line(
    data = df_raw,
    aes(x = Partition, y = Power, group = SimulationStudy),
    color = "grey70",
    alpha = 0.5,
    size = 0.4
  ) +

  # Vertical reference line
  geom_vline(xintercept = 2000, lty = 2, color = "black") +
  geom_ribbon(
    data = df_plot,
    aes(x = Partition, ymin = Lower, ymax = Upper),
    fill = "#D81B60",
    alpha = 0.4
  ) +

  # Mean line
  geom_line(
    data = df_plot,
    aes(x = Partition, y = Mean),
    color = "#C2185B",
    size = 0.8
  ) +
  # Horizontal lines (per scenario)
  geom_hline(
    data = hline_df,
    aes(yintercept = y),
    lty = 2,          # dashed
    color = "blue",   # color of horizontal lines
    size = 0.8,alpha=0.5
  ) +

  geom_hline(
    data = hline2_df,
    aes(yintercept = y),
    lty = 2,          # dashed
    color = "green",   # color of horizontal lines
    size = 0.8,alpha=0.8
  ) +

  # Facets
  facet_wrap(~Scenario, ncol = 5, scales = "free_y") +

  labs(
    x = "Number of Simulations",
    y = "Mean Power"
  ) +

  theme_bw()


#At 2000/10000 Simulations
filenamesdelta <- filenamesdelta_all[grep(paste0("_",2000,".rds"), filenamesdelta_all)]
Results <- readRDS(filenamesdelta)

MP <- matrix(NA,nrow=100,ncol=6)
for(i in 1:100){
  MP[i,] <- Results[[i]]$`Mean Power`
}
apply(MP,2,min)
apply(MP,2,max)
apply(MP,2,max)-apply(MP,2,min)
