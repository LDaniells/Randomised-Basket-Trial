#EXNEX-Liu
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/FinalCode/Calibration/EXNEXLiu_10000Cal")
files <- list.files(pattern="EXNEX_*.", full.names=TRUE)
file <- files[grep(paste0("EXNEX_LiuMyWay_Cal2_10000_Zeta",1,".rds"), files)]
Results <- readRDS(file)

Calfwer <- matrix(NA,nrow=8,ncol=6)
CalReject <- list()
for(z in 1:8){
file <- files[grep(paste0("EXNEX_LiuMyWay_Cal2_10000_Zeta",z,".rds"), files)]
Results <- readRDS(file)

PostProb <- Results[[1]]
quant.try <- seq(0.9,1,0.0001)
fwer <- c()
Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
for(j in 1:length(quant.try)){
  Q <- colQuantiles(PostProb,probs=quant.try[j])
  #0.947 needed to control FWER to 10%
  Q
  reject <- c()
  hypo <- matrix(NA,nrow=10000,ncol=5)
  for(i in 1:5){
    hypo[,i] <- PostProb[,i]>Q[i]
    reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
  }
  fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/10000 #FWER Under the Null
  Q.mat[j,] <- Q
}
fwer
fwer[which.min(abs(fwer-10))]
Q <- Q.mat[which.min(abs(fwer-10)),]
Q

fwer <- c()
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
Calfwer[z,] <- fwer
CalReject[[z]] <- Reject

}

Calfwer

#Mean Power
MP.mat <- matrix(NA,nrow=8,ncol=6)
for(z in 1:8){
MP <- c()
for(t in 1:6){
  if(t==1){
    MP[t] <- 0
  }else if(t==2){
    MP[t] <- CalReject[[z]][t,1]
  }else if(t==3){
    MP[t] <- mean(CalReject[[z]][t,1:2])
  }else if(t==4){
    MP[t] <- mean(CalReject[[z]][t,1:3])
  }else if(t==5){
    MP[t] <- mean(CalReject[[z]][t,1:4])
  }else{
    MP[t] <- mean(CalReject[[z]][t,1:5])
  }
}
MP.mat[z,] <- MP
}
MP.mat



#Combined
library(dplyr)
library(tidyr)
library(ggplot2)

zeta <- c(0.125,0.25,0.5,1,1.5,2.5,5,10)

# Power data
power_data <- data.frame(
  zeta = zeta,
  MP.mat
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "Mean Power"
  )

# FWER data
fwer_data <- data.frame(
  zeta = zeta,
  Calfwer
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "FWER"
  )

# Combine
plot_data <- bind_rows(power_data, fwer_data)

plot_data$zeta <- factor(
  plot_data$zeta,
  levels = c(0.125,0.25,0.5,1,1.5,2.5,5,10)
)

plot_data2 <- plot_data %>%
  filter(
    !(Scenario == "X1" & Measure == "Mean Power"),  # remove power under global null
    !(Scenario == "X6" & Measure == "FWER")         # remove FWER under scenario 5
  )

EXNEXLiu <- ggplot(plot_data2,
       aes(x = zeta,
           y = Value,
           colour = Measure,
           linetype = Measure,
           group = Measure)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_wrap(~Scenario,
             ncol = 6,
             scales = "fixed",
             labeller = labeller(
               Scenario = c(
                 X1 = "Global Null",
                 X2 = "Scenario 1",
                 X3 = "Scenario 2",
                 X4 = "Scenario 3",
                 X5 = "Scenario 4",
                 X6 = "Scenario 5"
               )
             )) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100,20)
  ) +
  
  geom_hline(
    data = subset(plot_data, Measure == "FWER"),
    aes(yintercept = 10),
    colour = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = expression(zeta),
    y = "Percentage (%)",
    colour = NULL,
    linetype = NULL,
    title='EXNEX Liu'
  ) +
  
  scale_colour_manual(
    values = c(
      "Mean Power" = "#C2185B",
      "FWER" = "#0072B2"
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )+theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )+theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

which.max(apply(MP.mat[,-1],1,prod)^(1/5))
which.max((apply(MP.mat[,-1],1,prod)^(1/5))-abs((10-apply(Calfwer[,-6],1,prod)^(1/5))))

#EXNEX Joint

#EXppNEX
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/FinalCode/Calibration/EXppNEX_Cal10000")
files <- list.files(pattern="EXppNEX_*.", full.names=TRUE)
file <- files[grep(paste0("EXppNEX_NonRobustified_Cal2_10000_",1,".rds"), files)]
Results <- readRDS(file)

Calfwer <- matrix(NA,nrow=8,ncol=6)
CalReject <- list()
for(z in 1:8){
  file <- files[grep(paste0("EXppNEX_NonRobustified_Cal2_10000_",z,".rds"), files)]
  Results <- readRDS(file)
  
  PostProb <- Results[[1]]
  quant.try <- seq(0.9,1,0.0001)
  fwer <- c()
  Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
  for(j in 1:length(quant.try)){
    Q <- colQuantiles(PostProb,probs=quant.try[j])
    #0.947 needed to control FWER to 10%
    Q
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/10000 #FWER Under the Null
    Q.mat[j,] <- Q
  }
  fwer
  fwer[which.min(abs(fwer-10))]
  Q <- Q.mat[which.min(abs(fwer-10)),]
  Q
  
  fwer <- c()
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
  Calfwer[z,] <- fwer
  CalReject[[z]] <- Reject
  
}

Calfwer

#Mean Power
MP.mat <- matrix(NA,nrow=8,ncol=6)
for(z in 1:8){
  MP <- c()
  for(t in 1:6){
    if(t==1){
      MP[t] <- 0
    }else if(t==2){
      MP[t] <- CalReject[[z]][t,1]
    }else if(t==3){
      MP[t] <- mean(CalReject[[z]][t,1:2])
    }else if(t==4){
      MP[t] <- mean(CalReject[[z]][t,1:3])
    }else if(t==5){
      MP[t] <- mean(CalReject[[z]][t,1:4])
    }else{
      MP[t] <- mean(CalReject[[z]][t,1:5])
    }
  }
  MP.mat[z,] <- MP
}
MP.mat



#Combined
library(dplyr)
library(tidyr)
library(ggplot2)

zeta <- c(0.125,0.25,0.5,1,1.5,2.5,5,10)

# Power data
power_data <- data.frame(
  zeta = zeta,
  MP.mat
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "Mean Power"
  )

# FWER data
fwer_data <- data.frame(
  zeta = zeta,
  Calfwer
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "FWER"
  )

# Combine
plot_data <- bind_rows(power_data, fwer_data)

plot_data$zeta <- factor(
  plot_data$zeta,
  levels = c(0.125,0.25,0.5,1,1.5,2.5,5,10)
)

plot_data2 <- plot_data %>%
  filter(
    !(Scenario == "X1" & Measure == "Mean Power"),  # remove power under global null
    !(Scenario == "X6" & Measure == "FWER")         # remove FWER under scenario 5
  )

EXppNEX <- ggplot(plot_data2,
       aes(x = zeta,
           y = Value,
           colour = Measure,
           linetype = Measure,
           group = Measure)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_wrap(~Scenario,
             ncol = 6,
             scales = "fixed",
             labeller = labeller(
               Scenario = c(
                 X1 = "Global Null",
                 X2 = "Scenario 1",
                 X3 = "Scenario 2",
                 X4 = "Scenario 3",
                 X5 = "Scenario 4",
                 X6 = "Scenario 5"
               )
             )) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100,20)
  ) +
  
  geom_hline(
    data = subset(plot_data, Measure == "FWER"),
    aes(yintercept = 10),
    colour = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = expression(zeta),
    y = "Percentage (%)",
    colour = NULL,
    linetype = NULL,
    title='EXppNEX'
  ) +
  
  scale_colour_manual(
    values = c(
      "Mean Power" = "#C2185B",
      "FWER" = "#0072B2"
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )+theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )+theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

which.max(apply(MP.mat[,-1],1,prod)^(1/5))
which.max((apply(MP.mat[,-1],1,prod)^(1/5))-abs((10-apply(Calfwer[,-6],1,prod)^(1/5))))


#EXppNEX Robust
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Calibration/EXppNEXR_Cal10000")
files <- list.files(pattern="EXppNEX_*.", full.names=TRUE)
file <- files[grep(paste0("EXppNEX_Robustified_Cal2_10000_",1,".rds"), files)]
Results <- readRDS(file)

Calfwer <- matrix(NA,nrow=8,ncol=6)
CalReject <- list()
for(z in 1:8){
  file <- files[grep(paste0("EXppNEX_Robustified_Cal2_10000_",z,".rds"), files)]
  Results <- readRDS(file)
  
  PostProb <- Results[[1]]
  quant.try <- seq(0.9,1,0.0001)
  fwer <- c()
  Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
  for(j in 1:length(quant.try)){
    Q <- colQuantiles(PostProb,probs=quant.try[j])
    #0.947 needed to control FWER to 10%
    Q
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/10000 #FWER Under the Null
    Q.mat[j,] <- Q
  }
  fwer
  fwer[which.min(abs(fwer-10))]
  Q <- Q.mat[which.min(abs(fwer-10)),]
  Q
  
  fwer <- c()
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
  Calfwer[z,] <- fwer
  CalReject[[z]] <- Reject
  
}

Calfwer

#Mean Power
MP.mat <- matrix(NA,nrow=8,ncol=6)
for(z in 1:8){
  MP <- c()
  for(t in 1:6){
    if(t==1){
      MP[t] <- 0
    }else if(t==2){
      MP[t] <- CalReject[[z]][t,1]
    }else if(t==3){
      MP[t] <- mean(CalReject[[z]][t,1:2])
    }else if(t==4){
      MP[t] <- mean(CalReject[[z]][t,1:3])
    }else if(t==5){
      MP[t] <- mean(CalReject[[z]][t,1:4])
    }else{
      MP[t] <- mean(CalReject[[z]][t,1:5])
    }
  }
  MP.mat[z,] <- MP
}
MP.mat



#Combined
library(dplyr)
library(tidyr)
library(ggplot2)

zeta <- c(0.125,0.25,0.5,1,1.5,2.5,5,10)

# Power data
power_data <- data.frame(
  zeta = zeta,
  MP.mat
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "Mean Power"
  )

# FWER data
fwer_data <- data.frame(
  zeta = zeta,
  Calfwer
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "FWER"
  )

# Combine
plot_data <- bind_rows(power_data, fwer_data)

plot_data$zeta <- factor(
  plot_data$zeta,
  levels = c(0.125,0.25,0.5,1,1.5,2.5,5,10)
)

plot_data2 <- plot_data %>%
  filter(
    !(Scenario == "X1" & Measure == "Mean Power"),  # remove power under global null
    !(Scenario == "X6" & Measure == "FWER")         # remove FWER under scenario 5
  )

EXppNEXR <- ggplot(plot_data2,
       aes(x = zeta,
           y = Value,
           colour = Measure,
           linetype = Measure,
           group = Measure)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_wrap(~Scenario,
             ncol = 6,
             scales = "fixed",
             labeller = labeller(
               Scenario = c(
                 X1 = "Global Null",
                 X2 = "Scenario 1",
                 X3 = "Scenario 2",
                 X4 = "Scenario 3",
                 X5 = "Scenario 4",
                 X6 = "Scenario 5"
               )
             )) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100,20)
  ) +
  
  geom_hline(
    data = subset(plot_data, Measure == "FWER"),
    aes(yintercept = 10),
    colour = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = expression(zeta),
    y = "Percentage (%)",
    colour = NULL,
    linetype = NULL,
    title='EXppNEX Robustified'
  ) +
  
  scale_colour_manual(
    values = c(
      "Mean Power" = "#C2185B",
      "FWER" = "#0072B2"
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )+theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )+theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

which.max(apply(MP.mat[,-1],1,prod)^(1/5))
Calfwer
which.max((apply(MP.mat[,-1],1,prod)^(1/5))-abs((10-apply(Calfwer[,-6],1,prod)^(1/5))))

#EXNEX2D
setwd("~/Library/CloudStorage/OneDrive-UniversityofCambridge,MRCBiostatisticsUnit/Research Work/Randomized Basket Trial/FinalCode/Calibration/EXNEX2D_Cal10000")
files <- list.files(pattern="EXNEX2D_*.", full.names=TRUE)
file <- files[grep(paste0("EXNEX2D_Liu_Cal2_Zeta",1,".rds"), files)]
Results <- readRDS(file)

Calfwer <- matrix(NA,nrow=8,ncol=6)
CalReject <- list()
for(z in 1:8){
  file <- files[grep(paste0("EXNEX2D_Liu_Cal2_Zeta",z,".rds"), files)]
  Results <- readRDS(file)
  
  PostProb <- Results[[1]]
  quant.try <- seq(0.9,1,0.0001)
  fwer <- c()
  Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
  for(j in 1:length(quant.try)){
    Q <- colQuantiles(PostProb,probs=quant.try[j])
    #0.947 needed to control FWER to 10%
    Q
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/10000 #FWER Under the Null
    Q.mat[j,] <- Q
  }
  fwer
  fwer[which.min(abs(fwer-10))]
  Q <- Q.mat[which.min(abs(fwer-10)),]
  Q
  
  fwer <- c()
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
  Calfwer[z,] <- fwer
  CalReject[[z]] <- Reject
  
}

Calfwer

#Mean Power
MP.mat <- matrix(NA,nrow=8,ncol=6)
for(z in 1:8){
  MP <- c()
  for(t in 1:6){
    if(t==1){
      MP[t] <- 0
    }else if(t==2){
      MP[t] <- CalReject[[z]][t,1]
    }else if(t==3){
      MP[t] <- mean(CalReject[[z]][t,1:2])
    }else if(t==4){
      MP[t] <- mean(CalReject[[z]][t,1:3])
    }else if(t==5){
      MP[t] <- mean(CalReject[[z]][t,1:4])
    }else{
      MP[t] <- mean(CalReject[[z]][t,1:5])
    }
  }
  MP.mat[z,] <- MP
}
MP.mat



#Combined
library(dplyr)
library(tidyr)
library(ggplot2)

zeta <- c(0.125,0.25,0.5,1,1.5,2.5,5,10)

# Power data
power_data <- data.frame(
  zeta = zeta,
  MP.mat
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "Mean Power"
  )

# FWER data
fwer_data <- data.frame(
  zeta = zeta,
  Calfwer
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "FWER"
  )

# Combine
plot_data <- bind_rows(power_data, fwer_data)

plot_data$zeta <- factor(
  plot_data$zeta,
  levels = c(0.125,0.25,0.5,1,1.5,2.5,5,10)
)

plot_data2 <- plot_data %>%
  filter(
    !(Scenario == "X1" & Measure == "Mean Power"),  # remove power under global null
    !(Scenario == "X6" & Measure == "FWER")         # remove FWER under scenario 5
  )

EXNEX2D <- ggplot(plot_data2,
       aes(x = zeta,
           y = Value,
           colour = Measure,
           linetype = Measure,
           group = Measure)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_wrap(~Scenario,
             ncol = 6,
             scales = "fixed",
             labeller = labeller(
               Scenario = c(
                 X1 = "Global Null",
                 X2 = "Scenario 1",
                 X3 = "Scenario 2",
                 X4 = "Scenario 3",
                 X5 = "Scenario 4",
                 X6 = "Scenario 5"
               )
             )) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100,20)
  ) +
  
  geom_hline(
    data = subset(plot_data, Measure == "FWER"),
    aes(yintercept = 10),
    colour = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = expression(zeta),
    y = "Percentage (%)",
    colour = NULL,
    linetype = NULL,
    title='EXNEX2D'
  ) +
  
  scale_colour_manual(
    values = c(
      "Mean Power" = "#C2185B",
      "FWER" = "#0072B2"
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )+theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )+theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )


which.max(apply(MP.mat[,-1],1,prod)^(1/5))
which.max((apply(MP.mat[,-1],1,prod)^(1/5))-abs((10-apply(Calfwer[,-6],1,prod)^(1/5))))



#EXNEX- Joint
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/Research Work/Randomized Basket Trial/FinalCode/Calibration/EXNEXJoint")
files <- list.files(pattern="EXNEX_*.", full.names=TRUE)
file <- files[grep(paste0("EXNEX_Joint_Cal2_10000_Zeta",1,".rds"), files)]
Results <- readRDS(file)

Calfwer <- matrix(NA,nrow=8,ncol=6)
CalReject <- list()
for(z in 1:8){
  file <- files[grep(paste0("EXNEX_Joint_Cal2_10000_Zeta",z,".rds"), files)]
  Results <- readRDS(file)
  
  PostProb <- Results[[1]]
  quant.try <- seq(0.9,1,0.0001)
  fwer <- c()
  Q.mat <- matrix(NA,nrow=length(quant.try),ncol=5)
  for(j in 1:length(quant.try)){
    Q <- colQuantiles(PostProb,probs=quant.try[j])
    #0.947 needed to control FWER to 10%
    Q
    reject <- c()
    hypo <- matrix(NA,nrow=10000,ncol=5)
    for(i in 1:5){
      hypo[,i] <- PostProb[,i]>Q[i]
      reject[i] <- 100*length(which(PostProb[,i]>Q[i]))/10000
    }
    fwer[j] <- 100*(length(which(apply(hypo,1,sum)!=0)))/10000 #FWER Under the Null
    Q.mat[j,] <- Q
  }
  fwer
  fwer[which.min(abs(fwer-10))]
  Q <- Q.mat[which.min(abs(fwer-10)),]
  Q
  
  fwer <- c()
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
  Calfwer[z,] <- fwer
  CalReject[[z]] <- Reject
  
}

Calfwer

#Mean Power
MP.mat <- matrix(NA,nrow=8,ncol=6)
for(z in 1:8){
  MP <- c()
  for(t in 1:6){
    if(t==1){
      MP[t] <- 0
    }else if(t==2){
      MP[t] <- CalReject[[z]][t,1]
    }else if(t==3){
      MP[t] <- mean(CalReject[[z]][t,1:2])
    }else if(t==4){
      MP[t] <- mean(CalReject[[z]][t,1:3])
    }else if(t==5){
      MP[t] <- mean(CalReject[[z]][t,1:4])
    }else{
      MP[t] <- mean(CalReject[[z]][t,1:5])
    }
  }
  MP.mat[z,] <- MP
}
MP.mat



#Combined
library(dplyr)
library(tidyr)
library(ggplot2)

zeta <- c(0.125,0.25,0.5,1,1.5,2.5,5,10)

# Power data
power_data <- data.frame(
  zeta = zeta,
  MP.mat
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "Mean Power"
  )

# FWER data
fwer_data <- data.frame(
  zeta = zeta,
  Calfwer
) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Measure = "FWER"
  )

# Combine
plot_data <- bind_rows(power_data, fwer_data)

plot_data$zeta <- factor(
  plot_data$zeta,
  levels = c(0.125,0.25,0.5,1,1.5,2.5,5,10)
)

plot_data2 <- plot_data %>%
  filter(
    !(Scenario == "X1" & Measure == "Mean Power"),  # remove power under global null
    !(Scenario == "X6" & Measure == "FWER")         # remove FWER under scenario 5
  )

EXNEXJoint <- ggplot(plot_data2,
                   aes(x = zeta,
                       y = Value,
                       colour = Measure,
                       linetype = Measure,
                       group = Measure)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_wrap(~Scenario,
             ncol = 6,
             scales = "fixed",
             labeller = labeller(
               Scenario = c(
                 X1 = "Global Null",
                 X2 = "Scenario 1",
                 X3 = "Scenario 2",
                 X4 = "Scenario 3",
                 X5 = "Scenario 4",
                 X6 = "Scenario 5"
               )
             )) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100,20)
  ) +
  
  geom_hline(
    data = subset(plot_data, Measure == "FWER"),
    aes(yintercept = 10),
    colour = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = expression(zeta),
    y = "Percentage (%)",
    colour = NULL,
    linetype = NULL,
    title='EXNEX Joint'
  ) +
  
  scale_colour_manual(
    values = c(
      "Mean Power" = "#C2185B",
      "FWER" = "#0072B2"
    )
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )+theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )+theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

which.max(apply(MP.mat[,-1],1,prod)^(1/5))
which.max((apply(MP.mat[,-1],1,prod)^(1/5))-abs((10-apply(Calfwer[,-6],1,prod)^(1/5))))


library(ggpubr)
ggarrange(
  EXNEXLiu,
  EXNEXJoint,
  EXNEX2D,
  EXppNEX,
  EXppNEXR,
  nrow = 5,
  ncol = 1,
  common.legend = TRUE,
  legend = "top"
)


library(patchwork)
library(ggpubr)

# Extract legend from one of the plots
shared_legend <- get_legend(EXNEXLiu + theme(legend.position = "top"))

# Convert legend to a patchwork-compatible grob
legend_plot <- as_ggplot(shared_legend)

# Stack: legend on top, then the 5 plots below with no individual legends
(
  legend_plot /
    no_x_axis(EXNEXLiu + theme(legend.position = "none")) /
    no_x_axis(EXNEXJoint + theme(legend.position = "none")) /
    no_x_axis(EXNEX2D + theme(legend.position = "none")) /
    no_x_axis(EXppNEX + theme(legend.position = "none")) /
    angle_x_axis(EXppNEXR + theme(legend.position = "none"))
) +
  plot_layout(
    ncol = 1,
    heights = c(0.1, 1, 1, 1, 1, 1)  # small height for legend row
  )
