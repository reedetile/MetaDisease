#Description-----------------------------------------
#A frequency dependent model of 2 patch 6 (possible) species metacommunities
#  09 Jan 2025
#RCS

#Initialize -----------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(betafunctions)
set.seed(1234)

#############################################################################################
##################Simulate over range of connectivities##################################
###########################################################################################
rm(list = ls())
main.wd <- "D:/gitrepos/MetaDisease"
graphs <- paste(getwd(),"/Graphs",sep = "")
setwd(main.wd)
source('Epimodel_funcs.R')
source('Metacomm_2patch.R')
#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)

birth <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
# if we change to a multi-season model, will need to add in some birth rate
# b <- c(0.6,0.5,0.4,0.3,0.2,0.1)/time #host birth rate
d <- c(0.024,
       NA,
       NA,
       0.0029,
       0.0144,
       0.04)
d[2] <- runif(n = 1, min = d[[4]], max = d[[1]])       
d[3] <- runif(n = 1, min = d[[4]], max = d[[2]])       
# Citations for death rate
# PREG = Jameson 1956
# ABOR = Pilliod 2010
# RCAT = Howell 2020
# RDRAY = Feller 2017
# TTOR = johnson 2013 + assumptions of dilution effects
# TGRAN =  johnson 2013 + assumptions of dilution effects
v <- c(NA,
       NA,
       NA,
       NA,
       0.012460223,
       NA)
# Citation for RCAT recovery = Daszak 2004
low_recovery <- runif(n = 4, min = 0, max = v[5])
v[1] <- min(low_recovery)
v[4] <- max(low_recovery)
mid_low_recovery <- low_recovery[low_recovery != v[1] & low_recovery != v[4]]
v[2] <- min(mid_low_recovery)
v[3] <- max(mid_low_recovery)
v[6] <- runif(n = 1, min = v[5], max = 0.1) #should find citation for max recovery rate
v

#Determining dispersal rate
psi <- c(0.83, 0.62, 0.14,NA, 0.12,0.64)
psi[4] <- runif(n = 1, min = psi[5], max = psi[3])

#Can't define N_meta here, so will need to do so within simulation
# N_meta <- colSums(N)
# phi <- ifelse(N_meta > 0, psi/N_meta, 0)
# phi
# alpha <- #disease specific mortality
# species_chara <- data.frame(Species = Species,
#                             birth = birth, 
#                             death = d, 
#                             recovery = v, 
#                             dispersal = phi)
### Transmission coefficient
# beta should be higher for intraspecific transmission than interspecific
#inter-specific should have higher rate from more competent species
# Most abundant species should have highest interspecific trans rates
# For now, keeping transmission from 1 spp -> all other spp the same
#Ex: P Regilla will have the same transmission to A Boreas, T. Taricha... R Draytonii
# We are using a beta distribution b/c that is the best for the probability scale
# we need to assign probabilities to each species
trans_rate <- function(n = 1,l = 0, u, alpha,beta){
  trans <- rBeta.4P(n = n, l = l, u = u, alpha = alpha, beta = beta)
  return(trans)
}


# PREG
intra_PREG <- trans_rate(l = 0, u = 1, alpha = 4,beta = 2)
inter_PREG <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.5)

# TGRAN
intra_TGRAN <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.25)
inter_TGRAN <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4, beta = 2.5)

#TTOR
intra_TTOR <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4,beta = 2)
inter_TTOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 3, beta = 2.75)


#ABOR
intra_ABOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 4,beta = 2)
inter_ABOR <- trans_rate(l = 0, u = intra_ABOR, alpha = 2.5, beta = 3.0)


#RCAT
intra_RCAT <- trans_rate(l = 0, u = intra_ABOR, alpha = 4,beta = 2)
inter_RCAT <- trans_rate(l = 0, u = intra_RCAT, alpha = 2.0, beta = 3.25)

#RDRAY
intra_RDRAY <- trans_rate(l = 0, u = intra_RCAT, alpha = 4,beta = 2)
inter_RDRAY <- trans_rate(l = 0, u = intra_RDRAY, alpha = 1.5, beta = 3.5)





beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- if(i == 1 & j == 1){
      intra_PREG}else if(i != 1 & j == 1){
        inter_PREG} else if(i == 2 & j == 2){
          intra_TGRAN} else if(i != 2 & j == 2){
            inter_TGRAN} else if(i == 3 & j == 3){
              intra_TTOR} else if(i != 3 & j == 3){
                inter_TTOR} else if(i == 4 & j == 4){
                  intra_ABOR} else if(i != 4 & j == 4){
                    inter_ABOR} else if(i == 5 & j == 5){
                      intra_RCAT} else if(i != 5 & j ==5){
                        inter_RCAT} else if(i == 6 & j == 6){
                          intra_RDRAY} else if(i != 6 & j == 6){
                            inter_RDRAY}
  }
}
beta <- beta/90


# Simulation over metacommunities
connect_sens <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 19))
colnames(connect_sens) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                       "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                       "BetaDiversity","Gamma_diversity",
                       "Beta_relative","LandscapeR0", "MetaCommID","Connectivity")
for (a in 1:length(meta_comm_list)) {
  stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
  go <- 1 - stay # probability individuals move. Should make sure this is always higher
  connect <- matrix(data = NA,
              nrow = num_patches, 
              ncol = num_patches)
  for(i in 1:ncol(connect)){
    for(j in 1:nrow(connect)){
      connect[i,j] <- ifelse(i == j, stay, go)
    }
  }
  S <- ceiling(meta_comm_list[[a]][,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #value of susceptibles
  I <- meta_comm_list[[a]][,1:6] - S #value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- ifelse(N_meta > 0, psi/N_meta, 0) #get dispersal rate
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc R0
  b <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc b
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-connect[p,q]*S[p,i] + connect[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-connect[p,q]*I[p,i] + connect[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
        }
      }
    }
  }
  # Calculate variables
  connect_sens[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    connect_sens[a,i] <- sum(N[,i-1])
    connect_sens[a,(i+num_spp)] <- connect_sens[a,i]/connect_sens[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'w')
  connect_sens[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  connect_sens[a,15] <- sum(ifelse(connect_sens[a,1:6] > 0, 1,0)) #gamma diversity
  
  connect_sens[a,16] <- connect_sens[a,7]/connect_sens[a,8] #beta diversity / total abundance
  
  #calculate landscape R0
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = connect,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  connect_sens[a,17] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  connect_sens[a,18] <- a #metacommunity ID
  connect_sens[a,19] <- go # probability leaving a patch
}



# lets look at the effect of connectivity
connect_plot <- ggplot(data = connect_sens, aes(x = Connectivity, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Connectivity")+
  ylab("Landscape R0")
connect_plot
setwd(graphs)
ggsave(filename = "connectplot.png", plot = connect_plot)


#############################################################################################
##################Simulate over range of transmission rate##################################
###########################################################################################
rm(list = ls())
main.wd <- "D:/gitrepos/MetaDisease"
graphs <- paste(main.wd,"/Graphs",sep = "")
setwd(main.wd)
source('Epimodel_funcs.R')
source('Metacomm_2patch.R')
#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)

birth <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
# if we change to a multi-season model, will need to add in some birth rate
# b <- c(0.6,0.5,0.4,0.3,0.2,0.1)/time #host birth rate
d <- c(0.024,
       NA,
       NA,
       0.0029,
       0.0144,
       0.04)
d[2] <- runif(n = 1, min = d[[4]], max = d[[1]])       
d[3] <- runif(n = 1, min = d[[4]], max = d[[2]])       
# Citations for death rate
# PREG = Jameson 1956
# ABOR = Pilliod 2010
# RCAT = Howell 2020
# RDRAY = Feller 2017
# TTOR = johnson 2013 + assumptions of dilution effects
# TGRAN =  johnson 2013 + assumptions of dilution effects
v <- c(NA,
       NA,
       NA,
       NA,
       0.012460223,
       NA)
# Citation for RCAT recovery = Daszak 2004
low_recovery <- runif(n = 4, min = 0, max = v[5])
v[1] <- min(low_recovery)
v[4] <- max(low_recovery)
mid_low_recovery <- low_recovery[low_recovery != v[1] & low_recovery != v[4]]
v[2] <- min(mid_low_recovery)
v[3] <- max(mid_low_recovery)
v[6] <- runif(n = 1, min = v[5], max = 0.1) #should find citation for max recovery rate
v

#Determining dispersal rate
psi <- c(0.83, 0.62, 0.14,NA, 0.12,0.64)
psi[4] <- runif(n = 1, min = psi[5], max = psi[3])
stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
go <- 1 - stay # probability individuals move. Should make sure this is always higher
connect <- matrix(data = NA,
                  nrow = num_patches, 
                  ncol = num_patches)
for(i in 1:ncol(connect)){
  for(j in 1:nrow(connect)){
    connect[i,j] <- ifelse(i == j, stay, go)
  }
}



trans_sens <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 19))
colnames(trans_sens) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                            "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                            "BetaDiversity","Gamma_diversity",
                            "Beta_relative","LandscapeR0", "MetaCommID","Max_trans")

# function for defining transmission rate
trans_rate <- function(n = 1,l = 0, u, alpha,beta){
  trans <- rBeta.4P(n = n, l = l, u = u, alpha = alpha, beta = beta)
  return(trans)
}

for (a in 1:length(meta_comm_list)) {
  # PREG
  intra_PREG <- trans_rate(l = 0, u = 1, alpha = 4,beta = 2)
  inter_PREG <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.5)
  
  # TGRAN
  intra_TGRAN <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.25)
  inter_TGRAN <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4, beta = 2.5)
  
  #TTOR
  intra_TTOR <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4,beta = 2)
  inter_TTOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 3, beta = 2.75)
  
  
  #ABOR
  intra_ABOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 4,beta = 2)
  inter_ABOR <- trans_rate(l = 0, u = intra_ABOR, alpha = 2.5, beta = 3.0)
  
  
  #RCAT
  intra_RCAT <- trans_rate(l = 0, u = intra_ABOR, alpha = 4,beta = 2)
  inter_RCAT <- trans_rate(l = 0, u = intra_RCAT, alpha = 2.0, beta = 3.25)
  
  #RDRAY
  intra_RDRAY <- trans_rate(l = 0, u = intra_RCAT, alpha = 4,beta = 2)
  inter_RDRAY <- trans_rate(l = 0, u = intra_RDRAY, alpha = 1.5, beta = 3.5)
  
  
  
  
  
  beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
  for (i in 1:nrow(beta)) {
    for (j in 1:ncol(beta)) {
      beta[i,j] <- if(i == 1 & j == 1){
        intra_PREG}else if(i != 1 & j == 1){
          inter_PREG} else if(i == 2 & j == 2){
            intra_TGRAN} else if(i != 2 & j == 2){
              inter_TGRAN} else if(i == 3 & j == 3){
                intra_TTOR} else if(i != 3 & j == 3){
                  inter_TTOR} else if(i == 4 & j == 4){
                    intra_ABOR} else if(i != 4 & j == 4){
                      inter_ABOR} else if(i == 5 & j == 5){
                        intra_RCAT} else if(i != 5 & j ==5){
                          inter_RCAT} else if(i == 6 & j == 6){
                            intra_RDRAY} else if(i != 6 & j == 6){
                              inter_RDRAY}
    }
  }
  beta <- beta/90
  S <- ceiling(meta_comm_list[[a]][,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #value of susceptibles
  I <- meta_comm_list[[a]][,1:6] - S #value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- ifelse(N_meta > 0, psi/N_meta, 0) #get dispersal rate
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc R0
  b <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc b
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-connect[p,q]*S[p,i] + connect[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-connect[p,q]*I[p,i] + connect[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
        }
      }
    }
  }
  # Calculate variables
  trans_sens[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    trans_sens[a,i] <- sum(N[,i-1])
    trans_sens[a,(i+num_spp)] <- trans_sens[a,i]/trans_sens[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'w')
  trans_sens[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  trans_sens[a,15] <- sum(ifelse(trans_sens[a,1:6] > 0, 1,0)) #gamma diversity
  
  trans_sens[a,16] <- trans_sens[a,7]/trans_sens[a,8] #beta diversity / total abundance
  
  #calculate landscape R0
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = connect,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  trans_sens[a,17] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  trans_sens[a,18] <- a #metacommunity ID
  trans_sens[a,19] <- intra_PREG # probability leaving a patch
}

trans_plot <- ggplot(data = trans_sens, aes(x = Max_trans, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Transmission rate")+
  ylab("Landscape R0")
trans_plot  
setwd(graphs)
ggsave(filename = "transplot.png", plot = trans_plot)


#############################################################################################
##################Simulate over range of death rates##################################
###########################################################################################
rm(list = ls())
main.wd <- "D:/gitrepos/MetaDisease"
graphs <- paste(main.wd,"/Graphs",sep = "")
setwd(main.wd)
source('Epimodel_funcs.R')
source('Metacomm_2patch.R')

#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)

birth <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
# if we change to a multi-season model, will need to add in some birth rate
# b <- c(0.6,0.5,0.4,0.3,0.2,0.1)/time #host birth rate

v <- c(NA,
       NA,
       NA,
       NA,
       0.012460223,
       NA)
# Citation for RCAT recovery = Daszak 2004
low_recovery <- runif(n = 4, min = 0, max = v[5])
v[1] <- min(low_recovery)
v[4] <- max(low_recovery)
mid_low_recovery <- low_recovery[low_recovery != v[1] & low_recovery != v[4]]
v[2] <- min(mid_low_recovery)
v[3] <- max(mid_low_recovery)
v[6] <- runif(n = 1, min = v[5], max = 0.1) #should find citation for max recovery rate
v

#Determining dispersal rate
psi <- c(0.83, 0.62, 0.14,NA, 0.12,0.64)
psi[4] <- runif(n = 1, min = psi[5], max = psi[3])

#Can't define N_meta here, so will need to do so within simulation
# N_meta <- colSums(N)
# phi <- ifelse(N_meta > 0, psi/N_meta, 0)
# phi
# alpha <- #disease specific mortality
# species_chara <- data.frame(Species = Species,
#                             birth = birth, 
#                             death = d, 
#                             recovery = v, 
#                             dispersal = phi)
### Transmission coefficient
# beta should be higher for intraspecific transmission than interspecific
#inter-specific should have higher rate from more competent species
# Most abundant species should have highest interspecific trans rates
# For now, keeping transmission from 1 spp -> all other spp the same
#Ex: P Regilla will have the same transmission to A Boreas, T. Taricha... R Draytonii
# We are using a beta distribution b/c that is the best for the probability scale
# we need to assign probabilities to each species
trans_rate <- function(n = 1,l = 0, u, alpha,beta){
  trans <- rBeta.4P(n = n, l = l, u = u, alpha = alpha, beta = beta)
  return(trans)
}


# PREG
intra_PREG <- trans_rate(l = 0, u = 1, alpha = 4,beta = 2)
inter_PREG <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.5)

# TGRAN
intra_TGRAN <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.25)
inter_TGRAN <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4, beta = 2.5)

#TTOR
intra_TTOR <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4,beta = 2)
inter_TTOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 3, beta = 2.75)


#ABOR
intra_ABOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 4,beta = 2)
inter_ABOR <- trans_rate(l = 0, u = intra_ABOR, alpha = 2.5, beta = 3.0)


#RCAT
intra_RCAT <- trans_rate(l = 0, u = intra_ABOR, alpha = 4,beta = 2)
inter_RCAT <- trans_rate(l = 0, u = intra_RCAT, alpha = 2.0, beta = 3.25)

#RDRAY
intra_RDRAY <- trans_rate(l = 0, u = intra_RCAT, alpha = 4,beta = 2)
inter_RDRAY <- trans_rate(l = 0, u = intra_RDRAY, alpha = 1.5, beta = 3.5)

# meta-community characteristics
#Connectivity of patches
stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
go <- 1 - stay # probability individuals move
connect <- matrix(data = NA,
                  nrow = num_patches, 
                  ncol = num_patches)
for(i in 1:ncol(connect)){
  for(j in 1:nrow(connect)){
    connect[i,j] <- ifelse(i == j, stay, go)
  }
}



beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- if(i == 1 & j == 1){
      intra_PREG}else if(i != 1 & j == 1){
        inter_PREG} else if(i == 2 & j == 2){
          intra_TGRAN} else if(i != 2 & j == 2){
            inter_TGRAN} else if(i == 3 & j == 3){
              intra_TTOR} else if(i != 3 & j == 3){
                inter_TTOR} else if(i == 4 & j == 4){
                  intra_ABOR} else if(i != 4 & j == 4){
                    inter_ABOR} else if(i == 5 & j == 5){
                      intra_RCAT} else if(i != 5 & j ==5){
                        inter_RCAT} else if(i == 6 & j == 6){
                          intra_RDRAY} else if(i != 6 & j == 6){
                            inter_RDRAY}
  }
}
beta <- beta/90

# run simulation #
death_sens <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 21))
colnames(death_sens) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                            "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                            "BetaDiversity","Gamma_diversity",
                            "Beta_relative","LandscapeR0", "MetaCommID","Death_Max","Death_Min","Death_range")
for (a in 1:length(meta_comm_list)) {
  d <- runif(6)
  S <- ceiling(meta_comm_list[[a]][,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #value of susceptibles
  I <- meta_comm_list[[a]][,1:6] - S #value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- ifelse(N_meta > 0, psi/N_meta, 0) #get dispersal rate
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc R0
  b <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc b
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-connect[p,q]*S[p,i] + connect[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-connect[p,q]*I[p,i] + connect[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
        }
      }
    }
  }
  # Calculate variables
  death_sens[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    death_sens[a,i] <- sum(N[,i-1])
    death_sens[a,(i+num_spp)] <- death_sens[a,i]/death_sens[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'w')
  death_sens[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  death_sens[a,15] <- sum(ifelse(death_sens[a,1:6] > 0, 1,0)) #gamma diversity
  
  death_sens[a,16] <- death_sens[a,7]/death_sens[a,8] #beta diversity / total abundance
  
  #calculate landscape R0
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = connect,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  death_sens[a,17] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  death_sens[a,18] <- a #metacommunity ID
  death_sens[a,19] <-  max(d)# probability leaving a patch
  death_sens[a,20] <- min(d)
  death_sens[a,21] <- max(d) - min(d)
}

death_max_plot <- ggplot(data = death_sens, mapping = aes(x = Death_Max, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Max death rate")+
  ylab("Landscape R0")+
  ylim(0,1)
death_max_plot  

death_min_plot <- ggplot(data = death_sens, mapping = aes(x = Death_Min, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'loess')+
  theme_classic()+
  xlab("Min death rate")+
  ylab("Landscape R0")+
  ylim(0,1)
death_min_plot  

death_range_plot <- ggplot(data = death_sens, mapping = aes(x = Death_range, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Death range")+
  ylab("Landscape R0")+
  ylim(0,1)
death_range_plot  

death_plots <- (death_max_plot + death_min_plot) / death_range_plot
death_plots
setwd(graphs)
ggsave(filename= "death_plots.png", plot = death_plots)


#############################################################################################
##################Simulate over range of recovery rates##################################
###########################################################################################
rm(list = ls())
main.wd <- "D:/gitrepos/MetaDisease"
graphs <- paste(main.wd,"/Graphs",sep = "")
setwd(main.wd)
source('Epimodel_funcs.R')
source('Metacomm_2patch.R')

#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)

birth <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
# if we change to a multi-season model, will need to add in some birth rate
# b <- c(0.6,0.5,0.4,0.3,0.2,0.1)/time #host birth rate
d <- c(0.024,
       NA,
       NA,
       0.0029,
       0.0144,
       0.04)
d[2] <- runif(n = 1, min = d[[4]], max = d[[1]])       
d[3] <- runif(n = 1, min = d[[4]], max = d[[2]]) 




#Determining dispersal rate
psi <- c(0.83, 0.62, 0.14,NA, 0.12,0.64)
psi[4] <- runif(n = 1, min = psi[5], max = psi[3])

#Can't define N_meta here, so will need to do so within simulation
# N_meta <- colSums(N)
# phi <- ifelse(N_meta > 0, psi/N_meta, 0)
# phi
# alpha <- #disease specific mortality
# species_chara <- data.frame(Species = Species,
#                             birth = birth, 
#                             death = d, 
#                             recovery = v, 
#                             dispersal = phi)
### Transmission coefficient
# beta should be higher for intraspecific transmission than interspecific
#inter-specific should have higher rate from more competent species
# Most abundant species should have highest interspecific trans rates
# For now, keeping transmission from 1 spp -> all other spp the same
#Ex: P Regilla will have the same transmission to A Boreas, T. Taricha... R Draytonii
# We are using a beta distribution b/c that is the best for the probability scale
# we need to assign probabilities to each species
trans_rate <- function(n = 1,l = 0, u, alpha,beta){
  trans <- rBeta.4P(n = n, l = l, u = u, alpha = alpha, beta = beta)
  return(trans)
}


# PREG
intra_PREG <- trans_rate(l = 0, u = 1, alpha = 4,beta = 2)
inter_PREG <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.5)

# TGRAN
intra_TGRAN <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.25)
inter_TGRAN <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4, beta = 2.5)

#TTOR
intra_TTOR <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4,beta = 2)
inter_TTOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 3, beta = 2.75)


#ABOR
intra_ABOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 4,beta = 2)
inter_ABOR <- trans_rate(l = 0, u = intra_ABOR, alpha = 2.5, beta = 3.0)


#RCAT
intra_RCAT <- trans_rate(l = 0, u = intra_ABOR, alpha = 4,beta = 2)
inter_RCAT <- trans_rate(l = 0, u = intra_RCAT, alpha = 2.0, beta = 3.25)

#RDRAY
intra_RDRAY <- trans_rate(l = 0, u = intra_RCAT, alpha = 4,beta = 2)
inter_RDRAY <- trans_rate(l = 0, u = intra_RDRAY, alpha = 1.5, beta = 3.5)

# meta-community characteristics
#Connectivity of patches
stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
go <- 1 - stay # probability individuals move
connect <- matrix(data = NA,
                  nrow = num_patches, 
                  ncol = num_patches)
for(i in 1:ncol(connect)){
  for(j in 1:nrow(connect)){
    connect[i,j] <- ifelse(i == j, stay, go)
  }
}



beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- if(i == 1 & j == 1){
      intra_PREG}else if(i != 1 & j == 1){
        inter_PREG} else if(i == 2 & j == 2){
          intra_TGRAN} else if(i != 2 & j == 2){
            inter_TGRAN} else if(i == 3 & j == 3){
              intra_TTOR} else if(i != 3 & j == 3){
                inter_TTOR} else if(i == 4 & j == 4){
                  intra_ABOR} else if(i != 4 & j == 4){
                    inter_ABOR} else if(i == 5 & j == 5){
                      intra_RCAT} else if(i != 5 & j ==5){
                        inter_RCAT} else if(i == 6 & j == 6){
                          intra_RDRAY} else if(i != 6 & j == 6){
                            inter_RDRAY}
  }
}
beta <- beta/90

# run simulation #
recover_sens <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 21))
colnames(recover_sens) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                          "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                          "BetaDiversity","Gamma_diversity",
                          "Beta_relative","LandscapeR0", "MetaCommID","recover_Max","recover_Min","recover_range")
for (a in 1:length(meta_comm_list)) {
  v <- runif(6)
  S <- ceiling(meta_comm_list[[a]][,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #value of susceptibles
  I <- meta_comm_list[[a]][,1:6] - S #value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- ifelse(N_meta > 0, psi/N_meta, 0) #get dispersal rate
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc R0
  b <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc b
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-connect[p,q]*S[p,i] + connect[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-connect[p,q]*I[p,i] + connect[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
        }
      }
    }
  }
  # Calculate variables
  recover_sens[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    recover_sens[a,i] <- sum(N[,i-1])
    recover_sens[a,(i+num_spp)] <- recover_sens[a,i]/recover_sens[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'w')
  recover_sens[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  recover_sens[a,15] <- sum(ifelse(recover_sens[a,1:6] > 0, 1,0)) #gamma diversity
  
  recover_sens[a,16] <- recover_sens[a,7]/recover_sens[a,8] #beta diversity / total abundance
  
  #calculate landscape R0
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = connect,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  recover_sens[a,17] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  recover_sens[a,18] <- a #metacommunity ID
  recover_sens[a,19] <-  max(v)# probability leaving of recovery
  recover_sens[a,20] <- min(v)
  recover_sens[a,21] <- max(v) - min(v)
}

recover_max_plot <- ggplot(data = recover_sens, mapping = aes(x = recover_Max, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Max recovery rate")+
  ylab("Landscape R0")+
  ylim(0,1)
recover_max_plot  

recover_min_plot <- ggplot(data = recover_sens, mapping = aes(x = recover_Min, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'loess')+
  theme_classic()+
  xlab("Min recovery rate")+
  ylab("Landscape R0")+
  ylim(0,1)
recover_min_plot  

recover_range_plot <- ggplot(data = recover_sens, mapping = aes(x = recover_range, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Recovery range")+
  ylab("Landscape R0")+
  ylim(0,1)
recover_range_plot  

recover_plots <- (recover_max_plot + recover_min_plot) / recover_range_plot
recover_plots

setwd(graphs)
ggsave(filename = "recoverygraphs.png", plot = recover_plots)



#############################################################################################
##################Simulate over range of dispersal rates##################################
###########################################################################################
rm(list = ls())
main.wd <- "D:/gitrepos/MetaDisease"
graphs <- paste(main.wd,"/Graphs",sep = "")
setwd(main.wd)
source('Epimodel_funcs.R')
source('Metacomm_2patch.R')

#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)

birth <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
# if we change to a multi-season model, will need to add in some birth rate
# b <- c(0.6,0.5,0.4,0.3,0.2,0.1)/time #host birth rate
d <- c(0.024,
       NA,
       NA,
       0.0029,
       0.0144,
       0.04)
d[2] <- runif(n = 1, min = d[[4]], max = d[[1]])       
d[3] <- runif(n = 1, min = d[[4]], max = d[[2]]) 


v <- c(NA,
       NA,
       NA,
       NA,
       0.012460223,
       NA)
# Citation for RCAT recovery = Daszak 2004
low_recovery <- runif(n = 4, min = 0, max = v[5])
v[1] <- min(low_recovery)
v[4] <- max(low_recovery)
mid_low_recovery <- low_recovery[low_recovery != v[1] & low_recovery != v[4]]
v[2] <- min(mid_low_recovery)
v[3] <- max(mid_low_recovery)
v[6] <- runif(n = 1, min = v[5], max = 0.1) #should find citation for max recovery rate
v

### Transmission coefficient
# beta should be higher for intraspecific transmission than interspecific
#inter-specific should have higher rate from more competent species
# Most abundant species should have highest interspecific trans rates
# For now, keeping transmission from 1 spp -> all other spp the same
#Ex: P Regilla will have the same transmission to A Boreas, T. Taricha... R Draytonii
# We are using a beta distribution b/c that is the best for the probability scale
# we need to assign probabilities to each species
trans_rate <- function(n = 1,l = 0, u, alpha,beta){
  trans <- rBeta.4P(n = n, l = l, u = u, alpha = alpha, beta = beta)
  return(trans)
}


# PREG
intra_PREG <- trans_rate(l = 0, u = 1, alpha = 4,beta = 2)
inter_PREG <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.5)

# TGRAN
intra_TGRAN <- trans_rate(l = 0, u = intra_PREG, alpha = 4, beta = 2.25)
inter_TGRAN <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4, beta = 2.5)

#TTOR
intra_TTOR <- trans_rate(l = 0, u = intra_TGRAN, alpha = 4,beta = 2)
inter_TTOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 3, beta = 2.75)


#ABOR
intra_ABOR <- trans_rate(l = 0, u = intra_TTOR, alpha = 4,beta = 2)
inter_ABOR <- trans_rate(l = 0, u = intra_ABOR, alpha = 2.5, beta = 3.0)


#RCAT
intra_RCAT <- trans_rate(l = 0, u = intra_ABOR, alpha = 4,beta = 2)
inter_RCAT <- trans_rate(l = 0, u = intra_RCAT, alpha = 2.0, beta = 3.25)

#RDRAY
intra_RDRAY <- trans_rate(l = 0, u = intra_RCAT, alpha = 4,beta = 2)
inter_RDRAY <- trans_rate(l = 0, u = intra_RDRAY, alpha = 1.5, beta = 3.5)

# meta-community characteristics
#Connectivity of patches
stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
go <- 1 - stay # probability individuals move
connect <- matrix(data = NA,
                  nrow = num_patches, 
                  ncol = num_patches)
for(i in 1:ncol(connect)){
  for(j in 1:nrow(connect)){
    connect[i,j] <- ifelse(i == j, stay, go)
  }
}



beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- if(i == 1 & j == 1){
      intra_PREG}else if(i != 1 & j == 1){
        inter_PREG} else if(i == 2 & j == 2){
          intra_TGRAN} else if(i != 2 & j == 2){
            inter_TGRAN} else if(i == 3 & j == 3){
              intra_TTOR} else if(i != 3 & j == 3){
                inter_TTOR} else if(i == 4 & j == 4){
                  intra_ABOR} else if(i != 4 & j == 4){
                    inter_ABOR} else if(i == 5 & j == 5){
                      intra_RCAT} else if(i != 5 & j ==5){
                        inter_RCAT} else if(i == 6 & j == 6){
                          intra_RDRAY} else if(i != 6 & j == 6){
                            inter_RDRAY}
  }
}
beta <- beta/90

# run simulation #
disperse_sens <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 21))
colnames(disperse_sens) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                            "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                            "BetaDiversity","Gamma_diversity",
                            "Beta_relative","LandscapeR0", "MetaCommID","disperse_Max","disperse_Min","disperse_range")
for (a in 1:length(meta_comm_list)) {
  S <- ceiling(meta_comm_list[[a]][,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #value of susceptibles
  I <- meta_comm_list[[a]][,1:6] - S #value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- runif(6) #get dispersal rate
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc R0
  b <- matrix(nrow = num_patches, ncol = num_spp) #create empty matrix, later will calc b
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-connect[p,q]*S[p,i] + connect[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-connect[p,q]*I[p,i] + connect[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
        }
      }
    }
  }
  # Calculate variables
  disperse_sens[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    disperse_sens[a,i] <- sum(N[,i-1])
    disperse_sens[a,(i+num_spp)] <- disperse_sens[a,i]/disperse_sens[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'w')
  disperse_sens[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  disperse_sens[a,15] <- sum(ifelse(disperse_sens[a,1:6] > 0, 1,0)) #gamma diversity
  
  disperse_sens[a,16] <- disperse_sens[a,7]/disperse_sens[a,8] #beta diversity / total abundance
  
  #calculate landscape R0
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = connect,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  disperse_sens[a,17] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  disperse_sens[a,18] <- a #metacommunity ID
  disperse_sens[a,19] <-  max(phi)# probability leaving of dispersey
  disperse_sens[a,20] <- min(phi)
  disperse_sens[a,21] <- max(phi) - min(phi)
}

disperse_max_plot <- ggplot(data = disperse_sens, mapping = aes(x = disperse_Max, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Max dispersal rate")+
  ylab("Landscape R0")+
  ylim(0,1)
disperse_max_plot  

disperse_min_plot <- ggplot(data = disperse_sens, mapping = aes(x = disperse_Min, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'loess')+
  theme_classic()+
  xlab("Min dispersal rate")+
  ylab("Landscape R0")+
  ylim(0,1)
disperse_min_plot  

disperse_range_plot <- ggplot(data = disperse_sens, mapping = aes(x = disperse_range, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Disperse range")+
  ylab("Landscape R0")+
  ylim(0,1)
disperse_range_plot  

disperse_plots <- (disperse_max_plot + disperse_min_plot)/disperse_range_plot
disperse_plots

setwd(graphs)
ggsave(filename = "disperseplots.png", plot = disperse_plots)
