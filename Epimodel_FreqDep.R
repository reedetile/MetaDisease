#Description-----------------------------------------
#description of script
#  09 Jan 2025
#RCS

#Initialize -----------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
source('Occu_abun_practice.R')
set.seed(1234)

main.wd <- getwd()
graphs <- paste(getwd(),"/Graphs",sep = "")
# Load functions--------------------------------------
source('Epimodel_funcs.R')

# Program Body------------------------------------------

########################################
### Frequency dependent time model ###
#######################################
# Parameters-------------------------------------
#starting values
S <- ceiling(meta_comm_list[[1]][,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #starting value of susceptibles
I <- meta_comm_list[[1]][,1:6] - S #starting value of infecteds
N <- S+I
time <- 90 #how many "days" do I want in the season

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
d <- c(0.006,
       NA,
       NA,
       0.0007,
       0.004,
       0.0099)
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
       0.003086911,
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
N_meta <- colSums(N)
phi <- ifelse(N_meta > 0, psi/N_meta, 0)
phi
# alpha <- #disease specific mortality
species_chara <- data.frame(Species = Species,
                            birth = birth, 
                            death = d, 
                            recovery = v, 
                            dispersal = phi)
### Transmission coefficient
# beta should be higher for intraspecific transmission than interspecific
#inter-specific should have higher rate from more competent species
# Most abundant species should have highest interspecific trans rates
# For now, keeping transmission from 1 spp -> all other spp the same
#Ex: P Regilla will have the same transmission to A Boreas, T. Taricha... R Draytonii
# We are using a beta distribution b/c that is the best for the probability scale
# we need to assign probabilities to each species
x <- seq(0,1,length = 100)
trans_rate <- function(n = 1,x = seq(0,1,length = 100),a,b){
  trans <- rbeta(n = n, shape1 = a, shape2 = b)
  plot(x = seq(0,1, length.out = 100), y = dbeta(x, shape1 = a, shape2 = b))
  return(trans)
}


# PREG
intra_PREG <- trans_rate(a = 4,b = 2)
inter_PREG <- trans_rate(a = 4, b = 2.5)

plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 4, shape2 = 2.5), col = "blue")


# TGRAN
intra_TGRAN <- trans_rate(a = 4, b = 2.25)
inter_TGRAN <- trans_rate(a = 4, b = 2.5)

plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.25), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 4, shape2 = 2.5), col = "blue")

#TTOR
intra_TTOR <- trans_rate(a = 3, b = 2.5)
inter_TTOR <- trans_rate(a = 3, b = 2.75)

plot(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.5), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 3, shape2 = 2.75), col = "blue")

#ABOR
intra_ABOR <- trans_rate(a = 2.5, b = 2.75)
inter_ABOR <- trans_rate(a = 2.5, b = 3.0)

plot(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 2.75), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 2.5, shape2 = 3.0), col = "blue")

#RCAT
intra_RCAT <- trans_rate(a = 2.0, b = 3.0)
inter_RCAT <- trans_rate(a = 2.0, b = 3.25)

plot(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.0), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 2.0, shape2 = 3.25), col = "blue")

#RDRAY
intra_RDRAY <- trans_rate(a = 1.5, b = 3.25)
inter_RDRAY <- trans_rate(a = 1.5, b = 3.5)

plot(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.25), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 1.5, shape2 = 3.5), col = "blue")

#let's plot all the intraspecific trans rate just to see
plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.0), ylab = "density", type = 'l', col = "black") #PREG
lines(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.25), col = "red") #TGRAN
lines(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.5), col = "blue") #TTOR
lines(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 2.75), col = "cyan4") #ABOR
lines(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.0), col = "blueviolet") #RCAT
lines(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.25), col = "deeppink") #RDRAY

#let's plot all the interspecific trans rate just to see
plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.5), ylab = "density", type = 'l', col = "black") #PREG
lines(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.5), col = "red") #TGRAN
lines(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.75), col = "blue") #TTOR
lines(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 3.0), col = "cyan4") #ABOR
lines(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.25), col = "blueviolet") #RCAT
lines(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.5), col = "deeppink") #RDRAY




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
beta <- beta/time #assumes that beta original is transmission over season. new beta is transmission
#over a single day (or time step)

# meta-community characteristics
#Connectivity of patches
stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
go <- rbeta(n = 1, shape1 = 1.5, shape2 = 3.5) # probability individuals move
c <- matrix(data = NA,
            nrow = num_patches, 
            ncol = num_patches)
for(i in 1:ncol(c)){
  for(j in 1:nrow(c)){
    c[i,j] <- ifelse(i == j, stay, go)
  }
}

#A <- rnorm(n = num_patches, mean = 0.5, sd = 0.1) #area ratios


# Program Body------------------------------------------
S <- as.matrix(S)
I <- as.matrix(I)
N <- as.matrix(N)
pop_list_Freq <- vector("list", length = time)
result <- data.frame(matrix(data = NA, nrow = time, ncol = 3))
colnames(result) <- c("BetaDiversity","LandscapeR0", "Time")

for (t in 1:time) {
  delta_s_matrix <- matrix(nrow = num_patches, ncol = num_spp)
  delta_I_matrix <- matrix(nrow = num_patches, ncol = num_spp)
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp)
  b <- matrix(nrow = num_patches, ncol = num_spp)
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-c[p,q]*S[p,i] + c[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-c[p,q]*I[p,i] + c[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          delta_s <- birth_rate - death_s - FI + connectivity_s
          delta_s_matrix[p,i] <- delta_s
          
          delta_I <- FI - loss_I + connectivity_I
          delta_I_matrix[p,i] <- delta_I
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
          
          # Below could be used for equilibrium calculation where beta_ss and b are unknown
          # Prev_prop <- ifelse(N[q,i] > 0 & N[p,i] > 0, 
          #                     (I[q,i]/N[q,i])/(I[p,i]/N[p,i]), 0)
          # r0_numerator <- ifelse(N[p,i] > 0,
          #                        1 + (phi[i]/b[p,i])*sum(c[p,q]-c[q,p]*Prev_prop*(N[q,i]/N[p,i])),
          #                        0)
          # r0_denominator <- ifelse(N[p,i] > 0 & sum(N[p,]) > 0 & N[p,j] > 0,
          #                          (1 - (I[p,i]/N[p,i]))*
          #                            (N[p,i]/sum(N[p,]))*
          #                            sum((beta[i,j]/beta[i,i])*
          #                                  (N[p,j]/N[p,i])*
          #                                  ((I[p,j]/N[p,j])/(I[p,i]/N[p,i]))),0)
          # r0_species_patch[p,i] <- ifelse(r0_denominator > 0, 
          #                                 r0_numerator/r0_denominator,
          #                                 0)
        }
      }
    }
  }
  # New population dynamics
  S <- S + delta_s_matrix
  I <- I + delta_I_matrix
  N <- S + I
  t = t
  pop <- list(Susceptible = S, Infectious = I, Total = N, Time = t)
  pop_list_Freq[[t]] <- pop
  beta_diversity <- betadiver(N, method = 'w')
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = c,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  result[t,1] <- mean(beta_diversity, na.rm =T)
  result[t,2] <- eigen(r0_landscape[[1]])$values[1]
  result[t,3] <- t
}

# R0s_spps = SxP matrix. Each entry is the species-specific R0 for species s in patch p
# Cmat = a PxP matrix. The colonization probabilities c_ij from patch j -> i
# phi = an array of length S. The dispersal rates for each species
# b = An SxP array. Relative loss rates of infecteds for each species in a patch
# pop_data_Freq <- lapply(pop_list_Freq, as.data.frame)
# View(pop_data_Freq[[1]]) #example of dataframe at time 1
# View(pop_data_Freq[[90]])  #example of dataframe at time 90

### Frequency dependent multiple metacommunities model ###
source('Occu_abun_practice.R')

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
x <- seq(0,1,length = 100)
trans_rate <- function(n = 1,x = seq(0,1,length = 100),a,b){
  trans <- rbeta(n = n, shape1 = a, shape2 = b)
  plot(x = seq(0,1, length.out = 100), y = dbeta(x, shape1 = a, shape2 = b))
  return(trans)
}


# PREG
intra_PREG <- trans_rate(a = 4,b = 2)
inter_PREG <- trans_rate(a = 4, b = 2.5)

plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 4, shape2 = 2.5), col = "blue")


# TGRAN
intra_TGRAN <- trans_rate(a = 4, b = 2.25)
inter_TGRAN <- trans_rate(a = 4, b = 2.5)

plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.25), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 4, shape2 = 2.5), col = "blue")

#TTOR
intra_TTOR <- trans_rate(a = 3, b = 2.5)
inter_TTOR <- trans_rate(a = 3, b = 2.75)

plot(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.5), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 3, shape2 = 2.75), col = "blue")

#ABOR
intra_ABOR <- trans_rate(a = 2.5, b = 2.75)
inter_ABOR <- trans_rate(a = 2.5, b = 3.0)

plot(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 2.75), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 2.5, shape2 = 3.0), col = "blue")

#RCAT
intra_RCAT <- trans_rate(a = 2.0, b = 3.0)
inter_RCAT <- trans_rate(a = 2.0, b = 3.25)

plot(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.0), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 2.0, shape2 = 3.25), col = "blue")

#RDRAY
intra_RDRAY <- trans_rate(a = 1.5, b = 3.25)
inter_RDRAY <- trans_rate(a = 1.5, b = 3.5)

plot(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.25), ylab = "density", type = 'l', col = "red")
lines(x=x, dbeta(x, shape1 = 1.5, shape2 = 3.5), col = "blue")

#let's plot all the intraspecific trans rate just to see
plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.0), ylab = "density", type = 'l', col = "black") #PREG
lines(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.25), col = "red") #TGRAN
lines(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.5), col = "blue") #TTOR
lines(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 2.75), col = "cyan4") #ABOR
lines(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.0), col = "blueviolet") #RCAT
lines(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.25), col = "deeppink") #RDRAY

#let's plot all the interspecific trans rate just to see
plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.5), ylab = "density", type = 'l', col = "black") #PREG
lines(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.5), col = "red") #TGRAN
lines(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.75), col = "blue") #TTOR
lines(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 3.0), col = "cyan4") #ABOR
lines(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.25), col = "blueviolet") #RCAT
lines(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.5), col = "deeppink") #RDRAY




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
# meta-community characteristics
#Connectivity of patches
stay <- rbeta(n = 1, shape1 = 4, shape2 = 2) #probability individuals stay in a patch?
go <- rbeta(n = 1, shape1 = 1.5, shape2 = 3.5) # probability individuals move
c <- matrix(data = NA,
            nrow = num_patches, 
            ncol = num_patches)
for(i in 1:ncol(c)){
  for(j in 1:nrow(c)){
    c[i,j] <- ifelse(i == j, stay, go)
  }
}

# Simulation over metacommunities
result2 <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 18))
colnames(result2) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                       "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                       "BetaDiversity","Gamma_diversity",
                       "Beta_relative","LandscapeR0", "MetaCommID")
for (a in 1:length(meta_comm_list)) {
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
          connectivity_s <- phi[i]*sum(-c[p,q]*S[p,i] + c[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-c[p,q]*I[p,i] + c[q,p]*I[q,i])
          birth_rate <- birth[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]+d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          
          #calculate R0 for each species at each patch
          b[p,i] <- v[i] + d[i] #gets loss rate
          r0_species_patch[p,i] <- beta[i,i]/b[p,i]
          
          # Below could be used for equilibrium calculation where beta_ss and b are unknown
          # Prev_prop <- ifelse(N[q,i] > 0 & N[p,i] > 0, 
          #                     (I[q,i]/N[q,i])/(I[p,i]/N[p,i]), 0)
          # r0_numerator <- ifelse(N[p,i] > 0,
          #                        1 + (phi[i]/b[p,i])*sum(c[p,q]-c[q,p]*Prev_prop*(N[q,i]/N[p,i])),
          #                        0)
          # r0_denominator <- ifelse(N[p,i] > 0 & sum(N[p,]) > 0 & N[p,j] > 0,
          #                          (1 - (I[p,i]/N[p,i]))*
          #                            (N[p,i]/sum(N[p,]))*
          #                            sum((beta[i,j]/beta[i,i])*
          #                                  (N[p,j]/N[p,i])*
          #                                  ((I[p,j]/N[p,j])/(I[p,i]/N[p,i]))),0)
          # r0_species_patch[p,i] <- ifelse(r0_denominator > 0, 
          #                                 r0_numerator/r0_denominator,
          #                                 0)
        }
      }
    }
  }
  # Calculate variables
  result2[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    result2[a,i] <- sum(N[,i-1])
    result2[a,(i+num_spp)] <- result2[a,i]/result2[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'w')
  result2[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  result2[a,15] <- sum(ifelse(result2[a,1:6] > 0, 1,0)) #gamma diversity
  
  result2[a,16] <- result2[a,7]/result2[a,8] #beta diversity / total abundance
  
  #calculate landscape R0
  r0_landscape <- landscape_R0(R0_spps = r0_species_patch,
                               Cmat = c,
                               b = b,
                               beta = beta,
                               phi = phi,
                               S = num_spp,
                               P = num_patches)
  
  result2[a,17] <- eigen(r0_landscape[[1]])$values[1] #landscape R0
  result2[a,18] <- a #metacommunity ID
}

beta_plot <- ggplot(data = result2, aes(x = BetaDiversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Beta Diversity")+
  ylab("Landscape R0")
beta_plot

abun_R0_plot <- ggplot(data = result2, aes(x = TotalAbundance, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Abundance")+
  ylab("Landscape R0")
abun_R0_plot

beta_relative_plot <- ggplot(data = result2, aes(x = Beta_relative, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Relative Beta Diversity")+
  ylab("Landscape R0")
beta_relative_plot

gamma_plot <- ggplot(data = result2, aes(x = Gamma_diversity, y = LandscapeR0))+
  geom_point(position = position_jitter(h=0.15,w=0.15))+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab("Gamma Diversity")+
  ylab("Landscape R0")
gamma_plot

meta_comm_effects <- (beta_plot + abun_R0_plot)/(beta_relative_plot + gamma_plot)
meta_comm_effects

#Species plots
#total abundance
Total_abund <- result2 %>% select(c(PREG:RDRAY, LandscapeR0))
Total_abund_longer <- pivot_longer(Total_abund, !LandscapeR0, names_to = "Species", values_to = "Abund")
Total_abund_plot <- ggplot(data = Total_abund_longer, mapping = aes(x = Abund, y = LandscapeR0, colour = Species))+
  geom_point()+
  geom_smooth(method = 'lm')
Total_abund_plot #all species

#let's do plots by species
#PREG
Total_PREG <- Total_abund_longer %>% filter(Species == 'PREG')
PREG_Total_plot <- ggplot(Total_PREG, mapping = aes(x = Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "PREG")

#ABOR
Total_ABOR <- Total_abund_longer %>% filter(Species == 'ABOR')
ABOR_Total_plot <- ggplot(Total_ABOR, mapping = aes(x = Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "ABOR")

#RCAT
Total_RCAT <- Total_abund_longer %>% filter(Species == 'RCAT')
RCAT_Total_plot <- ggplot(Total_RCAT, mapping = aes(x = Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "RCAT")

#RDRAY
Total_RDRAY <- Total_abund_longer %>% filter(Species == 'RDRAY')
RDRAY_Total_plot <- ggplot(Total_RDRAY, mapping = aes(x = Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "RDRAY")

#TGRAN
Total_TGRAN <- Total_abund_longer %>% filter(Species == 'TGRAN')
TGRAN_Total_plot <- ggplot(Total_TGRAN, mapping = aes(x = Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "TGRAN")

#TTOR
Total_TTOR <- Total_abund_longer %>% filter(Species == 'TTOR')
TTOR_Total_plot <- ggplot(Total_TTOR, mapping = aes(x = Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "TTOR")

Species_total_plots <- (ABOR_Total_plot + PREG_Total_plot + RCAT_Total_plot) / 
  (RDRAY_Total_plot + TGRAN_Total_plot + TTOR_Total_plot)
Species_total_plots



#relative abundance
Rel_abund <- result2 %>% select(c(PREG_rel:RDRAY_rel, LandscapeR0))
Rel_abund_longer <- pivot_longer(Rel_abund, !LandscapeR0, names_to = "Species", values_to = "Rel_Abund")
Rel_abund_plot <- ggplot(data = Rel_abund_longer, mapping = aes(x = Rel_Abund, y = LandscapeR0, colour = Species))+
  geom_point()+
  geom_smooth(method = 'lm')
Rel_abund_plot

#let's do plots by species
#PREG
Rel_PREG <- Rel_abund_longer %>% filter(Species == 'PREG_rel')
PREG_Rel_plot <- ggplot(Rel_PREG, mapping = aes(x = Rel_Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "PREG")

#ABOR
Rel_ABOR <- Rel_abund_longer %>% filter(Species == 'ABOR_rel')
ABOR_Rel_plot <- ggplot(Rel_ABOR, mapping = aes(x = Rel_Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "ABOR")

#RCAT
Rel_RCAT <- Rel_abund_longer %>% filter(Species == 'RCAT_rel')
RCAT_Rel_plot <- ggplot(Rel_RCAT, mapping = aes(x = Rel_Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "RCAT")

#RDRAY
Rel_RDRAY <- Rel_abund_longer %>% filter(Species == 'RDRAY_rel')
RDRAY_Rel_plot <- ggplot(Rel_RDRAY, mapping = aes(x = Rel_Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "RDRAY")

#TGRAN
Rel_TGRAN <- Rel_abund_longer %>% filter(Species == 'TGRAN_rel')
TGRAN_Rel_plot <- ggplot(Rel_TGRAN, mapping = aes(x = Rel_Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "TGRAN")

#TTOR
Rel_TTOR <- Rel_abund_longer %>% filter(Species == 'TTOR_rel')
TTOR_Rel_plot <- ggplot(Rel_TTOR, mapping = aes(x = Rel_Abund, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "TTOR")

Species_Rel_plots <- (ABOR_Rel_plot + PREG_Rel_plot + RCAT_Rel_plot) / 
  (RDRAY_Rel_plot + TGRAN_Rel_plot + TTOR_Rel_plot)
Species_Rel_plots

# Consolidate call for plots
beta_plot
abun_R0_plot
beta_relative_plot
gamma_plot
meta_comm_effects

#Lets consolidate the total + relative abundance plots for each species
ABOR_plots <- ABOR_Total_plot + ABOR_Rel_plot
PREG_plots <- PREG_Total_plot + PREG_Rel_plot
RCAT_plots <- RCAT_Total_plot + RCAT_Rel_plot
RDRAY_plots <- RDRAY_Total_plot + RDRAY_Rel_plot 
TTOR_plots <- TTOR_Total_plot + TTOR_Rel_plot
TGRAN_plots <- TGRAN_Total_plot + TGRAN_Rel_plot

Species_total_plots
Species_Rel_plots
#Save all your plots
setwd(graphs)
ggsave(filename = "beta_plot.png", plot = beta_plot)
ggsave(filename = "abun_R0_plot.png", plot = abun_R0_plot)
ggsave(filename = "beta_relative.png", plot = beta_relative_plot)
ggsave(filename = "gamma_plot.png", plot = gamma_plot)
ggsave(filename = "meta_comm_effects.png", plot = meta_comm_effects)
ggsave(filename = "ABOR_plots.png", plot = ABOR_plots)
ggsave(filename = "PREG_plots.png", plot = PREG_plots)
ggsave(filename = "RCAT_plots.png", plot = RCAT_plots)
ggsave(filename = "RDRAY_plots.png", plot = RDRAY_plots)
ggsave(filename = "TTOR_plots.png", plot = TTOR_plots)
ggsave(filename = "TGRAN_plots.png", plot = TGRAN_plots)
ggsave(filename = "Species_total_plots.png", Species_total_plots)
ggsave(filename = "Species_Rel_plots.png", Species_Rel_plots)
