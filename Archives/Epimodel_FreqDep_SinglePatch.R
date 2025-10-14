#Description-----------------------------------------
#Epidiomological model for 100 frequency dependent, single patch system
#  15 Jan 2025
#RCS

#Initialize -----------------------------------------
library(vegan)
library(dplyr)
library(ggplot2)
set.seed(1234)


# Program Body------------------------------------------

# Some notes
# This is simulation of 100 individual patches
# Goal: does alpha diversity correlate w. R_0P

# Step 1: Setup parameters
#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
num_spp <- 6
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


# Step 2: Form communities


S <- c(0.83, #PREG
       0.62, #TGRAN
       0.14, #TTOR
       NA, #ABOR
       0.12, #RCAT
       0.64) #RDRAY
S[4] <- runif(n = 1, min = S[5], max = S[3])
#an array of probability values for the occurence of each spp
#what if I over thought this, and I can just assign an occupancy probability?
K <- c(10,6,5,4,3,2) #this is just an example, but k is the abundance at each rank (i think?)

N <- 100 # number of communities
for(i in 1:length(S)){
  spp <- ifelse(rbinom(n = 100, size = 1, prob = S[i]) == 1,
                                     K[i],
                                     0)
  assign(paste0("spp_",i),spp)
}

comms <- data.frame(PREG = spp_1,
                    TGRAN = spp_2,
                    TTOR = spp_3,
                    ABOR = spp_4,
                    RCAT = spp_5,
                    RDRAY = spp_6)
# Cool, communities are made!
# Let's quickly check the nestedness of communities
nestedtemp(comm = comms)[7] # 25.33 which i believe indicates high nestedness

# Okay so let's start by measuring alpha, abundance and relative alpha
comms$alpha <- rowSums(ifelse(comms[,1:6] > 0, 1,0))
comms$abund <- rowSums(comms[,1:6])
comms <- comms %>% mutate(alpha_rel = case_when( abund > 0 ~ alpha/abund,
                                                 abund == 0 ~ 0))

# initializing for model
S <- ceiling(comms[,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #value of susceptibles
I <- comms[,1:6] - S #value of infecteds
N <- S+I #total pop of a patch
S <- as.matrix(S)
I <- as.matrix(I)
N <- as.matrix(N)
comms$R0 <- NA # adding in a R0 var which I will calc in the for loop

for(c in 1:nrow(comms)){
  G <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
  for(i in 1:num_spp){
    for(j in 1:num_spp){
      birth_rate <- birth[i]*N[c,i]
      death_s <- d[i] * S[c,i]
      loss_I <- (v[i]+d[i])*I[c,i]
      FI <- ifelse(N[c,j] > 0, beta[i,j]*S[c,i]*(I[c,j]/N[c,j]), 0)
      b <- v[i] + d[i]
      G[i,j] <- FI/b
    }
  }
  comms[c,10] <- max(abs(eigen(G)$values))
}

alpha_plot <- ggplot(data = comms, mapping = aes(x = alpha, y = R0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()
alpha_plot


alpha_rel_plot <- ggplot(data = comms, mapping = aes(x = alpha_rel, y = R0))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()
alpha_rel_plot
