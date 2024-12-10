#Description-----------------------------------------
#Rough draft script for a meta-community disease model
#  22 Nov 2024
#RCS

#Initialize -----------------------------------------
library(vegan)
source('Occu_abun_practice.R')

########################################
### Environmental Pool model
#######################################

# Parameters-------------------------------------
#starting values
#S <- ceiling(meta_comm1[,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #starting value of susceptibles
S <- meta_comm1[,1:6] #starting value of susceptibles
I <- meta_comm1[,1:6] - S #starting value of infecteds
N <- S+I
Z <- 1000 #starting amount of zoospores

#species characteristics
Species <- c("Spp1","Spp2","Spp3","Spp4","spp5","spp6")
b <- c(0.6,0.5,0.4,0.3,0.2,0.1) #host birth rate
# Justification = ...? NEED CITATION(S)
d <- c(0.06,0.05,0.04,0.03,0.02,0.01) #host death rate
# Justification = ...? NEED CITATION(S)
beta <- c(0.00013,0.00012,0.00011,0.00010,0.00009,0.00008) #transmission per species
# Justification = ...? NEED CITATION(S)
v <- c(0.4,0.5,0.6,0.7,0.8,0.9) #recovery rate
# Justification = ...? NEED CITATION(S)
phi <- c(0.09,0.08,0.07,0.06,0.05,0.04) #dispersal rate
# Justification = ...? NEED CITATION(S)
lambda <- c(3.0,2.5,2.0,1.5,1.0,0.5) #log(bd) load. losely based off figure 2 from wilber 2020
# Justification = Wilber 2020. Could use more citations
species_chara <- data.frame(Species = Species,
                            birth = b, 
                            death = d, 
                            trans = beta, 
                            recovery = v, 
                            dispersal = phi, 
                            shedding = lambda)

# meta-community characteristics
gamma <- 0.001 #zoospore decay rate
c <- matrix(data = rnorm(n = num_patches^2, mean = 0.5, sd = 0.1),
            nrow = num_patches, 
            ncol = num_patches)
#Connectivity of patches
A <- rnorm(n = num_patches, mean = 0.5, sd = 0.1) #area ratios
time <- 90 #how many "days" do I want in the season

# Program Body------------------------------------------
pop_list_EP <- vector("list", length = time)
for (t in 1:time) {
  delta_s <- b*N - d*S - beta*S*Z + v*I + phi*sum(-c*S + c*S) #for now have excluded area of patches
  #may want to add that back in though
  delta_I <- beta*S*Z-(v+d)*I+phi*sum(-c*I + c*I)#for now have excluded area of patches
  #may want to add that back in though
  delta_Z <- sum(lambda*I - gamma*Z)
  S <- S+delta_s
  I <- I+delta_I
  z <- Z+delta_Z
  N <- S + I
  t = t
  pop <- list(Susceptible = S, Infectious = I, Zoospores = Z, Total = N, Time = t)
  pop_list_EP[[t]] <- pop
}
  
pop_data_EP <- do.call(rbind, pop_list_EP)
View(pop_data_EP)  

########################################
### Frequency dependent model ###
#######################################
rm(list = ls())
source('Occu_abun_practice.R')

# Parameters-------------------------------------
#starting values
S <- ceiling(meta_comm1[,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #starting value of susceptibles
I <- meta_comm1[,1:6] - S #starting value of infecteds
N <- S+I
time <- 90 #how many "days" do I want in the season

#species characteristics
Species <- c("PREG","ABOR","RCAT","RDRAY","TTOR","TGRAN")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)
b <- c(0.6,0.5,0.4,0.3,0.2,0.1)/time #host birth rate
d <- c(0.06,0.05,0.04,0.03,0.02,0.01)/time #host death rate
v <- c(0.4,0.5,0.6,0.7,0.8,0.9)/time #recovery rate
phi <- c(0.09,0.08,0.07,0.06,0.05,0.04)/time #dispersal rate
species_chara <- data.frame(Species = Species,
                            birth = b, 
                            death = d, 
                            recovery = v, 
                            dispersal = phi)
### Transmission coefficient
beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- ifelse(i == j, rbeta(n = 1, shape1 = 2, shape2 = 10),rbeta(n = 1, shape1 = 1, shape2 = 1)) #need to make this more realistic
  }
}


# meta-community characteristics
c <- matrix(data = rnorm(n = num_patches^2, mean = 0.5, sd = 0.1),
            nrow = num_patches, 
            ncol = num_patches)
#Connectivity of patches
A <- rnorm(n = num_patches, mean = 0.5, sd = 0.1) #area ratios


# Program Body------------------------------------------
S <- as.matrix(S)
I <- as.matrix(I)
N <- as.matrix(N)
pop_list_Freq <- vector("list", length = time)

for (t in 1:time) {
  delta_s_matrix <- matrix(nrow = num_patches, ncol = num_spp)
  delta_I_matrix <- matrix(nrow = num_patches, ncol = num_spp)
  for (p in 1:num_patches) {
    for(q in 1:num_patches){
      for (i in 1:num_spp) {
        for (j in 1:num_spp) {
          #establish parameters for time t of species s
          connectivity_s <- phi[i]*sum(-c[p,q]*S[p,i] + c[q,p]*S[q,i])
          connectivity_I <- phi[i]*sum(-c[p,q]*I[p,i] + c[q,p]*I[q,i])
          birth <- b[i]*N[p,i]
          death_s <- d[i]*S[p,i]
          loss_I <- (v[i]*d[i])*I[p,i]
          FI <- ifelse(N[p,j] > 0, beta[i,j]*S[p,i]*(I[p,j]/N[p,j]), 0)
          
          #change in rate of susceptible + infectious individuals
          delta_s <- birth - death_s - FI + connectivity_s
          delta_s_matrix[p,i] <- delta_s
          
          delta_I <- FI - loss_I + connectivity_I
          delta_I_matrix[p,i] <- delta_I
        }
      }
    }
  }
   # New population dynamics
  S <- S + delta_s_matrix
  I <- I + delta_I_matrix
  N <- S + I
  # t = t
  pop <- list(Susceptible = S, Infectious = I, Total = N, Time = t)
  pop_list_Freq[[t]] <- pop
}
pop_data_Freq <- lapply(pop_list_Freq, as.data.frame)
View(pop_data_Freq[[1]]) #example of dataframe at time 1
View(pop_data_Freq[[90]])  #example of dataframe at time 90

