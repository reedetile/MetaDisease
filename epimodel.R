#Description-----------------------------------------
#Rough draft script for a meta-community disease model
#  22 Nov 2024
#RCS

#Initialize -----------------------------------------
library(vegan)
library(matlib)
library(epimdr2)
source('Occu_abun_practice.R')
set.seed(1234)

# Define Fucntions------------------------------------
#######################################################
# FUNCTION: build_R_with_data
# Purpose:
#input:R0s, b, beta, S, P
#R0s = An SxP Matrix with the species species specific R0 at each patch
# b = An SxP matrix. The relative loss rate of infecteds for each species in a patch
# beta = an SxS matrix. The transmission coefficient between species i and species j
# S = number of species
# P = number of patches
#output: fullR.an S*P by S*P matrix. 
#-----------------------------------------------------
build_R_with_data <- function(R0s,b,beta,S,P) {
  fullR <- matrix(nrow = s*P, ncol = S*P)
  for (p in 1:P) {
    tR0 <- R0s[,p]
    Rmat <- matrix(data = rep(tR0,S), nrow = S, ncol = S)
    tb <- b[,p]
    bmat <- matrix(data = rep(tb,S), nrow = S, ncol = S)
    # I have removed the following code (in .py) but should ask mark if I need to add something for
    # Force of infection
    # tλ = λs[:, p]
    # λmat = np.repeat(tλ, S).reshape(S, S).T
    # λ_ratios = λmat / λmat.T    
    Rmat <- Rmat * bmat
    Rmat[is.nan(Rmat)] = 0
    Rmat[is.infinte(Rmat)] = 0
    start = p*S
    stop = start + S
    fullR[start:stop, start:stop] = Rmat
  }
  return(fullR)
}

#######################################################
# FUNCTION: build_B_with_data
# Purpose: Cmat, As, psi, b, S, P
# Cmat = 
#input:
#output:
#-----------------------------------------------------
build_B_with_data <- function(Cmat, As = 1, psi, b, S, P) {
  message("testing...build_B_with_data")
}


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

b <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
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
       0.03,
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
#A <- rnorm(n = num_patches, mean = 0.5, sd = 0.1) #area ratios


# Program Body------------------------------------------
S <- as.matrix(S)
I <- as.matrix(I)
N <- as.matrix(N)
pop_list_Freq <- vector("list", length = time)
dilute_effect <- data.frame(matrix(data = NA, nrow = time, ncol = 2))


#params for next generation matrix model
istates <- "I"
flist <- quote(beta * S * I/N)
Vlist <- quote(

for (t in 1:time) {
  delta_s_matrix <- matrix(nrow = num_patches, ncol = num_spp)
  delta_I_matrix <- matrix(nrow = num_patches, ncol = num_spp)
  r0_species_patch <- matrix(nrow = num_patches, ncol = num_spp)
  r0_patch <- matrix(nrow = num_patches, ncol = num_patches)
  B <- matrix(nrow = num_patches, ncol = num_patches)
  
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
          
<<<<<<< HEAD
          #parameters to calc R0. Need help from mark to get this to run
          r0_species_patch[p,j] <-  ifelse(N[p,j] > 0, FI/loss_I, 0)
          r0_patch[p,q] <- ifelse(p == q, mean(r0_species_patch, na.rm = T),0) 
          DP <- matrix(data = ifelse(i == j, -(v[i]*d[i]) - phi[i], 0),
                       nrow = num_spp, ncol = num_spp)
          EIP <- matrix(data = ifelse(i == j,(A[p]/A[q])*c[p,q]*psi[i],0),
                        nrow = num_spp, ncol = num_spp)
          B[p,q] <- ifelse(p == q, DP, EIP)
=======
          #parameters to calc R0
          # r0_species_patch[p,j] <-  ifelse(N[p,j] > 0, FI/loss_I, 0)
          # DP <- matrix(data = ifelse(i == j, -(v[i]*d[i]) - phi[i], 0),
          #              nrow = num_spp, ncol = num_spp)
          # EIP <- matrix(data = ifelse(i == j,c[p,q]*psi[i],0),
          #               nrow = num_spp, ncol = num_spp)
          # B[p,q] <- ifelse(p == q, DP, EIP)
          F_mat <- matrix(data = c(0,beta[i,j]/v[i],0,0), nrow = 2, ncol = 2)
          v_mat <- matrix(data = c(v[i],0,0,v[i]), nrow = 2, ncol = 2)
>>>>>>> bd4399cb88c855e0557883052ca933dac4489bf0
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
<<<<<<< HEAD
  beta_diversity <- betadiver(N, method = 'w')
  r0_landscape <- eigen(r0_patch * (-B^(-1)))[1]
  dilute_effect[i,1] <- mean(beta_diversity, na.rm =T)
=======
  betadiver <- betadiver(N, method = 'w')
  #r0_landscape <- eigen(r0_species_patch *(-inv(B)))[1]
  F_mat <- matrix
  dilute_effect[i,1] <- betadiver
>>>>>>> bd4399cb88c855e0557883052ca933dac4489bf0
  dilute_effect[i,2] <- r0_landscape
}

pop_data_Freq <- lapply(pop_list_Freq, as.data.frame)
View(pop_data_Freq[[1]]) #example of dataframe at time 1
View(pop_data_Freq[[90]])  #example of dataframe at time 90

