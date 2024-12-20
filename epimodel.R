#Description-----------------------------------------
#Rough draft script for a meta-community disease model
#  22 Nov 2024
#RCS

#Initialize -----------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)
source('Occu_abun_practice.R')
set.seed(1234)

# Define Fucntions------------------------------------
#######################################################
# FUNCTION: build_R_with_data
# Purpose: to calculate R matrix with data
#input:R0s, b, beta, S, P
#R0s = An Pxs Matrix with the species species specific R0 at each patch
# b = An PxS matrix. The relative loss rate of infecteds for each species in a patch
# beta = an SxS matrix. The transmission coefficient between species i and species j
# S = number of species
# P = number of patches
#output: fullR.an S*P by S*P matrix. 
#-----------------------------------------------------

build_R_with_data <- function(R0s,b,beta,S,P) {
  fullR <- matrix(data = 0, nrow = P*S, ncol = P*S)
  for (p in 1:P) {
    tR0 <- R0s[p,]
    Rmat <- matrix(data = rep(unlist(tR0),S), nrow = S, ncol = S)
    tb <- b[p,]
    bmat <- matrix(data = rep(unlist(tb),S), nrow = S, ncol = S)
    # I have removed the following code (in .py) but should ask mark if I need to add something for
    # Force of infection
    # tλ = λs[:, p]
    # λmat = np.repeat(tλ, S).reshape(S, S).T
    # λ_ratios = λmat / λmat.T    
    Rmat <- Rmat * bmat
    Rmat[is.nan(Rmat)] = 0
    Rmat[is.infinite(Rmat)] = 0
    Rmat <- as.matrix(Rmat)
    start <- ifelse(p == 1, p,(S*p)-(S-1))
    stop <- start + S - 1
    fullR[start:stop, start:stop] <- Rmat
  }
  
  return(fullR)
}

#######################################################
# FUNCTION: build_B_with_data
# Purpose: To calculate B matrix with data
# input:Cmat, As, psi, b, S, P
# Cmat = A pXp matrix. Colonization probability c_ij from patch j -> i
# As = relative patch areas (set at 1)
# phi = dispersal rate for each species.
# b = a SxP matrix. Relative loss rates of infected for each species
# S = int. number of species
# P = int. number of patches
#
#output:
#-----------------------------------------------------
build_B_with_data <- function(Cmat, phi, b, S, P) {
  Cmat_axis <- colSums(Cmat)
  diag_list <- vector("list",P^2)
  diag_list_names <- array(dim = c(sqrt(length(diag_list)), sqrt(length(diag_list))))
  for (p in 1:P) { #loop over columns
    for (j in 1:P) { #loop over rows
      diag_list_names[j,p] <- paste(p,j,sep="_")
    }
  }
  diag_list_names <- as.vector(diag_list_names)
  names(diag_list) <- diag_list_names
  for (p in 1:P) { #loop over columns
    for (j in 1:P) { #loop over rows
      tZ <- matrix(0,nrow = S, ncol = S)
      new_diag <- array(dim = S)
      new_diag <- if(p == j){
        (-1)*b[p,] - phi*Cmat_axis[p]
        } else{
          phi*Cmat[j,p]
          }
      diag(tZ) <- new_diag
      x <- paste(p,j,sep="_")
      diag_list[[x]] <- tZ 
    }
  }
  diag_matrix <- matrix(diag_list, nrow = P, ncol = P)
  tB <- list()
  for (j in 1:P) {
    tB[[j]] <- matrix(data = unlist(diag_matrix[j,]), nrow = S*P, ncol = S, byrow = T) #looks like this worked!
  }
  
  B <- do.call(cbind, tB[1:P]) #need to determine if this is what B should look like
  return(B) 
}

####################################################################
# Function: landscape_R0
# Purpose: Calculate the landscape level R0 from the data
#Input: R0s, Cmat, phi, b, S, P
# R0s_spps = SxP matrix. Each entry is the species-specific R0 for species s in patch p
# Cmat = a PxP matrix. The colonization probabilities c_ij from patch j -> i
# phi = an array of length S. The dispersal rates for each species
# b = An SxP array. Relative loss rates of infecteds for each species in a patch
# beta = a transmission coefficient between species
# S = int. Number of species
# P = int. Number of patches
#Output: K. an S*P x S*P matrix. The matrix should be ordered by patch then by species. max eigenvalue of K = landscape R0.
# R, a component of K, and B, a component of K
landscape_R0 <- function(R0_spps,Cmat,b, beta, phi,S,P){
  R <- build_R_with_data(R0s = R0_spps,
                         b = b,
                         beta = beta,
                         S = S,
                         P = P)
  B <- build_B_with_data(Cmat = Cmat,
                         phi = phi,
                         b = b,
                         S = S,
                         P = P)
  K <- R %*% solve(-1*B)
  return(list(K = K, R = R, B = B))
}

# Example of using the 3 functions
S <- 6 #number of species
P <- 5 #number of patches
R0s <- matrix(data = rnorm(n = S*P, mean = 1, sd = 0.1), nrow = P, ncol = S)
b <- matrix(data = rnorm(n = S*P, mean = 1, sd = 0.1), nrow = P, ncol = S)
phi <- rnorm(n = 6, mean = 0.5, sd = 0.1)
#transmission
beta <- matrix(data = NA, nrow = S, ncol = S)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- ifelse(i == j, rbeta(n = 1, shape1 = 2, shape2 = 10),rbeta(n = 1, shape1 = 1, shape2 = 1)) 
  }
}
#connectivity
c <- matrix(data = rnorm(n = P^2, mean = 0.5, sd = 0.1),
                 nrow = P, 
                 ncol = P)


r <- build_R_with_data(R0s = R0s,
                       b = b,
                       beta = beta,
                       S = S,
                       P = P)
B <-build_B_with_data(Cmat = c,
                      phi = phi,
                      b = b,
                      S = S,
                      P = P) 
K <- landscape_R0(R0_spps = R0s,
                  Cmat = c,
                  b = b, 
                  beta = beta, 
                  phi = phi,
                  S = S,
                  P = P)
eigen(K[[1]])$values[1]
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


# meta-community characteristics
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
#Connectivity of patches
#A <- rnorm(n = num_patches, mean = 0.5, sd = 0.1) #area ratios


# Program Body------------------------------------------
S <- as.matrix(S)
I <- as.matrix(I)
N <- as.matrix(N)
pop_list_Freq <- vector("list", length = time)
dilute_effect <- data.frame(matrix(data = NA, nrow = time, ncol = 2))

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
          
          #parameters to calc R0. Need help from mark to get this to run
          r0_species_patch[p,j] <-  ifelse(N[p,j] > 0, FI/loss_I, 0)
          r0_patch[p,q] <- ifelse(p == q, mean(r0_species_patch, na.rm = T),0) 
          DP <- matrix(data = ifelse(i == j, -(v[i]*d[i]) - phi[i], 0),
                       nrow = num_spp, ncol = num_spp)
          EIP <- matrix(data = ifelse(i == j,(A[p]/A[q])*c[p,q]*psi[i],0),
                        nrow = num_spp, ncol = num_spp)
          B[p,q] <- ifelse(p == q, DP, EIP)
          #parameters to calc R0
          # r0_species_patch[p,j] <-  ifelse(N[p,j] > 0, FI/loss_I, 0)
          # DP <- matrix(data = ifelse(i == j, -(v[i]*d[i]) - phi[i], 0),
          #              nrow = num_spp, ncol = num_spp)
          # EIP <- matrix(data = ifelse(i == j,c[p,q]*psi[i],0),
          #               nrow = num_spp, ncol = num_spp)
          # B[p,q] <- ifelse(p == q, DP, EIP)
          F_mat <- matrix(data = c(0,beta[i,j]/v[i],0,0), nrow = 2, ncol = 2)
          v_mat <- matrix(data = c(v[i],0,0,v[i]), nrow = 2, ncol = 2)
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
  beta_diversity <- betadiver(N, method = 'w')
  r0_landscape <- eigen(r0_patch * (-B^(-1)))[1]
  dilute_effect[i,1] <- mean(beta_diversity, na.rm =T)
  betadiver <- betadiver(N, method = 'w')
  #r0_landscape <- eigen(r0_species_patch *(-inv(B)))[1]
  F_mat <- matrix
  dilute_effect[i,1] <- betadiver
  dilute_effect[i,2] <- r0_landscape
}

pop_data_Freq <- lapply(pop_list_Freq, as.data.frame)
View(pop_data_Freq[[1]]) #example of dataframe at time 1
View(pop_data_Freq[[90]])  #example of dataframe at time 90

