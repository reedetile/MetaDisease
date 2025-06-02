#Description-----------------------------------------
#Build functions for EpiModel
#  09 Jan 2025
#RCS

#Initialize -----------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load functions--------------------------------------


# Global Variables-------------------------------------


# Program Body------------------------------------------
# Define Fucntions------------------------------------
#######################################################
# FUNCTION: build_R_with_data_env
# Purpose: to calculate R matrix with data
#input:R0s, b, beta, S, P
#R0s = An Pxs Matrix with the species species specific R0 at each patch
# b = An PxS matrix. The relative loss rate of infecteds for each species in a patch
# beta = an SxS matrix. The transmission coefficient between species i and species j
# S = number of species
# P = number of patches
#output: fullR.an S*P by S*P matrix. 
#-----------------------------------------------------

build_R_with_data_env <- function(R0s,b,beta,S,P) {
  fullR <- matrix(data = 0, nrow = P*S, ncol = P*S)
  for (p in 1:P) {
    tR0 <- R0s[p,]
    Rmat <- matrix(data = rep(unlist(tR0),S), nrow = S, ncol = S)
    tb <- b[p,]
    bmat <- matrix(data = rep(unlist(tb),S), nrow = S, ncol = S)
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
# FUNCTION: build_R_with_data_freq
# Purpose: to calculate R matrix with data
#input:beta, I, N_TP, S, P
# beta = an SxS matrix. The transmission coefficient between species i and species s
# I = an PxS matrix. The number of infectious individuals of each species in each patch
# N = A P*S matrix. The total abundance of each species in each patch 
# S = number of species
# P = number of patches
#output: fullR.an S*P by S*P matrix. 
#-----------------------------------------------------

build_R_with_data_freq <- function(beta,I,N,S,P) {
  fullR <- matrix(data = 0, nrow = P*S, ncol = P*S)
  for (p in 1:P) {
    Rmat <- matrix(data = 0, nrow = S, ncol = S)
    for (s in 1:S){
      for (i in 1:S){
        Rmat[s,i] <- (beta[s,i]*I[p,i])/sum(N[p,])
      }
    }
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
# b = an array of length S. The relative loss rate for each spp
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
        #build the diagonal matrix
        (-1)*b - phi*Cmat_axis[p]
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
# Function: landscape_R0_env
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
landscape_R0_env <- function(R0_spps,Cmat,b, beta, phi,S,P){
  R <- build_R_with_data_env(R0s = R0_spps,
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


####################################################################
# Function: landscape_R0_freq
# Purpose: Calculate the landscape level R0 from the data for a frequency dependent system
#Input: beta, I, N_TP, Cmat, phi, b, S, P
# beta: an SxS matrix.
# Cmat = a PxP matrix. The colonization probabilities c_ij from patch j -> i
# phi = an array of length S. The dispersal rates for each species
# b = An SxP array. Relative loss rates of infecteds for each species in a patch
# beta = a transmission coefficient between species
# S = int. Number of species
# P = int. Number of patches
#Output: K. an S*P x S*P matrix. The matrix should be ordered by patch then by species. max eigenvalue of K = landscape R0.
# R, a component of K, and B, a component of K
landscape_R0_freq <- function(beta,I,N,Cmat,b, phi,S,P){
  R <- build_R_with_data_freq(beta = beta,
                         I = I,
                         N = N,
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


# # Examples of using the functions 
Spp <- 2 #number of species
Patches <- 2 #number of patches
R0s <- matrix(data = rnorm(n = Spp*Patches, mean = 1, sd = 0.1), nrow = Patches, ncol = Spp)
b <- rnorm(n = Spp, mean = 1, sd = 0.1)
phi <- rnorm(n = 2, mean = 0.5, sd = 0.1)
#transmission
beta <- matrix(data = NA, nrow = Spp, ncol = Spp)
for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    beta[i,j] <- ifelse(i == j, rbeta(n = 1, shape1 = 2, shape2 = 10),rbeta(n = 1, shape1 = 1, shape2 = 1))
  }
}
#connectivity
Cmat <- matrix(data = rnorm(n = Patches^2, mean = 0.5, sd = 0.1),
            nrow = Patches,
            ncol = Patches)
# 
# ### Environmental zoospore model
# r <- build_R_with_data_env(R0s = R0s,
#                        b = b,
#                        beta = beta,
#                        S = Spp,
#                        P = Patches)
# B <- build_B_with_data(Cmat = Cmat,
#                        phi = phi,
#                        b = b,
#                        S = Species,
#                        P = Patches)
# K2 <- landscape_R0_env(R0_spps = R0s,
#                    Cmat = Cmat,
#                    b = b,
#                    beta = beta,
#                    phi = phi,
#                    S = Spp,
#                    P = Patches)
# max(abs(eigen(K2[[1]])$values))
# eigen(K2[[1]])$values[1]
# View(K2[[1]])
# 
# ### Freq dependent model
# N <- matrix(rep(NA,4), nrow = 2, ncol = 2)
# N[,1] <- runif(2, min = 10, max = 20)
# for(i in 1:2){
#   N[i,2] <- runif(1, min = 1, max = N[i,1])
# }
# S <- N*c(0.8,0.95)
# I <- N - S
# r_freq <- build_R_with_data_freq(beta = beta,
#                                  I = I,
#                                  N = N,
#                                  S = Spp,
#                                  P = Patches)
# 
# # Not currently working. Need to probably adjust B matrix
# K_freq <- landscape_R0_freq(beta = beta,
#                             I = I,
#                             N = N,
#                             Cmat = Cmat,
#                             b = b,
#                             phi = phi,
#                             S = Spp,
#                             P = Patches)
# max(abs(eigen(K_freq[[1]])$values))
# eigen(K_freq[[1]])$values[1]
# View(K_freq[[1]])
