#Description-----------------------------------------
#Build functions for EpiModel. Patch lvl
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

patch_build_R <- function(beta,sus,N,S,P) {
  R_list <- vector("list", length = P)
  for (p in 1:P) {
    Rmat <- matrix(data = 0, nrow = S, ncol = S)
    for (s in 1:S){
      for (i in 1:S){
        Rmat[s,i] <- (beta[s,i]*sus[p,i])/sum(N[p,])
      }
    }
    R_list[[p]] <- Rmat
  }
  return(R_list)
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
patch_build_B_ISO <- function(b, S, P) {
  B_list <- vector("list",length = P)
  for(p in 1:P){
    Bmat <- matrix(0,nrow = S, ncol = S)
    for(i in 1:S){
      for(j in 1:S){
        Bmat[i,j] <- ifelse(i == j,-b[[i]],0)
      }
    }
    B_list[[p]] <- Bmat
  }
  return(B_list)
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
patch_build_B_CON <- function(b, S, P, phi, c) {
  B_list <- vector("list",length = P)
  Cmat_nodiag <- c
  diag(Cmat_nodiag) <- 0
  Cmat_axis <- colSums(Cmat_nodiag)
  for(p in 1:P){
    Bmat <- matrix(0,nrow = S, ncol = S)
    for(i in 1:S){
      for(j in 1:S){
        Bmat[i,j] <- ifelse(i == j,-b[[i]]-phi*Cmat_axis[p],0)
      }
    }
    B_list[[p]] <- Bmat
  }
  return(B_list)
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
patch_R0_freq_ISO <- function(beta,sus,N,b,S,P){
  K_list <- vector("list", length = P)
  R <- patch_build_R(beta = beta,
                         sus = sus,
                         N = N,
                         S = S,
                         P = P)
  B <- patch_build_B_ISO(b = b,
                     S = S,
                     P = P)
  for (i in 1:P) {
    K_list[[i]] <- -R[[i]] %*% solve(-1*-B[[i]])
  }
  return(list(K_list = K_list, R_list = R, B_list = B))
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
patch_R0_freq_CON <- function(beta,sus,N,b,S,P,phi,c){
  K_list <- vector("list", length = P)
  R <- patch_build_R(beta = beta,
                     sus = sus,
                     N = N,
                     S = S,
                     P = P)
  B <- patch_build_B_CON(b = b,
                         S = S,
                         P = P,
                         phi = phi,
                         c = c)
  for (i in 1:P) {
    K_list[[i]] <- -R[[i]] %*% solve(-1*-B[[i]])
  }
  return(list(K_list = K_list, R_list = R, B_list = B))
}


# # # Examples of using the functions 
# Spp <- 2 #number of species
# Patches <- 2 #number of patches
# R0s <- matrix(data = rnorm(n = Spp*Patches, mean = 1, sd = 0.1), nrow = Patches, ncol = Spp)
# b <- rnorm(n = Spp, mean = 1, sd = 0.1)
# phi <- rnorm(n = 2, mean = 0.5, sd = 0.1)
# #transmission
# beta <- matrix(data = NA, nrow = Spp, ncol = Spp)
# for (i in 1:nrow(beta)) {
#   for (j in 1:ncol(beta)) {
#     beta[i,j] <- ifelse(i == j, rbeta(n = 1, shape1 = 2, shape2 = 10),rbeta(n = 1, shape1 = 1, shape2 = 1))
#   }
# }
# 
# 
# ### Freq dependent model
# N <- matrix(rep(NA,4), nrow = 2, ncol = 2)
# N[,1] <- runif(2, min = 10, max = 20)
# for(i in 1:2){
#   N[i,2] <- runif(1, min = 1, max = N[i,1])
# }
# Sus <- N*c(0.8,0.95)
# I <- N - Sus
# r_freq <- patch_build_R(beta = beta,
#                         sus = Sus,
#                         N = N,
#                         S = Spp,
#                         P = Patches)
# 
# b_freq <- patch_build_B(b = b,
#                         S = Spp,
#                         P = Patches)
# 
# # Not currently working. Need to probably adjust B matrix
# K_freq <- patch_R0_freq(beta = beta,
#                         sus = Sus,
#                         N = N,
#                         b = b,
#                         S = Spp,
#                         P = Patches)
# max(abs(eigen(K_freq[[1]][[1]])$values))
# eigen(K_freq[[1]][[1]])$values[1]
# View(K_freq[[1]][[1]])
