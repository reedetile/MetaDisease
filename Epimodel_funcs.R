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
  K <- R * solve(-1*B)
  return(list(K = K, R = R, B = B))
}

# Example of using the 3 functions
# S <- 6 #number of species
# P <- 5 #number of patches
# R0s <- matrix(data = rnorm(n = S*P, mean = 1, sd = 0.1), nrow = P, ncol = S)
# b <- matrix(data = rnorm(n = S*P, mean = 1, sd = 0.1), nrow = P, ncol = S)
# phi <- rnorm(n = 6, mean = 0.5, sd = 0.1)
# #transmission
# beta <- matrix(data = NA, nrow = S, ncol = S)
# for (i in 1:nrow(beta)) {
#   for (j in 1:ncol(beta)) {
#     beta[i,j] <- ifelse(i == j, rbeta(n = 1, shape1 = 2, shape2 = 10),rbeta(n = 1, shape1 = 1, shape2 = 1)) 
#   }
# }
# #connectivity
# c <- matrix(data = rnorm(n = P^2, mean = 0.5, sd = 0.1),
#             nrow = P, 
#             ncol = P)
# 
# 
# r <- build_R_with_data(R0s = R0s,
#                        b = b,
#                        beta = beta,
#                        S = S,
#                        P = P)
# B <- build_B_with_data(Cmat = c,
#                        phi = phi,
#                        b = b,
#                        S = S,
#                        P = P) 
# K2 <- landscape_R0(R0_spps = R0s,
#                    Cmat = c,
#                    b = b, 
#                    beta = beta, 
#                    phi = phi,
#                    S = S,
#                    P = P)
# eigen(K2[[1]])$values[1]
# View(K2[[1]])