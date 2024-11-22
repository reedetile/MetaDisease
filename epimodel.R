#Description-----------------------------------------
#Rough draft script for a meta-community disease model
#  22 Nov 2024
#RCS

#Initialize -----------------------------------------
library(vegan)

# Load functions--------------------------------------


# Parameters-------------------------------------
p <- #number of communities
sp <- #number of species POSSIBLE in metacommunity

#how to determine rank abundance of species
M <- 3
Y0 <- 10
q <- 0.1 #constant. Used q instead of z b/c z represents zoospores in this work
P <- 1:8 #a vector. The rank for each spp, which should correspond to the assigned abundance
rank_abun <- Y0*exp(-(q*P-M)^2)
plot(rank_abun)

for(i in length(P)){
  

N <- rnorm(
  #A matrix of p * sp. Total number of individuals in patch p of species sp. Should be S + I

S <- #A matrix of p * sp. number of susceptible individuals for each species at each patch
I <- #A matrix of p * sp. number of infectious individuals for each species at each patch




# Program Body------------------------------------------
delta_s <-b*N - d*S - beta*S*Z+v*I+phi*sum(-c*S+c*S*A[j]/A[p])
delta_I <- beta*S*Z-(v+d)*I+phi*sum(-c*I + c*I*A[j]/A[p])
delta_Z <- sum(lambda*I - gamma*Z)


