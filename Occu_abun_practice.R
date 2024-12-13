#Description-----------------------------------------
#Practicing preston's law in R
#  22 Nov 2024
#RCS

# Load packages---------------------------------
library(vegan)
# Parameters-------------------------------------
#setting up meta-community parameters
set.seed(1234)
num_patches <- 5 #number of patches in metacommunity
num_spp <- 6 #number of POSSIBLE spp in metacommunity

meta_comm1 <- data.frame(matrix(NA, nrow = num_patches, ncol = num_spp))
S <- c(0.95,0.9,0.75,0.50,0.20,0.1) #an array of probability values for the occurence of each spp
#what if I over thought this, and I can just assign an occupancy probability?
K <- c(10,6,5,4,3,2) #this is just an example, but k is the abundance at each rank (i think?)
#this could later become some sort of rnorm(). EX rnorm(n = 1, mean = 10, sd = 1). This provides a starting abdunance for each spp if present

### Creating a meta-community###
#For right now I have it set where the presence / absence of a spp is a binomial variable with some
# Probability S[i]
#IF the species is present, then it is present at an abundance of K[i]. Right now this abundance is fixed.
# One idea for future development is to make abundance a distribution (could even just turn it into a poisson distribution?)
#Other idea/issue: if abundance of a "prior" spp is 0 than abundance of all following spp is 0. This is based on the idea that ecological communities are nested within one-another. However, instead of guaranteeing that abundance of next spp = 0, I wonder if it makes more sense to just make the probability of occ / abundance lower (but not set to zero.

# I think there's a lot that could be done here. A lot could be adapted. But this set of for loops is a decent start


# #determine occurence of spp 1
# for(i in 1:nrow(meta_comm1)){
#   for(j in 1:ncol(meta_comm1)){
#     meta_comm1[i,1] <- rbinom(n = 1, size = 1, prob = S[1])
#     meta_comm1[i,1] <- ifelse(meta_comm1[i,1] == 1, K[1],0)
#   }
# }
# 
# #determine occurrence of spp 2
# for(i in 1:nrow(meta_comm1)){
#   for(j in 1:ncol(meta_comm1)){
#     meta_comm1[i,2] <- ifelse(meta_comm1[i,1] > 0, rbinom(n = 1, size = 1, prob = S[2]),0)
#     meta_comm1[i,2] <- ifelse(meta_comm1[i,2] == 1, K[2],0)
#   }
# }
# 
# #determine occurrence of spp 3
# for(i in 1:nrow(meta_comm1)){
#   for(j in 1:ncol(meta_comm1)){
#     meta_comm1[i,3] <- ifelse(meta_comm1[i,2] > 0, rbinom(n = 1, size = 1, prob = S[3]),0)
#     meta_comm1[i,3] <- ifelse(meta_comm1[i,3] == 1, K[3],0)
#   }
# }
# 
# #determin occurrence of spp 4
# for(i in 1:nrow(meta_comm1)){
#   for(j in 1:ncol(meta_comm1)){
#     meta_comm1[i,4] <- ifelse(meta_comm1[i,3] > 0, rbinom(n = 1, size = 1, prob = S[4]),0)
#     meta_comm1[i,4] <- ifelse(meta_comm1[i,4] == 1, K[4],0)
#   }
# }
# 
# #determine occupancy of spp 5
# for(i in 1:nrow(meta_comm1)){
#   for(j in 1:ncol(meta_comm1)){
#     meta_comm1[i,5] <- ifelse(meta_comm1[i,4] > 0, rbinom(n = 1, size = 1, prob = S[5]),0)
#     meta_comm1[i,5] <- ifelse(meta_comm1[i,5] == 1, K[5],0)
#   }
# }
# 
# #determine occupancy of spp 6
# for(i in 1:nrow(meta_comm1)){
#   for(j in 1:ncol(meta_comm1)){
#     meta_comm1[i,6] <- ifelse(meta_comm1[i,5] > 0, rbinom(n = 1, size = 1, prob = S[6]),0)
#     meta_comm1[i,6] <- ifelse(meta_comm1[i,6] == 1, K[6],0)
#   }
# }

# I THINK this is a better way to determine abundance
for(i in 1:nrow(meta_comm1)){
  for(j in 1:ncol(meta_comm1)){
    #determine occurrence and abundance of spp 1
    meta_comm1[i,1] <- rbinom(n = 1, size = 1, prob = S[1])
    meta_comm1[i,1] <- ifelse(meta_comm1[i,1] == 1, K[1],0)
    
    #determine occurrence and abundance of spp2
    meta_comm1[i,2] <- ifelse(meta_comm1[i,1] > 0, rbinom(n = 1, size = 1, prob = S[2]),0)
    meta_comm1[i,2] <- ifelse(meta_comm1[i,2] == 1, K[2],0)
    
    #determine occurrence and abundance of spp 3
    meta_comm1[i,3] <- ifelse(meta_comm1[i,2] > 0, rbinom(n = 1, size = 1, prob = S[3]),0)
    meta_comm1[i,3] <- ifelse(meta_comm1[i,3] == 1, K[3],0)
    
    #determine occurrence and abundance of spp 4
    meta_comm1[i,4] <- ifelse(meta_comm1[i,3] > 0, rbinom(n = 1, size = 1, prob = S[4]),0)
    meta_comm1[i,4] <- ifelse(meta_comm1[i,4] == 1, K[4],0)
    
    #determine occurrence and abundance of spp 5
    meta_comm1[i,5] <- ifelse(meta_comm1[i,4] > 0, rbinom(n = 1, size = 1, prob = S[5]),0)
    meta_comm1[i,5] <- ifelse(meta_comm1[i,5] == 1, K[5],0)
    
    #determine occurrence and abundance of spp 6
    meta_comm1[i,6] <- ifelse(meta_comm1[i,5] > 0, rbinom(n = 1, size = 1, prob = S[6]),0)
    meta_comm1[i,6] <- ifelse(meta_comm1[i,6] == 1, K[6],0)
  }
}

meta_comm1$Patch <- c('Patch1','Patch2','Patch3','Patch4','Patch5')
colnames(meta_comm1) <- c('Spp1','Spp2','Spp3','spp4','spp5','Spp6','PatchID')
beta_diversity <- betadiver(meta_comm1[,1:6], method = 'w')
#plot(beta_diversity)


### Below is use of Preston's law. Not sure we'll actually be doing this, so I have commented it out for
##now
# M <- 3
# Y0 <- 10
# z <- 0.1 #constant.
# P <- 1:num_spp #a vector. Preston's rankThe rank for each spp, which should correspond to the assigned abundance
# K <- c(80,70,60,50,40,30,20,10) #this is just an example, but k is the abundance at each rank (i think?)
# #this could later become some sort of rnorm(). EX rnorm(n = 1, mean = 10, sd = 1). This provides a starting abdunance for each spp if present
# rank_abun <- Y0*exp(-(z*P-M)^2)
# plot(rank_abun)
# 
# abun <- data.frame(K, rank_abun)
# plot(abun$rank_abun, abun$K)  
