#Description-----------------------------------------
#Practicing preston's law in R
#  22 Nov 2024
#RCS

# Parameters-------------------------------------
#setting up meta-community parameters
num_patches <- 5 #number of patches in metacommunity
num_spp <- 6 #number of POSSIBLE spp in metacommunity

spp_occ <- data.frame(matrix(NA, nrow = num_patches, ncol = num_spp))
S <- c(0.95,0.9,0.75,0.50,0.20,0.1) #an array of probability values for the occurence of each spp

#determine occurence of spp 1
for(i in 1:nrow(spp_occ)){
  for(j in 1:ncol(spp_occ)){
    spp_occ[i,1] <- rbinom(n = 1, size = 1, prob = S[1])
  }
}

#determine occurrence of spp 2
for(i in 1:nrow(spp_occ)){
  for(j in 1:ncol(spp_occ)){
    spp_occ[i,2] <- ifelse(spp_occ[i,1] == 1, rbinom(n = 1, size = 1, prob = S[2]),0)
  }
}

#determine occurrence of spp 3
for(i in 1:nrow(spp_occ)){
  for(j in 1:ncol(spp_occ)){
    spp_occ[i,3] <- ifelse(spp_occ[i,2] == 1, rbinom(n = 1, size = 1, prob = S[3]),0)
  }
}

#determin occurrence of spp 4
for(i in 1:nrow(spp_occ)){
  for(j in 1:ncol(spp_occ)){
    spp_occ[i,4] <- ifelse(spp_occ[i,3] == 1, rbinom(n = 1, size = 1, prob = S[4]),0)
  }
}

#determine occupancy of spp 5
for(i in 1:nrow(spp_occ)){
  for(j in 1:ncol(spp_occ)){
    spp_occ[i,5] <- ifelse(spp_occ[i,4] == 1, rbinom(n = 1, size = 1, prob = S[5]),0)
  }
}

#determine occupancy of spp 6
for(i in 1:nrow(spp_occ)){
  for(j in 1:ncol(spp_occ)){
    spp_occ[i,6] <- ifelse(spp_occ[i,5] == 1, rbinom(n = 1, size = 1, prob = S[6]),0)
  }
}

M <- 3
Y0 <- 10
z <- 0.1 #constant.
P <- 1:num_spp #a vector. Preston's rankThe rank for each spp, which should correspond to the assigned abundance
K <- c(80,70,60,50,40,30,20,10) #this is just an example, but k is the abundance at each rank (i think?)
#this could later become some sort of rnorm(). EX rnorm(n = 1, mean = 10, sd = 1). This provides a starting abdunance for each spp if present
rank_abun <- Y0*exp(-(z*P-M)^2)
plot(rank_abun)

abun <- data.frame(K, rank_abun)
plot(abun$rank_abun, abun$K)  
