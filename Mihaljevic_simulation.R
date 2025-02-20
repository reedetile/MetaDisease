#Description-----------------------------------------
#Attempt to replicate results of 
#  18 Feb 2025
#RCS

#Initialize -----------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(betafunctions)
library(ggplot2)
# Load functions--------------------------------------


# Global Variables-------------------------------------
num_spp <- 6
K <- c(10,6,5,4,3,2) #rank abundance of each species

# psi
psi <- c(0.83, #PREG
       0.62, #TGRAN
       0.14, #TTOR
       NA, #ABOR
       0.12, #RCAT
       0.64) #RDRAY
psi[4] <- runif(n = 1, min = psi[5], max = psi[3])

comms_additive <- vector("list", length = 100)
for(c in 1:length(comms_additive)){
  N <- rep(NA,6)
  for(i in 1:length(N)){
    N[i] <- ifelse(rbinom(n=1,size=1,prob=psi[i]),K[i],0)
  }
  comms_additive[[c]] <- N
}

##skipping this for now
# comms_compens <- vector("list", length = 100)
# for(c in 1:length(comms)){
#   N <- rep(NA,6)
#   for(i in 1:length(N)){
#     N[i] <- ifelse(rbinom(n=1,size=1,prob=psi[i]),K[i],0)
#   }
#   comms[[c]] <- N
# }

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


### Running simulations

# Density dependent, additive models #
den_dep_add_r0 <- vector("list", length = length(comms_additive))
# density dependent model
for(c in 1:length(comms_additive)){
  N <- comms_additive[[c]]
  p <- N
  G <- matrix(NA, nrow = num_spp, ncol = num_spp)
  for(i in 1:num_spp){
    for (j in 1:num_spp) {
      G[i,j] <- (beta[i,j]*p[j])/(v[j]+d[j])
    }
  }
  R0 <- eigen(G)$values[1]
  den_dep_add_r0[[c]] <- R0
}

comms_df <- data.frame(matrix(unlist(comms_additive), nrow = length(comms_additive), ncol = num_spp, byrow = T))
den_dep_add_df <- cbind(comms_df, unlist(den_dep_add_r0))
den_dep_add_df$alpha <- rowSums(den_dep_add_df[,1:6] > 0)
den_dep_add_df$abund <- rowSums(den_dep_add_df[,1:6])
names(den_dep_add_df)[1:7] <- c("Spp1","Spp2","Spp3","Spp4","Spp5","Spp6","R0")
den_dep_add_Plot <- ggplot(data = den_dep_add_df,mapping = aes(x = alpha, y = R0))+
  geom_point()+
  geom_smooth(method = "lm")
den_dep_add_Plot
# Okay so this is a density dependent additive model
# It makes sense that in this case we see R0 go up with amplification

#Next step: Density dependent compensatory model






