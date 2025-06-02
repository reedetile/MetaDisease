#Description-----------------------------------------
#description of script
#  25 Apr 2025
#RCS

#Initialize -----------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(truncdist)
set.seed(1234)

main.wd <- getwd()
graphs <- paste(getwd(),"/Graphs",sep = "")
# Load functions--------------------------------------
source('Epimodel_funcs.R')

# Global variables----------------------------------
meta_comm_list <- readRDS("metacomm_2Patch.RDS")
num_patches <- 2
num_spp <- 6
time <- 90 #how many "days" do I want in the season

### Species Characteristics
#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
# PREG = Pseudacris Regilla (Pacific tree frog)
# ABOR = Anaxyrus boreas (western toad)
# RCAT = Rana catesbeiana (American bullfrog)
# RDRAY = Rana draytonii (Califronia red legged frog)
# TTOR = Taricha torosa (California newt)
# TGRAN = Taricha granulosa (rough-skinned newt)

birth <- rep(0, times = 6) #currently have birth rate set to 0 b/c this is a within season model
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
       0.003086911,
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
phi <- runif(6) # need to check simulation for determining if diserpsal matters
# alpha <- #disease specific mortality
mass <- c(0.00035, #0.35 g, 0.00035 kg
          0.0103, # 10.3 g, 0.0103 kg
          0.0150, # 15.0 g, 0.0150 kg 
          0.0230, # 23.0 g, 0.0230 kg
          0.0380, # 38.0 g, 0.0380 kg
          0.0320) # 32.0 g, 0.0320 kg

alpha <- c(0.004320053, # PREG McMahon 2023
           NA,
           NA,
           0.007167381, # ABOR Carey 2006 
           0, # RCAT Daszak 2004
           0) # RDRAY Padgett Flohr 2008

alpha[2] <- runif(1, alpha[1], alpha[4])
alpha[3] <- runif(1, alpha[2], alpha[4])



R0 <- rtrunc(n = 6, spec = "gamma", a = 0, b = 2, shape = 2)
R0 <- sort(R0, decreasing = T)
## To plot what the distribution looks like
# x <- seq(0,2, by = 0.1)
# plot(x,dtrunc(x, spec = "gamma", a = 0, b = 2, shape =2))
K <- c(20,16,15,7,4,2) #this is just an example, but k is the abundance at each rank (i think?)
#Prop that will be susceptible
Infectious <- c(0.67, # Reeder 2012
         0.41, # Jost 2025
         0.4, # Jost 2025
         0.4, # Peralta Garcia 2018
         0.192, # Huss 2019
         0.84) #Adams 2022
Susceptible <- 1 - sus #prop of pop that will be infected

species_chara <- data.frame(Species = Species,
                            birth = birth, 
                            death = d, 
                            recovery = v, 
                            dispersal = phi,
                            mass = mass,
                            alpha = alpha,
                            R0 = R0,
                            K = K,
                            Susceptible = Susceptible,
                            Infectious = Infectious)
### Transmission coefficient
# beta should be higher for intraspecific transmission than interspecific
#inter-specific should have higher rate from more competent species
# Most abundant species should have highest interspecific trans rates
# For now, keeping transmission from 1 spp -> all other spp the same
#Ex: P Regilla will have the same transmission to A Boreas, T. Taricha... R Draytonii
# We are using a beta distribution b/c that is the best for the probability scale
# we need to assign probabilities to each species
species_chara <- species_chara %>% mutate(beta_intra = (R0*K*(death+alpha))/((Susceptible*K)*(Infectious*K)))
#species_chara <- species_chara %>% mutate(beta_intra = (R0*(death+alpha)))

closeness <- 0.5
beta <- matrix(nrow = num_spp, ncol = num_spp)

for (s in 1:num_spp) {
  for (i in 1:num_spp) {
    beta[s,i] <- if(s == i){
      species_chara[s,12]} else{
        closeness * ((species_chara[s,12]*species_chara[i,12])/2)}
  }
  
}

# x <- seq(0,1,length = 100)
# trans_rate <- function(n = 1,x = seq(0,1,length = 100),a,b){
#   trans <- rbeta(n = n, shape1 = a, shape2 = b)
#   plot(x = seq(0,1, length.out = 100), y = dbeta(x, shape1 = a, shape2 = b))
#   return(trans)
# }
# 
# 
# # PREG
# intra_PREG <- trans_rate(a = 4,b = 2)
# inter_PREG <- intra_PREG*0.85 #trans_rate(a = 4, b = 2.5)
# 
# plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2), ylab = "density", type = 'l', col = "red")
# lines(x=x, dbeta(x, shape1 = 4, shape2 = 2.5), col = "blue")
# 
# 
# # TGRAN
# intra_TGRAN <- trans_rate(a = 4, b = 2.25)
# inter_TGRAN <- intra_TGRAN*0.85 #trans_rate(a = 4, b = 2.5)
# 
# plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.25), ylab = "density", type = 'l', col = "red")
# lines(x=x, dbeta(x, shape1 = 4, shape2 = 2.5), col = "blue")
# 
# #TTOR
# intra_TTOR <- trans_rate(a = 3, b = 2.5)
# inter_TTOR <- intra_TTOR*0.85 #trans_rate(a = 3, b = 2.75)
# 
# plot(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.5), ylab = "density", type = 'l', col = "red")
# lines(x=x, dbeta(x, shape1 = 3, shape2 = 2.75), col = "blue")
# 
# #ABOR
# intra_ABOR <- trans_rate(a = 2.5, b = 2.75)
# inter_ABOR <- intra_ABOR*0.85 #trans_rate(a = 2.5, b = 3.0)
# 
# plot(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 2.75), ylab = "density", type = 'l', col = "red")
# lines(x=x, dbeta(x, shape1 = 2.5, shape2 = 3.0), col = "blue")
# 
# #RCAT
# intra_RCAT <- trans_rate(a = 2.0, b = 3.0)
# inter_RCAT <- intra_RCAT*0.85 #trans_rate(a = 2.0, b = 3.25)
# 
# plot(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.0), ylab = "density", type = 'l', col = "red")
# lines(x=x, dbeta(x, shape1 = 2.0, shape2 = 3.25), col = "blue")
# 
# #RDRAY
# intra_RDRAY <- trans_rate(a = 1.5, b = 3.25)
# inter_RDRAY <- intra_RDRAY*0.85#trans_rate(a = 1.5, b = 3.5)
# 
# plot(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.25), ylab = "density", type = 'l', col = "red")
# lines(x=x, dbeta(x, shape1 = 1.5, shape2 = 3.5), col = "blue")
# 
# #let's plot all the intraspecific trans rate just to see
# plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.0), ylab = "density", type = 'l', col = "black") #PREG
# lines(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.25), col = "red") #TGRAN
# lines(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.5), col = "blue") #TTOR
# lines(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 2.75), col = "cyan4") #ABOR
# lines(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.0), col = "blueviolet") #RCAT
# lines(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.25), col = "deeppink") #RDRAY
# 
# #let's plot all the interspecific trans rate just to see
# plot(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.5), ylab = "density", type = 'l', col = "black") #PREG
# lines(x = x,dbeta(x = x,shape1 = 4,shape2 = 2.5), col = "red") #TGRAN
# lines(x = x,dbeta(x = x,shape1 = 3,shape2 = 2.75), col = "blue") #TTOR
# lines(x = x,dbeta(x = x,shape1 = 2.5,shape2 = 3.0), col = "cyan4") #ABOR
# lines(x = x,dbeta(x = x,shape1 = 2.0,shape2 = 3.25), col = "blueviolet") #RCAT
# lines(x = x,dbeta(x = x,shape1 = 1.5,shape2 = 3.5), col = "deeppink") #RDRAY
# 
# 
# 
# 
# beta <- matrix(data = NA, nrow = num_spp, ncol = num_spp)
# for (i in 1:nrow(beta)) {
#   for (j in 1:ncol(beta)) {
#     beta[i,j] <- if(i == 1 & j == 1){
#       intra_PREG}else if(i != 1 & j == 1){
#         inter_PREG} else if(i == 2 & j == 2){
#           intra_TGRAN} else if(i != 2 & j == 2){
#             inter_TGRAN} else if(i == 3 & j == 3){
#               intra_TTOR} else if(i != 3 & j == 3){
#                 inter_TTOR} else if(i == 4 & j == 4){
#                   intra_ABOR} else if(i != 4 & j == 4){
#                     inter_ABOR} else if(i == 5 & j == 5){
#                       intra_RCAT} else if(i != 5 & j ==5){
#                         inter_RCAT} else if(i == 6 & j == 6){
#                           intra_RDRAY} else if(i != 6 & j == 6){
#                             inter_RDRAY}
#   }
# }
# beta #assumes that beta original is transmission over season. new beta is transmission
#over a single day (or time step)
#other way to measure beta
# beta_intra <- 2.47*10^(-2)*((species_chara$alpha/species_chara$death)
#                              +species_chara$recovery/species_chara$death)*species_chara$mass^(0.260)

#### meta-community characteristics
#Connectivity of patches
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

result <- data.frame(matrix(data = NA, nrow = length(meta_comm_list), ncol = 17))
colnames(result) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                       "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                       "BetaDiversity","Gamma_diversity","LandscapeR0", "MetaCommID")

for (a in 1:length(meta_comm_list)) {
  S <- meta_comm_list[[a]][,1:6]*species_chara$Susceptible #value of susceptibles
  I <- meta_comm_list[[a]][,1:6]*species_chara$Infectious#value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- phi #get dispersal rate
  b <- v + d  # total loss rate
  alpha_div <- rowSums(N)
  
  # Calculate variables
  result[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    result[a,i] <- sum(N[,i-1])
    result[a,(i+num_spp)] <- result[a,i]/result[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'sor')
  result[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  result[a,15] <- sum(ifelse(result[a,1:6] > 0, 1,0)) #gamma diversity
  
  #calculate landscape R0
  r0_landscape <- landscape_R0_freq(beta = beta,
                                    I = I,
                                    N = N,
                                    Cmat = c,
                                    b = b,
                                    phi = phi,
                                    S = num_spp,
                                    P = num_patches)
  
  result[a,16] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  result[a,17] <- a #metacommunity ID
}


# plot beta X R0
ggplot(data = result, mapping = aes(x = BetaDiversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth()+
  theme_classic()
ggplot(data = result, mapping = aes(x = Gamma_diversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth()+
  theme_classic()
ggplot(data = result, mapping = aes(x = TotalAbundance, y = LandscapeR0))+
  geom_point()+
  geom_smooth()+
  theme_classic()

#want to make sure this is still following a sat curve
ggplot(data = result, mapping = aes(x = Gamma_diversity, y = TotalAbundance))+
  geom_point()+
  geom_smooth()+
  theme_classic()


