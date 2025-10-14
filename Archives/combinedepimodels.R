#Description-----------------------------------------
#Script combining 2 and 5 patch epi models
#  03 Jun 2025
#RCS

#Initialize -----------------------------------------
library(vegan)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(truncdist)
library(mgcv)
library(patchwork)
set.seed(1234)

main.wd <- getwd()
graphs <- paste(getwd(),"/Graphs",sep = "")
# Load functions--------------------------------------
source('Epimodel_funcs.R')

#Let's create species characteristics first
#species characteristics
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
num_spp <- 6 #set num of possible spp

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
Size <- 1000 # size in m^2
K_rel <- c(10, 7.5, 7, 8, 4, 3)
K <- Size*K_rel

#Prop that will be susceptible
Infectious <- c(0.67, # Reeder 2012
                0.41, # Jost 2025
                0.4, # Jost 2025
                0.4, # Peralta Garcia 2018
                0.192, # Huss 2019
                0.84) #Adams 2022
Susceptible <- 1 - Infectious #prop of pop that will be infected

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

species_chara <- species_chara %>% mutate(beta_intra = (R0*(death+recovery))/(1-alpha))
#species_chara <- species_chara %>% mutate(beta_intra = (R0*(death+alpha)))

closeness <- 0.05 #based on mihaljevic 2014
beta <- matrix(nrow = num_spp, ncol = num_spp)

for (s in 1:num_spp) {
  for (i in 1:num_spp) {
    beta[s,i] <- if(s == i){
      species_chara[s,12]} else{
        closeness * ((species_chara[s,12]*species_chara[i,12])/2)}
  }
  
}


### 2 Patch System
# Set metacommunity params
meta_comm_list_2patch <- readRDS("metacomm_2Patch.RDS")
num_patches <- 2

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

result_2patch <- data.frame(matrix(data = NA, nrow = length(meta_comm_list_2patch), ncol = 17))
colnames(result_2patch) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                      "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                      "BetaDiversity","Gamma_diversity","LandscapeR0", "MetaCommID")

for (a in 1:length(meta_comm_list_2patch)) {
  S <- meta_comm_list_2patch[[a]][,1:6]*species_chara$Susceptible #value of susceptibles
  I <- meta_comm_list_2patch[[a]][,1:6]*species_chara$Infectious#value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- phi #get dispersal rate
  b <- v + d  # total loss rate
  alpha_div <- rowSums(N)
  
  # Calculate variables
  result_2patch[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    result_2patch[a,i] <- sum(N[,i-1])
    result_2patch[a,(i+num_spp)] <- result_2patch[a,i]/result_2patch[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'sor')
  result_2patch[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  result_2patch[a,15] <- sum(ifelse(result_2patch[a,1:6] > 0, 1,0)) #gamma diversity
  
  #calculate landscape R0
  r0_landscape <- landscape_R0_freq(beta = beta,
                                    I = I,
                                    N = N,
                                    Cmat = c,
                                    b = b,
                                    phi = phi,
                                    S = num_spp,
                                    P = num_patches)
  
  result_2patch[a,16] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  result_2patch[a,17] <- a #metacommunity ID
}

# Run GAMs
beta_2patch_GAM <- gam(LandscapeR0 ~ s(BetaDiversity, bs = "cr", k = 3), data = result_2patch)
gamma_2patch_GAM <- gam(LandscapeR0 ~ s(Gamma_diversity, bs = "cr", k = 3), data = result_2patch)
summary(beta_2patch_GAM)
summary(gamma_2patch_GAM)

AIC(beta_2patch_GAM, gamma_2patch_GAM)

# plot beta X R0
betaplot_2patch<- ggplot(data = result_2patch, mapping = aes(x = BetaDiversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cr", k = 3), colour = "black")+
  ylab("Landscape R0")+
  xlab("Beta Diversity")+
  theme_classic()
betaplot_2patch
gammaplot_2patch <- ggplot(data = result_2patch, mapping = aes(x = Gamma_diversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 3))+
  theme_classic()
gammaplot_2patch
#want to make sure this is still following a sat curve
ggplot(data = result_2patch, mapping = aes(x = Gamma_diversity, y = TotalAbundance))+
  geom_point()+
  geom_smooth(method ="loess")+
  theme_classic()

### 5 patch system
meta_comm_list_5patch <- readRDS("metacomm_5Patch.RDS")
num_patches <- 5
#Connectivity of patches
c <- matrix(data = NA,
            nrow = num_patches, 
            ncol = num_patches)
for(i in 1:ncol(c)){
  for(j in 1:nrow(c)){
    c[i,j] <- ifelse(i == j, stay, go)
  }
}

result_5patch <- data.frame(matrix(data = NA, nrow = length(meta_comm_list_5patch), ncol = 17))
colnames(result_5patch) <- c("TotalAbundance","PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY",
                      "PREG_rel","TGRAN_rel","TTOR_rel","ABOR_rel","RCAT_rel","RDRAY_rel",
                      "BetaDiversity","Gamma_diversity","LandscapeR0", "MetaCommID")

for (a in 1:length(meta_comm_list_5patch)) {
  S <- meta_comm_list_5patch[[a]][,1:6]*species_chara$Susceptible #value of susceptibles
  I <- meta_comm_list_5patch[[a]][,1:6]*species_chara$Infectious#value of infecteds
  N <- S+I #total pop of a patch
  S <- as.matrix(S)
  I <- as.matrix(I)
  N <- as.matrix(N)
  N_meta <- colSums(N)
  phi <- phi #get dispersal rate
  b <- v + d  # total loss rate
  alpha_div <- rowSums(N)
  
  # Calculate variables
  result_5patch[a,1] <- sum(N) #total abundance
  #abdunance of each species
  for (i in 2:(1+num_spp)) {
    result_5patch[a,i] <- sum(N[,i-1])
    result_5patch[a,(i+num_spp)] <- result_5patch[a,i]/result_5patch[a,1]
  }
  
  #add beta diversity
  beta_diversity <- betadiver(N, method = 'sor')
  result_5patch[a,14] <- mean(beta_diversity, na.rm =T) #beta diversity
  result_5patch[a,15] <- sum(ifelse(result_5patch[a,1:6] > 0, 1,0)) #gamma diversity
  
  #calculate landscape R0
  r0_landscape <- landscape_R0_freq(beta = beta,
                                    I = I,
                                    N = N,
                                    Cmat = c,
                                    b = b,
                                    phi = phi,
                                    S = num_spp,
                                    P = num_patches)
  
  result_5patch[a,16] <- max(abs(eigen(r0_landscape[[1]])$values)) #landscape R0
  result_5patch[a,17] <- a #metacommunity ID
}


# Run GAMs
beta_5patch_GAM <- gam(LandscapeR0 ~ s(BetaDiversity, bs = "cr", k = 3), data = result_5patch)
gamma_5patch_GAM <- gam(LandscapeR0 ~ s(Gamma_diversity, bs = "cr", k = 3), data = result_5patch)

summary(beta_5patch_GAM)
AIC(beta_5patch_GAM,gamma_5patch_GAM)



# Plot
betaplot_5patch <- ggplot(data = result_5patch, mapping = aes(x = BetaDiversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 3), colour = "black")+
  ylab("Landscape R0")+
  xlab("Beta Diversity")+
  theme_classic()
betaplot_5patch
gammaplot_5patch <- ggplot(data = result_5patch, mapping = aes(x = Gamma_diversity, y = LandscapeR0))+
  geom_point()+
  geom_smooth(method = "lm", colour = "black")+
  theme_classic()

# ggplot(data = result_5patch, mapping = aes(x = TotalAbundance, y = LandscapeR0))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   theme_classic()

# combine 2 patch and 5 patch plots
beta_plots <- betaplot_2patch + betaplot_5patch + plot_annotation(tag_levels = "a")
beta_plots
setwd(graphs)
ggsave(filename = "beta_plots.png", beta_plots)



# Lets create AIC tables
summary(gamma_2patch_GAM) #p < 2e-16, Rsq = 0.342, Deviance explained = 35.3%
summary(beta_2patch_GAM) # p < 2e-16, Rsq = 0.155, Deviance explained = 17.1%
summary(beta_5patch_GAM) # p < 2e-16, Rsq = 0.381, Deviance explauned = 39.2%
