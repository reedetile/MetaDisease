#Description-----------------------------------------
# Creating highly controlled 2 patch ex
#  22 October 2025
#RCS

# Load packages---------------------------------
library(vegan)
library(ggplot2)
library(tidyr)
library(patchwork)
library(truncdist)

#setup -------------------------------------------
repo <- "C:/Users/rcscott/MetaDisease"
graphs <- paste(repo,"/Graphs",sep="")

# Load function ----------------------------------
setwd(repo)
source('Epimodel_funcs.R')
source("Epimodel_funcs_patch_lvl.R")

# Parameters-------------------------------------
set.seed(1234)
#setup functions for simulations
vary_disp_R0L <- function(meta_comm_df,species_chara,beta,num_spp,num_patches){
  set.seed(1234)
  stay <- 0.01
  go <- 1-stay
  c <- matrix(data = NA,
              nrow = num_patches, 
              ncol = num_patches)
  for(i in 1:ncol(c)){
    for(m in 1:nrow(c)){
      c[i,m] <- ifelse(i == m, stay, go)
    }
  }
  max_disp <- seq(0,0.1,by = 0.001)
  result_disp <- data.frame(matrix(data = NA, nrow = 0, ncol = 2))
  
  for(z in 1:length(max_disp)){
    phi <- runif(6, max = max_disp[[z]])
    result <- data.frame(matrix(data = NA, nrow = 1, ncol = 2))
    colnames(result) <- c("LandscapeR0","max_disp")
    
    Sus <- meta_comm_df*species_chara$Susceptible #value of susceptibles
    I <- meta_comm_df*species_chara$Infectious#value of infecteds
    N <- Sus+I #total pop of a patch
    Sus <- as.matrix(Sus)
    I <- as.matrix(I)
    N <- as.matrix(N)
    N_meta <- colSums(N)
    b <- v + d  # total loss rate
    alpha_div <- rowSums(N)
    
    # Calculate variables
    
    #calculate landscape R0
    r0_landscape <- landscape_R0_freq(beta = beta,
                                      sus = Sus,
                                      N = N,
                                      Cmat = c,
                                      b = b,
                                      phi = phi,
                                      S = num_spp,
                                      P = num_patches)
    
    result[1,1] <- max(abs(eigen(r0_landscape[[1]])$values))
    result[1,2] <- max(phi)
    result_disp <- rbind(result_disp, result)
  }
  return(result_disp)
}

vary_dispR0P_iso <- function(meta_comm_df,species_chara,beta,num_spp,num_patches){
  set.seed(1234)
  stay <- 0.01
  go <- 1-stay
  c <- matrix(data = NA,
              nrow = num_patches, 
              ncol = num_patches)
  for(i in 1:ncol(c)){
    for(m in 1:nrow(c)){
      c[i,m] <- ifelse(i == m, stay, go)
    }
  }
  max_disp <- seq(0,0.1,by = 0.001)
  result_disp <- data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  
  for(z in 1:length(max_disp)){
    phi <- runif(6, max = max_disp[[z]])
    result <- data.frame(matrix(data = NA, nrow = 1, ncol = 3))
    colnames(result) <- c("PatchR0","Dispersal","Patch")
    
    Sus <- meta_comm_df*species_chara$Susceptible #value of susceptibles
    I <- meta_comm_df*species_chara$Infectious#value of infecteds
    N <- Sus+I #total pop of a patch
    Sus <- as.matrix(Sus)
    I <- as.matrix(I)
    N <- as.matrix(N)
    N_meta <- colSums(N)
    b <- species_chara$death + 
      species_chara$recovery + 
      species_chara$alpha  # total loss rate
    alpha_div <- rowSums(N)
    
    # Calculate variables
    
    #calculate landscape R0
    r0_patch <- patch_R0_freq_ISO(beta = beta,
                              sus = Sus,
                              N = N,
                              b = b,
                              S = num_spp,
                              P = num_patches)
    for (p in 1:num_patches) {
      result[p,1] <- max(abs(eigen(r0_patch[[1]][[p]])$values))
      result[p,2] <- max(phi)
      result[p,3] <- p
    }
    result_disp <- rbind(result_disp, result)
  }
  result_disp$Patch <- as.factor(result_disp$Patch)
  return(result_disp)
}


vary_dispR0P_connected <- function(meta_comm_df,species_chara,beta,num_spp,num_patches){
  set.seed(1234)
  stay <- 0.01
  go <- 1-stay
  c <- matrix(data = NA,
              nrow = num_patches, 
              ncol = num_patches)
  for(i in 1:ncol(c)){
    for(m in 1:nrow(c)){
      c[i,m] <- ifelse(i == m, stay, go)
    }
  }
  max_disp <- seq(0,0.1,by = 0.001)
  result_disp <- data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  
  for(z in 1:length(max_disp)){
    phi <- runif(6, max = max_disp[[z]])
    result <- data.frame(matrix(data = NA, nrow = 1, ncol = 3))
    colnames(result) <- c("PatchR0","Dispersal","Patch")
    
    Sus <- meta_comm_df*species_chara$Susceptible #value of susceptibles
    I <- meta_comm_df*species_chara$Infectious#value of infecteds
    N <- Sus+I #total pop of a patch
    Sus <- as.matrix(Sus)
    I <- as.matrix(I)
    N <- as.matrix(N)
    N_meta <- colSums(N)
    b <- species_chara$death + 
      species_chara$recovery + 
      species_chara$alpha  # total loss rate
    alpha_div <- rowSums(N)
    
    # Calculate variables
    
    #calculate landscape R0
    r0_patch <- patch_R0_freq_CON(beta = beta,
                              sus = Sus,
                              N = N,
                              b = b,
                              S = num_spp,
                              P = num_patches,
                              phi = phi,
                              c = c)
    for (p in 1:num_patches) {
      result[p,1] <- max(abs(eigen(r0_patch[[1]][[p]])$values))
      result[p,2] <- max(phi)
      result[p,3] <- p
    }
    result_disp <- rbind(result_disp, result)
  }
  result_disp$Patch <- as.factor(result_disp$Patch)
  return(result_disp)
}

vary_conR0P_ISO <- function(meta_comm_df,species_chara,beta,num_spp,num_patches){
  phi <- rep(1, 6)
  result_con <- data.frame(matrix(data = NA, nrow = 0, ncol = 2))
  max_con <- seq(0,0.1,by = 0.001)
  for (z in 1:length(max_con)) {
    go <- max_con[[z]]
    stay <- 1-go
    c <- matrix(data = NA,
                nrow = num_patches, 
                ncol = num_patches)
    for(i in 1:ncol(c)){
      for(m in 1:nrow(c)){
        c[i,m] <- ifelse(i == m, stay, go)
      }
    }
    result <- data.frame(matrix(data = NA, nrow = 1, ncol = 3))
    colnames(result) <- c("PatchR0","connectivity","Patch")
    Sus <- meta_comm_df*species_chara$Susceptible #value of susceptibles
    I <- meta_comm_df*species_chara$Infectious#value of infecteds
    N <- Sus+I #total pop of a patch
    Sus <- as.matrix(Sus)
    I <- as.matrix(I)
    N <- as.matrix(N)
    N_meta <- colSums(N)
    b <- species_chara$death + 
      species_chara$recovery + 
      species_chara$alpha  # total loss rate
    alpha_div <- rowSums(N)
      
    # Calculate variables
      
    #calculate landscape R0
      r0_patch <- patch_R0_freq_ISO(beta = beta,
                                  sus = Sus,
                                  N = N,
                                  b = b,
                                  S = num_spp,
                                  P = num_patches)
      for (p in 1:num_patches) {
        result[p,1] <- max(abs(eigen(r0_patch[[1]][[p]])$values))
        result[p,2] <- go
        result[p,3] <- p
      }
        result_con <- rbind(result_con, result)
  }
  result_con$Patch <- as.factor(result_con$Patch)
  return(result_con)
}

vary_conR0P_connected <- function(meta_comm_df,species_chara,beta,num_spp,num_patches){
  set.seed(1234)
  phi <- rep(1, 6)
  result_con <- data.frame(matrix(data = NA, nrow = 0, ncol = 2))
  max_con <- seq(0,0.1,by = 0.001)
  result_con <- data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  for(z in 1:length(max_con)){
    go <- max_con[[z]]
    stay <- 1-go
    c <- matrix(data = NA,
                nrow = num_patches, 
                ncol = num_patches)
    for(i in 1:ncol(c)){
      for(m in 1:nrow(c)){
        c[i,m] <- ifelse(i == m, stay, go)
      }
    }
    result <- data.frame(matrix(data = NA, nrow = 1, ncol = 3))
    colnames(result) <- c("PatchR0","connectivity","Patch")
    
    Sus <- meta_comm_df*species_chara$Susceptible #value of susceptibles
    I <- meta_comm_df*species_chara$Infectious#value of infecteds
    N <- Sus+I #total pop of a patch
    Sus <- as.matrix(Sus)
    I <- as.matrix(I)
    N <- as.matrix(N)
    N_meta <- colSums(N)
    b <- species_chara$death + 
      species_chara$recovery + 
      species_chara$alpha  # total loss rate
    alpha_div <- rowSums(N)
    
    # Calculate variables
    
    #calculate landscape R0
    r0_patch <- patch_R0_freq_CON(beta = beta,
                                  sus = Sus,
                                  N = N,
                                  b = b,
                                  S = num_spp,
                                  P = num_patches,
                                  phi = phi,
                                  c = c)
    for (p in 1:num_patches) {
      result[p,1] <- max(abs(eigen(r0_patch[[1]][[p]])$values))
      result[p,2] <- go
      result[p,3] <- p
    }
    result_con <- rbind(result_con, result)
  }
  result_con$Patch <- as.factor(result_con$Patch)
  return(result_con)
}

vary_con <- function(meta_comm_df,species_chara,beta,num_spp,num_patches){
  phi <- rep(1, 6)
  result_con <- data.frame(matrix(data = NA, nrow = 0, ncol = 2))
  max_con <- seq(0,0.1,by = 0.001)
  for (z in 1:length(max_con)) {
    go <- max_con[[z]]
    stay <- 1-go
    c <- matrix(data = NA,
                nrow = num_patches, 
                ncol = num_patches)
    for(i in 1:ncol(c)){
      for(m in 1:nrow(c)){
        c[i,m] <- ifelse(i == m, stay, go)
      }
    } 
    result <- data.frame(matrix(data = NA, nrow = 1, ncol = 2))
    colnames(result) <- c("LandscapeR0","con")
    
    Sus <- meta_comm_df*species_chara$Susceptible #value of susceptibles
    I <- meta_comm_df*species_chara$Infectious#value of infecteds
    N <- Sus+I #total pop of a patch
    Sus <- as.matrix(Sus)
    I <- as.matrix(I)
    N <- as.matrix(N)
    N_meta <- colSums(N)
    b <- v + d  # total loss rate
    alpha_div <- rowSums(N)
    
    # Calculate variables
    
    #calculate landscape R0
    r0_landscape <- landscape_R0_freq(beta = beta,
                                      sus = Sus,
                                      N = N,
                                      Cmat = c,
                                      b = b,
                                      phi = phi,
                                      S = num_spp,
                                      P = num_patches)
    
    result[1,1] <- max(abs(eigen(r0_landscape[[1]])$values))
    result[1,2] <- go
    result_con <- rbind(result_con, result)
  }
  return(result_con)
}
metacomm_func <- function(alpha, patches = 2, species = 6){
  #setting up meta-community parameters
  set.seed(1234)
  num_patches <- patches #number of patches in metacommunity
  num_spp <- species #number of POSSIBLE spp in metacommunity
  Size <- 1000 # size in m^2
  K_rel <- c(10, 7.5, 7, 8, 4, 3)
  K <- Size*K_rel
  a <- 2000
  b <- 400
  max_abund <- sum(K)/4 # a somewhat arbitrarily decided upon max abundance
  meta_comm <- vector("list", length = num_patches)
  for(c in 1:num_patches){
    alpha_div <- alpha[[c]]
    #need to create a relationship between species richness and abundance
    a <- 2000
    b <- 400
    error <- rnorm(length(alpha_div), mean = 0, sd = 1)
    abundance <- a*log(b*alpha_div)+error
    R <- alpha_div
    KCOM <- max_abund/(1+300*exp(-0.30*(R+50)))
    KS <- KCOM/abundance
    meta_comm[[c]] <- if(abundance >= KCOM){
      c(K[1:alpha_div]*KS,rep(0, num_spp - alpha_div))} else{
        c(K[1:alpha_div],rep(0, num_spp - alpha_div))
      }
  }
  meta_comm_df <- as.data.frame(matrix(unlist(meta_comm),nrow = num_patches, ncol = num_spp, byrow = T))
  return(meta_comm_df)
}




#setup spp parameters
Species <- c("PREG","TGRAN","TTOR","ABOR","RCAT","RDRAY")
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
# Size <- 1000 # size in m^2
# K_rel <- c(10, 7.5, 7, 8, 4, 3)
# K <- Size*K_rel

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
                            #mass = mass,
                            alpha = alpha,
                            R0 = R0,
                            #K = K,
                            Susceptible = Susceptible,
                            Infectious = Infectious)

closeness <- 0.05 #based on mihaljevic 2014
beta <- matrix(nrow = length(Species), ncol = length(Species))

for (s in 1:length(Species)) {
  for (i in 1:length(Species)) {
    beta[s,i] <- if(s == i){
      species_chara[s,ncol(species_chara)]} else{
        closeness * ((species_chara[s,ncol(species_chara)]*species_chara[i,ncol(species_chara)])/2)}
  }
  
}
# 2 patches. patch 1 has 1 spp patch 2 has 6 spp
#setting up meta-community parameters
num_patches <- 2
num_spp <- 6
alpha <- c(1,6)
sim1 <- metacomm_func(alpha = alpha)

# test <- patch_R0_freq(beta = beta,
#                       sus = sim1*species_chara$Susceptible,
#                       N = sim1,
#                       b = rowSums(species_chara[,c(3,4,6)]),
#                       S = num_spp,
#                       P = num_patches)
# max(abs(eigen(test[[1]][[2]])$values))

# Doing actual simulation
#Extreme connectivity
#Patches treated as isolated
sim1_disp_iso <- vary_dispR0P_iso(meta_comm_df = sim1,
                                  species_chara = species_chara,
                                  beta = beta,
                                  num_spp = num_spp,
                                  num_patches = num_patches)
sim1_disp_iso_plot <- ggplot(data = sim1_disp_iso, mapping = aes(x = Dispersal, y = PatchR0))+
  geom_point()+
  labs(x = "Dispersal", y = expression("Patch"~R[0]))+
  theme_classic()

sim1_disp_iso_plot # behaves how expected!!!

# Patches treated as connected
sim1_disp_CON <- vary_dispR0P_connected(meta_comm_df = sim1,
                                        species_chara = species_chara,
                                        beta = beta,
                                        num_spp = num_spp,
                                        num_patches = num_patches)
sim1_disp_con_plot <- ggplot(data = sim1_disp_CON, mapping = aes(x = Dispersal, y = PatchR0,colour = Patch))+
  geom_point()+
  geom_line(data = sim1_disp_iso)+
  labs(x = "Dispersal", y = expression("Patch"~R[0]))+
  theme_classic()

sim1_disp_con_plot # behaves how expected!!!

sim1_disp_plot <- ggplot(data = sim1_disp_CON, mapping = aes(x = Dispersal, y = PatchR0, colour = Patch))+
  scale_color_grey()+
  geom_point(size = 2)+
  geom_line(data = sim1_disp_iso, linewidth = 2)+
  ylim(0,35)+
  labs(x = "Dispersal", y = expression("Patch"~R[0]))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(size = 0, hjust = 0.5),
        axis.text = element_text(size = 0),
        axis.line = element_line(linewidth = 1.5))

sim1_disp_plot
# # R0L
# disp_results_R0L <- vary_disp_R0L(meta_comm_df = sim1,
#                           species_chara = species_chara,
#                           beta = beta,
#                           num_spp = num_spp,
#                           num_patches = num_patches)
# disp_R0L_plot <- ggplot(data = disp_results_R0L, mapping = aes(x = max_disp, y = LandscapeR0))+
#   geom_point(size = 2)+
#   geom_smooth(method = "glm", formula = y~x,
#               method.args = list(family = gaussian(link = 'log')),
#               colour = "black")+
#   labs(x = "Dispersal", y = expression(R[0]))+
#   theme_classic()+
#   theme(axis.title = element_text(size = 25),
#         plot.title = element_text(size = 40, hjust = 0.5),
#         axis.text = element_text(size = 0),
#         axis.line = element_line(linewidth = 1.5))
# disp_R0L_plot
#Extreme dispersal
# Patches are treated as isolated
sim1_con_iso <- vary_conR0P_ISO(meta_comm_df = sim1,
                                species_chara = species_chara,
                                beta = beta,
                                num_spp = num_spp,
                                num_patches = num_patches)
sim1_con_CON <- vary_conR0P_connected(meta_comm_df = sim1,
                                      species_chara = species_chara,
                                      beta = beta,
                                      num_spp = num_spp,
                                      num_patches = num_patches)

sim1_con_plot <- ggplot(data = sim1_con_CON, aes(x = connectivity, y = PatchR0, colour = Patch))+
  scale_color_grey()+
  geom_point(size = 2)+
  geom_line(data = sim1_con_iso,linewidth = 2)+
  labs(x = "Connectivity", y = expression("Patch"~R[0]))+
  ylim(0,35)+
  theme_classic()+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(size = 0, hjust = 0.5),
        axis.text = element_text(size = 0),
        axis.line = element_line(linewidth = 1.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
sim1_con_plot

# con_results <- vary_con(meta_comm_df = sim1, 
#          species_chara = species_chara, 
#          beta = beta,
#          num_spp = num_spp, 
#          num_patches = num_patches)
# con_plot <- ggplot(data = con_results, mapping = aes(x = con, y = LandscapeR0))+
#   geom_point(size = 2)+
#   geom_smooth(method = "glm", formula = y~x,
#               method.args = list(family = gaussian(link = 'log')),
#               colour = "black")+
#   # ggtitle("Connectivity (sim1)")+
#   labs(x = "Connectivity", y = expression(R[0]))+
#   theme_classic()+
#   theme(axis.title = element_text(size = 25),
#         plot.title = element_text(size = 40, hjust = 0.5),
#         axis.text = element_text(size = 0),
#         axis.line = element_line(linewidth = 1.5))
# 
# con_plot 

# 2 patches. both have 1 spp
#setting up meta-community parameters
num_patches <- 2
num_spp <- 6
alpha <- c(1,1)
sim2 <- metacomm_func(alpha = alpha)

# Doing actual simulation
#Extreme connectivity

sim2_disp_iso <- vary_dispR0P_iso(meta_comm_df = sim2,
                                  species_chara = species_chara,
                                  beta = beta,
                                  num_spp = num_spp,
                                  num_patches = num_patches)
sim2_disp_con <- vary_dispR0P_connected(meta_comm_df = sim2,
                                  species_chara = species_chara,
                                  beta = beta,
                                  num_spp = num_spp,
                                  num_patches = num_patches)

sim2_disp_plot <- ggplot(data = sim2_disp_con, mapping = aes(x = Dispersal, y = PatchR0, colour = Patch))+
  scale_color_grey()+
  geom_point(size = 3)+
  geom_line(data = sim2_disp_iso, linewidth = 2)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "Dispersal", y = expression("Patch"~R[0]))+
  ylim(0,35)+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(size = 0, hjust = 0.5),
        axis.text = element_text(size = 0),
        axis.line = element_line(linewidth = 1.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
sim2_disp_plot

# disp_results_2 <- vary_disp(meta_comm_df = sim2,
#                           species_chara = species_chara,
#                           beta = beta,
#                           num_spp = num_spp,
#                           num_patches = num_patches)
# disp_plot_2 <- ggplot(data = disp_results_2, mapping = aes(x = max_disp, y = LandscapeR0))+
#   geom_point(size = 2)+
#   geom_smooth(method = "glm", formula = y~x,
#               method.args = list(family = gaussian(link = 'log')),
#               colour = "black")+
#   labs(x = "Dispersal", y = expression(R[0]))+
#   theme_classic()+
#   theme(axis.title = element_text(size = 25),
#         plot.title = element_text(size = 40, hjust = 0.5),
#         axis.text = element_text(size = 0),
#         axis.line = element_line(linewidth = 1.5))
# 
# disp_plot_2
#Extreme dispersal
sim2_con_iso <- vary_conR0P_ISO(meta_comm_df = sim2,
                                species_chara = species_chara,
                                beta = beta,
                                num_spp = num_spp,
                                num_patches = num_patches)

sim2_con_CON <- vary_conR0P_connected(meta_comm_df = sim2,
                                species_chara = species_chara,
                                beta = beta,
                                num_spp = num_spp,
                                num_patches = num_patches)
sim2_con_plot <- ggplot(data = sim2_con_CON, aes(x = connectivity, y = PatchR0, colour = Patch))+
  scale_color_grey()+
  geom_point(size = 2)+
  geom_line(data = sim2_con_iso, linewidth = 2)+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,35)+
  labs(x = "Connectivity", y = expression("Patch"~R[0]))+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(size = 0, hjust = 0.5),
        axis.text = element_text(size = 0),
        axis.line = element_line(linewidth = 1.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
sim2_con_plot


# con_results_2 <- vary_con(meta_comm_df = sim2, 
#                         species_chara = species_chara, 
#                         beta = beta,num_spp = num_spp, 
#                         num_patches = num_patches)
# con_plot_2 <- ggplot(data = con_results_2, mapping = aes(x = con, y = LandscapeR0))+
#   geom_point(size = 2)+
#   geom_smooth(method = "glm", formula = y~x,
#               method.args = list(family = gaussian(link = 'log')),
#               colour = "black")+
#   # ggtitle("Connectivity (sim1)")+
#   labs(x = "Connectivity", y = expression(R[0]))+
#   theme_classic()+
#   theme(axis.title = element_text(size = 25),
#         plot.title = element_text(size = 40, hjust = 0.5),
#         axis.text = element_text(size = 0),
#         axis.line = element_line(linewidth = 1.5))

# con_plot_2 


# 2 patches, both have 6 species

#setting up meta-community parameters
num_patches <- 2
num_spp <- 6
alpha <- c(6,6)
sim3 <- metacomm_func(alpha = alpha)

# Doing actual simulation
#Extreme connectivity

sim3_disp_iso <- vary_dispR0P_iso(meta_comm_df = sim3,
                                  species_chara = species_chara,
                                  beta = beta,
                                  num_spp = num_spp,
                                  num_patches = num_patches)
sim3_disp_con <- vary_dispR0P_connected(meta_comm_df = sim3,
                                        species_chara = species_chara,
                                        beta = beta,
                                        num_spp = num_spp,
                                        num_patches = num_patches)

sim3_disp_plot <- ggplot(data = sim3_disp_con, mapping = aes(x = Dispersal, y = PatchR0, colour = Patch))+
  scale_color_grey()+
  geom_point(size = 2)+
  geom_line(data = sim3_disp_iso, linewidth = 2)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "Dispersal", y = expression("Patch"~R[0]))+
  ylim(0,35)+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(size = 0, hjust = 0.5),
        axis.text = element_text(size = 0),
        axis.line = element_line(linewidth = 1.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
sim3_disp_plot


# disp_results_3 <- vary_disp(meta_comm_df = sim3,
#                             species_chara = species_chara,
#                             beta = beta,
#                             num_spp = num_spp,
#                             num_patches = num_patches)
# disp_plot_3 <- ggplot(data = disp_results_3, mapping = aes(x = max_disp, y = LandscapeR0))+
#   geom_point(size = 2)+
#   geom_smooth(method = "glm", formula = y~x,
#               method.args = list(family = gaussian(link = 'log')),
#               colour = "black")+
#   labs(x = "Dispersal", y = expression(R[0]))+
#   theme_classic()+
#   theme(axis.title = element_text(size = 25),
#         plot.title = element_text(size = 40, hjust = 0.5),
#         axis.text = element_text(size = 0),
#         axis.line = element_line(linewidth = 1.5))
# disp_plot_3
#Extreme dispersal

sim3_con_iso <- vary_conR0P_ISO(meta_comm_df = sim3,
                                species_chara = species_chara,
                                beta = beta,
                                num_spp = num_spp,
                                num_patches = num_patches)

sim3_con_CON <- vary_conR0P_connected(meta_comm_df = sim3,
                                      species_chara = species_chara,
                                      beta = beta,
                                      num_spp = num_spp,
                                      num_patches = num_patches)
sim3_con_plot <- ggplot(data = sim3_con_CON, aes(x = connectivity, y = PatchR0, colour = Patch))+
  scale_color_grey()+
  geom_point(size = 2)+
  geom_line(data = sim3_con_iso, linewidth = 2)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "Connectivity", y = expression("Patch"~R[0]))+
  ylim(0,35)+
  theme(axis.title = element_text(size = 25),
        plot.title = element_text(size = 0, hjust = 0.5),
        axis.text = element_text(size = 0),
        axis.line = element_line(linewidth = 1.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
sim3_con_plot



# con_results_3 <- vary_con(meta_comm_df = sim3, 
#                           species_chara = species_chara, 
#                           beta = beta,num_spp = num_spp, 
#                           num_patches = num_patches)
# con_plot_3 <- ggplot(data = con_results_3, mapping = aes(x = con, y = LandscapeR0))+
#   geom_point(size = 2)+
#   geom_smooth(method = "glm", formula = y~x,
#               method.args = list(family = gaussian(link = 'log')),
#               colour = "black")+
#   # ggtitle("Connectivity (sim1)")+
#   labs(x = "Connectivity", y = expression(R[0]))+
#   theme_classic()+
#   theme(axis.title = element_text(size = 25),
#         plot.title = element_text(size = 40, hjust = 0.5),
#         axis.text = element_text(size = 0),
#         axis.line = element_line(linewidth = 1.5))
# 
# con_plot_3 


# view disp plots together
Sim1_plot <- (sim1_disp_plot + sim1_con_plot)
Sim1_plot # Is actually scenario 2 in the figure in manuscript
Sim2_plot <- (sim2_disp_plot + sim2_con_plot)
Sim2_plot # Is actually scenario 1 in the figure in manusript
Sim3_plot <- (sim3_disp_plot + sim3_con_plot)
Sim3_plot
