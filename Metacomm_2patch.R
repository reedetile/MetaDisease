#Description-----------------------------------------
#PCreating 100 metacommunities: metacommunities consist of 2 patches with 6 species
#  22 Nov 2024
#RCS

# Load packages---------------------------------
library(vegan)
library(ggplot2)
library(tidyr)
library(patchwork)

#setup -------------------------------------------
repo <- "D:/gitrepos/MetaDisease"
graphs <- paste(repo,"/Graphs",sep="")
# Parameters-------------------------------------
#setting up meta-community parameters
set.seed(1234)
num_patches <- 2 #number of patches in metacommunity
num_spp <- 6 #number of POSSIBLE spp in metacommunity
N <- 100 #number of metacommunity simulations to run

meta_comm1 <- data.frame(matrix(NA, nrow = num_patches, ncol = num_spp))
S <- c(0.83, #PREG
       0.62, #TGRAN
       0.14, #TTOR
       NA, #ABOR
       0.12, #RCAT
       0.64) #RDRAY
S[4] <- runif(n = 1, min = S[5], max = S[3])
#an array of probability values for the occurence of each spp
#what if I over thought this, and I can just assign an occupancy probability?
K <- c(20,16,15,7,4,2) #this is just an example, but k is the abundance at each rank (i think?)
max_abund <- 65 # a somewhat arbitrarily decided upon max abundance
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

meta_comm_list <- vector("list",N)
beta_list <- vector("list", N)
nestedness_list <- vector("list", N)
for(n in 1:N){
  meta_comm <- data.frame(matrix(NA, nrow = num_patches, ncol = num_spp))
  for(i in 1:nrow(meta_comm)){
    for(j in 1:ncol(meta_comm)){
      #determine occurrence and abundance of spp 1
      meta_comm[i,1] <- rbinom(n = 1, size = 1, prob = S[1])
      meta_comm[i,1] <- ifelse(meta_comm[i,1] == 1, K[1],0)
      
      #determine occurrence and abundance of spp2
      meta_comm[i,2] <- rbinom(n = 1, size = 1, prob = S[2])
      meta_comm[i,2] <- ifelse(meta_comm[i,2] == 1, K[2],0)
      
      #determine occurrence and abundance of spp 3
      meta_comm[i,3] <- rbinom(n = 1, size = 1, prob = S[3])
      meta_comm[i,3] <- ifelse(meta_comm[i,3] == 1, K[3],0)
      
      #determine occurrence and abundance of spp 4
      meta_comm[i,4] <- rbinom(n = 1, size = 1, prob = S[4])
      meta_comm[i,4] <- ifelse(meta_comm[i,4] == 1, K[4],0)
      
      #determine occurrence and abundance of spp 5
      meta_comm[i,5] <- rbinom(n = 1, size = 1, prob = S[5])
      meta_comm[i,5] <- ifelse(meta_comm[i,5] == 1, K[5],0)
      
      #determine occurrence and abundance of spp 6
      meta_comm[i,6] <- rbinom(n = 1, size = 1, prob = S[6])
      meta_comm[i,6] <- ifelse(meta_comm[i,6] == 1, K[6],0)
    }
  }
  meta_comm$Patch <- c('Patch1','Patch2')
  colnames(meta_comm) <- c('Spp1','Spp2','Spp3','spp4','spp5','Spp6','PatchID')
  beta_list[[n]] <- mean(betadiver(meta_comm[,1:6], method = 'w'))
  nestedness_list[[n]] <- nestedtemp(comm = meta_comm[,1:6])[7]
  meta_comm_list[[n]] <- meta_comm
}

# nestedness <- unlist(x = nestedness_list)
# beta <- unlist(x = beta_list)
# meta_comm_vars <- data.frame(Diversity = beta, Nest = nestedness) 
# 
# 
# meta_comm_list[[2]]
# beta_diversity <- betadiver(meta_comm1[,1:6], method = 'w')
#plot(beta_diversity)

# 03/31/2025: changing from an additive model to a saturatured model
meta_comm_list <- vector("list",N)
beta_list <- vector("list", N)
nestedness_list <- vector("list", N)

for(n in 1:N){
  meta_comm <- vector("list", length = num_patches)
                      #data.frame(matrix(NA, nrow = num_patches, ncol = num_spp)) #create metacom
  for(c in 1:num_patches){
    alpha <- as.numeric(sample(1:6, size = 1, replace = T))
    #need to create a relationship between species richness and abundance
    a <- 10
    b <- 4
    error <- rnorm(length(alpha), mean = 0, sd = 1)
    abundance <- a*log(b*alpha)+error
    R <- alpha
    KCOM <- max_abund/(1+3*exp(-0.05*(R)))
    KS <- KCOM/abundance
    meta_comm[[c]] <- if(abundance >= KCOM){
      c(K[1:alpha]*KS,rep(0, num_spp - alpha))} else{
        c(K[1:alpha],rep(0, num_spp - alpha))
      }
  }
  meta_comm_list[[n]] <- meta_comm
}

# okay so I've made my 2 patch metacomms. Now I want to show that this still follows a saturated model
matrix(data = unlist(meta_comm_list[[1]]), nrow = num_patches, ncol = num_spp, byrow = T)
alpha_df <- data.frame(matrix(unlist(meta_comm_list), nrow = num_patches*N, ncol = num_spp, byrow = T))
alpha_df$abund <- rowSums(alpha_df[1:6])
max(alpha_df$abund)
alpha_df$richness <- rowSums(alpha_df[1:6] > 0)
plot(alpha_df$richness, alpha_df$abund) # alpha diversity is a saturated curve!
alpha_sat_plot <- ggplot(data = alpha_df, aes(x = richness, y = abund))+
  geom_point()+ 
  xlab("Species Richness")+
  ylab("Species Abundance")+
#  ggtitle("Alpha diversity relationship with abundance")+
  theme_classic()
alpha_sat_plot


gamma_df <- data.frame(richness = rep(NA, length = N), abund = rep(NA, length = N))
for(n in 1:N){
  meta_comm <- data.frame(matrix(unlist(meta_comm_list[[n]]), nrow = num_patches, ncol = num_spp, byrow=T))
  gamma_df[n,2] <- sum(colSums(meta_comm))
  gamma_df[n,1]<- sum(colSums(meta_comm) > 0)
}
plot(gamma_df$richness, gamma_df$abund)
gamma_sat_plot <- ggplot(data = gamma_df, aes(x = richness, y = abund))+
  geom_point()+ 
  xlab("Species Richness")+
  ylab("Species Abundance")+
#  ggtitle("Gamma diversity relationship with abundance")+
  theme_classic()
gamma_sat_plot
# let's save these
sat_plots <- alpha_sat_plot / gamma_sat_plot + plot_annotation(tag_levels = "A")
setwd(graphs)
ggsave(filename = 'sat_plots.png', plot = sat_plots)

# save meta comm list
setwd(repo)
meta_comm_list <- meta_comm_list %>% lapply(unlist) %>% lapply(matrix,nrow = num_patches, ncol = num_spp, byrow = T) %>% lapply(data.frame)


saveRDS(meta_comm_list, "metacomm_2Patch.RDS")


# this is also a saturated curve!

# Next step: I need to show how communities change, or stay the same over time
Time <- 90 #assume a 90 day breeding season
phi <- runif(6)# Need to establish dispersal metric. May need to determine more realistic values (see notes
# from meeting with mark + brittany on 03/28/25)
meta_comm_change <- data.frame(meta_com = 1:length(meta_comm_list), time = NA)
for(c in 1:nrow(meta_comm_change)){
  meta_comm <- meta_comm_list[[c]]
  meta_comm_df <- data.frame(matrix(unlist(meta_comm), nrow = num_patches, ncol = num_spp, byrow = T))
  for(t in 1:Time) {
    deltaP <- data.frame(matrix(data = NA, nrow = num_patches, ncol = num_spp))
    for (i in 1:num_patches) {
      for(j in 1:num_spp) {
        other_patch <- as.numeric(meta_comm_df[-i,])
        deltaP[i,j] <- phi[[j]]*other_patch[[j]] - phi[[j]]*meta_comm_df[i,j]
      }
    }
    meta_comm_change[c,2] <- if(sum(abs(deltaP)) > 0.01){ 
      #if deltaP is > 0 then the system is not yet at equil.
      t} else{
        next}
    new_abund <- meta_comm_df + deltaP
    meta_comm_df <- new_abund
  }
}

  
# let's plot meta comm 1 as an example to show B + M
# setup initial conditions
meta_comm <- data.frame(matrix(unlist(meta_comm_list[[1]]), nrow = num_patches, ncol = num_spp, byrow = T))
timeXchange <- vector("list", length = Time)
timeXchange[[1]] <- meta_comm
deltaP_list  <- vector("list", length = Time-1)

for(t in 2:Time){
  meta_comm <- timeXchange[[t-1]]
  deltaP <- data.frame(matrix(data = NA, nrow = num_patches, ncol = num_spp))
  for (i in 1:num_patches) {
    for(j in 1:num_spp) {
      other_patch <- as.numeric(meta_comm[-i,])
      deltaP[i,j] <- phi[[j]]*other_patch[[j]] - phi[[j]]*meta_comm[i,j]
    }
  }
  meta_comm <- meta_comm + deltaP
  timeXchange[[t]] <- meta_comm
  deltaP_list[[t-1]] <- sum(abs(deltaP))
}

# add in time covariate
for (t in 1:Time) {
  timeXchange[[t]]$Time <- rep(t,2)
}


timeXchange_df <- do.call(rbind, timeXchange)

deltaP_df <- data.frame(deltaP = unlist(deltaP_list), Time = 2:90) 

timeXchange_df$Patch <- rep(1:2, times = 90)
colnames(timeXchange_df) <- c("Spp1","Spp2","Spp3","Spp4","Spp5","Spp6","Time","Patch")
timeXchange_df <- timeXchange_df %>% pivot_longer(cols = Spp1:Spp6, names_to = "Species")

timeXchange_plot <- ggplot(data = timeXchange_df, 
                           mapping = aes(x = Time, y =  value, shape = as.factor(Patch), colour = Species))+
  geom_point()+
  geom_line()+
  theme_classic()
timeXchange_plot

deltaP_plot <- ggplot(data = deltaP_df, mapping = aes(x = Time, y = deltaP))+
  geom_point()+
  geom_line()+
  theme_classic()
deltaP_plot


setwd(graphs)
ggsave(filename = "equil_plot.png", plot = timeXchange_plot)
ggsave(filename = "deltaP_plot.png", plot = deltaP_plot)


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
