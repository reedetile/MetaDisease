#Description-----------------------------------------
#Script to build sets of 5 patch metacommunities
#  29 Apr 2025
#RCS

# Load packages---------------------------------
library(vegan)
library(ggplot2)
library(tidyr)
library(patchwork)
# Parameters------------------------------------
repo <- "C:/Users/rcscott/MetaDisease" # laptop repo
# repo <- "D:/gitrepos/MetaDisease" #PC repo
graphs <- paste(repo,"/Graphs",sep="")
#setting up meta-community parameters
set.seed(1234)
num_patches <- 5 #number of patches in metacommunity
num_spp <- 6 #number of POSSIBLE spp in metacommunity

N <- 100 #number of metacommunity simulations to run

Size <- 1000 # size in m^2
K_rel <- c(10, 7.5, 7, 8, 4, 3)
K <- Size*K_rel

max_abund <- sum(K)/4 # a somewhat arbitrarily decided upon max abundance

#Building the metacommunities
meta_comm_list <- vector("list",N)
beta_list <- vector("list", N)
nestedness_list <- vector("list", N)

for(n in 1:N){
  meta_comm <- vector("list", length = num_patches)
  #data.frame(matrix(NA, nrow = num_patches, ncol = num_spp)) #create metacom
  for(c in 1:num_patches){
    alpha <- as.numeric(sample(1:6, size = 1, replace = T))
    #need to create a relationship between species richness and abundance
    a <- 2000
    b <- 400
    error <- rnorm(length(alpha), mean = 0, sd = 1)
    abundance <- a*log(b*alpha)+error
    R <- alpha
    KCOM <- max_abund/(1+300*exp(-0.30*(R+50)))
    KS <- KCOM/abundance
    meta_comm[[c]] <- if(abundance >= KCOM){
      c(K[1:alpha]*KS,rep(0, num_spp - alpha))} else{
        c(K[1:alpha],rep(0, num_spp - alpha))
      }
  }
  meta_comm_list[[n]] <- meta_comm
}

# okay so I've made my 5 patch metacomms. Now I want to show that this still follows a saturated model
matrix(data = unlist(meta_comm_list[[1]]), nrow = num_patches, ncol = num_spp, byrow = T)
alpha_df <- data.frame(matrix(unlist(meta_comm_list), nrow = num_patches*N, ncol = num_spp, byrow = T))
alpha_df$abund <- rowSums(alpha_df[1:6])
max(alpha_df$abund)
alpha_df$richness <- rowSums(alpha_df[1:6] > 0)
plot(alpha_df$richness, alpha_df$abund) # alpha diversity is a saturated curve!
alpha_sat_plot <- ggplot(data = alpha_df, aes(x = richness, y = abund))+
  geom_jitter(width = 0.5, height = 1000)+
  geom_point()+ 
  geom_line()+
  xlab("Species Richness")+
  ylab("Species Abundance")+
  #  ggtitle("Alpha diversity relationship with abundance")+
  theme_classic()
alpha_sat_plot
setwd(graphs)
saveRDS(alpha_sat_plot, file = "alpha_sat_plot5.RDS")
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

# save meta comm list
setwd(repo)
meta_comm_list <- meta_comm_list %>% lapply(unlist) %>% lapply(matrix,nrow = num_patches, ncol = num_spp, byrow = T) %>% lapply(data.frame)
saveRDS(meta_comm_list, "metacomm_5Patch.RDS")








