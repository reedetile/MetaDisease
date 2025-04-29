#Description-----------------------------------------
#Script to build sets of 5 patch metacommunities
#  29 Apr 2025
#RCS

# Load packages---------------------------------
library(vegan)
# Parameters------------------------------------
repo <- "D:/gitrepos/MetaDisease"
graphs <- paste(repo,"/Graphs",sep="")
#setting up meta-community parameters
set.seed(1234)
num_patches <- 5 #number of patches in metacommunity
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

# okay so I've made my 5 patch metacomms. Now I want to show that this still follows a saturated model
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

# save meta comm list
setwd(repo)
meta_comm_list <- meta_comm_list %>% lapply(unlist) %>% lapply(matrix,nrow = num_patches, ncol = num_spp, byrow = T) %>% lapply(data.frame)
saveRDS(meta_comm_list, "metacomm_5Patch.RDS")








