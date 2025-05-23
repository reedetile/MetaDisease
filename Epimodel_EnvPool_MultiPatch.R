#Description-----------------------------------------
#Epidiomological model for an environmental pool disease
#  09 Jan 2025
#RCS

#Initialize -----------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
source('Occu_abun_practice.R')
source('Epimodel_funcs.R')
# Load functions--------------------------------------


# Global Variables-------------------------------------


# Program Body------------------------------------------

# Parameters-------------------------------------
#starting values
#S <- ceiling(meta_comm1[,1:6]*c(0.8,0.85,0.9,0.93,0.95,.99)) #starting value of susceptibles
S <- meta_comm1[,1:6] #starting value of susceptibles
I <- meta_comm1[,1:6] - S #starting value of infecteds
N <- S+I
Z <- 1000 #starting amount of zoospores

#species characteristics
Species <- c("Spp1","Spp2","Spp3","Spp4","spp5","spp6")
b <- c(0.6,0.5,0.4,0.3,0.2,0.1) #host birth rate
# Justification = ...? NEED CITATION(S)
d <- c(0.06,0.05,0.04,0.03,0.02,0.01) #host death rate
# Justification = ...? NEED CITATION(S)
beta <- c(0.00013,0.00012,0.00011,0.00010,0.00009,0.00008) #transmission per species
# Justification = ...? NEED CITATION(S)
v <- c(0.4,0.5,0.6,0.7,0.8,0.9) #recovery rate
# Justification = ...? NEED CITATION(S)
phi <- c(0.09,0.08,0.07,0.06,0.05,0.04) #dispersal rate
# Justification = ...? NEED CITATION(S)
lambda <- c(3.0,2.5,2.0,1.5,1.0,0.5) #log(bd) load. losely based off figure 2 from wilber 2020
# Justification = Wilber 2020. Could use more citations
species_chara <- data.frame(Species = Species,
                            birth = b, 
                            death = d, 
                            trans = beta, 
                            recovery = v, 
                            dispersal = phi, 
                            shedding = lambda)

# meta-community characteristics
gamma <- 0.001 #zoospore decay rate
c <- matrix(data = rnorm(n = num_patches^2, mean = 0.5, sd = 0.1),
            nrow = num_patches, 
            ncol = num_patches)
#Connectivity of patches
A <- rnorm(n = num_patches, mean = 0.5, sd = 0.1) #area ratios
time <- 90 #how many "days" do I want in the season

# Program Body------------------------------------------
pop_list_EP <- vector("list", length = time)
for (t in 1:time) {
  delta_s <- b*N - d*S - beta*S*Z + v*I + phi*sum(-c*S + c*S) #for now have excluded area of patches
  #may want to add that back in though
  delta_I <- beta*S*Z-(v+d)*I+phi*sum(-c*I + c*I)#for now have excluded area of patches
  #may want to add that back in though
  delta_Z <- sum(lambda*I - gamma*Z)
  S <- S+delta_s
  I <- I+delta_I
  z <- Z+delta_Z
  N <- S + I
  t = t
  pop <- list(Susceptible = S, Infectious = I, Zoospores = Z, Total = N, Time = t)
  pop_list_EP[[t]] <- pop
}

pop_data_EP <- do.call(rbind, pop_list_EP)
View(pop_data_EP)  