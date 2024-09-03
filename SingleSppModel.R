#Description-----------------------------------------
#This script provides an example of how to build a single species metapopulation model
#  03 Sep 2024
#RCS

#Initialize -----------------------------------------

#### The Simple SIR Model ####

#Parameters
S <- 99 # Number of susceptible individuals in population
I <- 1 # Number of infected individuals in population
R <- 0 # Number of recovered individuals in population
N <- S+I+R
B <- 0.1 #birth rate
d <- 0.5 # death rate
beta <- 0.4 #transmission rate
gamma <- 0.9 #recovery rate

#Diff. equations#
delta_S <- B*N-beta*((S*I)/N)-d*S
delta_I <- beta*((S*I)/N)-gamma*I-d*I
delta_R <- gamma*I-d*R


#### Coupling A 2 Patch Model ####
S <- c(99,99)
I <- c(1,1)
R <- c(0,0)
N <- S+I+R

B <- 0.1 #birth rate
d <- 0.5 # death rate
beta <- 0.4 #transmission rate
r <- 0.9 #recovery rate

#sigma = rate of coupling between two populations
p <- c(0.2,0.3) #rate of commuting
tau <- c(0.9,0.9) #rate of returning
gamma <- p/tau #ratio of commuting to return rates


#Equations for sigma are based on Equations 17.3 + 17.4 and Hanski + Gaggioti 2004
sigma <- matrix(data = NA, nro = 2, ncol = 2)
sigma[1,1] <- ((1-gamma[1])^2)/((1+gamma[1])*N[1]+gamma[2]*N[2])+
  (gamma[1]^2)/(gamma[1]*N[1]+(1-gamma[2])*N[2])
sigma[1,2] <- (((1-gamma[1])*gamma[1])/((1-gamma[1])*N[1]+gamma[2]*N[2]))+
  ((gamma[1]*(1-gamma[2]))/(gamma[1]*N[1]+(1-gamma[1])*N[2]))

sigma[2,1] <- sigma[1,2]

delta_Si <- b*N[1]-beta*S[1]*(sigma[]*I[1]+sigma[]*I[2])-d*S[1]


