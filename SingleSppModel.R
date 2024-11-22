# #Description-----------------------------------------
# #This script provides an example of how to build a single species metapopulation model
# #  03 Sep 2024
# #RCS
# 
# #Initialize -----------------------------------------
# 
# #### The Simple SIR Model ####
# 
# #Parameters
# S <- 99 # Number of susceptible individuals in population
# I <- 1 # Number of infected individuals in population
# R <- 0 # Number of recovered individuals in population
# N <- S+I+R
# B <- 0.1 #birth rate
# d <- 0.5 # death rate
# beta <- 0.4 #transmission rate
# gamma <- 0.9 #recovery rate
# 
# #Diff. equations#
# delta_S <- B*N-beta*((S*I)/N)-d*S
# delta_I <- beta*((S*I)/N)-gamma*I-d*I
# delta_R <- gamma*I-d*R
# 
# 
# #### Coupling A 2 Patch Model ####
# S <- c(99,99)
# I <- c(1,1)
# R <- c(0,0)
# N <- S+I+R
# 
# B <- 0.1 #birth rate
# d <- 0.5 # death rate
# beta <- 0.4 #transmission rate
# r <- 0.9 #recovery rate
# 
# #sigma = rate of coupling between two populations
# p <- c(0.2,0.3) #rate of commuting
# tau <- c(0.9,0.9) #rate of returning
# gamma <- p/tau #ratio of commuting to return rates
# 
# 
# #Equations for sigma are based on Equations 17.3 + 17.4 and Hanski + Gaggioti 2004
# sigma <- matrix(data = NA, nro = 2, ncol = 2)
# sigma[1,1] <- ((1-gamma[1])^2)/((1+gamma[1])*N[1]+gamma[2]*N[2])+
#   (gamma[1]^2)/(gamma[1]*N[1]+(1-gamma[2])*N[2])
# sigma[1,2] <- (((1-gamma[1])*gamma[1])/((1-gamma[1])*N[1]+gamma[2]*N[2]))+
#   ((gamma[1]*(1-gamma[2]))/(gamma[1]*N[1]+(1-gamma[1])*N[2]))
# 
# sigma[2,1] <- sigma[1,2]
# 
# delta_Si <- b*N[1]-beta*S[1]*(sigma[]*I[1]+sigma[]*I[2])-d*S[1]


#### Example take from  https://rpubs.com/bhraynor/MetapopulationModel based on Hanski + Gaggioti 2004 ####

#Clean environment
rm(list = ls())

#Load libraries
library(dplyr)
library(ggplot2)
library(deSolve)

# Functions
  #Initial state condition
    FormatInit <- function(init.data.frame){
      
      #Create flexible naming str for vectors
      l <- nrow(init.data.frame)
      state.names <- as.vector(colnames(init.data.frame))
      varnames <- NULL
      for(i in 1:length(state.names)){
        varnames <- cbind(varnames,paste(state.names[i],seq(1:l),sep=""))
      }
      as.vector(varnames)
      
      #vectorize matrix format
      init.vector <- as.vector(as.matrix(init.data.frame))
      
      names(init.vector) <- varnames
      return(init.vector)
    }
  
  #Parameters
    FormatParams <- function(param.df){
      
      #Create flexible naming str for vectors
      l <- nrow(param.df)
      param.names <- as.vector(colnames(param.df))
      varnames <- NULL
      for( i in 1:length(param.names)){
        varnames <- cbind(varnames,paste(param.names[i],seq(1:l),sep=""))
      }
      as.vector(varnames)
      
      #vectorize matrix input
      param.vector <- as.vector(as.matrix(param.df))
      
      names(param.vector) <- varnames
      return(param.vector)
    }

  # contact matrix
    formatContact <- function(df.contact,beta.contact_constant){
      l <- length(df.contact)
      matrix.contact <- as.matrix(df.contact, ncol = 4)
      return(beta.contact_constant*matrix.contact)
    }
    
# Specify the model
    MODEL <- function(time,state,parameters,beta.contact){
      with(as.list(c(state,parameters)),{
        
        #number of cells
        l1 <- length(state)
        l2 <- length(unique(substr(names(state),1,1)))
        l <- l1/l2
        
        #define initial conditions
        S <- matrix(state[1:l],ncol=1)
        I <- matrix(state[(l+1):(2*l)],ncol=1)
        R = matrix(state[(2*l+1):(3*l)], ncol=1)
        
        #Define params
        beta <- matrix(parameters[paste0("beta", 1:l)], ncol=1) #transmission coefficient
        gamma <- matrix(parameters[paste0("gamma", 1:l)], ncol=1) #recovery rate
        alpha <- matrix(parameters[paste0("alpha", 1:l)], ncol=1)  #pathogen induced mortality
        mu <- matrix(parameters[paste0("mu", 1:l)], ncol=1) # mortality rate
        
        #System of eqns
        #bSI = (beta) %*% t(S) %*% I
        #bSI = beta*S*I #density dependent transmission
        bSI = (beta*S*I)/(S+I+R)
        bSI.contact = beta.contact%*%I*S
        birth = mu*(S+I+R) + alpha*I
        
        
        
        dS <- birth -bSI - bSI.contact-mu*S
        dI <-  bSI + bSI.contact - gamma*I -alpha*I - mu*I
        dR <-  gamma*I - mu*R
        
        #Output
        return(list(c(dS, dI, dR)))
      })
    }
    
    
#Initial conditions + parameter specification
init.df <- data.frame(S=c(100,100,100),
                   I= c(10,10,10),
                   R= c(0,0,0))
init <- FormatInit(init.data.frame = init.df)


#parameters
params.df <- data.frame(beta = c(3.98*10^-2,3.98*10^-2,3.98*10^-2), #beta based on Rachowicz + Briggs 2007.
                        #Beta is constant right now but could be adjusted to vary by pop 
                        #(look at range of values in R+B 2007)
                        gamma = c(0.386,0.386,0.386), #gamma based on Briggs et al 2010
                        #Again, gamma is constant but could be varied
                        mu = c(1/1000,1/1000,1/1000),
                        alpha = c(0.2,0.3,0.4)) %>%
  FormatParams()
#contact matrix
beta.contact_constant <- 0.001

Time = 100
dt = 1 #time step size
times <- seq(0,Time,by = dt)


#Example 1: all patches interact with each otherequally
#setup contact matrix for example
df.contact <- data.frame(q1 = c(0,1,1),
                         q2 = c(1,0,1),
                         q3 = c(1,1,0))
beta.contact <- formatContact(df.contact,beta.contact_constant)

#run simulation
output <- ode(y = init, 
              times = times, 
              func = MODEL, 
              parms = params.df,
              beta.contact = beta.contact)

#format output
df <- reshape2::melt(as.data.frame(as.matrix(output)),id = 'time') %>%
  mutate(state = substring(variable,1,1),
         patch = substring(variable,2,2))

#Visualize
Pallete1 <- c('S' = 'gold',
              'I' = 'red3',
              'R' = 'dodgerblue3')
fig1 <- ggplot()+
  theme_classic()+
  geom_line(data=df, aes(x=time,y=value,group=variable,color=state),alpha=0.5,size=2)+
  scale_color_manual(values = Pallete1,
                     name = "Disease State", breaks = c("S","I","R"))+
  ggtitle("Example1")+
  xlab("time")+
  ylab("Number ind. in state")
fig1



### Need to continue working on this model. Need to think about parameter estimations for alpha +
### mu. But that also requires me to understand what those parameters actually are...
