#Description-----------------------------------------
#Part 1: A model simulation for the number of zoospores 
# and sporangium in a single population

#Part 2: A model simulation for the number of zoospores
# and sporangium in a meta-population
#  05 Sep 2024
#RCS



# Part 1: A model simulation for the number of zoospore + sporangium in a SINGLE population
rm(list=ls()) #clear environment

### Parameter values ###
N <- 100 #Frog population size
Time <- 365 #365 days in a year
n <- 17.5 #production rate of zoospores from sporangium
V_big <- 1 #volume of pond
v_little <- 0.5 #fraction of zoospores that succ.infect frog skin
f <- 0.5 #fraction of zoospores that immediately reinfect same host
sigma <- 0.2 #loss rate of sporangia
gamma <- 10 # Rate of encounter between zoospores + frogs. Briggs did not provide a helpful est. here
mu <- 1 #Loss rate of zoospores in pool
delta_t <- 1
beta <- V_big * v_little
# Initial Values
N <- 100
Si <- 1
Zi <- 1

c1 <- n*v_little*f-sigma  #Number of new sporangia
c2 <- n*(1-f)*Si
c3 <- (gamma/V_big)*N+mu


disease_df <- data.frame(Sporangia = integer(),
                         Zoospores = integer())
disease_list = vector("list", length = Time)
for(i in 1:Time){
  if(i == 1){
    disease_new <- data.frame(Sporangia = Si,Zoospores = Zi)
    disease_list[[i]] <- disease_new 
  }
  else(
  Si_new <- (c1*exp(c1*delta_t)*Si+(beta/V_big)*(exp(c1*delta_t)-1)*Zi)/(c1)
  Zi_new <- (exp(-c3*delta_t)*(c2*(exp(c3*delta_t)-1)+c3*Zi))/(c3)
  disease_new <- data.frame(Sporangia = Si_new,Zoospores = Zi_new)
  disease_list[[i]] <- disease_new 
  Si <- Si_new
  Zi <- Zi_new
}

disease_df <- do.call(rbind,disease_list) 
View(disease_df)
