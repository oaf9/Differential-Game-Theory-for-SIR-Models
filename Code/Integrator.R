
# set a time sequence ... were doing this over 
# the length of our original interval
t = 22

times <- seq(0, t, by = 1)

S_0 <- S[1:1,]
I_0 <- I[1:1,]
R_0 <- R[1:1,]

#initial values and paramaters must be vectorized
V_0 <- c(as.numeric(S_0), as.numeric(I_0), as.numeric(R_0))

#deSolve is very sensitive to the types of the parmamaters.
paramaters <- c(beta = as.numeric(beta), 
                gamma = as.numeric(3*gamma), 
                N = as.numeric(N))

#this function describes the dynamic of the system w.r.t. t
SIR <- function(t, state, paramaters){
  with (as.list(c(state, paramaters)),{
    
      dS <- (-beta/N)*I*S
      dI <- (beta/N)*I*S - gamma*I
      dR = gamma*I
      
      #R will automatically return the last last line in the function body
      #so no return keyword is required
      list(c(dS,dI,dR))
    })
}

#We Use deSolve to solve the PDE System
library(deSolve)
require(deSolve)

#The equations must be integrated through time. 
predictions <- ode(y = V_0, times = times, func = SIR, parms = paramaters)
colnames(predictions) <- c("time", "S", "I", "R")      
        
#plot the results 
library(ggplot2)
require(ggplot2)

SIR_pred <- ggplot(data = predictions, aes(x = times))+
  geom_line(aes(y = R), 
            size = 1,
            color = "#D5C36C")+
  geom_line(aes(y = S), 
            size = 1,
            color = "#C0E4ED")+ 
  geom_line(aes(y = I), 
            size = 1,
            color = "#397689")+
  labs(title = 'SIR Predictions', 
       subtitle = 'Population = 10.2M',
       x = "Date",
       y = "Number of People")

#display results
SIR_pred

