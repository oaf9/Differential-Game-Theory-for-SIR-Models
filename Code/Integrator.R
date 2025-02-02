
# set a time sequence ... were doing this over 
# the length of our original interval
t = 22
times <- seq(1, t, by = 1)

#get initial values from the data
S_0 <- as.numeric(S$positive[1])
I_0 <- as.numeric(I$recovered[1])
R_0 <- as.numeric(R$recovered[1])

#initial values and paramaters must be vectorized
state <- c(S = S_0, I = I_0, R = R_0)

# deSolve is very sensitive to the types of the parmamaters.
# since beta is estimated in terms of days, it has to be scaled since we are 
# predicting in terms of weeks
paramaters <- c(b = as.numeric(7*beta), 
                g = as.numeric(gamma), 
                n = as.numeric(N))


#this function describes the dynamic of the system w.r.t. t
SIR <- function(t, state, paramaters){
  with (as.list(c(state, paramaters)),{
    
      dS <- (-b/n)*I*S
      dI <- (b/n)*I*S - g*I
      dR = g*I
    
      #R will automatically return the last last line in the function body
      #so no return keyword is required
      return(list(c(dS, dI, dR)))
    })
}

#We Use deSolve to solve the PDE System
library(deSolve)
require(deSolve)

#The equations must be integrated through time. 
predictions <- ode(y = state, times = times, func = SIR, parms = paramaters)
colnames(predictions) <- c("time", "S", "I", "R")      

#Creating a datetime vector for plotting the results
start_time <- as.Date("2020-10-06")
date_times <- start_time + (7*times)

#plot the results 
library(ggplot2)
require(ggplot2)

#
SIR_pred <- ggplot(data = predictions, aes(x = date_times))+
  geom_line(aes(y = R, color = "Recovered"), 
            size = 2)+
  geom_line(aes(y = S, color = "Susceptible"), 
            size = 2)+ 
  geom_line(aes(y = I, color = "Infected"), 
            size = 2)+
  labs(title = 'SIR Predictions', 
       subtitle = 'Population = 10.2M',
       x = "Date",
       y = "Number of People", 
       color = "Category")+ 
   scale_color_manual(values = c("Recovered" = "#D5C36C", 
                                 "Susceptible" = "#C0E4ED", 
                                 "Infected" = "#397689"))

#display results
SIR_pred


#plot the difference between actual and predicted values: 
#lets create a New data frame with real and predicted values: 
output <- cbind(I, predictions)
#claculate RSS

ggplot(data = output, aes(x = date_times))+
  geom_line(aes(y = I, color = "Predicted Counts"), 
  size = 2, linetype  = 'dashed')+
  geom_line(aes(y = recovered, color = "Actual Infection Counts"), 
            size = 2)+ 
  scale_color_manual(values = c("Predicted Counts" = "red", 
                                "Actual Infection Counts" = "black"))+ 
  labs(title = "Predicted vs. Actual Infections", 
       subtitle = TeX("Paramaters: $ \\beta = .784, \\gamma = .107, N = 10.2M, 
                      I_0 = 44,357, S_0 = 8,761,220, R_0 = 995$"),
       x = "Month",
       y = "Number of Infections")

