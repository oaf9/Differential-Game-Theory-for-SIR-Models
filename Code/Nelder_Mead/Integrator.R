
# set a time sequence ... were doing this over 
# the length of our original interval
t_0 = 22
times_0 <- seq(1, t_0, by = 1)

#get initial values from the data
S_0 <- as.numeric(S$positive[1])
I_0 <- as.numeric(I$recovered[1])
R_0 <- as.numeric(R$recovered[1])

N <- N

#initial values and paramaters must be vectorized
V_0 <- c(S = S_0, I = I_0, R = R_0)


# deSolve is very sensitive to the types of the parmamaters.
# since beta is estimated in terms of days, it has to be scaled since we are 
# predicting in terms of weeks

#params_0 <- c(b = as.numeric(beta), g = as.numeric(gamma))
params_0 <- c(b =.588, g = .222)



#this function describes the dynamic of the system w.r.t. t
SIR <- function(times, state, paramaters){
  with (as.list(c(state, paramaters)),{
    
      dS <- (-b/N)*I*S
      dI <- (b/N)*I*S - g*I
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

predictions <- ode(y = V_0, 
                   times = times_0, 
                   func = SIR, 
                   parms = params_0)


colnames(predictions) <- c("time", "S", "I", "R")      

#Creating a datetime vector for plotting the results
start_time <- as.Date("2020-10-06")
date_times <- start_time + (7*times_0)

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
#calculate RSS

I_predictions_v_reality <- ggplot(data = output, aes(x = date_times))+
  geom_line(aes(y = I, color = "Predicted Counts"), 
  size = 2, linetype  = 'dashed')+
  geom_line(aes(y = recovered, color = "Actual Infection Counts"), 
            size = 2)+ 
  scale_color_manual(values = c("Predicted Counts" = "red", 
                                "Actual Infection Counts" = "black"))+ 
  labs(title = "Predicted vs. Actual Infections", 
       subtitle = TeX("Paramaters: $ \\beta = .784, \\gamma = .107, N = 10.2M, 
                      I_0 = 44,357, S_0 = 8,761,220, R_0 = 995,210$"),
       x = "Month",
       y = "Number of Infections")


#grid search for optimal beta and gamma: 

#first, we need to create a function that 
#evaluates RSS for a given value of beta and gamma: 

RSS <- function(y_pred, y){
  
  sum((y-y_pred)*(y-y_pred))
  
}

SIR_loss <- function(paramaters, state, times, I_true){

  model_out <- ode(y = state, 
                   times = times, 
                   func = SIR, 
                   parms = paramaters)
 
  I_pred <- model_out[,3]
  
  return (RSS(I_pred, I_true))
}

loss_wrt_params <- function(paramaters){
  
  state <-  c(S = 8761220, I = 443570, R = 995210)
  N <- 10200000
  t = 22
  times <- seq(1, t, by = 1)
  
  I_true <- I[,1]
  
  SIR_loss(paramaters = paramaters, 
           state = state, 
           times = times, 
           I_true = I_true)
}



#get the fitted values from Nelder Mead
x_optimal <- optim(par = params_0, 
                  fn = loss_wrt_params,
                  method = "Nelder-Mead")
x_optimal <- x_optimal$par


#plot the fitted output
fitted_predictions <- ode(y = V_0, 
                      times = times_0, 
                      func = SIR, 
                      parms = x_optimal)

fitted_output <- cbind(I, fitted_predictions)



fitted_predictions_plot <- ggplot(data = fitted_output, aes(x = date_times))+
  geom_line(aes(y = I, color = "Predicted Counts"), 
            size = 2, linetype  = 'dashed')+
  geom_line(aes(y = recovered, color = "Actual Infection Counts"), 
            size = 2)+ 
  scale_color_manual(values = c("Predicted Counts" = "red", 
                                "Actual Infection Counts" = "black"))+ 
  labs(title = "Predicted vs. Actual Infections", 
       subtitle = TeX("Paramaters: $ \\beta = .588, \\gamma = .222, N = 10.2M, 
                      I_0 = 44,357, S_0 = 8,761,220, R_0 = 995,210$"),
       x = "Month",
       y = "Number of Infections")

fitted_predictions_plot



















