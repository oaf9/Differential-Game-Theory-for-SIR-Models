library(ggplot2)
library(latex2exp)

path <- "/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/michigan-covid-data.csv"

michigan_data <- read.csv(path)
michigan_data <- michigan_data[order(michigan_data$date, decreasing = FALSE),]
print('DATA LOADED Successfully')


# the data is given in units of millions. Hence everything must be scaled by 10
# since the MI population is ~10.2
michigan_data <- michigan_data[c(seq(220, nrow(michigan_data), 7)),]
N <- 10200000 #population estimate in 2021 
R <- 10*michigan_data["recovered"]
S <- N - michigan_data["positive"]*10
I <- N - R - S

#adding the columns back into the data frame ... this makes plotting a little
#easier 
michigan_data["Susceptible"] = S
michigan_data["Infected"] = I
michigan_data["Recovered"] = R
michigan_data['date'] <- as.Date(michigan_data$date)
michigan_data["log_infected"] <- log(I)

#plotting the SIR data
sir_data_plot <- ggplot(data = michigan_data, aes(x = date))+
  geom_line(aes(y = Recovered), 
            size = 1,
            color = "#D5C36C")+
  geom_line(aes(y = Susceptible), 
            size = 1,
            color = "#C0E4ED")+ 
  geom_line(aes(y = Infected), 
            size = 1,
            color = "#397689")+
  scale_x_date(date_breaks = "1 month")+
  labs(title = 'SIR Data For The State of Michigan', 
       subtitle = 'Population = 10.2M',
       x = "Date",
       y = "Number of People")

sir_data_plot

#gamma is estimated over a period of one week. 
gamma <- (R[2,] - R[1,])/I[1,]

#plotting the regression for the logarithm of the data
ggplot(data= michigan_data[1:10,], aes(x = date, y = log_infected))+
  geom_point(size = 3,
             color = "#C0E6ED")+
  geom_smooth(method = "lm",
              size = 1,
              color = "#397689", 
              se = 0)+
  labs(title = TeX("Regression Line for Estimating $\\beta - \\gamma$"), 
       subtitle = "Michigan Data",
       x = "Date",
       y = TeX("$log(I(t)$"))
  
  
#solving for Beta
beta_mod <- lm(log_infected~date, data = michigan_data)
beta_minus_gamma <- beta_mod$coefficients[2]
beta <- gamma + beta_minus_gamma




