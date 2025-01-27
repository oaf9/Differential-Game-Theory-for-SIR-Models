library(ggplot2)


michigan_data <- read.csv("/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/michigan-covid-data.csv")
michigan_data <- michigan_data[order(michigan_data$date, decreasing = FALSE),]
print('DATA LOADED Successfully')

michigan_data <- michigan_data[220:373,]
N <- 10500000 #population estimate in 2021 
R <- 10*michigan_data["recovered"]
S <- N - michigan_data["positive"]*10
I <- N - R - S

michigan_data["Susceptible"] = S
michigan_data["Infected"] = I
michigan_data["Recovered"] = R
michigan_data['date'] <- as.Date(michigan_data$date)





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
       subtitle = 'Population = 10.5M',
       x = "Date",
       y = "Number of People")












