load(ggplot)

sf_data <- read.csv("/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/san_francisco.csv")
sf_data <- sf_data[631:731,]

N = 810000 #population size in 2021
#percentage of infected individuals in 2021
 
I = N*sf_data['pct']














