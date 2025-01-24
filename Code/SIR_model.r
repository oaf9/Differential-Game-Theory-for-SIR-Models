sf_data <- read.csv("""/Users/omarafifi/MyFolders/Differential-Game-Theory-for-
                    SIR-Models/Data/san_francisco.csv""")[631:731,]


N = 810000 #population size in 2021
 #percentage of infected individuals in 2021
percent_infected = sf_data['pos']/sf_data['tests']

print(percent_infected)













