#A quick script to reformat the data nicely

import pandas as pd

data = pd.read_csv("/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/michigan-covid-data.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.set_index('date').sort_index()
data = data[220: ].iloc[::7]


N = 10200000 #population estimate in 2021 

R = 10*data["recovered"]
S = N - data["positive"]*10
I = N - R - S

data["S"]=  S
data["I"]=  I
data["R"]=  R

data = data[["S", "I", "R"]]

data.to_csv("/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/formatted_data.csv")
