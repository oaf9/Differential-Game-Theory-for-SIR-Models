import pandas as pd
import numpy as np
from scipy import interpolate

data = pd.read_csv("/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/michigan-covid-data.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.set_index('date').sort_index()
data = data[220: ].iloc[::7]


N = 10200000 #population estimate in 2021 

R = 10*data["recovered"]
S = N - data["positive"]*10
I = N - R - S

data["S"]=  S/N
data["I"]=  I/N
data["R"]=  R/N

f_S = interpolate.interp1d(np.arange(len(data)), data['S'], kind='quadratic', axis=-1)
f_I = interpolate.interp1d(np.arange(len(data)), data['I'], kind='quadratic', axis=-1)
f_R = interpolate.interp1d(np.arange(len(data)), data['R'], kind='quadratic', axis=-1)

S = [f_S(t) for t in np.linspace(0, len(data)-1, 152)]
I = [f_I(t) for t in np.linspace(0, len(data)-1, 152)]
R = [f_R(t) for t in np.linspace(0, len(data)-1, 152)]

output = pd.DataFrame()
output["S"] = S
output["I"] = I
output["R"] = R


data =  output

data.to_csv("/Users/omarafifi/MyFolders/Differential-Game-Theory-for-SIR-Models/Data/_normalized.csv")
