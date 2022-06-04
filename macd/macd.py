import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib import patches

plt.close("all")

def calculateEMA(data, period):
    alpha = 2/(period+1)
    factor = 1 - alpha
    ema = []
    for i in range(0, len(data)):
        numerator = 0
        denominator = 0
        for j in range (0, period+1):
            if (i - j >= 0):
                numerator += data.iloc[i-j, 1] * pow(factor, j)
                denominator += pow(factor, j)
        ema.append(numerator/denominator)
    return ema

def calculateMACD(data):
    # calculate EMA
    ema = calculateEMA(data, 12)
    data['ema12'] = ema
    ema = calculateEMA(data, 26)
    data['ema26'] = ema
    macd = []
    for i in range(0, len(data)):
        macd.append(data.iloc[i, 2] - data.iloc[i, 3])

    del data["ema12"]
    del data["ema26"]
    data['MACD'] = macd
    return data

def calculateSignal2(data):
    ema = calculateEMA(data.drop(columns=["Otwarcie"]), 9)
    data['SIGNAL'] = ema
    return data

def buy(units, money, price):
    units += money//price 
    money %= price
    return units, money

def sell(units, money, price):    
    money += price * units
    units = 0   
    return units, money

def invest(data, units):
    money = 0
    macdOverSignal = True
    macdOverZero = True
    startMoney = data.iat[0,1]

    if data.iat[1,2] < data.iat[1,3]: # macd lower than signal
        macdOverSignal = False
    if data.iat[1,2] < 0: # macd at start lower than 0
        macdOverZero = False
    
    for i in range(2, len(data)):
        if data.iat[i,2] > data.iat[i,3] and not macdOverSignal: # macd crosses signal from bottom
            if not units:
                units, money = buy(units, money, data.iat[i,1])
            macdOverSignal = True
        elif data.iat[i,2] < data.iat[i,3] and macdOverSignal: # macd crosses signal from top
            if units:
                units, money = sell(units, money, data.iat[i,1])
            macdOverSignal = False
        elif data.iat[i,2] < 0 and macdOverZero: # macd goes under 0
            units, money = sell(units, money, data.iat[i,1])
            macdOverZero = False
        elif data.iat[i,2] > 0 and not macdOverZero: # macd goes above 0
            units, money = buy(units, money, data.iat[i,1])
            macdOverZero = True
    print("Różnica pieniężna: " ,((units*data.iat[len(data)-1, 1]+money) - (1000*startMoney)), "\nStosunek: " ,(units*data.iat[len(data)-1, 1]+money)/(1000*startMoney))
    return money

data = pd.read_csv("^ipc_d.csv")
data = data.drop(columns=["Najwyzszy","Najnizszy","Zamkniecie","Wolumen"])
data = calculateMACD(data)
data = calculateSignal2(data)

invest(data, 1000)

fig,ax1= plt.subplots(2, sharex=True)

date = []
for i in range(len(data)):
    date.append(datetime.strptime(data.to_numpy()[i,0], '%Y-%m-%d'))

green_patch = patches.Patch(color='green', label='Otwarcie')
red_patch = patches.Patch(color='red', label='MACD')
blue_patch = patches.Patch(color='blue', label='SIGNAL')
ax1[0].legend(handles=[green_patch,blue_patch,red_patch])

ax1[0].plot(date, data.to_numpy()[:,1], 'g-', label='Otwarcie')
ax1[1].plot(date, data.to_numpy()[:,2], 'r-', label='MACD')
ax1[1].plot(date, data.to_numpy()[:,3], 'b-', label='SIGNAL')
ax1[1].set_xlabel("Data")
ax1[1].set_ylabel("Wartosc")
ax1[0].grid("both")
ax1[1].grid("both")
plt.show()