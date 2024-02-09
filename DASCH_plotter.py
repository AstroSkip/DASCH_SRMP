import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#data =  np.loadtxt('/Users/conorransome/TCrB.txt', delimiter = ',')

## OR we can use


data = pd.read_csv('/Users/conorransome/TCrB.txt')

# separate out the bits of data we want

times =  data['Year']
mags = data['mag']
mag_err = data['mag_err']

# print(mags)
# print(len(mags))

## We use PYPLOT to plot. First, start a figure:

plt.figure() 

# There are lots of differnet arguments that can go in a plt figure.

#Now we can plot things! Plotting goes as plot(x,y)

plt.scatter(times, mags)
plt.errorbar(times, mags, yerr = mag_err)

# wait, we are astronomers, we use magnitues, so we need to invert the y-axis gca = grab current axis

plt.gca().invert_yaxis()

# while we are at it, we should add some axis labels!

plt.xlabel('Year')
plt.ylabel('Apparent magnitude')
#now show!

plt.show()





