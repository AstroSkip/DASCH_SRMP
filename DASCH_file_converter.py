import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = '/Users/conorransome/Downloads/short_15-59-30.1622_25-55-12.613_APASS_J155930.1+255512_0001.txt'

data = pd.read_csv(data, delimiter = '	', skiprows = [0,2])

print(data['year'])

year = data['year']
mag = data['magcal_magdep']
mag_err = data['magcal_local_rms']

# We can save a new file if we want... we can use pandas to concatenate these arrays

new_file = pd.concat([year, mag, mag_err], axis = 1, keys = ['year', 'mag', 'mag_err'])

new_file.to_csv('/Users/conorransome/Downloads/TCrB.csv', index = False)

# But we don't REALLY need to do that -- could be nice for later though, and we could add whatever DASCH
# columns we want

# now we have these things, we want to plot them, WITH error bars. But the errors are a little all over the place
# so we want to be able to filter out bad data

# In photometry this is quite easy to have a 'first look' at -- if the error is more than 0.4 mag or so, it's probably
# bad data -- this is also the cut DASCH uses

# We can use LOOPS to filter our data to fifure out what we want to plot -- we will start with a simple cut
# filtering out points with large uncertainties

# lets start our figure

plt.figure(figsize = (5,5))

# in a loop we use an index as well -- choose the index and it loops throught the data

for i in range(len(year)):
	if mag_err[i] < 0.4:
		plt.scatter(year[i], mag[i], color = 'crimson', alpha  = 0.5, s = 25, edgecolor = 'black')
		plt.errorbar(year[i], mag[i], color = 'black', alpha = 0.5, yerr = mag_err[i])

plt.xlabel('Year')
plt.ylabel('Apparent magnitude')

plt.gca().invert_yaxis()
plt.show()











