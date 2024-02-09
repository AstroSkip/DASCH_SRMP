# Putting everything together!

# Import statements 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.visualization import simple_norm

# First for the lightcurve:

data = pd.read_csv('/Users/conorransome/TCrB.txt')

# separate out the bits of data we want

times =  data['Year']
mags = data['mag']
mag_err = data['mag_err']

# Now the image:

filepath = '/Users/conorransome/Downloads/cutout_rings.v3.skycell.1903.048.stk.r.unconv.fits'

# now, we can open our fits file -- the main argument of the fits.open function is the 
# file directory, so we add in the filepath variable we made

tcrb = fits.open(filepath)

# Now we can access the different level of our fits file.
# First let's see what's in the main header

header = tcrb[0].header

# how do we get to the image data?

data = tcrb[0].data

# This is now just a data array with the shape and size of the image, with entries which 
# are just the pixel values


# Let's make our figure, we can add a few arguments to resize our figure

plt.figure(figsize = (8,4))

# lets define some 'axes' which are just each plot in our subplot

ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)

# ax1 will be our light curve

ax1.scatter(times, mags, color = 'crimson', alpha  = 0.5, s = 25, edgecolor = 'black')
#ax1.errorbar(times, mags, yerr = mag_err)

# wait, we are astronomers, we use magnitues, so we need to invert the y-axis gca = grab current axis

ax1.invert_yaxis()

#  while we are at it, we should add some axis labels!

ax1.set_xlabel('Year')
ax1.set_ylabel('Apparent magnitude')

# ax2 will be our image!

# to scale, we can try some extra arguments in imshow

norm = simple_norm(data ,'log', invalid = 0,min_percent = 10, max_percent = 99)

ax2.imshow(data, cmap='magma', origin='lower' , norm = norm)

ax2.tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False)
plt.tight_layout()

plt.savefig('/Users/conorransome/TCrB.pdf')

plt.show()












