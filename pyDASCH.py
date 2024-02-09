import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm
import ps1getter as ps1
from astropy.table import Table
import requests
import time
from io import StringIO

def DASCHconvert(path, file, filename):

	directory = path + file
	data = pd.read_csv(directory, delimiter = '	', skiprows = [0,2])

	year = data['year']
	mag = data['magcal_magdep']
	mag_err = data['magcal_local_rms']

	new_file = pd.concat([year, mag, mag_err], axis = 1, keys = ['year', 'mag', 'mag_err'])

	new_file.to_csv(path+filename, index = False)
	print('new file made! file',filename, 'made from DASCH file', file)



#DASCHconvert(path = '/Users/conorransome/Downloads/', file = 'short_15-59-30.1622_25-55-12.613_APASS_J155930.1+255512_0001.txt', filename = 'TCrBtest.csv')

path = '/Users/conorransome/Downloads/'
file = 'short_15-59-30.1622_25-55-12.613_APASS_J155930.1+255512_0001.txt'
filename = 'TCrBtest.csv'

#DASCHconvert(path = path, file = file, filename = filename)

def DASCHplot(path, file, filename, save = True, show = True):

	directory = path + file
	data = pd.read_csv(directory, delimiter = '	', skiprows = [0,2])

	year = data['year']
	mag = data['magcal_magdep']
	mag_err = data['magcal_local_rms']

	plt.figure(figsize = (5,5))

# in a loop we use an index as well -- choose the index and it loops throught the data

	for i in range(len(year)):
		if mag_err[i] < 0.4:
			plt.scatter(year[i], mag[i], color = 'crimson', alpha  = 0.5, s = 25, edgecolor = 'black')
			plt.errorbar(year[i], mag[i], color = 'black', alpha = 0.5, yerr = mag_err[i], ls = 'none')

	plt.xlabel('Year')
	plt.ylabel('Apparent magnitude')

	plt.gca().invert_yaxis()
	if save == True:
		plt.savefig(path+filename)
		print('saved figure as', filename, '!')
	if show == True:
		plt.show()

#DASCHplot(path = '/Users/conorransome/Downloads/', file = 'short_15-59-30.1622_25-55-12.613_APASS_J155930.1+255512_0001.txt', filename = 'TCrBtest.pdf', save = False, show = True)

def DASCHimplot(path, file, imfile, filename,  pixmax, pixmin, colormap, save = True, show = True):
	directory = path + file
	lc_data = pd.read_csv(directory, delimiter = '	', skiprows = [0,2])

	year = lc_data['year']
	mag = lc_data['magcal_magdep']
	mag_err = lc_data['magcal_local_rms']

	tcrb = fits.open(path+imfile)

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

	for i in range(len(year)):
		if mag_err[i] < 0.4:
			#plt.scatter(year[i], mag[i], color = 'crimson', alpha  = 0.5, s = 25, edgecolor = 'black')
			#plt.errorbar(year[i], mag[i], color = 'black', alpha = 0.5, yerr = mag_err[i])

	# ax1 will be our light curve

			ax1.scatter(year[i], mag[i], color = 'crimson', alpha  = 0.5, s = 25, edgecolor = 'black')
			ax1.errorbar(year[i], mag[i], yerr = mag_err[i], alpha = 0.5, color = 'black')

	# wait, we are astronomers, we use magnitues, so we need to invert the y-axis gca = grab current axis

	ax1.invert_yaxis()

	#  while we are at it, we should add some axis labels!

	ax1.set_xlabel('Year')
	ax1.set_ylabel('Apparent magnitude')

	# ax2 will be our image!

	# to scale, we can try some extra arguments in imshow

	norm = simple_norm(data ,'log', invalid = 0,min_percent = pixmin, max_percent = pixmax)

	ax2.imshow(data, cmap=colormap, origin='lower' , norm = norm)

	ax2.tick_params(left = False, right = False , labelleft = False , 
	                labelbottom = False, bottom = False)
	plt.tight_layout()

	if save == True:
		plt.savefig(path+filename)
		print('new plot saved!', filename, 'saved in', path)
	if show == True:
		plt.show()

#DASCHimplot('/Users/conorransome/Downloads/', 'short_15-59-30.1622_25-55-12.613_APASS_J155930.1+255512_0001.txt', 'cutout_rings.v3.skycell.1903.048.stk.r.unconv.fits', 'TCrB_implot.pdf',  99, 5, save = True, show = True, colormap = 'magma')

def DASCH_plotter(path, file, filename, imfile, tra, tdec, filt, pixmax, pixmin, colormap, save = True, show = True):
	t0 = time.time()

	directory = path + file
	lc_data = pd.read_csv(directory, delimiter = '	', skiprows = [0,2])

	year = lc_data['year']
	mag = lc_data['magcal_magdep']
	mag_err = lc_data['magcal_local_rms']

	table = ps1.getimages(tra,tdec,filters=filt)
	print("{:.1f} s: got list of {} images for {} positions".format(time.time()-t0,len(table),len(tra)))
	table.sort(['projcell','subcell','filter'])

	# extract cutout for each position/filter combination
	for row in table:
	    ra = row['ra']
	    dec = row['dec']
	    projcell = row['projcell']
	    subcell = row['subcell']
	    filter = row['filter']

	    # create a name for the image -- could also include the projection cell or other info
	    fname = "t{:08.4f}{:+07.4f}.{}.fits".format(ra,dec,filter)

	    url = row["url"]
	    print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
	    r = requests.get(url)
	    open(fname,"wb").write(r.content)
	print("{:.1f} s: retrieved {} FITS files for {} positions".format(time.time()-t0,len(table),len(tra)))

	img = fits.open(fname)
	header = img[0].header
	data = img[0].data

	plt.figure(figsize = (8,4))

	ax1 = plt.subplot(1,2,1)
	ax2 = plt.subplot(1,2,2)

	for i in range(len(year)):
		if mag_err[i] < 0.4:

			ax1.scatter(year[i], mag[i], color = 'crimson', alpha  = 0.5, s = 25, edgecolor = 'black')
			ax1.errorbar(year[i], mag[i], yerr = mag_err[i], alpha = 0.5, color = 'black')

	ax1.invert_yaxis()

	ax1.set_xlabel('Year')
	ax1.set_ylabel('Apparent magnitude')

	norm = simple_norm(data ,'log', invalid = 0,min_percent = pixmin, max_percent = pixmax)
	ax2.imshow(data, cmap=colormap, origin='lower' , norm = norm)
	ax2.tick_params(left = False, right = False , labelleft = False , 
	                labelbottom = False, bottom = False)
	plt.tight_layout()

	if save == True:
		plt.savefig(path+imfile)
		print('new plot saved!', imfile, 'saved in', path)
	if show == True:
		plt.show()


DASCH_plotter(path = path, file = file, filename = filename, imfile = 'TCrBPS1.pdf', tra= [239.8750], tdec =[25.9201389], filt = 'r', pixmax = 99, pixmin = 10, colormap = 'magma', save = True, show = True)    




print(ps1.getimages([120.0], [45.0], size=120, filters="r", format="fits", imagetypes="stack"))

t0 = time.time()
 
# create a test set of image positions
tdec = np.append(np.arange(31)*3.95 - 29.1, 88.0)
tra = np.append(np.arange(31)*12., 0.0)

# get the PS1 info for those positions
table = ps1.getimages(tra,tdec,filters="ri")
print("{:.1f} s: got list of {} images for {} positions".format(time.time()-t0,len(table),len(tra)))

# if you are extracting images that are close together on the sky,
# sorting by skycell and filter will improve the performance because it takes
# advantage of file system caching on the server
table.sort(['projcell','subcell','filter'])

# extract cutout for each position/filter combination
for row in table:
    ra = row['ra']
    dec = row['dec']
    projcell = row['projcell']
    subcell = row['subcell']
    filter = row['filter']

    # create a name for the image -- could also include the projection cell or other info
    fname = "t{:08.4f}{:+07.4f}.{}.fits".format(ra,dec,filter)

    url = row["url"]
    print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
    r = requests.get(url)
    open(fname,"wb").write(r.content)
print("{:.1f} s: retrieved {} FITS files for {} positions".format(time.time()-t0,len(table),len(tra)))
