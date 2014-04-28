#Testing fits I/O with astropy
#Documentation can be found here:
#http://docs.astropy.org/en/stable/io/fits/index.html
#4/27/2014
from astropy.io import fits
import scipy
import ds9
import subprocess
import numpy as np
import pylab as pl
import sklearn.cluster as sk
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

np.set_printoptions(threshold=np.nan)

hdulist = fits.open("Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits")
#print hdulist.info()
#print "Observer: " + hdulist[0].header['OBSRVR']
scidata = hdulist[0].data
#print "Shape: ", scidata.shape
#print "Number of frames: ", len(scidata)


#print "First pixel in first frame: ", frame_one[0][0]
#subprocess.call(["ds9", "-zoom","8", \
#					"Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits"])

frame_one = scidata[0]

#returns index of brightest pixel
def brightest_pixel(frame):
	max_ind = (0, 0)
	brightest = 0
	for i in range(len(frame)):
		for j in range(len(frame[i])):
			if frame[i][j] > brightest:
				brightest = frame[i][j]
				max_ind = (i, j)
	return max_ind

#returns data with items outside circle (of given radius) masked
def mk_circle(frame, center_tuple, radius):
	(x, y) = center_tuple
	new_frame = np.array([[0]*32]*32)
	h_mid = 16
	counter = 0
	inc = radius/32
	for i in range(32 - radius, y):
		new_frame[i,(x - counter) : (x + counter)] = frame[i, (x - counter) : (x + counter)]
		counter += 1
	counter -= 1
	for i in range(32 - radius, 31):
		new_frame[i,(x - counter) : (x + counter)] = frame[i, (x - counter) : (x + counter)]
		counter -= 1
	return new_frame

print mk_circle(frame_one, (16,16), 10)
		

		

	
	


