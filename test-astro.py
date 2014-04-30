#Aperture production with Astropy, Bresenham circle algorithm 
#Documentation can be found here:
#http://docs.astropy.org/en/stable/io/fits/index.html
#Started 4/27/2014
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

NUM_ROWS = 32
NUM_COLS = 32

#Converts frames into a manageable unit: MegaJanskys per Arc
def flux_convert(hdu_list, frame_num):
	pixel_rows = hdu_list[0].header['pxscal1']  # arcseconds per pixel (1D)
	pixel_columns = hdu_list[0].header['pxscal2'] # arcseconds per pixel (1D)
	pixel = pixel_rows*pixel_columns  # arcseconds^2 per pixel (2D)
	converted_frame = []
	for row in hdu_list[0].data[frame_num]:
		n_row = []
		for element in row:
			MJy_per_arc2 = element*2.35045E-11
			MJy_per_pix = MJy_per_arc2*pixel
			n_row.append(MJy_per_pix)
		converted_frame.append(n_row)
	return np.array(converted_frame)

#Takes numpy array, returns index of brightest pixel
def brightest_pixel(frame):
	max_ind = (0, 0)
	brightest = 0
	for i in range(len(frame)):
		for j in range(len(frame[i])):
			if frame[i][j] > brightest:
				brightest = frame[i][j]
				max_ind = (i, j)
	return max_ind

def frame_neighbors(frame, index_tuple):
	i = index_tuple[0]
	j = index_tuple[1]
	if i != 0 and j != 0 and i != (NUM_ROWS - 1) and j != (NUM_COLS - 1):
		A = frame[i - 1, j - 1]
		B = frame[i, j - 1]
		C = frame[i + 1, j - 1]
		D = frame[i - 1, j]
		E = frame[i + 1, j]
		F = frame[i - 1, j + 1]
		G = frame[i, j + 1]
		H = frame[i + 1, j + 1]
		return [A, B, C, D, E, F, G, H]
	elif i == 0 and j != 0 and i != (NUM_ROWS - 1) and j != (NUM_COLS - 1):
		B = frame[i, j - 1]
		C = frame[i + 1, j - 1]
		D = frame[i + 1, j]
		E = frame[i, j + 1]
		F = frame[i + 1, j + 1]
		return [B,C,D,E,F]
	elif i != 0 and j == 0 and i != (NUM_ROWS - 1) and j != (NUM_COLS - 1):
		D = frame[i - 1, j]
		E = frame[i + 1, j]
		F = frame[i - 1, j + 1]
		G = frame[i, j + 1]
		H = frame[i + 1, j + 1]
		return [D, E, F, G, H]
	elif i == 0 and j == 0 and i != (NUM_ROWS - 1) and j != (NUM_COLS - 1):
		D = frame[i + 1, j]
		E = frame[i, j + 1]
		F = frame[i + 1, j + 1]
		return [D, E, F]
	elif i == (NUM_ROWS - 1) and j != (NUM_COLS - 1):
		A = frame[i - 1, j - 1]
		B = frame[i, j - 1]
		C = frame[i - 1, j]
		D = frame[i - 1, j + 1]
		E = frame[i, j + 1]
		return [A, B, C, D, E]
	elif i != (NUM_ROWS - 1) and j == (NUM_COLS - 1):
		A = frame[i - 1, j - 1]
		B = frame[i, j - 1]
		C = frame[i + 1, j - 1]
		D = frame[i - 1, j]
		E = frame[i + 1, j]
		return [A, B, C, D, E]
	elif i == (NUM_ROWS - 1) and j == (NUM_COLS - 1):
		A = frame[i - 1, j - 1]
		B = frame[i, j - 1]
		D = frame[i - 1, j]
		return [A, B, D]

#Takes numpy array, index of pixel in question, and integer number of standard deviations
#If pixel in question and its immediate neighbors are n SD outside of the mean, the pixel is considered 'significant'
def is_significant(frame, index_tuple, deviations = 1, mean = -1, SD = -1):
	if mean == -1 or SD == -1:
		frame_mean = np.mean(frame)
		frame_SD = np.std(frame)
	else:
		frame_mean = mean
		frame_SD = SD
	lower_bound = frame_mean + deviations * frame_SD
	if frame[index_tuple] > lower_bound:
		neighbors = frame_neighbors(frame, index_tuple)
		for nbr in frame_neighbors:
			if frame[nbr] <= lower_bound:
				return False
		return True
	return False

#Returns data with items outside circle (of given radius) masked by zero
#Uses the Bresenham Circle Algorithm to draw nice looking circles in taxicab geometry
def mk_circle(frame, center_tuple, radius):
	(x0, y0) = center_tuple
	new_frame = np.array([[0]*32]*32)
	x = radius
	y = 0
	radiusError = 1 - x
	while x >= y:
		#print new_frame[0,1]
		new_frame[x+x0,y+y0] = frame[x+x0, y+y0]
		new_frame[y+x0,x+y0] = frame[y+x0, x+y0]
		new_frame[-x+x0,y+y0] = frame[-x+x0, y+y0]
		new_frame[-y+x0,x+y0] = frame[-y+x0, x+y0]
		new_frame[-x+x0,-y+y0] = frame[-x+x0, -y+y0]
		new_frame[-y+x0,-x+y0] = frame[-y+x0, -x+y0]
		new_frame[x+x0,-y+y0] = frame[x+x0, -y+y0]
		new_frame[y+x0,-x+y0] = frame[y+x0, -x+y0]
		y += 1
		if radiusError < 0:
			radiusError += 2*y + 1
		else:
			x -= 1
			radiusError += 2*(y - x + 1)
	return new_frame

#Adaptation of Bresenham Circle Algorithm for disks
def mk_disk(frame, center_tuple, radius):
	new_frame = []
	frame_list_form = np.ndarray.tolist(frame)
	circle_frame = np.ndarray.tolist(mk_circle(frame, center_tuple, radius))
	for row_num in range(len(circle_frame)):
		indA, indB = 16, 16
		for i in range(16):
			if circle_frame[row_num][i] != 0:
				indA = i
				break
		for j in range(31,16,-1):
			if circle_frame[row_num][j] != 0:
				indB = j
				break
		n_row = circle_frame[row_num][0:indA] + frame_list_form[row_num][indA:indB] + circle_frame[row_num][indB:]
		new_frame.append(n_row)
	return np.array(new_frame)

#Print a fits frame as a pretty array of numbers
def print_frame(frame, spacing = 1, num_sigfigs = 0):
	for row in np.ndarray.tolist(frame):		
		for element in row:
			if element != 0:
				a = len(str(round(element,num_sigfigs)))
				print str(round(element,num_sigfigs)) + " " * (spacing-a),
			else:
				print "0", " " * (spacing-1),
		print ""

#Returns numerical measurement describing quality of mask
def aperture_metric(frame, mask):
	#Qualities we want in a mask:
		#Conserves valuable information.
			#What information is valuable?
			#Information in a moderately sized cluster of similar values;
			#High pixel spread; high brightness
			#Possibility: find a bright pixel whose neighbors have luminosities above the average
			#This can be recursive up to some number of levels (must be finite)
			#Such a pixel is a good candidate for being non-noise, depending on recursion level
	pass

def test_aperture(frame, center_tuple, radius):
	pass

def main():
	fits_file = "Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits"
	hdulist = fits.open(fits_file)
	frame_list = hdulist[0].data
	frame_one = frame_list[0]
	
	#Prints a version of the star with center at brightest pixel, and radius r
	max_pixel = brightest_pixel(frame_one)
	test_disk = mk_disk(frame_one, max_pixel, 3)
	print_frame(test_disk, spacing = 3)
	
	#subprocess.call(["ds9", "-zoom","8", fits_file])

if __name__ == "__main__":
	main()