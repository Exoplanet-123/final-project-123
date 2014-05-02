#Aperture production with Astropy, Bresenham circle algorithm 
#Documentation can be found here:
#http://docs.astropy.org/en/stable/io/fits/index.html
#Started 4/27/2014
from astropy.io import fits
import scipy
import ds9
import subprocess
import numpy as np
import matplotlib.pyplot as plt
#import sklearn.cluster as sk
#from sklearn import metrics
#from sklearn.cluster import KMeans
#from sklearn.datasets import load_digits
#from sklearn.decomposition import PCA
#from sklearn.preprocessing import scale

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

#Find 2x2 area with greatest total luminosity
#Return index of brightest pixel in the 2x2 area
def brightest_region(frame):
	N = NUM_ROWS
	M = NUM_COLS
	# Since we're looking for 2x2 regions using the top left index, 
	# we will only search N-1 rows and M - 1 columns.
	max_lum = 0
	best_ind = (0, 0) 
	for i in range(N - 1):
		for j in range(M - 1):
			pix_A = frame[i, j]
			pix_B = frame[i + 1, j]
			pix_C = frame[i, j + 1]
			pix_D = frame[i + 1, j + 1]
			z = [pix_A, pix_B, pix_D, pix_D]
			total_lum = sum(z)
			if total_lum > max_lum:
				max_lum = total_lum
				if pix_A == max(z):
					best_ind = (i, j)
				if pix_B == max(z):
					best_ind = (i + 1, j)
				if pix_C == max(z):
					best_ind = (i, j + 1)
				else:
					best_ind = (i + 1, j + 1)
	return best_ind

#Returns up to 8 neighbors of a pixel at radius 1 
#Handles edge cases silently
#Gives neighbors in a square, NOT a circle
def frame_neighbors(frame, index_tuple):
	i = index_tuple[0]
	j = index_tuple[1]
	min_bd = 0
	level = 1
	left_bd = NUM_ROWS - 1
	right_bd = NUM_COLS - 1
	#Handling the twelve standard edge cases for a grid
	if i != min_bd and j != min_bd and i != left_bd and j != right_bd:
		A = frame[i - level, j - level]
		B = frame[i, j - level]
		C = frame[i + level, j - level]
		D = frame[i - level, j]
		E = frame[i + level, j]
		F = frame[i - level, j + level]
		G = frame[i, j + level]
		H = frame[i + level, j + level]
		return [A, B, C, D, E, F, G, H]
	elif i == min_bd and j != min_bd and i != left_bd and j != right_bd:
		B = frame[i, j - level]
		C = frame[i + level, j - level]
		D = frame[i + level, j]
		E = frame[i, j + level]
		F = frame[i + level, j + level]
		return [B,C,D,E,F]
	elif i != min_bd and j == min_bd and i != left_bd and j != right_bd:
		D = frame[i - level, j]
		E = frame[i + level, j]
		F = frame[i - level, j + 1]
		G = frame[i, j + level]
		H = frame[i + level, j + level]
		return [D, E, F, G, H]
	elif i == min_bd and j == min_bd and i != left_bd and j != right_bd:
		D = frame[i + level, j]
		E = frame[i, j + level]
		F = frame[i + level, j + level]
		return [D, E, F]
	elif i != min_bd and j != min_bd and i == left_bd and j != right_bd:
		A = frame[i - level, j - level]
		B = frame[i, j - level]
		C = frame[i - level, j]
		D = frame[i - level, j + level]
		E = frame[i, j + level]
		return [A, B, C, D, E]
	elif i != min_bd and j != min_bd and i != left_bd and j == right_bd:
		A = frame[i - level, j - level]
		B = frame[i, j - level]
		C = frame[i + level, j - level]
		D = frame[i - level, j]
		E = frame[i + level, j]
		return [A, B, C, D, E]
	elif i != min_bd and j != min_bd and i == left_bd and j == right_bd:
		A = frame[i - level, j - level]
		B = frame[i, j - level]
		C = frame[i - level, j]
		return [A, B, C]
	elif i == min_bd and j != min_bd and i != left_bd and j == right_bd:
		A = frame[i, j - level]
		B = frame[i + level, j - level]
		D = frame[i + level, j]
		return [A, B, D]
	elif i != min_bd and j == min_bd and i == left_bd and j != right_bd:
		A = frame[i - level, j]
		B = frame[i - level, j + level]
		C = frame[i, j + level]
		return [A, B, C]
	else:
		print "Error in frame neighbor computer."

#Returns neighbors of a pixel at arbitrary radius
#Neighbors are returned in a circle rather than a square
def frame_neighbors_n(frame, index_tuple, level=1):
	if level == 0:
		return [frame[index_tuple]]
	if level == 1:
		return frame_neighbors(frame, index_tuple)
	masked_frame = mk_circle(frame, index_tuple, level)
	rv = masked_frame[np.nonzero(masked_frame)]
	return rv

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
		for nbr in neighbors:
			if nbr <= lower_bound:
				return False
		return True
	return False

#Returns data with items not lying on circle (of given radius) masked by zero
#Uses the Bresenham Circle Algorithm to draw nice looking circles in taxicab geometry
def mk_circle(frame, center_tuple, radius):
	(x0, y0) = center_tuple
	new_frame = np.array([[0.0]*32]*32)
	x = float(radius)
	y = 0.0
	radiusError = 1.0 - x
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
			radiusError += 2.0*y + 1.0
		else:
			x -= 1
			radiusError += 2.0*(y - x + 1.0)
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
				print "0", " " * (spacing - 1),
		print ""

#Plots light curve with respect to a given pixel as the origin
#Saves image in images folder
def light_curve(frame, index_tuple, max_radius = 5):
	x_points = []
	y_points = []
	x = 0
	for rad in range(max_radius + 1):
		z = frame_neighbors_n(frame, index_tuple, level = rad)
		print z
		avg = (sum(z))/(float(len(z)))
		x_points.append(x)
		y_points.append(avg)
		x += 1
	plt.plot(x_points, y_points, 'r--')
	plt.title('Light Curve from pixel ' + str(index_tuple))
	plt.xlabel('Radius (pixels)')
	plt.ylabel('Average luminosity')
	plt.savefig("images/light_curve" + str(index_tuple[0]) + "-" + str(index_tuple[0]) + ".png")
	plt.show()

#Should return numerical measurement describing quality of mask (in progress)
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
	
	#Prints a version of the star with center at top left of brightest region, and radius r
	max_pixel = brightest_region(frame_one)
	test_disk = mk_circle(frame_one, (0,0), radius = 2)
	#print_frame(test_disk, spacing = 2, num_sigfigs = 0)
	light_curve(frame_one, max_pixel)
	#print "Max pixel: ",max_pixel
	#print is_significant(frame_one, max_pixel, deviations = 1)
	
	#subprocess.call(["ds9", "-zoom","8", fits_file])

if __name__ == "__main__":
	main()