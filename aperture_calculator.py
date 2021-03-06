#################################################################
# Aperture calculation with MPI, Numpy, and Astropy 
# This script reads in FITS files and outputs apertures using a PSF approximation
# Astropy documentation can be found here: http://docs.astropy.org/en/stable/io/fits/index.html
# Project started 4/27/2014
# Team: Hannah-Diamond Lowe, Zakir Gowani
################################################################
from astropy.io import fits
import scipy
import subprocess
import numpy as np
import math
import matplotlib.pyplot as plt
try:
	import Image
except:
	print "Image library load failed. Contact RCC or load different Python version"
import random

NUM_ROWS = 32					# The dimensions of a FITS image
NUM_COLS = 32

#Spitzer parameters:
FOCAL_LENGTH = 10.2				# In meters
DIAMETER = 0.85					# In meters
F_NUM = FOCAL_LENGTH/DIAMETER
WAVELENGTH = 3.6				# In microns

####################################################################################
# Interpreting and processing raw FITS data #
####################################################################################

def brightest_region(frame):
	"""Find 2x2 area with greatest total luminosity
	Return index of brightest pixel in the 2x2 area
	Finding the brightest pixel this way eliminates the concern for a false 
	positive brightest pixel, which could occur with a spike in detector noise
	"""
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
			z = [pix_A, pix_B, pix_C, pix_D]
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


def avg_brightest_region(frame_list):
	"""Returns the index of the brightest region across a list of frames"""
	brightest_indices = []
	#avg_brightest_ind = (0, 0)
	for frame in frame_list:
		brightest_indices.append(brightest_region(frame))
	brightest_ind_total = [sum(x) for x in zip(*brightest_indices)]
	avg_brightest_ind = (brightest_ind_total[0]/len(frame_list), brightest_ind_total[1]/len(frame_list))
	return avg_brightest_ind

def mk_circle(frame, center_tuple, radius):	
	"""Returns data with items not lying on circle (of given radius) masked by zero
	Uses the Bresenham Circle Algorithm to draw nice looking circles in taxicab geometry
	"""
	(x0, y0) = center_tuple
	new_frame = np.array([[0.0]*NUM_COLS] * NUM_ROWS)
	if radius == 1:
		i = x0
		j = y0
		min_bd = 0
		level = 1
		left_bd = NUM_ROWS - 1
		right_bd = NUM_COLS - 1
		z = [(i - level, j - level),(i, j - level),(i + level, j - level),(i - level, j),(i + level, j),\
		(i - level, j + level), (i, j + level),(i + level, j + level)]
		for tuple in z:
			try:
				new_frame[tuple] = frame[tuple]
			except:
				pass
		return new_frame
	x = float(radius)
	y = 0.0
	radiusError = 1.0 - x
	while x >= y:
		#print new_frame[0,1]
		new_frame[x + x0, y + y0] = frame[x + x0, y + y0]
		new_frame[y + x0, x + y0] = frame[y + x0, x + y0]
		new_frame[-x + x0, y + y0] = frame[-x + x0, y + y0]
		new_frame[-y + x0, x + y0] = frame[-y + x0, x + y0]
		new_frame[-x + x0, -y + y0] = frame[-x + x0, -y + y0]
		new_frame[-y + x0, -x + y0] = frame[-y + x0, -x + y0]
		new_frame[x + x0, -y + y0] = frame[x + x0, -y + y0]
		new_frame[y + x0, -x + y0] = frame[y + x0, -x + y0]
		y += 1
		if radiusError < 0:
			radiusError += 2.0*y + 1.0
		else:
			x -= 1
			radiusError += 2.0*(y - x + 1.0)
	return new_frame

def mk_disk(frame, center_tuple, radius):
	"""Adaptation of Bresenham Circle Algorithm for disks"""	
	new_frame = []
	circle_frame = mk_circle(frame, center_tuple, radius)
	for row_num in range(NUM_COLS):
		indA, indB = NUM_ROWS/2, NUM_COLS/2 
		for i in range(NUM_ROWS/2):
			if circle_frame[row_num, i] != 0:
				indA = i
				break
		for j in range(NUM_COLS - 1, NUM_COLS/2, -1):
			if circle_frame[row_num, j] != 0:
				indB = j
				break
		n_row = np.concatenate([circle_frame[row_num, 0:indA], frame[row_num, indA:indB], circle_frame[row_num, indB:]])
		new_frame.append(np.atleast_2d(n_row))
	return np.vstack(tuple(new_frame))

def punch_hole(frame, center_tuple, radius):
	"""Inversion of the mk_disk function (zeros in the center)"""
	new_frame = []
	#Circle_frame has radius one pixel larger than what we want
	circle_frame = mk_circle(frame, center_tuple, radius + 1)
	for row_num in range(NUM_COLS):
		indA, indB = NUM_ROWS/2, NUM_COLS/2
		for i in range(NUM_ROWS/2):
			if circle_frame[row_num, i] != 0:
				indA = i
				break
		for j in range(NUM_COLS - 1, NUM_COLS/2, -1):
			if circle_frame[row_num, j] != 0:
				indB = j
				break
		n_row = np.concatenate([frame[row_num, 0:indA], circle_frame[row_num, indA:indB], frame[row_num, indB:]])
		new_frame.append(np.atleast_2d(n_row))
	return np.vstack(tuple(new_frame))

def frame_merge(frame1, frame2):
	"""Merges frames, replacing zero elements in one frame with nonzero elements from the other
	This function assumes no overlapping elements"""
	frame1[frame1 == 0] = frame2[frame1 == 0]
	return frame1	

def thick_ring(frame, center_tuple, inner_rad, outer_rad):
	"""Returns frame completely masked (by zeros) except for a ring with specified inner and outer radii"""
	new_frame = []
	#Combine the punch_hole and mk_disk functions
	outer_circle = mk_circle(frame, center_tuple, outer_rad)
	inner_circle = mk_circle(frame, center_tuple, inner_rad)
	cf = frame_merge(inner_circle, outer_circle) #cf stands for circle_frame
	for row_num in range(NUM_ROWS):
		indA, indB = NUM_ROWS/2, NUM_COLS/2
		for i in range(NUM_ROWS - 1):
			if cf[row_num,i] != 0:
				indA = i
				break
		for j in range(NUM_COLS - 1, NUM_COLS/2, -1):
			if cf[row_num,j] != 0:
				indB = j
				break
		indC, indD = indA, indB
		for i in range(indA + 2, NUM_ROWS-1):
			if cf[row_num,i] != 0:
				indC = i
				break
		for j in range(indB - 2, 0, -1):
			if cf[row_num,j] != 0:
				indD = j
				break
		if indD < indC:
			n_row = np.concatenate([cf[row_num,0:indA], frame[row_num,indA:indB], cf[row_num,indB:]])
			new_frame.append(n_row)
			continue
		n_row = np.concatenate([cf[row_num,0:indA], frame[row_num,indA:indC], cf[row_num,indC:indD], frame[row_num,indD:indB], cf[row_num,indB:]])
		new_frame.append(np.atleast_2d(n_row))
	return np.vstack(tuple(new_frame))

def print_frame(frame, spacing = 1, num_sigfigs = 0):
	"""Print a fits frame as a pretty array of numbers"""
	for row in frame:		
		for element in row:
			if element != 0:
				a = len(str(round(element,num_sigfigs)))
				print str(round(element,num_sigfigs)) + " " * (spacing-a),
			else:
				print "0", " " * (spacing - 1),
		print ""

def numpy2image(new_filename, frame):
	"""Writes frame to as an image to a png (ideally) file"""
	#Normalize frame luminosities to list in range 0 to 255
	frame1 = frame/np.max(np.abs(frame))
	frame1 *= 255
	im = Image.fromarray(np.uint8(frame1))
	im.save(new_filename)
	new_img = Image.open(new_filename)
	new_img = new_img.resize((288, 288), Image.BILINEAR)
	new_img.save(new_filename)

##################################################################################
# Working with a FITS frame to get the best aperture size #
##################################################################################

def mk_avg_frame(frame_list):
	"""Returns a frame (32x32 array) that averages all 64 data frames in each fits file together """
	avg_frame = np.zeros((NUM_ROWS, NUM_COLS), dtype=float)
	for frame in frame_list:
		for i in range(NUM_ROWS):
			for j in range(NUM_COLS):
				avg_frame[i][j] += frame[i][j]
	return avg_frame/len(frame_list)

def calculate_frame_flux(frame, avg_max_pixel, radius):
	"""Returns the total flux in a disk created from the frame argument"""
	#disk = mk_disk(frame, avg_max_pixel, radius)
	total_signal_flux = 0
	for annulus in range(radius):
		weight = weighting_function(annulus)
		annulus_sum = np.sum(thick_ring(frame, avg_max_pixel, annulus, annulus+1))
		total_signal_flux += weight*annulus_sum
	return total_signal_flux

def weighting_function(radius):
	"""Returns a weight by radius, based on a Gaussian function approximating a point spread function (equation from Wikipedia)
	Lower radii (closer in to the star) get a higher weight than those farther away"""
	return math.exp( -(radius * radius) / (2 * 0.45 * WAVELENGTH * F_NUM) )

def penalty_function(radius):
	"""Returns a penalty for using larger radii"""
	return math.exp(radius - WAVELENGTH)

def calculate_noise_flux(frame, avg_max_pixel, radius):
	"""Returns the total noise in a frame (anything in a frame that is not the signal)
	Total flux is calculated by taking into account annuli of increasing radius and weighting the flux therein accordingly"""
	total_noise_flux = 0
	for annulus in range(radius, NUM_COLS/2):
		weight = weighting_function(annulus)
		annulus_sum = np.sum(thick_ring(frame, avg_max_pixel, annulus, annulus+1))
		total_noise_flux += weight*annulus_sum
	return total_noise_flux

def frame_aperture(frame, max_pixel):
	"""Returns two best aperture sizes. 
	Given the way the signal to noise and penalty functions work, these two apertures should be consecutive - they are the range in which we can find the ideal aperture"""
	aperture_1 = 0
	aperture_2 = 0
	signal_to_noise_1 = 0	
	signal_to_noise_1 = 0
	for aperture_radius in range(1, NUM_ROWS/2):
		noise = calculate_noise_flux(frame, max_pixel, aperture_radius)
		signal = calculate_frame_flux(frame, max_pixel, aperture_radius)
		signal_to_noise = signal/noise - penalty_function(aperture_radius)
		#print signal/noise
		#print signal_to_noise
		if signal_to_noise > signal_to_noise_1:
			signal_to_noise_2 = signal_to_noise_1
			aperture_2 = aperture_1
			signal_to_noise_1 = signal_to_noise
			aperture_1 = aperture_radius
		elif signal_to_noise > signal_to_noise_2:
			signal_to_noise_2 = signal_to_noise
			aperture_2 = aperture_radius
	return (aperture_1, aperture_2)

def avg_aperture(frame_list):
	apertures = []
	for frame in frame_list:
		max_pixel = brightest_region(frame)
		apertures.append(frame_aperture(frame, max_pixel))
	return (sum(ap[0] for ap in apertures)/len(apertures), sum(ap[1] for ap in apertures)/len(apertures))


###################################################################################
# Methods of producing the best aperture range #
###################################################################################

def best_ap(fits_file):
	"""A very encapsulated aperture finding function
	Input: filename. Output: aperture. Easy, one-step."""
	hdulist = fits.open(fits_file)
	frame_list = hdulist[0].data
	rv = avg_aperture(frame_list)
	return rv

def best_ap_avg_frame(fits_file):
	"""Alternate version of best_ap(); 
	instead calculates a single average frame from the 64."""
	hdulist = fits.open(fits_file)
	frame_list = hdulist[0].data
	avg_frame = mk_avg_frame(frame_list)
	max_pixel = brightest_region(avg_frame)
	rv = frame_aperture(avg_frame, max_pixel)
	return rv

def best_ap_random(fits_file):
	"""Calculates best aperture using 4 random representative frames
	chosen from distinct 16-frame sectors.
	For more detailed description of implementation, 
	see docstring for Method D in parallelizer.py.
	"""
	hdulist = fits.open(fits_file)
	frame_list = hdulist[0].data
	index_choices = []
	for d in range(0,3):
		#16-frame sector
		z = range(16*d, 16*(d+1))
		choice = random.choice(z)
		index_choices.append(choice)
	new_frame_list = [frame_list[i] for i in index_choices]
	return avg_aperture(new_frame_list)




def main():
	fits_file = "prototype_data/SPITZER_I1_41629440_0000_0000_1_bcd.fits"
	#hdulist = fits.open(fits_file)
	#frame_list = hdulist[0].data
	#frame_one = frame_list[0]
	#br = brightest_region(frame_one)
	#test_disk = thick_ring(frame_one, br, 4, 6)
	#print_frame(test_disk, spacing = 3)
	
	aperture = best_ap_avg_frame(fits_file)
	print "best aperture for frame list:", aperture
	#numpy2image("images/frame_one.png", frame_one)
	#numpy2image("images/frame_one_mask_rad_2.png", test_disk)
	
	#subprocess.call(["ds9", "-zoom","8", fits_file])


if __name__ == "__main__":
	main()
