#Aperture production with Astropy, Bresenham circle algorithm 
#Documentation can be found here:
#http://docs.astropy.org/en/stable/io/fits/index.html
#Started 4/27/2014
from astropy.io import fits
import scipy
import ds9
import subprocess
import numpy as np
import math
import matplotlib.pyplot as plt
import PIL
import Image
#import sklearn.cluster as sk
#from sklearn import metrics
#from sklearn.cluster import KMeans
#from sklearn.datasets import load_digits
#from sklearn.decomposition import PCA
#from sklearn.preprocessing import scale

NUM_ROWS = 32
NUM_COLS = 32

#Spitzer params:
FOCAL_LENGTH = 10.2		#meters
DIAMETER = 0.85				#meters
F_NUM = FOCAL_LENGTH/DIAMETER
WAVELENGTH = 3.6			#microns


#Find 2x2 area with greatest total luminosity
#Return index of brightest pixel in the 2x2 area
#Finding the brightest pixel this way eliminates the concern for a false 
#positive brightest pixel, which could occur with a spike in detector noise
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

#Returns the index of the brightest region across a list of frames
def avg_brightest_region(frame_list):
	brightest_indices = []
	#avg_brightest_ind = (0, 0)
	for frame in frame_list:
		brightest_indices.append(brightest_region(frame))
	brightest_ind_total = [sum(x) for x in zip(*brightest_indices)]
	avg_brightest_ind = (brightest_ind_total[0]/len(frame_list), brightest_ind_total[1]/len(frame_list))
	return avg_brightest_ind


#Returns data with items not lying on circle (of given radius) masked by zero
#Uses the Bresenham Circle Algorithm to draw nice looking circles in taxicab geometry
def mk_circle(frame, center_tuple, radius):	
	(x0, y0) = center_tuple
	new_frame = np.array([[0.0]*NUM_COLS]*NUM_ROWS)
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


#Gives inversion of the mk_disk function (zeros in the center)
def punch_hole(frame, center_tuple, radius):
	new_frame = []
	frame_list_form = np.ndarray.tolist(frame)
	#Circle_frame has radius one pixel larger than what we want
	circle_frame = np.ndarray.tolist(mk_circle(frame, center_tuple, radius + 1))
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
		n_row = frame_list_form[row_num][0:indA] + circle_frame[row_num][indA:indB] + frame_list_form[row_num][indB:]
		new_frame.append(n_row)
	return np.array(new_frame)


#Merges frames, replacing zero elements in one frame with nonzero elements from the other
#This function assumes no overlapping elements
def frame_merge(frame1, frame2):
	new_frame = frame1
	for i in range(NUM_ROWS):
		for j in range(NUM_COLS):
			if frame2[i][j] != 0:
				new_frame[i][j] = frame2[i][j]
	return new_frame
			

#Returns frame completely masked (by zeros) except for a ring with specified inner and outer radii
def thick_ring(frame, center_tuple, inner_rad, outer_rad):
	new_frame = []
	#flf stands for frame_list_form
	flf = np.ndarray.tolist(frame)
	#Combine the punch_hole and mk_disk functions
	outer_circle = np.ndarray.tolist(mk_circle(frame, center_tuple, outer_rad))
	inner_circle = np.ndarray.tolist(mk_circle(frame, center_tuple, inner_rad))
	cf = frame_merge(inner_circle, outer_circle) #cf stands for circle_frame
	for row_num in range(len(cf)):
		indA, indB = 16, 16
		for i in range(31):
			if cf[row_num][i] != 0:
				indA = i
				break
		for j in range(31,16,-1):
			if cf[row_num][j] != 0:
				indB = j
				break
		indC, indD = indA, indB
		for i in range(indA+2, 31):
			if cf[row_num][i] != 0:
				indC = i
				break
		for j in range(indB-2,0,-1):
			if cf[row_num][j] != 0:
				indD = j
				break
		if indD < indC:
			n_row = cf[row_num][0:indA] + flf[row_num][indA:indB] + cf[row_num][indB:]
			new_frame.append(n_row)
			continue
		n_row = cf[row_num][0:indA] + flf[row_num][indA:indC] + cf[row_num][indC:indD] + flf[row_num][indD:indB] + cf[row_num][indB:]
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


def numpy2image(new_filename, frame):
	#Normalize frame luminosities to list in range 0 to 255
	frame1 = frame/np.max(np.abs(frame))
	frame1 *= 255
	im = Image.fromarray(np.uint8(frame1))
	im.save(new_filename)
	new_img = Image.open(new_filename)
	new_img = new_img.resize((288, 288), Image.BILINEAR)
	new_img.save(new_filename)


#Plots light curve with respect to a given pixel as the origin
#Saves image in images folder
#Returns a list of the average fluxes plotted
def light_curve(frame, index_tuple, max_radius = 5):
	x_points = []
	y_points = []
	x = 0
	for rad in range(max_radius + 1):
		z = get_annulus_frames(frame, index_tuple, rad) # frame_neighbors_n(frame, index_tuple, level = rad)
		#print z
		avg = (sum(z))/(float(len(z)))
		x_points.append(x)
		y_points.append(avg)
		x += 1
	print "flux by anulus", y_points
	plt.plot(x_points, y_points, 'r--')
	plt.title('Light Curve from pixel ' + str(index_tuple))
	plt.xlabel('Radius (pixels)')
	plt.ylabel('Average luminosity')
	plt.savefig("images/light_curve" + str(index_tuple[0]) + "-" + str(index_tuple[0]) + ".png")
	#plt.show()
	return y_points


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


def calculate_avg_flux(frame_list, avg_max_pixel, radius, is_signal):
	avg_frame_flux = 0
	for frame in frame_list:
		if is_signal:
			disk = mk_disk(frame, avg_max_pixel, radius)
		else:
			disk = punch_hole(frame, avg_max_pixel, radius)
		frame_flux = 0
		for row in disk:
			for flux in row:
				if flux > 0:
					frame_flux += flux
		avg_frame_flux += frame_flux
	avg_frame_flux = avg_frame_flux/len(frame_list)
	return avg_frame_flux


#Returns the total flux in a disk created from the frame argument
def calculate_frame_flux(frame, avg_max_pixel, radius):
	disk = mk_disk(frame, avg_max_pixel, radius)
	return sum(sum(disk))


#Returns a weight by radius, based on a Gaussian function approximating a point spread function
#Lower radii (closer in to the star) get a higher weight than those farther away
def weighting_function(max_flux, radius):
	return math.exp(-(radius*radius)/(2*0.45*WAVELENGTH*F_NUM))


#Returns a penalty for using more radii (still refining best function to use)
def penalty_function(radius):
	return math.exp(radius - 5.0)


#Returns the total noise in a frame (anything in a frame that is not the signal)
#Total flux is calcuated by taking into account annuli of increasing radius and weighting the flux therein accordingly
def calculate_noise_flux(frame, avg_max_pixel, radius):
	total_noise_flux = 0
	for annulus in range(radius, len(frame[0])/2):
		weight = weighting_function(frame[avg_max_pixel[0]][avg_max_pixel[1]], annulus)
		annulus_sum = sum(sum(thick_ring(frame, avg_max_pixel, annulus, annulus+1)))
		total_noise_flux += weight*annulus_sum
	return total_noise_flux


def calculate_avg_noise(frame_list, avg_max_pixel, radius):
	for frame in frame_list:
		noise_flux = 0
	pass

def test_aperture(frame, center_tuple, radius):
	pass



def main():

	fits_file = "Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits"
	hdulist = fits.open(fits_file)
	frame_list = hdulist[0].data
	frame_one = frame_list[0]
	
	#noise_file = "Exoplanet123_Prototype/SPITZER_I1_41608704_0000_1_C8734032_sdark.fits"
	#noiselist = fits.open(noise_file)
	#noise_list = noiselist[0].data
	#noise_one = noise_list[0]

	#Prints a version of the star with center at brightest region, and radius r
	avg_max_pixel = avg_brightest_region(frame_list)
	print "average maximum pixel location", avg_max_pixel
	#test_disk = mk_disk(frame_one, max_pixel, radius = 4)
	#test_noise_disk = mk_disk(noise_one, max_pixel, radius = 4)

	#avg_noise = calculate_avg_flux(noise_list, avg_max_pixel, radius = 1)
	#noise_frame = punch_hole(frame_one, avg_max_pixel, 4)
	#print_frame(noise_frame, spacing = 2)

	radius = 1

	noise_frame = calculate_noise_flux(frame_one, avg_max_pixel, radius)
	print "weighted noise in frame one", noise_frame
	flux_frame = calculate_frame_flux(frame_one, avg_max_pixel, radius)
	print "flux from frame on", flux_frame
	#flux = calculate_avg_flux(frame_list, avg_max_pixel, 2, True)
	#noise = calculate_noise_flux(frame_list, avg_max_pixel, 2, False)
	#print "average noise", avg_noise
	#print "average flux", avg_flux

	sig_noise_r4 = flux_frame/noise_frame - penalty_function(radius)
	print "signal-to-noise ratio for a radius of " + str(radius) + " in frame one:", sig_noise_r4

	#test_disk2 = punch_hole(frame_one, max_pixel, radius = 4)
	
	#thick = thick_ring(frame_one, avg_max_pixel, 4,7)
	#print_frame(thick, spacing = 2)
	
	#print_frame(test_disk, spacing = 2)
	#print_frame(test_noise_disk, spacing = 2)
	#Testing punch_hole
	#print "-"*50
	#print_frame(test_disk2, spacing = 2)
	
	
	#test_plot = mk_plot
	#significance = is_significant(frame_one, (15, 15), 1)
	
	#plot_flux_vs_pos(frame_one, max_pixel, 3)
	#y_points = light_curve(frame_one, max_pixel)
	#slopes = plot_slopes(y_points)	
	#ratios = flux_ratios(y_points)
	#print "slopes", slopes
	#print "ratios", ratios
	#print "sinal-to-noise ratios by annulus", signal_to_noise(frame_one, inverted_frame, y_points)
	#numpy2image("images/frame_one.png", frame_one)
	#numpy2image("images/frame_one_mask_rad_2.png", test_disk)
	
	#annulus_frames = get_annulus_frames(frame_one, max_pixel, 2)
	
	#subprocess.call(["ds9", "-zoom","8", fits_file])


if __name__ == "__main__":
	main()
