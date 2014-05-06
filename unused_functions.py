def rms_xy(frame, max_pixel, radius):
	i = max_pixel[0]
	j = max_pixel[1]
	mean = frame[i][j]

	rms_x = 0
	rms_y = 0

	row_pixels = []
	col_pixels = []

	for pixel in frame[i][j-radius:j+radius+1]:
		row_pixels.append(pixel*pixel)
	rms_x = np.mean(row_pixels)

	for pixel in frame.T[j-radius:j+radius+1]:
		col_pixels.append(pixel*pixel)
	rms_y = np.mean(col_pixels)

	return (rms_x, rms_y)

def weighting_function(max_flux, radius)
	weight = max_flux

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

def get_annulus_frames(frame, center_tuple, annulus_num):
	i = center_tuple[0]
	j = center_tuple[1]
	annulus_frames = []
	x_shift = range(i - annulus_num, i + annulus_num+1)
	y_shift = range(j - annulus_num+1, j + annulus_num)
	for x in x_shift:
		if frame[x, j - annulus_num] != 0:
			annulus_frames.append(frame[x, j - annulus_num])
		if frame[x, j + annulus_num] != 0:
			annulus_frames.append(frame[x, j + annulus_num])
	for y in y_shift:
		if frame [i - annulus_num, y] != 0:
			annulus_frames.append(frame[i - annulus_num, y])
		if frame [i + annulus_num, y] != 0:
			annulus_frames.append(frame[i + annulus_num, y])
	return annulus_frames

#Takes numpy array, index of pixel in question, and integer number of standard deviations
#If pixel in question and its immediate neighbors are n SD outside of the mean, the pixel is considered 'significant'
def is_significant(frame, index_tuple, deviations = 1, mean = -1, SD = -1):
	if mean == -1 or SD == -1:
		frame_mean = np.mean(frame)
		#print frame_mean
		frame_SD = np.std(frame)
		#print frame_SD
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

#Returns a list of the slopes seen in the plot
def plot_slopes(y_points):
	slopes = []
	for y in range(len(y_points)):
		if y == len(y_points)-1:
			break
		else:
			slopes.append(y_points[y+1] - y_points[y])
	return slopes


#returns a list of ratios of different annuli fluxes
def flux_ratios(fluxes):
	ratios = []
	for flux in range(len(fluxes)):
		if flux == len(fluxes)-1:
			break
		else:
			ratios.append(fluxes[flux]/fluxes[flux+1])
	return ratios


def signal_to_noise(frame, inverted_frame, annulus_flux):
	return
	avg_noise = 0
	length = 0
	for noise in inverted_frame:
		#Some flux values for the background are negative. 
		#This is not a realistic measurement and is more likely a dead pixel. 
		#We do not want to include these in our average background noise
		if noise > 0:
			avg_noise += noise
			length += 1
	avg_noise += 1

	ratios = []
	for flux in annulus_flux:
		ratios.append(flux/avg_noise)

	return ratios


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


def calculate_avg_noise(frame_list, avg_max_pixel, radius):
	for frame in frame_list:
		noise_flux = 0
	pass
