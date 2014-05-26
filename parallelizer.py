#################################################################
# Aperture calculation with MPI, Numpy, and Astropy  
# This script parallelizes the action of the aperture_calculator script
# Project started 4/27/2014
# Team: Hannah-Diamond Lowe, Zakir Gowani
################################################################
import aperture_calculator as AC
import numpy as np
from mpi4py import MPI
import os
import sys
from astropy.io import fits

NUM_FITS_IMAGES = 64	# The number of images in a single fits file

comm = MPI.COMM_WORLD
size = comm.Get_size()

def get_fits(directory):
	"""This function takes a directory and returns a list of files to analyze
		Condition: the filename ends with 'bcd.fits'
	"""
	all_files = [os.path.join(directory, file) for file in os.listdir(directory)]
	fits_images = []
	for filename in all_files:
		ext = filename.split('.')[-1]
		if ext == "fits" and filename[-8:-5] == "bcd":
			fits_images.append(filename)
	return fits_images

def parallelize_single_file(filename, output_file):
	"""This function analyzes a single file in parallel.
	64 FITS images are spread across our available ranks
	This implementation would be distinct from the parallelize_file_set function
	This breaks a single file's analysis into steps, but parallelize_file_set function spreads a set of files
	Beware: this function will fail on a system that has less than 2 threads
	"""
	print '(' + str(comm.Get_rank()) + ')', "Entered parallelizer."
	rank = comm.Get_rank()
	hdulist = fits.open(filename)
	frame_list = hdulist[0].data
	# Responsibilities for process 0, the file distributor
	if rank == 0:
		if size > NUM_FITS_IMAGES:
			print "You have more threads than images in a FITS file."
			print "Suggestion: revise your functions to take advantage of all threads."
		remainder = (NUM_FITS_IMAGES) % size
		imgs_per_rank = NUM_FITS_IMAGES / size
		print "(0) Assigned imgs ", 0, " through ", imgs_per_rank, " to process ", rank
		for process in range(1, size):
			#Distribute image numbers to each process
			if process != (size - 1):
				upper_bound = (process+1) * imgs_per_rank
			else:
				upper_bound = (process+1) * imgs_per_rank + remainder
			comm.send(range(process*imgs_per_rank, upper_bound), dest=process, tag=process)
			print "(0) Assigned imgs ", process*imgs_per_rank, " through ", upper_bound, " to process ", process
	# Responsibilities for all processes
	if rank != 0:
		my_range = comm.recv(source=0, tag=rank)
	if rank == 0:
		my_range = range(0, imgs_per_rank)
	sum = 0
	for i in my_range:
		max_pixel = AC.brightest_region(frame_list[i])
		sum += AC.frame_aperture(frame_list[i], max_pixel)
	print '(' + str(comm.Get_rank()) + ')', "Intermediate sum: ", sum
	if rank != 0:
		comm.send(sum, dest=0, tag=rank)
	if rank == 0:
		total = sum
		for process in range(1, size):
			z = comm.recv(source=process, tag=process)
			total += z
		avg_aperture = float(total) / (float(NUM_FITS_IMAGES))
		f = open(output_file, 'a')
		f.write(filename + '\t' + str(avg_aperture))
		print "Final aperture val:", filename + '\t' + str(avg_aperture)
		f.close()
	return
		
def prll_dir_mthd_A(directory, output_file):
	"""Uses the parallel_single_file function to process each file in a directory.
	This approach is very different from the method B function, prll_dir_mthd_B. 
	In this approach, each thread contributes to processing all files (cooperative work).
	In method B, each thread gets its own files (individual work).
	Timing of this function's operation on 5 files yields: 1m13.984s (real), which is faster than mthd B
	"""
	input_files = get_fits(directory)
	for file in input_files:
		parallelize_single_file(file, output_file)
	return

def prll_dir_mthd_B(directory, output_file):
	"""This function analyzes all the appropriate files in an input directory
	Individual fits files are *not* spread among the ranks, unlike the parallelize_single_file function.
	Instead, a set of files is distributed to each thread to process on its own.
	Timing of this function on a set of 5 files yields: 1m55.745s(real)
	"""
	print '(' + str(comm.Get_rank()) + ')', "Entered parallelizer."
	rank = comm.Get_rank()
	
	# Responsibilities for process 0, the file distributor
	if rank == 0:
		input_files = get_fits(directory)
		# How evenly can the files be distributed among the ranks?
		remainder = len(input_files) % size
		files_per_rank = len(input_files) / size
		for process in range(size):
			# Distribute files to each process
			if process != (size - 1):
				upper_bound = (process+1) * files_per_rank
			else:
				upper_bound = (process+1) * files_per_rank + remainder
			comm.send(input_files[process*files_per_rank : upper_bound], dest=process, tag=process)
			print "(0) Sent files ", process*files_per_rank, " through ", upper_bound, " to process ", process 
	# Responsibilities for all processes
	input_files = comm.recv(source=0, tag=rank)
	aperture_dict = {}
	for file in input_files:
		#aperture = 2
		aperture = AC.best_ap(file)
		aperture_dict[file] = aperture
	# The amount of space each process needs in the file is determined by the size of aperture_dict
	# Each rank will tell one other rank how much space it needs
	tot_bytes = 0
	for key in aperture_dict.keys():
		tot_bytes += sys.getsizeof(key + '\t' + str(aperture_dict[key]) + '\n')
	comm.send(tot_bytes, dest=0)
	offset_dict = {}
	if rank == 0:
		for process in range(size):
			byte_offset = comm.recv(source=process)
			offset_dict[process] = byte_offset
		for process in range(1, size):
			comm.send(offset_dict, dest=process, tag=process)
	elif rank != 0:
		offset_dict = comm.recv(source=0, tag=rank)
	my_offset = 0
	for key in offset_dict.keys():
		if key < rank:
			my_offset += offset_dict[key]
	print '(' + str(rank) + ')' + " Byte offset is", my_offset
	for process in range(size):
		if rank == process:
			f = open(output_file, 'a')
			f.seek(my_offset)
			for key in aperture_dict.keys():
				f.write(key + '\t' + str(aperture_dict[key]) + '\n')
				print key + '\t' + str(aperture_dict[key])
			f.close()
		else:
			continue
	return

def main():
	prll_dir_mthd_A("prototype_data", "combined_outputB.txt")
	#prll_dir_mthd_B("prototype_data", "combined_outputB.txt")
	#parallelize_single_file("prototype_data/abcd.fits", "abcd_output.txt")

if __name__ == "__main__":
	main()
