import aperture_calculator as AC
import numpy as np
from mpi4py import MPI
import os

#The number of images in a single fits file
NUM_FITS_IMAGES = 64

comm = MPI.COMM_WORLD
size = comm.Get_size()

#This function takes a directory and returns a list of files to analyze
#Condition: the filename ends with "bcd.fits"
def get_fits(directory):
	all_files = [os.path.join(directory, file) for file in os.listdir(directory)]
	fits_images = []
	for filename in all_files:
		ext = filename.split('.')[-1]
		if ext == "fits" and filename[-8:-5] == "bcd":
			fits_images.append(filename)
	return fits_images

#This function analyzes a single file in parallel
#64 FITS images are spread across our available ranks
#This implementation would be distinct from the parallelize_file_set function
#This breaks a single file's analysis into steps, but parallelize_file_set function spreads a set of files
#Function still in progress
def parallelize_single_file(filename):
	pass

#This function analyzes the appropriate files in an input directory
#Individual fits files are *not* spread among the ranks, unlike the parallelize_single_file
#This is an alternate approach to parallelization
#This function is working, it's just waiting on the AC.best_ap function.
def parallelize_file_set(directory):
	print '('+str(comm.Get_rank())+')',"Entered parallelizer."
	rank = comm.Get_rank()
	
	#Responsibilities for process 0, the file distributor
	if rank == 0:
		input_files = get_fits(directory)
		#How evenly can the files be distributed among the ranks?
		remainder = len(input_files) % size
		files_per_rank = len(input_files) / size
		for process in range(size):
			#Distribute files to each process
			if process != (size - 1):
				upper_bound = (process + 1) * files_per_rank
			else:
				upper_bound = (process + 1) * files_per_rank + remainder
			comm.send(input_files[process * files_per_rank : upper_bound], dest=process, tag=process)
			print "Thread (0) sent files ", process * files_per_rank, " through ", upper_bound, " to process ", process 
	
	#Responsibilities for all processes
	input_files = comm.recv(source=0, tag=rank)
	print "Thread ", rank, " got files ", [(str(file) + " ") for file in input_files]
	aperture_dictionary = {}
	for file in input_files:
		print '('+str(comm.Get_rank())+')',"Analyzing " + file
		aperture = AC.best_ap(file)
		aperture_dictionary[file] = aperture
	f = open('process_'+str(rank)+'_output', 'w')
	for key in aperture_dictionary:
		f.write(key + '\t' + str(aperture_dictionary[key]) + '\n')
	f.close()
	return


def main():
	parallelize_file_set("prototype_data")


if __name__ == "__main__":
	main()
