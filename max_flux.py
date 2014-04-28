from astropy.io import fits
import numpy as np

hdu_list = fits.open("Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits")

def find_max_flux(hdu_list):

  total_max_flux = 0
  total_max_x_pix = 0
  total_max_y_pix = 0

  for frame in range(len(hdu_list[0].data)):    

    max_flux = 0
    max_x_pix = 0
    max_y_pix = 0


    for x in range(len(hdu_list[0].data[frame])):
      for y in range(len(hdu_list[0].data[frame][x])):

        if hdu_list[0].data[frame][x][y] > max_flux:
          max_flux = hdu_list[0].data[frame][x][y]
          max_x_pix = x
          max_y_pix = y

    total_max_flux = total_max_flux + max_flux
    total_max_x_pix = total_max_x_pix + max_x_pix
    total_max_y_pix = total_max_y_pix + max_y_pix

  avg_max_flux = total_max_flux/float(len(hdu_list[0].data))
  avg_max_x_pix = total_max_x_pix/float(len(hdu_list[0].data))
  avg_max_y_pix = total_max_y_pix/float(len(hdu_list[0].data))

  print avg_max_flux
  print avg_max_x_pix
  print avg_max_x_pix

find_max_flux(hdu_list)
