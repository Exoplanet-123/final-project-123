from astropy.io import fits
import numpy as np

hdu_list = fits.open("Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits")

pixel_rows = hdu_list[0].header['pxscal1']  # arcseconds per pixel (1D)
pixel_columns = hdu_list[0].header['pxscal2'] # arcseconds per pixel (1D)

pixel = pixel_rows*pixel_columns  # arcseconds^2 per pixel (2D)

data_MJy_per_pix = []

for i in hdu_list[0].data[0][0]:
  MJy_per_arc2 = i*2.35045E-11
  MJy_per_pix = MJy_per_arc2*pixel
  data_MJy_per_pix.append(MJy_per_pix)

data = np.array(data_MJy_per_pix)

