#Testing fits I/O with astropy
#Documentation can be found here:
#http://docs.astropy.org/en/stable/io/fits/index.html
#4/27/2014
from astropy.io import fits

hdu_list = fits.open("Exoplanet123_Prototype/SPITZER_I1_41629440_0000_0000_1_bcd.fits")
print hdu_list.info()
print "Observer: " + hdu_list[0].header['OBSRVR']
hdu_list.close()