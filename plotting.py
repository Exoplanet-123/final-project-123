import math
import matplotlib.pyplot as plt
import aperture_calculator as ac
from astropy.io import fits

FOCAL_LENGTH = 10.2				# In meters
DIAMETER = 0.85					# In meters
F_NUM = FOCAL_LENGTH/DIAMETER
WAVELENGTH = 3.6				# In microns

NUM_ROWS = 32					# The dimensions of a FITS image
NUM_COLS = 32

fits_file = "prototype_data/SPITZER_I1_41629440_0000_0000_1_bcd.fits"
hdulist = fits.open(fits_file)
frame_list = hdulist[0].data
avg_frame = ac.mk_avg_frame(frame_list)
max_pixel = ac.brightest_region(avg_frame)
rv = ac.frame_aperture(avg_frame, max_pixel)

s_n = [0.154547861818, 0.403526244264, 0.772498520473, 1.33493102866, 2.11287033914, 3.44958935702, 5.5838167695, 9.28302033016, 15.7997153381, 27.0342961138, 49.8872234858, 95.2443001898, 199.353286902, 454.247471923, 1366.3357413]
#for i in range(len(s_n)):
#	s_n[i] = s_n[i]/s_n[-1]

an_sum = [18560.9447522, 21859.9067966, 22793.8054371, 23751.4771974, 22396.5774943, 25454.3094973, 25498.6133756, 26717.4479049, 27129.5026967, 26564.0617414, 29080.7097892, 28853.3001293, 30392.3182606, 29923.7592573, 31420.1164373, 33063.7496585]

an_weight = [18560.9447522, 21304.8352603, 20565.3608689, 18843.4098908, 14840.8798195, 13381.6042331, 10101.6436409, 7576.32888009, 5230.62754907, 3307.61010714, 2221.22445384, 1284.13386742, 748.627449789, 387.493955085, 203.171683728, 101.40838172]
#for i in range(len(an_weight))


x = (range(NUM_ROWS/2))
y = []
p = []
flux = []

for radius in range(NUM_ROWS/2):
	y.append(math.exp( -(radius * radius) / (2 * 0.45 * WAVELENGTH * F_NUM) ))
	p.append(math.exp(radius - 3.6))
	flux.append(ac.calculate_noise_flux(avg_frame, (15, 15), radius))

for i in range(len(an_weight)):
	print an_weight[i] - p[i]


#for i in range(len(p)):
#	p[i] = p[i]/p[-1]

plt.figure(figsize=(9, 6))
plt.axes([.11, .08, .85, .85]) #[left, bottom, width, height]
#plt.plot(x, y, "k-", linewidth=2, label="Gaussian approximation")
#plt.xlabel("aperture radius")
#plt.ylabel("normalized Gaussian")
#plt.title("weighting function v. aperture radius")
#plt.legend()
#plt.ylim(-0.01, 1.01)
#plt.xlim(0.0, 15.0)

plt.plot(x, an_weight, "b-", linewidth=2, label="annulus flux")
plt.plot(x, p, "b--", linewidth=2, label="penalty")
plt.xlabel("aperture radius")
plt.ylabel("flux")
plt.title("flux v. aperture radius")
plt.ylim(-500, 60000)
plt.xlim(0.0, 15.0)
plt.legend(loc=2)

plt.show()
