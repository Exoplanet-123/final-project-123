#Optimal Aperture Calculation
####Parallelized Analysis of the Spitzer Telescope Dataset
Team: Hannah Diamond-Lowe, Zakir Gowani

The Spitzer telescope exoplanet observation dataset provides 32x32 images of stars over a time period comparable to the orbit of the (potential) orbiting exoplanet. The time-dependence of brightness variations describes qualities of the exoplanet. However, light noise is a problem. Current methods rely on manual aperture adjustment to throw out unwanted data. Our project is to automate this process and extend the accuracy by accounting for a finer degree of time variation. 

Our immediate purpose is finding intelligent, time-dependent aperture sizes for star observations using the Spitzer telescope which reduce noise but conserve valuable star image information. The ultimate goal is to to produce the optimal aperture size masks for a considerable number of stars (each star is ~5 GB worth of images). A stretch goal is to compare models existing in the literature to a model output based on our aperture calculations. 



Packages to be used will include:
* [Hadoopy](http://www.hadoopy.com/en/latest/), a Python wrapper for Hadoop written in Cython.
* [Astropy](https://astropy.readthedocs.org/en/stable/overview.html), a popular astronomy library for Python which deals with .fits files, the standard data format for our Spitzer star observation data.
* [Python Imaging Library](http://www.pythonware.com/products/pil/)


#####Timeline:
* By May 1: Prototype will demonstrate aperture calculation algorithm on a single frame for a single star
* By May 18: Expand algorithm's operation to consider time variation over frames
* By May 25: Move algorithm into Hadoop framework, linked by Hadoopy
* By June 3: Final implementation of project which compares models, accompanied by visualizations and a writeup
