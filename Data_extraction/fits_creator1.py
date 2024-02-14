from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from Supporting_func.stokes_coefficients import path_to_data


fits_image_filename = '20240212_122749_sun+0_out.fits'
hdul = fits.open(fits_image_filename)
print(hdul.info())
hdr = hdul[0].header
print(f'hdr[0] = {hdr[5]}')
a = hdul[0].data
b = hdul[1].data
c = hdr.cards
a1 = a[:, 0, :]
print(c)
hdul.close()
pass
