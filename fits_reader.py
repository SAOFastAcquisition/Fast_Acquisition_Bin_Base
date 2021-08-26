import astropy
from astropy.io import fits
import matplotlib.pyplot as plt


hdulist = fits.open('NICMOSn4hk12010_mos.fits.txt')
print(hdulist)
hdulist.info()
hdu = hdulist[2]
print(hdu.header)
a = hdu.header
b = hdu.data
print(hdu.data)
m = hdu.data.shape
print(m)
plt.imshow(hdu.data[:,:], origin='lower')
plt.show()