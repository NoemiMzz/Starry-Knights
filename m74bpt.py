# -*- coding: utf-8 -*-
"""m74bpt.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1prJoAPWlrkJnKDRnqyRpPUwW_lZeYzLY
"""

!pip install photutils

import os
from google.colab import drive
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
path = '/content/drive/MyDrive/'

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import photutils.aperture as pha
from tqdm import tqdm
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry

filename = path+'m74datacube.fits'
with fits.open(filename) as hdul:
    # Esamina le informazioni delle estensioni nel file
    hdul.info()

    # Estrai il data cube dalla prima estensione (solitamente 0 è l'header primario)
    data_cube = hdul[0].data
    header = hdul[0].header

    imHa = hdul['HA6562_FLUX'].data
    imSii= hdul['SII6716_FLUX'].data
    imOii= hdul['OIII4958_FLUX'].data
    imHb= hdul['HB4861_FLUX'].data

    #errors
    imHa_err = hdul['HA6562_FLUX_ERR'].data
    imSii_err= hdul['SII6716_FLUX_ERR'].data
    imOii_err= hdul['OIII4958_FLUX_ERR'].data
    imHb_err= hdul['HB4861_FLUX_ERR'].data
im = [imHa, imSii, imOii, imHb]
im_err = [imHa_err, imSii_err, imOii_err, imHb_err]


# Esamina i dati
#print(f"Dimensioni del data cube: {data_cube.shape}")
print(f"Header: {header}")

plt.imshow(imHa, cmap='viridis', clim = [0, 1000])
plt.colorbar()
plt.show()

aperture = CircularAperture((700, 700), r = 700)
for i in range(4):
  plt.imshow(im[i], cmap='viridis', clim = [0, 1000])
  aperture.plot(color = 'red', lw = 2)
  plt.colorbar()
  plt.show()

#S/N

for i in range(4):
  plt.imshow(im[i]/im_err[i], cmap='viridis', clim = [0, 50])
  #plt.plot(aperture, color = 'red')
  plt.colorbar()
  plt.show()

SN = im[i]/im_err[i]
mask = SN < 5
temp = im
temp[0][mask] = np.nan
temp[1][mask] = np.nan
temp[2][mask] = np.nan
temp[3][mask] = np.nan
imHa_sn = temp[0]
imSii_sn = temp[1]
imOii_sn = temp[2]
imHb_sn = temp[3]

plt.imshow(imHa_sn, cmap='viridis', clim = [0, 1000])
plt.colorbar()
plt.show()

flux_ha = aperture_photometry(imHa_sn, aperture)
flux_sii = aperture_photometry(imSii_sn, aperture)
flux_oii = aperture_photometry(imOii_sn, aperture)
flux_hb = aperture_photometry(imHb_sn, aperture)

f_ha = flux_ha['aperture_sum']
f_sii = flux_sii['aperture_sum']
f_oii = flux_oii['aperture_sum']
f_hb = flux_hb['aperture_sum']

