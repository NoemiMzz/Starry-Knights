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
   # imNii = hdul['NII6548_FLUX'].data

    #errors
    imHa_err = hdul['HA6562_FLUX_ERR'].data
    imSii_err= hdul['SII6716_FLUX_ERR'].data
    imOii_err= hdul['OIII4958_FLUX_ERR'].data
    imHb_err= hdul['HB4861_FLUX_ERR'].data
    #imNii_err = hdul['NII6548_FLUX_ERR'].data
im = [imHa, imSii, imOii, imHb] #, imNii
im_err = [imHa_err, imSii_err, imOii_err, imHb_err ] #imNii_err


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
#temp[4][mask] = np.nan
imHa_sn = temp[0]
imSii_sn = temp[1]
imOii_sn = temp[2]
imHb_sn = temp[3]
#imNii_sn = temp[4]
sn_array = [imHa_sn, imSii_sn, imOii_sn, imHb_sn] #imNii_sn

plt.imshow(imHa_sn, cmap='viridis', clim = [0, 1000])
plt.colorbar() #848 , 1131
plt.show()

aperture_array = []

# Parametri di input
radius = 25
x = np.arange(0, 1401, 50)
y = np.arange(0, 1401, 50)

# Array per memorizzare le aperture
aperture_array = [CircularAperture((xi, yi), r=radius) for xi in x for yi in y]

# Array per memorizzare gli indici delle aperture valide
valid_aperture_idx = []

# Filtra le aperture valide
for i in range(len(sn_array)):
    for j, aperture in enumerate(aperture_array):
        # Coordinate dei pixel coperte dall'apertura
        mask = aperture.to_mask(method='center') #maschera che considera solo i dati nell'apertura
        cutout = mask.to_image(sn_array[i].shape) #adatta la maschera alla dimensione della immagine, in particolare 1 se pixel è coperto un valore intermedio se parzialmente coperto e 0 se no

        # Verifica se l'apertura contiene solo valori validi
        if not np.isnan(sn_array[i][cutout.astype(bool)]).all():
            valid_aperture_idx.append(j)

# Plot dell'immagine con le aperture valide
for i in range(len(sn_array)):
    plt.imshow(sn_array[i], cmap='viridis', clim=[0, 1000])
    for j in valid_aperture_idx:
        aperture_array[j].plot(color='red', lw=1)
    plt.colorbar()
    plt.show()
aperture_array = [aperture_array[i] for i in valid_aperture_idx]



"""
imHa_sn = np.where(np.isnan(imHa_sn), 0, imHa_sn)
imOii_sn = np.where(np.isnan(imOii_sn), 0, imOii_sn)
imHb_sn = np.where(np.isnan(imHb_sn), 0, imHb_sn)
imSii_sn = np.where(np.isnan(imSii_sn), 0, imSii_sn)"""
flux_ha = []
flux_sii = []
flux_oii = []
flux_hb = []
#flux_nii = []
f_ha = []
f_sii = []
f_oii = []
f_hb = []
#f_nii = []
imHa_sn = np.where(np.isnan(imHa_sn), 0, imHa_sn)
imHb_sn = np.where(np.isnan(imHb_sn), 0, imHb_sn)
imOii_sn = np.where(np.isnan(imOii_sn), 0, imOii_sn)
imSii_sn = np.where(np.isnan(imSii_sn), 0, imSii_sn)
for i in range(len(aperture_array)):

  flux_ha.append(aperture_photometry(imHa_sn, aperture_array[i])['aperture_sum'])
  flux_sii.append(aperture_photometry(imSii_sn, aperture_array[i])['aperture_sum'])
  flux_oii.append(aperture_photometry(imOii_sn, aperture_array[i])['aperture_sum'])
  flux_hb.append(aperture_photometry(imHb_sn, aperture_array[i])['aperture_sum'])
  #flux_nii.append(aperture_photometry(imNii_sn, aperture_array[i])['aperture_sum'])
  #print(flux_ha[i])
for i in range(len(flux_ha)):

  f_ha.append(flux_ha[i][0])
  f_sii.append(flux_sii[i][0])
  f_oii.append(flux_oii[i][0])
  f_hb.append(flux_hb[i][0])
  #f_nii.append(flux_nii[i][0])
ratio_x = []
ratio_y = []
print(len(f_ha))
#ratio_xN= []
for i in range(len(f_ha)):
  ratio_x.append(f_sii[i]/f_ha[i])
  ratio_y.append(f_oii[i]/f_hb[i])
  #ratio_xN.append(f_nii[i]/f_ha[i])
print(len(ratio_x))
print(ratio_x)
"""
print(len(f_ha))

print(f_ha)
print(f_sii)
print(f_oii)
print(f_hb)
"""

plt.scatter(np.log10(ratio_x), np.log10(ratio_y), s= 10)
plt.xlabel('log(Sii/Ha)')
plt.ylabel('log(Oii/Hb)')
plt.title('BPT Diagram')

plt.show()

"""plt.scatter(np.log10(ratio_xN), np.log10(ratio_y), s= 10)
plt.xlabel('log(Nii/Ha)')
plt.ylabel('log(Oii/Hb)')
plt.title('BPT Diagram')
#plt.ylim(-2,0)
plt.show()"""