import numpy as np
import matplotlib.pyplot as plt

### PARAMETERS #################################################################################################

F0 = 3631   #reference flux in Junsky

ron = 9   #number of electrons
dark = 0.25   #electrons/sec at T=0°C

size = 10.2 * 60   #size of the target in arcsec
px = 0.44   #arcsec per pixel
radius = size / px / 2   #galaxy radius in pixels
area = np.pi * radius**2   #area of the galaxy in pixels


### broad filters ###
g = [10.92, 11.3]   #magnetude
lambda_g = 4770   #central wavelenght of the filter in A

r = [10.4, 10.66]   #magnetude
lambda_r = 6580   #central wavelenght of the filter in A

i = [10.1, 10.29]   #magnetude
lambda_i = 8060   #central wavelenght of the filter in A


### narrow filters ###
lambda_Hb = 4860   #wavelenght in A

lambda_OIII = 5000   #wavelenght in A

lambda_Ha = 6560   #wavelenght in A

lambda_SII = 6720   #wavelenght in A


#%%
### FUNCTIONS ##################################################################################################

def m_to_flux(mag, lambda_c):
    f0 = F0 * 10**(-23)   #reference flux in erg/s/cm2/Hz
    f = f0 * 10**(-mag/2.5)   #flux in erg/s/cm2/Hz
    return f * 3*10**(18) / lambda_c**2

def flux_to_rate(flux, band):
    if band=='g':
        return flux * 10**(17.61)
    if band=='r':
        return flux * 10**(18.12)
    if band=='i':
        return flux * 10**(18.08)
    if band=='Hbeta':
        return flux * 10**(16.06)
    if band=='OIII':
        return flux * 10**(16.00)
    if band=='Halpha':
        return flux * 10**(16.58)
    if band=='SII':
        return flux * 10**(16.52)
    
def e_sky(band):
    if band=='g':
        return px**2 * 10**(-17.47)
    if band=='r':
        return px**2 * 10**(-17.98)
    if band=='i':
        return px**2 * 10**(-18.88)
    if band=='Hbeta':
        return px**2 * 10**(-17.62)
    if band=='OIII':
        return px**2 * 10**(-17.54)
    if band=='Halpha':
        return px**2 * 10**(-18.21)
    if band=='SII':
        return px**2 * 10**(-18.26)

def SN(flux, t, band):
    ne_sky = e_sky(band)
    ne_dark = dark * t
    ne_s = flux_to_rate(flux, band) * t
    return ne_s / np.sqrt(ne_s + area * (ne_sky + ne_dark + ron**2))


#%%
### FLUXES #####################################################################################################

### broad filters ###
F_g = m_to_flux(np.mean(g), lambda_g)
F_r = m_to_flux(np.mean(r), lambda_r)
F_i = m_to_flux(np.mean(i), lambda_i)


### narrow filters ###
#the flux is multiplied by some factor to take into account both continuum and line
F_Hb = m_to_flux(np.mean(g), lambda_Hb) * 1.5
F_OIII = m_to_flux(np.mean(g), lambda_OIII) * 1.5
F_Ha = m_to_flux(np.mean(r), lambda_Ha) * 2
F_SII = m_to_flux(np.mean(r), lambda_SII) * 2


#%%
### S/N RATIO ##################################################################################################

t = np.linspace(0.2, 30*60, 100)   #time of exposure

SN_g = SN(F_g, t, 'g')   #signal to noise computed for each band
SN_r = SN(F_r, t, 'r')
SN_i = SN(F_i, t, 'i')
SN_Hb = SN(F_Hb, t, 'Hbeta')
SN_OIII = SN(F_OIII, t, 'OIII')
SN_Ha = SN(F_Ha, t, 'Halpha')
SN_SII = SN(F_SII, t, 'SII')

plt.figure()
plt.plot(t, SN_g, color='g', label='g')
plt.plot(t, SN_r, color='r', label='r')
plt.plot(t, SN_i, color='b', label='i')
plt.title('Signal to Noise ratio - broad filters')
plt.xlabel('exposure time (s)')
plt.ylabel('S / N')
plt.legend(title='Band:')
plt.show()

plt.figure()
plt.plot(t, SN_Hb, color='tomato', label='Hb')
plt.plot(t, SN_OIII, color='dodgerblue', label='OIII')
plt.plot(t, SN_Ha, color='orange', label='Ha')
plt.plot(t, SN_SII, color='mediumseagreen', label='SII')
plt.title('Signal to Noise ratio - narrow filters')
plt.xlabel('exposure time (s)')
plt.ylabel('S / N')
plt.legend(title='Band:')
plt.show()