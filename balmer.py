import numpy as np
import os
from glob import glob
from astropy.io import fits
from pathlib import Path
import matplotlib.pyplot as plt

def makeSpectrograph(fluxFile,noiseFile,redshift):

    #Extracting flux and noise data from the FITS file
    fluxArray=fits.open(fluxFile)[0].data
    noiseArray=fits.open(noiseFile)[0].data
    
    #Creating a wavelength array
    wavelengthInitial=fits.open(fluxFile)[0].header['CRVAL1']
    wavelengthIncrement=fits.open(fluxFile)[0].header['CDELT1']
    waveArray=10**(np.arange(len(fluxArray))*wavelengthIncrement+wavelengthInitial)
    waveArray=waveArray/(1+redshift)

    #Trimming the arrays
    fluxArray=fluxArray[waveArray >= 6557.79]
    noiseArray=noiseArray[waveArray >= 6557.79]
    waveArray=waveArray[waveArray >= 6557.79]
    fluxArray=fluxArray[waveArray <= 6567.79]
    noiseArray=noiseArray[waveArray <= 6567.79]
    waveArray=waveArray[waveArray <= 6567.79]

    #Converting flux units
    fluxArray=(3.336*(10**(-19)))*((waveArray)**2)*fluxArray*(10**(-28))
    noiseArray=(3.336*(10**(-19)))*((waveArray)**2)*noiseArray*(10**(-28))

    #fluxArray=np.where((wavelengthArray>6559.3) & (wavelengthArray<6569.3),fluxArray,0)
    #maskedNoiseArray=np.where((wavelengthArray>6559.3) & (wavelengthArray<6569.3),noiseArray,0)
   
    plt.plot(waveArray,fluxArray,drawstyle='steps-mid',color='black',label='Flux')
    plt.plot(waveArray,noiseArray,drawstyle='steps-mid',color='red',label='Flux St. Dev.')

    titleInput=input('Name the plot: ')
    plt.title(titleInput)

    plt.xlabel('Rest Wavelength (Å)')
    plt.ylabel('Flux (10$^{-28}$ erg/cm$^2$/s/Hz)')

    plt.xlim(6557.79,6567.79)

    plt.legend(loc='upper right')

    plt.show()

fluxFileList=glob(os.path.join('.','Data/Spectra/FIRE Spectra/Flux/*.fits'))
noiseFileList=glob(os.path.join('.','Data/Spectra/FIRE Spectra/Error/*.fits'))
#sixth and eighth redshifts are approximate; redshift .txt doc unclear
redshiftList=[0,2.37015,2.36993,2.37067,2.37045,2.37,2.37054,2.37,2.37042]

for (fluxFile,noiseFile,redshift) in zip(fluxFileList,noiseFileList,redshiftList):
    makeSpectrograph(fluxFile,noiseFile,redshift)
