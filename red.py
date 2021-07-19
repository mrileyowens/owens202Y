import pandas as pd
import os
from glob import glob
import numpy as np
from pathlib import Path
from astropy.io import fits
import scipy.interpolate
import matplotlib.pyplot as plt

mageFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Data/Spectra/Lyman Spectra/planckarc_pos1-comb1_MWdr.txt'
dataCubeFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/PSZ311_CUBE_wcscorr.fits'
#testRedTableFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/testRedTable.fits'


def func(mageFile,dataCubeFile):
# Converting MagE .txt file into Pandas data frame 
    # and removing data collection and labeling information outlined in .txt header
    dataFrame=pd.read_csv(mageFile,delimiter='\t',header=0,skiprows=13)

    # Dropping NaN values and erroneously large values in the data frame
    dataFrame=dataFrame.apply(pd.to_numeric,errors='coerce').dropna()
    dataFrame=dataFrame[dataFrame.iloc[:,1]<1e-20]

    # Dropping empty data frame rows to reset the index so that
    # the first row with an entry is at the zero index
    dataFrame=dataFrame.reset_index(drop=True)

    # Extracting data series for the wavelength, flux and error from the MagE data frame
    mageWaveSeries=dataFrame.iloc[:,0]
    mageFluxSeries=dataFrame.iloc[:,1]
    mageErrorSeries=dataFrame.iloc[:,2]
    redCorrSeries=dataFrame.iloc[:,-1]

    # Converting MagE wavelength/flux/error data series into Numpy arrays
    mageFluxArray=mageFluxSeries.to_numpy()
    mageWaveArray=mageWaveSeries.to_numpy()
    mageErrorArray=mageErrorSeries.to_numpy()
    redCorrArray=redCorrSeries.to_numpy()

    startWave=fits.open(dataCubeFile)[1].header['CRVAL3']
    deltaWave=fits.open(dataCubeFile)[1].header['CD3_3']
    layers=fits.open(dataCubeFile)[1].header['NAXIS3']

    vltWaveArray=np.arange(layers)*deltaWave+startWave

    vltWaveArray=vltWaveArray[vltWaveArray <= np.max(mageWaveArray)]
#    redCorrArray=redCorrArray[mageWaveArray >= np.min(vltWaveArray)]
#    mageWaveArray=mageWaveArray[mageWaveArray >= np.min(vltWaveArray)]
#    redCorrArray=redCorrArray[mageWaveArray <= np.max(vltWaveArray)]
#    mageWaveArray=mageWaveArray[mageWaveArray <= np.max(vltWaveArray)]

    intFuncRed=scipy.interpolate.interp1d(mageWaveArray, redCorrArray,kind='linear',axis=0)

    redCorrArrayRebin=intFuncRed(vltWaveArray)

#    plt.plot(vltWaveArray,redCorrArrayRebin,color='orange',label='VLT Reddening',drawstyle='steps-mid',linewidth=0.5)
#    plt.plot(mageWaveArray, redCorrArray,color='blue',label='MagE Reddening',drawstyle='steps-mid',linewidth=0.5)
#    plt.legend(loc='upper right')
#    plt.xlim(xmin=np.min(vltWaveArray),xmax=np.max(vltWaveArray))
#    plt.ylim(1.2,1.4)
#    plt.show()

    c1=fits.Column(name='Wavelength (A)',array=vltWaveArray,format='E')
    c2=fits.Column(name='Reddening Correction',array=redCorrArrayRebin,format='E')
    cols=fits.ColDefs([c1,c2])
    t=fits.BinTableHDU.from_columns(cols)
    t.writeto('slit1RedCorrMUSE.fits')
    
func(mageFile,dataCubeFile)
