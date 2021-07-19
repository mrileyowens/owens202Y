import pandas as pd
import matplotlib.pyplot as plt
import os
from glob import glob
import numpy as np
from pathlib import Path
from astropy.io import fits
from matplotlib.offsetbox import AnchoredText
import scipy.interpolate
import scipy.ndimage
from scipy.optimize import curve_fit

def plotSpectrographs(mageFile,slitMaskFile,arcMaskFile,dataCubeFile):

    # ------------------------------------------------------------------------
    # Adjusting MagE data into Numpy arrays from .txt format and fixing errors
    # ------------------------------------------------------------------------

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

    # -----------------------------------------------------------------------
    # Adjusting VLT data into Numpy array from .FITS format and fixing errors
    # -----------------------------------------------------------------------

    # Opening VLT .FITS datacube and retrieving data array containing
    # measured flux and pixel location at each wavelength
    dataCubeFlux=fits.open(dataCubeFile)[1].data
    dataCubeError=fits.open(dataCubeFile)[2].data

    # Reassigning NaN values to 0 in the VLT data
    dataCubeFlux[np.where(np.isnan(dataCubeFlux))]=0.0
    dataCubeError[np.where(np.isnan(dataCubeError))]=0.0

    # ----------------------------------------------------------------------
    # Converting arc mask from .FITS to array and adjusting to binary values
    # ----------------------------------------------------------------------

    # Opening .FITS mask of the arc and retrieving data array,
    # which specifies where the arc is and isn't on the sky
    arcMask=fits.open(arcMaskFile)[0].data

    # To simplify, we reassign positive values of the mask
    # to be 1 and zero or negative values to be 0
    arcMask[arcMask > 0]=1.0
    arcMask[arcMask <= 0]=0

    # -----------------------------------------------------------------------
    # Converting slit mask from .FITS to array and adjusting to binary values
    # -----------------------------------------------------------------------

    # Opening .FITS mask of a slit, which specifies where in the sky
    # a telescope slit did and did not collect light
    slitMask=fits.getdata(slitMaskFile)

    # To simplify, we reassign positive values of the mask
    # to be 1 and zero or negative values to be 0
    slitMask[slitMask > 0]=1.0
    slitMask[slitMask <= 0]=0

    # -------------------------------------------------------------------
    # Calculating the intersection of the arc and slit masks and creating
    # datacube mask
    # -------------------------------------------------------------------

    # Adding the two masks together into a new array
    maskSum=np.add(arcMask,slitMask)

    # The intersection of the two masks is where maskSum
    # is 2, since we adjusted the slit and arc masks to have
    # binary values, so we make a new array which reassigns the
    # array locations which don't meet this condition to zero
    maskInt=np.where(maskSum==2,maskSum,0)

    # To simplify maskInt, its values are reassigned to be binary
    maskInt[maskInt > 0]=1.0

    # maskInt is only 2D, but we have to apply it to a 3D datacube,
    # so we broadcast it to the datacube's shape
    maskInt3D=np.broadcast_to(maskInt,dataCubeFlux.shape)

    # Multiplying the 3D mask array with the datacube array reduces
    # non-intersecting values to zero but preserves the datacube's
    # original values which intersect with the mask, since the
    # mask's values are binary
    maskedFlux=maskInt3D*dataCubeFlux
    maskedError=maskInt3D*dataCubeError

    #fits.writeto('C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Test Directory/slitH5pxlSrc.fits',maskInt,header=fits.open(dataCubeFile)[1].header)
    #fits.writeto('C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Test Directory/slitH5pxlSrc.fits',maskedFlux,header=fits.open(dataCubeFile)[1].header)

    # ---------------------------------------------------------------
    # Creating wavelength, flux and error arrays from masked datacube
    # ---------------------------------------------------------------

    # Extracting the wavelength range and interval from the .FITS
    # datacube header
    startWavelength=fits.open(dataCubeFile)[1].header['CRVAL3']
    deltaWavelength=fits.open(dataCubeFile)[1].header['CD3_3']
    layers=fits.open(dataCubeFile)[1].header['NAXIS3']

    # Creating a wavelength array for the VLT data from the
    # wavelength range and interval
    vltWaveArray=np.arange(layers)*deltaWavelength+startWavelength

    # Summing the flux and error of the masked datacube at each sampled wavelength
    vltFluxArray=np.sum(maskedFlux,axis=(1,2))
    vltErrorArray=np.sqrt(np.sum(np.square(maskedError),axis=(1,2)))

    # Reassigning negative flux values to 0
    vltFluxArray[vltFluxArray <= 0]=0

    # The datacube's flux is measured in terms of fnu, but the MagE data
    # uses flambda, so we adjust the units here
    vltFluxArray=(3.336*(10**(-19)))*((vltWaveArray)**2)*vltFluxArray*(10**(-20))
    vltErrorArray=(3.336*(10**(-19)))*((vltWaveArray)**2)*vltErrorArray*(10**(-20))

    # -------------------------------------------
    # Trimming MagE and VLT data before rebinning
    # -------------------------------------------

    # The two wavelength arrays from MagE and VLT do not have the same
    # range, so we trim them to only their intersecting wavelength range
    mageFluxArray=mageFluxArray[mageWaveArray >= np.min(vltWaveArray)]
    mageErrorArray=mageErrorArray[mageWaveArray >= np.min(vltWaveArray)]
    redCorrArray=redCorrArray[mageWaveArray >= np.min(vltWaveArray)]
    mageWaveArray=mageWaveArray[mageWaveArray >= np.min(vltWaveArray)]

    mageFluxArray=mageFluxArray[mageWaveArray <= 7000]
    mageErrorArray=mageErrorArray[mageWaveArray <= 7000]
    redCorrArray=redCorrArray[mageWaveArray <= 7000]
    mageWaveArray=mageWaveArray[mageWaveArray <= 7000]
    mageErrorArray=mageErrorArray[mageWaveArray <= 7000]
    vltFluxArray=vltFluxArray[vltWaveArray <= 7000]
    vltErrorArray=vltErrorArray[vltWaveArray <= 7000]
    vltWaveArray=vltWaveArray[vltWaveArray <= 7000]

    # ---------------------------------------------------------------------------
    # Rebinning the MagE flux array to the dimensions of the VLT wavelength array
    # ---------------------------------------------------------------------------

    # If the source array is not already composed of floats, this is adjusted here
    if not mageFluxArray.dtype in [np.float64, np.float32]:
        mageFluxArray=np.cast[float](mageFluxArray)

    # This quantity determines the resampling factor
    #minusOne=np.cast[int](False)

    # ofs adjusts the interpolation points to the center or front of the bins
    #ofs=np.cast[int](True)*0.5

    # Creates an array storing the shape of the MagE flux array before
    # interpolation
    #oldShapeArray=np.array(mageFluxArray.shape)

    # Specifies the number of dimensions of the MagE flux array and raises an
    # error if this is not consistent with the dimensions of the VLT
    # wavelength array
    nDims=len(mageFluxArray.shape)
    if len(vltWaveArray.shape) != nDims:
        print("Dimensions error: must rebin to the same number of dimensions.")
        return None

    # Stores the desired dimensions as an array of integers
    #newDimArray=np.asarray(vltWaveArray.shape,dtype=int)
    #dimList=[]

    method=input('Enter rebinning method (neighbor, linear, or spline): ')

    if method == 'neighbor':
        for i in range(nDims):
            base=np.indices(newDimArray,dtype=float)[i]
            dimList.append(((oldShapeArray[i]-minusOne)/(newDimArray[i]-minusOne))*(base+ofs)-ofs)
        cd=np.array(dimList).round().astype(int)
        mageFluxArrayRebin=mageFluxArray[tuple(list(cd))]

    elif method == 'linear':

        # This returns an interpolated function which approximates the MagE flux and error
        intFuncFlux=scipy.interpolate.interp1d(mageWaveArray,mageFluxArray,kind='linear',axis=0)
        intFuncError=scipy.interpolate.interp1d(mageWaveArray,mageErrorArray,kind='linear',axis=0)
        intFuncRed=scipy.interpolate.interp1d(mageWaveArray,redCorrArray,kind='linear',axis=0)

        # The interpolated function strictly only covers the domain between the minimum and
        # maximum wavelengths of the MagE wavelength array, so we have to drop any points of
        # the VLT arrays which do not fall between those boundaries; these would require
        # extrapolation
        while np.min(vltWaveArray)<=np.min(mageWaveArray):
            vltWaveArray=np.delete(vltWaveArray,0)
            vltFluxArray=np.delete(vltFluxArray,0)
            vltErrorArray=np.delete(vltErrorArray,0)

        while np.max(vltWaveArray)>=np.max(mageWaveArray):
            vltWaveArray=np.delete(vltWaveArray,-1)
            vltFluxArray=np.delete(vltFluxArray,-1)
            vltErrorArray=np.delete(vltErrorArray,-1)

        # Calling the interpolation function at the VLT wavelength samples, which uses the
        # interpolated function to approximate the MagE flux at those points in a new array,
        # with the
        mageFluxArrayRebin=intFuncFlux(vltWaveArray)
        mageErrorArrayRebin=intFuncError(vltWaveArray)
        redCorrArrayRebin=intFuncRed(vltWaveArray)

    elif method == 'spline':
        oSlices=[slice(0,j) for j in oldShapeArray]
        oldCoords=np.ogrid[oSlices]
        nSlices=[slice(0,j) for j in list(newDimArray)]
        newCoords=np.mgrid[nSlices]
        newCoordsDims=list(range(np.ndim(newCoords)))
        newCoordsDims.append(newCoordsDims.pop(0))
        newCoordsTr=newCoords.transpose(newCoordsDims)
        #newCoordsTr+=ofs
        #deltas=(np.asarray(oldShapeArray)-minusOne)/(newDimArray-minusOne)
        #newCoordsTr*=deltas
        #newCoordsTr-=ofs
        mageFluxArrayRebin=scipy.ndimage.map_coordinates(mageFluxArray,newCoords)

    else:
        print("Error: unrecognized interpolation type.")

    redInput=input('Would you like to write the reddening corrections? (Y/N): ')
    if input==Y:
        writeArray=np.stack(vltWaveArray,redCorrArrayRebin)
        fits.writeto(Path(mageFile).stem+'MUSEredCorr.fits',writeArray,header=None)
    else:
        return None

    # Applying the de-reddening correction to the VLT spectrum
    vltFluxArray=redCorrArrayRebin*vltFluxArray

    # The two spectra have much different flux intensities, so to better
    # compare them, we normalize the flux and error
    normMageFluxArray=mageFluxArrayRebin/np.median(mageFluxArray[(5700 <= mageWaveArray) & (mageWaveArray <= 6200)])
    normVltFluxArray=vltFluxArray/np.median(vltFluxArray[(5700 <= vltWaveArray) & (vltWaveArray <= 6200)])
    normMageErrorArray=(normMageFluxArray/mageFluxArrayRebin)*mageErrorArrayRebin
    normVltErrorArray=(normVltFluxArray/vltFluxArray)*vltErrorArray


    mageErrorWeightsArray=1/np.square(normMageErrorArray)
    vltErrorWeightsArray=1/np.square(normVltErrorArray)
    residualFluxArray=((normMageFluxArray-normVltFluxArray)/normMageFluxArray)
    residualErrorArray=np.sqrt(np.square(normMageErrorArray)+np.square(normVltErrorArray))
    residualSigmaArray=residualFluxArray/residualErrorArray

    ratioArray=(normVltFluxArray/normMageFluxArray)

    def funcLin(x,a,b):
        return a*x+b

    linFit=scipy.optimize.curve_fit(funcLin,vltWaveArray,ratioArray)

    linParams=linFit[0]

    correctedMageFluxArray=np.multiply(normMageFluxArray,funcLin(vltWaveArray,linParams[0],linParams[1]))

    avgFluxArray=((mageErrorWeightsArray*normMageFluxArray)+(vltErrorWeightsArray*normVltFluxArray))/(mageErrorWeightsArray+vltErrorWeightsArray)
    avgErrorArray=np.sqrt(((mageErrorWeightsArray*normMageErrorArray)**2+(vltErrorWeightsArray*normVltErrorArray))/((mageErrorWeightsArray**2)+(vltErrorWeightsArray**2)))

    #Titling the plot
    titleInput=input('Name the plot: ')

    fig, axs=plt.subplots(4,sharex=True)
    fig.suptitle(titleInput)
    #axs[3].plot(vltWaveArray,ratioArray-1,drawstyle='steps-mid',color='red',linewidth=0.5,label='Flux Ratio')
    #axs[3].plot(vltWaveArray,funcLin(vltWaveArray,linParams[0],(linParams[1]-1)),color='black',linestyle='--',linewidth=0.5,label='Linear Fit, '+"{:.2e}".format(linParams[0])
#+'x'+"{:.2f}".format((linParams[1]-1)))
    #axs[3].plot(vltWaveArray,funcPow(vltWaveArray,[0],parameters[1]),color='black',linestyle='--',linewidth=0.5,label='Linear Fit, '+"{:.2e}".format(parameters[0])+'x+'+"{:.2f}".format(parameters[1]))
    #axs[3].set_ylabel('Flux Ratio with Linear Fit \n [VLT/MagE]')
    #axs[3].legend(loc='upper right')

    axs[1].plot(vltWaveArray,normMageFluxArray,drawstyle='steps-mid',color='blue',linewidth=0.5,label='MagE')
    axs[1].plot(vltWaveArray,normVltFluxArray,drawstyle='steps-mid',color='orange',linewidth=0.5,label='VLT')
    axs[1].set_ylabel('Normalized MagE/VLT Flux')
    axs[1].set_ylim(ymin=0)
    axs[1].legend(loc='upper right')

    # Plot the residual between the MagE and MUSE flux in terms of the uncertainty
    axs[2].plot(vltWaveArray, residualSigmaArray,drawstyle='steps-mid',color='red',linewidth=0.5)
    axs[2].axhline(y=0, color='black', linestyle='--',linewidth=0.5)
    axs[2].set_ylabel('Residual ($\sigma$) \n [(MagE - VLT) / MagE]')

    # Plot the averaged normalized flux of MagE and MUSE, and the associated error
    axs[0].plot(vltWaveArray,avgFluxArray,drawstyle='steps-mid',color='black',linewidth=0.5,label='Average Flux')
    axs[0].plot(vltWaveArray,avgErrorArray,drawstyle='steps-mid',color='red',linewidth=0.5,label='Uncertainty')
    axs[0].set_ylabel('Variance Weighted Flux\nand Noise')
    axs[0].legend(loc='upper right')
    axs[0].set_ylim(ymin=0,ymax=3)

    #Plot the corrected MUSE flux against the MagE flux
    axs[3].plot(vltWaveArray,mageFluxArrayRebin,drawstyle='steps-mid',color='blue',label='MagE',linewidth=0.5)
    axs[3].plot(vltWaveArray,vltFluxArray,drawstyle='steps-mid',color='orange',label='VLT',linewidth=0.5)
    axs[3].legend(loc='upper right')
    axs[3].set_ylabel('MagE vs. Corrected MUSE Flux')

    plt.xlim(np.min(vltWaveArray),7000)
    plt.xlabel('Wavelength (Å)')

    plt.show()

mageFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Data/Spectra/Lyman Spectra/psz-arcslit-f-comb1_MWdr.txt'
slitMaskFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/Slit Masks (.fits)/VLT Masks/slitH5VLT.fits'
arcMaskFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/psz1550_1909_v2_mask.fits'
dataCubeFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/PSZ311_CUBE_wcscorr.fits'

plotSpectrographs(mageFile,slitMaskFile,arcMaskFile,dataCubeFile)
