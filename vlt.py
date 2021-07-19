import pandas as pd
import regions
from regions import read_ds9,write_ds9,ds9_objects_to_string,write_fits_region, fits_region_objects_to_table
import matplotlib.pyplot as plt
import os
from glob import glob
import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy import wcs
from matplotlib.colors import LogNorm
from pathlib import Path
from astropy.visualization import astropy_mpl_style

hubbleImageFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/V2.0_PSZ1G311.65-18.48_F606W_0.03g0.8_cr2.5_0.7_drc_sci.fits'
slitRegionFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/Slit Regions (.reg)/mageSlit1.reg'

regionFilePath=os.path.join('.','Region Images/Slit Regions (.reg)/*.reg')
regionFileList=glob(regionFilePath)


def makeMask(regionFile,fitsFile):

    # 
    fitsData,fitsHeader=fits.getdata(fitsFile,header=True)
    w=wcs.WCS(fitsHeader)

    region=read_ds9(regionFile)
    pixelRegion=region[0].to_pixel(w)
    maskRegion=pixelRegion.to_mask(mode='center').to_image((fitsHeader['NAXIS1'],fitsHeader['NAXIS2']))

    fits.writeto(Path(regionFile).stem+'.fits',maskRegion.data,header=fitsHeader)

for regionFile in regionFileList:
    makeMask(regionFile,hubbleImageFile)

def extractWCS(fitsFile):
    fitsData=fits.open(fitsFile,memap=True)
    w=wcs.WCS(fitsData[0].header)
    return w
        
def viewData(fitsFile):
    imageData=fits.getdata(fitsFile)
    plt.imshow(imageData,cmap='gray',norm=LogNorm(),vmax=10,vmin=1,origin='lower')
    plt.colorbar()
    plt.show()

def makeMask(regionFile,fitsFile):
    region=read_ds9(regionFile)
    pixelRegion=region[0].to_pixel(extractWCS(fitsFile))
    maskRegion=pixelRegion.to_mask(mode='center').to_image((319,315))
    return maskRegion

def plotMask(regionFile,fitsFile):
    plt.imshow(makeMask(regionFile,fitsFile).to_image((319,315)),cmap='gray',origin='lower')
    plt.show()

def writeMask(regionFile,fitsFile):
    hdu=fits.PrimaryHDU(makeMask(regionFile,fitsFile).data)
    hdu.writeto(Path(regionFile).stem+'VLT.fits')

def writeMaskInt(regionFile1,regionFile2,fitsFile):
    maskRegion1=makeMask(regionFile1,fitsFile)
    maskRegion2=makeMask(regionFile2,fitsFile)
    maskSum=np.add(maskRegion1.data,maskRegion2.data)
    maskInt=np.where(maskSum == 2,maskSum,0)
    hdu=fits.PrimaryHDU(maskInt)
    hdu.writeto(Path(regionFile1).stem+Path(regionFile2).stem+'Int-1500A.fits')

def writeMaskDiff(regionFile1,regionFile2,fitsFile):
    maskRegion1=makeMask(regionFile1,fitsFile)
    maskRegion2=makeMask(regionFile2,fitsFile)
    maskSum=np.add(maskRegion1.data,maskRegion2.data)
    maskDiff=np.where(maskSum != 2,maskSum,0)
    hdu=fits.PrimaryHDU(maskDiff)
    hdu.writeto(Path(regionFile1).stem+Path(regionFile2).stem+'Diff-1500A.fits')

regionFilePath=os.path.join('.','Region Images/Slit Regions (.reg)/*.reg')
regionFileList=glob(regionFilePath)

slitFilePath=os.path.join('.','Region Images/Slit Masks (.fits)/VLT Masks/*.fits')
slitFileList=glob(slitFilePath)
        
hubbleImageFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/Region Images/V2.0_PSZ1G311.65-18.48_F606W_0.03g0.8_cr2.5_0.7_drc_sci.fits'
#slitFile='C:/Users/15136/OneDrive/Documents/Research/Region Images/Slit Masks (.fits)/VLT Masks/slitH9VLT.fits'
slitNameList=['MagE Slit 1','Slit F0','Slit F1','Slit F2','Slit F3','Slit F4', 'Slit F5','Slit F6','Slit F7','Slit F8','Slit F9','Slit H1a','Slit H1','SLit H2','Slit H3','Slit H4','Slit H5','Slit H6','Slit H9']

#for (slitName,slitFile) in zip(slitNameList,slitFileList):
#    writeSpectra(dataCubeFile,arcMaskFile,slitFile,slitName)
