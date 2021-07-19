import pandas as pd
import matplotlib.pyplot as plt
import os
from glob import glob
import astropy.units as u
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from lmfit.models import SkewedGaussianModel
from scipy import stats

spectraFile='C://Users/15136/OneDrive - University of Cincinnati/Documents/Research/Data/Spectra/MagE Spectra/psz-arcslit-h9-comb1_MWdr.txt'
redshift=2.37054

def extractData(spectra):

    #Converting .txt file into a data frame
    dataFrame=pd.read_csv(spectra,delimiter='\t',header=0,skiprows=13)
	
    #Discarding missing data
    dataFrame=dataFrame.apply(pd.to_numeric,errors='coerce').dropna()
	
    #Discarding extreme outliers
    dataFrame=dataFrame[dataFrame.iloc[:,1]<1e-20]
	
    #Resetting the index
    dataFrame=dataFrame.reset_index(drop=True)

    # Extracting wavelength, flux, and noise
    wArr=dataFrame.iloc[:,0].to_numpy()
    fArr=dataFrame.iloc[:,1].to_numpy()
    nArr=dataFrame.iloc[:,2].to_numpy()

    return wArr,fArr,nArr

def gaussFit(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fitNormPeak(vArr,fArr,sampleMin,sampleMax,params,bounds):

    fArr=fArr[vArr > sampleMin]
    fArr=fArr[vArr < sampleMax]
    vArr=vArr[vArr > sampleMin]
    vArr=vArr[vArr < sampleMax]

    popt,pcov=curve_fit(gaussFit,vArr,fArr,p0=params,bounds=bounds)

    return popt,pcov

def fitSkewPeak(vArr,fArr,sampleMin,sampleMax,amplitude,center,sigma,gamma,color):

    fArr=fArr[vArr >= sampleMin]
    vArr=vArr[vArr >= sampleMin]
    fArr=fArr[vArr <= sampleMax]
    vArr=vArr[vArr <= sampleMax]

    #ae,loce,scalee=stats.skewnorm.fit(vArr,fArr)

    model=SkewedGaussianModel()

    guessParams=model.make_params(amplitude=amplitude,center=center,sigma=sigma,gamma=gamma)
    result=model.fit(fArr,x=vArr)
    print(result.fit_report())

    plt.plot(vArr,result.best_fit,color=color)

    return result

def centroid(spectra,z,cPeakParams):

    wArr,fArr,nArr=extractData(spectra)
    fArr=1e27*fArr
    vArr=(3e5)*(((wArr/(1+z))/1215.67)-1)

    brfArr=fArr-gaussFit(vArr,*cPeakParams)

    fArrB=brfArr[vArr < 15.0]
    vArrB=vArr[vArr < 15.0]
    fArrB=fArrB[vArrB > -400.0]
    vArrB=vArrB[vArrB > -400.0]

    fArrR=brfArr[vArr < 750.0]
    vArrR=vArr[vArr < 750.0]
    fArrR=fArrR[vArrR > 95.0]
    vArrR=vArrR[vArrR > 95.0]


    #fitB=fitSkewPeak(vArr,fArr,-346.0,12.0,1.0,-204.0,60.0,-3.0,'blue')
    #fitR=fitSkewPeak(vArr,fArr,97.0,500.0,1.39,200.0,50.0,4.0,'red')

    plt.plot(vArr,brfArr,drawstyle='steps-mid',color='black')
    plt.plot(vArrB,fArrB,drawstyle='steps-mid',color='blue',label=r'Blueshifted Peak')
    plt.plot(vArrR,fArrR,drawstyle='steps-mid',color='red',label=r'Redshifted Peak')
    plt.plot(vArr,gaussFit(vArr,*cPeakParams),drawstyle='steps-mid',linestyle='-.',color='black',label='Central Peak')

    plt.xlim(-1000,1000)
    plt.ylim(ymin=0)
    plt.xlabel('Velocity (km s$^{-1}$)')
    plt.ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    plt.legend(loc='upper right',fontsize=20)

    plt.show()

#centroid(spectraFile,redshift,[3.00,57.35,43.56])

def linFit(x,m,b):
    return m*x+b

def fitCont(wArr,fArr,leftBound1,leftBound2,rightBound1,rightBound2):
    
    fArr=fArr[((wArr >= leftBound1) & (wArr <= leftBound2)) | ((wArr >= rightBound1) & (wArr <= rightBound2))]
    wArr=wArr[((wArr >= leftBound1) & (wArr <= leftBound2)) | ((wArr >= rightBound1) & (wArr <= rightBound2))]

    popt,pcov=curve_fit(linFit,wArr,fArr)

    return popt,pcov

def int(vArr,fArr,bound1,bound2):

    fArr=fArr[(vArr >= bound1) & (vArr <= bound2)]
    vArr=vArr[(vArr >= bound1) & (vArr <= bound2)]

    areaTot=np.trapz(fArr,vArr)

    for v in vArr:
        fArr=fArr[(vArr >= v) & (vArr <= bound2)]
        vArr=vArr[(vArr >= v) & (vArr <= bound2)]

        percentile=np.trapz(fArr,vArr)/areaTot

        print(v, percentile)

wArr,fArr,nArr=extractData(spectraFile)
vArr=(3e5)*(((wArr/(1+redshift))/1393.76)-1)
fArr=1e27*fArr
popt,pcov=fitCont(wArr,fArr,4540.0,4630.0,4750.0,4925.0)
fit=linFit(wArr,popt[0],popt[1])

int(vArr,fit-fArr,-663.0,52.0)

plt.plot(wArr,fArr,drawstyle='steps-mid',color='black',linewidth=0.5)
plt.plot(wArr,fit,color='red',linestyle='dashed')
plt.xlim(wArr[0],wArr[-1])
plt.ylim(ymin=0.0)
plt.show()

plt.plot(vArr,fArr,drawstyle='steps-mid',color='black')
plt.plot(vArr,fit,drawstyle='steps-mid',color='red',linestyle='dashed')
plt.xlim(-1000.0,1000.0)
plt.ylim(0.0,0.175)
plt.show()

plt.plot(vArr,fit-fArr,drawstyle='steps-mid',color='black')
plt.xlim(-1000.0,1000.0)
plt.ylim(0.0,0.1)
plt.show()

def showVelocityPlot(spectraFile,redshift,lineLambda):

    #Converting .txt file into a data frame
    dataFrame=pd.read_csv(spectraFile,delimiter='\t',header=0,skiprows=13)
	
    #Discarding missing data
    dataFrame=dataFrame.apply(pd.to_numeric,errors='coerce').dropna()
	
    #Discarding extreme outliers
    dataFrame=dataFrame[dataFrame.iloc[:,1]<1e-20]
	
    #Resetting the index
    dataFrame=dataFrame.reset_index(drop=True)

    # Extracting wavelength, flux, and noise
    waveSeries=dataFrame.iloc[:,0]
    fluxSeries=1e27*dataFrame.iloc[:,1]
    noiseSeries=dataFrame.iloc[:,2]

    # Converting wavelength to velocity for several lines
    velocitySeries=(3e5)*(((waveSeries/(1+redshift))/lineLambda)-1)
    velocitySeriesCIV1548=(3e5)*(((waveSeries/(1+redshift))/1548.2)-1)
    velocitySeriesCIV1550=(3e5)*(((waveSeries/(1+redshift))/1550.77)-1)
    velocitySeriesSiIV1393=(3e5)*(((waveSeries/(1+redshift))/1393.76)-1)
    velocitySeriesSiIV1402=(3e5)*(((waveSeries/(1+redshift))/1402.77)-1)
    velocitySeriesSiII1259=(3e5)*(((waveSeries/(1+redshift))/1259.52)-1)
    velocitySeriesSiII1260=(3e5)*(((waveSeries/(1+redshift))/1260.42)-1)
    velocitySeriesOI1302=(3e5)*(((waveSeries/(1+redshift))/1302.0)-1)
    velocitySeriesOIII=(3e5)*(((waveSeries/(1+redshift))/5007.0)-1)
    velocitySeriesOII=(3e5)*(((waveSeries/(1+redshift))/3727.0)-1)

    fig,axs=plt.subplots(2,sharex=True)

    axs[0].plot(velocitySeriesOII,fluxSeries,drawstyle='steps-mid',color='black')
    axs[0].set_ylim(ymin=0)
    axs[1].plot(velocitySeriesOIII,fluxSeries,drawstyle='steps-mid',color='black')
    axs[1].set_ylim(ymin=0)

    plt.xlim(-1000,1000)

    plt.show()
    
    velocityDataFrame=pd.DataFrame({'Velocity':velocitySeries,'Flux':fluxSeries})

    velocityDataFrameShort=pd.DataFrame.reset_index(velocityDataFrame.loc[(velocitySeries>13.0) & (velocitySeries<98.0)])

    def gaussFit(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    popt,pcov=curve_fit(gaussFit,velocityDataFrameShort.iloc[:,1].to_numpy(),velocityDataFrameShort.iloc[:,2].to_numpy(),p0=[2.96,56.85,46.99])
    twoPeaks=velocityDataFrame.iloc[:,1].to_numpy()-gaussFit(velocityDataFrame.iloc[:,0],*popt)

    fig,axs=plt.subplots(2,sharex=True)
    plt.plot(velocityDataFrame.iloc[:,0],velocityDataFrame.iloc[:,1],drawstyle='steps-mid',color='black',label=r'Ly$\alpha$')
    plt.plot(velocityDataFrame.iloc[:,0],gaussFit(velocityDataFrame.iloc[:,0],*popt),drawstyle='steps-mid',color='red',linestyle='--',label='Gaussian Fit, '+'$a$='+"{:.2f}".format(popt[0])+', $x_0$='+"{:.2f}".format(popt[1])+', $\sigma$='+"{:.2f}".format(popt[2]))
    plt.xlim(-1000,1000)
    plt.ylim(ymin=0.0)
    plt.legend(loc='upper right')
    plt.xlabel('Velocity (km s$^{-1}$)')
    plt.ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    plt.title(r'Slit H9: Ly$\alpha$ and Gaussian Fit')
    
    plt.show()

    fig,axs=plt.subplots(3,sharex=True)
    fig.suptitle(r'Slit H1: Ly$\alpha$ and Absorption Profiles')
    axs[0].plot(velocityDataFrame.iloc[:,0],twoPeaks,drawstyle='steps-mid',color='black',label=r'Subt. Ly$\alpha$')
    axs[0].set_ylim(0,0.5)
    axs[0].set_ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    axs[0].legend(loc='upper right')

    axs[1].plot(velocitySeriesSiII1259,fluxSeries,drawstyle='steps-mid',color='orange',linestyle='dashed',label='SiII 1259')
    axs[1].plot(velocitySeriesSiII1260,fluxSeries,drawstyle='steps-mid',color='blue',linestyle='dashed',label='SiII 1260')
    axs[1].plot(velocitySeriesOI1302,fluxSeries,drawstyle='steps-mid',color='green',linestyle='dashed',label='OI 1302')
    axs[1].set_ylim(0.0,0.17)
    axs[1].set_ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    axs[1].legend(loc='upper right')

    axs[2].plot(velocitySeriesSiIV1393,fluxSeries,drawstyle='steps-mid',color='orange',linestyle='dashed',label='SiIV 1393')
    axs[2].plot(velocitySeriesSiIV1402,fluxSeries,drawstyle='steps-mid',color='blue',linestyle='dashed',label='SiIV 1402')
    axs[2].set_ylim(0.0,0.15)
    axs[2].set_ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    axs[2].legend(loc='upper right')

    plt.xlim(-1000,1000)
    plt.xlabel('Velocity (km s$^{-1}$)')

    plt.show()

    fig,axs=plt.subplots(2,sharex=True)
    fig.suptitle(r'Slit H1: Ly$\alpha$ and CIV Absorption Profile')

    axs[0].plot(velocityDataFrame.iloc[:,0],twoPeaks,drawstyle='steps-mid',color='black',label=r'Subt. Ly$\alpha$')
    axs[0].set_ylim(0,0.5)
    axs[0].set_ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
    axs[0].legend(loc='upper right')

    axs[1].plot(velocitySeriesCIV1548,fluxSeries,drawstyle='steps-mid',color='orange',linestyle='dashed',label='CIV 1548')
    axs[1].plot(velocitySeriesCIV1550,fluxSeries,drawstyle='steps-mid',color='blue',linestyle='dashed',label='CIV 1550')
    axs[1].set_ylim(0,0.15)
    axs[1].set_ylabel('Flux (10$^{-27}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')

    axs[1].legend(loc='upper right')

    plt.xlim(-1000,1000)
    plt.xlabel('Velocity (km s$^{-1}$)')

    plt.show()

def convertData(file):
	
	#Converting .txt file into a data frame
	dataFrameRaw=pd.read_csv(file,delimiter='\t',header=0,skiprows=13)
	
        #Discarding missing data
	dataFrameOutliers=dataFrameRaw.apply(pd.to_numeric,errors='coerce').dropna()
	
	#Discarding extreme outliers
	dataFrameSkewedIndex=dataFrameOutliers[dataFrameOutliers.iloc[:,1]<1e-20]
	
	#Resetting the index
	dataFrame=dataFrameSkewedIndex.reset_index(drop=True)
	
	return dataFrame
	
def showMovingFrameSpectralPlots(fileList):
	
    for file in fileList:
    
        fileName=fileNameList[fileList.index(file)]
    
        #Converting the raw data into a usable form
        dataFrame=convertData(file)
    
        #Splitting the data columns into series
        waveseries=dataframe.iloc[:,0]
        fluxseries=dataframe.iloc[:,1]*1e27
        noiseseries=dataframe.iloc[:,2]*1e27
        
        #finding the maximum and minimum flux
        fluxseriesmin=fluxseries.loc[fluxseries.idxmin()]
        fluxseriesmax=fluxseries.loc[fluxseries.idxmax()]

        #plotting flux vs. wavelength and associated noise
        plt.plot(waveseries,fluxseries,linewidth=0.2,color='black',label='observed flux')
        plt.plot(waveseries,noiseseries,linewidth=0.5,color='red',label='noise')
        plt.legend(loc='upper right')
	
        #titling plot
        plt.title('spectrum of '+filename+' (moving frame)')
	
        #labeling x and y axes
        plt.xlabel('observed wavelength ('+u'\u212B'.encode('utf-8')+')')
        plt.ylabel('flux (erg/s/cm^2/hz)*10^27')
	
        #scaling x and y axes
        plt.xlim([waveseries.loc[0],waveseries.loc[dataframe.index[-1]]])
        #plt.ylim([fluxseriesmin-1e-27,fluxseriesmax+1e-27])
	
        plt.show()
	
def showrestframespectralplots(filelist):
    for file in filelist:
    
        filename=filenamelist[filelist.index(file)]
        
        #converting the raw data into a usable form
        dataframe=convertdata(file)
        
        #splitting the data columns into series
        waveseries=dataframe.iloc[:,0]
        fluxseries=dataframe.iloc[:,1]*1e27
        noiseseries=dataframe.iloc[:,2]*1e27
        
        #finding the maximum and minimum flux
        fluxseriesmin=fluxseries.loc[fluxseries.idxmin()]
        fluxseriesmax=fluxseries.loc[fluxseries.idxmax()]
        
        #plotting flux vs. wavelength and associated noise
        plt.plot(waveseries,fluxseries,linewidth=0.2,color='black',label='flux')
        plt.plot(waveseries,noiseseries,linewidth=0.5,color='red',label='noise')
        plt.legend(loc='upper right')
        
        #titling plot
        plt.title('spectrum of '+filename+' (rest frame)')
        
        #labeling x and y axes
        plt.xlabel('observed wavelength ('+u'\u212B'.encode('utf-8')+')')
        plt.ylabel('flux (erg/s/cm^2/hz)*10^27')
        
        #scaling x and y axes
        plt.xlim([waveseries.loc[0],waveseries.loc[dataframe.index[-1]]])
        #plt.ylim([fluxseriesmin-1e-27,fluxseriesmax+1e-27])
        
        plt.show()
	
def findredshift(file,filelist):
	
	redshiftseries=pd.series([2.37067,2.37031,2.37067,2.37042,2.37045,2.36993,2.37015,2.37054])
	
	redshift=redshiftseries.loc[filelist.index(file)]
	
	return redshift

def makedataframe(line):

    velocityseries=(3e5)*(((waveseries/(1+findredshift(file,filelist)))/line)-1)
    
    velocitydataframeraw=pd.dataframe({'velocity':velocityseries,'flux':fluxseries})
    velocitydataframe=velocitydataframeraw.loc[(velocityseries>-600)&(velocityseries<600)]
    
    return velocitydataframe
    
def showvelocityplots(filelist):

    for file in filelist:
    
        filename=filenamelist[filelist.index(file)]
    
        dataframe=convertdata(file)
        
        waveseries=dataframe.iloc[:,0]
        fluxseries=dataframe.iloc[:,1]*1e27
        noiseseries=dataframe.iloc[:,2]*1e27
        
        for (line,linename) in zip(linelist,linenamelist):


            lymanvelocityseries=(3e5)*(((waveseries/(1+findredshift(file,filelist)))/1215.67)-1)
    
            lymanvelocitydataframeraw=pd.dataframe({'velocity':lymanvelocityseries,'flux':fluxseries})
            lymanvelocitydataframe=lymanvelocitydataframeraw.loc[(lymanvelocityseries>-600)&(lymanvelocityseries<600)]

            velocityseries=(3e5)*(((waveseries/(1+findredshift(file,filelist)))/line)-1)
    
            velocitydataframeraw=pd.dataframe({'velocity':velocityseries,'flux':fluxseries})
            velocitydataframe=velocitydataframeraw.loc[(velocityseries>-600)&(velocityseries<600)]
            
            plt.plot(lymanvelocitydataframe.iloc[:,0],lymanvelocitydataframe.iloc[:,1],linewidth=1,color='black',label='lyman-a',linestyle='dashed')
            
            plt.plot(velocitydataframe.iloc[:,0],velocitydataframe.iloc[:,1]*5,linewidth=0.5,color='black',label='5'+linename)
            
            plt.title('flux vs. velocity of '+linename+' line vs. lyman-alpha line for '+filename)
            
            plt.xlabel('velocity (km/s)')
            plt.ylabel('flux (erg/s/cm^2/hz)*10^27')
            
            plt.xlim([-600,600])
            
            plt.legend(loc='upper left')
            
            plt.show()

def showlymanalphavelocityplots(filelist):
    
    for file in filelist:
    
        filename=filenamelist[filelist.index(file)]
        
        if filename=='arc slit h3':
            linestyle='dashed'
        elif filename=='arc slit h1':
            linestyle='dashed'
        else:
            linestyle='solid'
         
        color=colorlist[filelist.index(file)]
    
        dataframe=convertdata(file)
        
        waveseries=dataframe.iloc[:,0]
        fluxseries=dataframe.iloc[:,1]
        noiseSeries=dataFrame.iloc[:,2]
        
        velocitySeries=(3e5)*(((waveSeries/(1+findRedshift(file,fileList)))/1215.67)-1)
        
        velocityDataFrameRaw=pd.DataFrame({'Velocity':velocitySeries,'Flux':fluxSeries})
        velocityDataFrame=velocityDataFrameRaw.loc[(velocitySeries>-600)&(velocitySeries<600)]
        
        normalizedFluxSeries=(velocityDataFrame.iloc[:,1]+velocityDataFrame.iloc[:,1].min())*(1/((velocityDataFrame.iloc[:,1]+velocityDataFrame.iloc[:,1].min()).max()))
            
        plt.plot(velocityDataFrame.iloc[:,0],normalizedFluxSeries,linewidth=1,color=color,label=fileName,linestyle=lineStyle)
        
        plt.title('Normalized Flux vs. Velocity of Lyman-a Lines')
        
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('Normalized Flux')
        
    plt.legend(loc='upper left')
    
    plt.show()

def convertBalmerData(fFile,eFile):

    #Creating a flux series
    fFitsData=fits.open(fFile,memmap=True)
    fDataFrame=pd.DataFrame(fFitsData[0].data)
    fSeries=fDataFrame.iloc[:,0]
    
    #Creating an error series
    eFitsData=fits.open(eFile,memmap=True)
    eDataFrame=pd.DataFrame(eFitsData[0].data)
    eSeries=eDataFrame.iloc[:,0]
    
    #Identifying the wavelength increment and the initial wavelength
    wInitial=10**fFitsData[0].header['CRVAL1']
    wIncrement=10**fFitsData[0].header['CDELT1']
    
    #Creating a wavelength series with a dummy array
    dArray=np.arange(len(fSeries))
    dSeries=pd.Series(dArray) 
    wSeries=dSeries*wIncrement+wInitial
   
    #Creating a data frame for the wavelength, flux, and error
    dataFrame=pd.DataFrame({'Wavelength':wSeries,'Flux':fSeries,'Error':eSeries})
    
    return dataFrame

def showBalmerSpectralPlots(balmerFFileList,balmerEFileList):
    
    for (fFile,eFile) in zip(balmerFFileList,balmerEFileList):
        
        dataFrame=convertBalmerData(fFile,eFile)
        
        wSeries=dataFrame.iloc[:,0]
        fSeries=dataFrame.iloc[:,1]
        eSeries=dataFrame.iloc[:,2]
        
        plt.plot(wSeries,fSeries,linewidth=0.5,color='black',label='Observed Flux')
        plt.plot(wSeries,eSeries,linewidth=0.5,color='red',label='Noise')
        plt.legend(loc='upper left')
        
        plt.xlabel('Observed Wavelength ')
        
        plt.show()
    
    
lymanGlobpath=os.path.join('.','Data/Spectra/Lyman Spectra/*MWdr.txt')
lymanFileList=glob(lymanGlobpath)

lymanFileNameList=['Arc Slit P1','Arc Slit F','Arc Slit H1','Arc Slit H2','Arc Slit H3','Arc Slit H4','Arc Slit H6','Arc Slit H9']

balmerFGlobpath=os.path.join('.','Data/Spectra/Balmer Spectra/Observed/*.fits')
balmerFFileList=glob(balmerFGlobpath)

balmerEGlobpath=os.path.join('.','Data/Spectra/Balmer Spectra/Error/*.fits')
balmerEFileList=glob(balmerEGlobpath)

balmerFileNameList=['','','','','','','','','']

lineList=[1215.67,1548.2,1550.77,1260.42,1393.76,1242.8,1862.79,1670.79]
lineNameList=['Lyman-α','Carbon-IVa','Carbon-IVe','Silicon-II','Silicon-IV','Nitrogen-V','Aluminum-III','Aluminum-II']

colorList=['blue','green','red','cyan','magenta','yellow','black','saddlebrown']

input1=input('Display graphs (y/n)?')
if input1=='y':
    input2=input('Show spectral plots (s), velocity plots (v), or both (b)?')
    if input2=='s':
        input3=input('For spectral plots, show moving frame (m), rest frame (r), or both frames (b)?')
        if input3=='m':
            showMovingFrameSpectralPlots(lymanFileList)
        elif input3=='r':
            showRestFrameSpectralPlots(lymanFileList)
        else:
            showMovingFrameSpectralPlots(lymanFileList)
            showRestFrameSpectralPlots(lymanFileList)
    elif input2=='v':
        input5=input('Show Lyman-α profile compilation (l), Lyman-α vs. other lines (o), or both (b)?')
        if input5=='l':
            showLymanAlphaVelocityPlots(lymanFileList)
        elif input5=='o':
            showVelocityPlots(lymanFileList)
        else:
            showLymanAlphaVelocityPlots(lymanFileList)
            showVelocityPlots(lymanFileList)
    else:
        input4=input('For spectral plots, show moving frame (m), rest frame (r), or both frames (b)?')
        if input4=='m':
            showMovingFrameSpectralPlots(lymanFileList)
        elif input4=='r':
            showRestFrameSpectralPlots(lymanFileList)
        else:
            showMovingFrameSpectralPlots(lymanFileList)
            showRestFrameSpectralPlots(lymanFileList)
        showVelocityPlots(lymanFileList)
else:
    print('cww')
