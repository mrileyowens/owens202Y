import sys

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=100
plt.rcParams['lines.linewidth']=0.5

home='C://Users/15136/OneDrive - University of Cincinnati/Documents/Research/sunburstarc'
data=home+'/data'
figs=home+'/figs'

file=data+'/stats.txt'

df=pd.read_csv(file,delimiter=' ',header=0,skiprows=3)

ew=df.iloc[:,1].to_numpy()
psep=np.delete(df.iloc[:,4].to_numpy(),[1,3])
fwhmr=df.iloc[:,8].to_numpy()
fwhmb=np.delete(df.iloc[:,10].to_numpy(),[1,3])
lyc=df.iloc[:,11].to_numpy()

plt.close('all')
fig,ax=plt.subplots(1,1)

plt.scatter(ew,fwhmr,c='black',marker='.')
plt.scatter(np.delete(ew,lyc==False),np.delete(fwhmr,lyc==False),c='red',marker='.',label='LyC Leakers')
plt.scatter(np.mean(ew),np.mean(fwhmr),c='black',marker='.')
plt.errorbar(np.mean(ew),np.mean(fwhmr),xerr=np.std(ew)/np.sqrt(8),yerr=np.std(fwhmr)/np.sqrt(8),c='black')
ax.set_xlabel('Equivalent Width (Å)')
ax.set_ylabel('FWHM (Red Peak) (km s$^{-1}$)')
ax.legend(loc='upper right')

plt.savefig(figs+'/ewfwhmr.png',dpi=100,bbox_inches='tight',overwrite=True)
plt.show()

plt.close('all')
fig,ax=plt.subplots(1,1)

plt.scatter(np.delete(ew,[1,3]),fwhmb,c='black',marker='.')
plt.scatter(np.delete(np.delete(ew,[1,3]),np.delete(lyc,[1,3])==False),np.delete(fwhmb,np.delete(lyc,[1,3])==False),c='red',marker='.',label='LyC Leakers')
plt.scatter(np.mean(np.delete(ew,[1,3])),np.mean(fwhmb),c='black',marker='.')
plt.errorbar(np.mean(np.delete(ew,[1,3])),np.mean(fwhmb),xerr=np.std(np.delete(ew,[1,3]))/np.sqrt(6),yerr=np.std(fwhmb)/np.sqrt(6),c='black')
ax.set_xlabel('Equivalent Width (Å)')
ax.set_ylabel('FWHM (Blue Peak) (km s$^{-1}$)')
ax.legend(loc='upper right')

plt.savefig(figs+'/ewfwhmb.png',dpi=100,bbox_inches='tight',overwrite=True)
plt.show()

plt.close('all')
fig,ax=plt.subplots(1,1)

plt.scatter(np.delete(ew,[1,3]),psep,c='black',marker='.')
plt.scatter(np.delete(np.delete(ew,[1,3]),np.delete(lyc,[1,3])==False),np.delete(psep,np.delete(lyc,[1,3])==False),c='red',marker='.',label='LyC Leakers')
plt.scatter(np.mean(np.delete(ew,[1,3])),np.mean(psep),c='black',marker='.')
plt.errorbar(np.mean(np.delete(ew,[1,3])),np.mean(psep),xerr=np.std(np.delete(ew,[1,3]))/np.sqrt(6),yerr=np.std(np.delete(psep,[1,3]))/np.sqrt(6),c='black')
ax.set_xlabel('Equivalent Width (Å)')
ax.set_ylabel('Peak Separation (km s$^{-1}$)')
ax.legend(loc='upper right')

plt.savefig(figs+'/ewpsep.png',dpi=100,bbox_inches='tight',overwrite=True)
plt.show()
