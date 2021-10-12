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
psep=df.iloc[:,4].to_numpy()
fwhmr=df.iloc[:,5].to_numpy()
fwhmb=df.iloc[:,7].to_numpy()

plt.close('all')
fig,ax=plt.subplots(1,1)

plt.scatter(ew,fwhmr,c='black',marker='.')
ax.set_xlabel('Equivalent Width (Å)')
ax.set_ylabel('FWHM (Red Peak) (Å)')

plt.savefig(figs+'/ewfwhmr.png',dpi=100,bbox_inches='tight',overwrite=True)
plt.show()

plt.close('all')
fig,ax=plt.subplots(1,1)

plt.scatter(np.delete(ew,[1,3]),np.delete(fwhmb,[1,3]),c='black',marker='.')
ax.set_xlabel('Equivalent Width (Å)')
ax.set_ylabel('FWHM (Blue Peak) (Å)')

plt.savefig(figs+'/ewfwhmb.png',dpi=100,bbox_inches='tight',overwrite=True)
plt.show()

plt.close('all')
fig,ax=plt.subplots(1,1)

plt.scatter(np.delete(ew,[1,3]),np.delete(psep,[1,3]),c='black',marker='.')
ax.set_xlabel('Equivalent Width (Å)')
ax.set_ylabel('Peak Separation (km s$^{-1}$)')

plt.savefig(figs+'/ewpsep.png',dpi=100,bbox_inches='tight',overwrite=True)
plt.show()
