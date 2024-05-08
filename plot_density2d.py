'''  Density 2D-Plot,
     written by : Krai dated Aug11,2023
     system name : RCO3 and Glucose  '''

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR, XTC
from MDAnalysis.analysis import density
import xdrlib
from gridData import Grid

import numpy as np
import matplotlib.pyplot as plt

#creating sysytem universe 
#change here according to your system files
u = mda.Universe('../../lipid_WI.psf','../final.dcd')
from MDAnalysis.analysis.density import DensityAnalysis
#make selction for calculating density
glu = u.select_atoms("resname AGLC")
D = DensityAnalysis(glu, delta=1)
D.run()
D.density.export("Glu.dx", type="double")
grid = D.results.density.grid
#it iwll give you shape of grid matrix
grid.shape
avg = grid.mean(axis=1)  # Modify axis here to plot along Y-axis
avg.shape
fig, ax = plt.subplots()
colormap = 'cividis'
im = ax.imshow(avg.T, origin='lower',cmap=colormap)
cbar = plt.colorbar(im)
cbar.set_label('Mean Density ($\AA$$^{-3}$)',fontsize=14)
#plt.ylim([10,150])
#plt.subplots_adjust(left=0., right=0.9, top=1, bottom=0.1)
plt.xlabel('X-axis ($\AA$)',fontsize=20)
plt.ylabel('Z-axis ($\AA$)',fontsize=20)  # Change label to Z-axis
plt.tight_layout
plt.savefig("densityheatmap.png", dpi=300)


g = Grid("Glu.dx")
data1=g.grid
data2=g.midpoints
avgz=np.mean(data1,axis=(0,1))
np.savetxt('glu_z.dat',np.vstack([data2[2],avgz]).T)
file = np.loadtxt('glu_z.dat',unpack=True)
fig, ax1 = plt.subplots( )
ax1.plot(file[0,:],file[1,:]*10000,linewidth=3.0,color='red')
#plt.xlim([-75,75])
plt.xlabel('Z-axis ($\AA$)',fontsize=20)
plt.ylabel('Mean Density($x10^{-4}$ $\AA$$^{-3}$)',fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig("lineardensity.png",dpi=700)
