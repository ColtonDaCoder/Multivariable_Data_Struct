#!/usr/bin/env python
# coding: utf-8

#6-25-21: LH: gammadion
#6-28-21: RH gammadion
#6-29-21: Racemic gammadion
#6-30-21: LH spiral
#7-1-21: RH spiral
#7-2-21: Racemic spiral

# Standard library imports
import os
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import interpolate
import math
#import scipy

from matplotlib.ticker import FormatStrFormatter 
from matplotlib.ticker import FormatStrFormatter 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy.linalg import logm, expm

# Import custom MM functions
from MM_functions import extract_jcmwave, Cloude_Decomposition, sum_MMs

# Load LH gammadion data into resultsLH
os.chdir(r"C:\Users\kmcpe\Documents\Leite_etal_bianisotropic\jcmwave\MM_Spectrum_Al_TS_RH_gammadion_per_AOI20_pitch840_arm120_t40") #Define path where 
os.getcwd()
azi = []
azi = np.arange(0,62,2)
resultsLH = []
azimuthLH = []

for i in range(len(azi)):
    angle = azi[i]
    azimuthLH.append(angle)
    txt = "LH_gammadion_Azi" + str(angle) + ".txt"
    MM = extract_jcmwave(txt,210,400) #load Muelle Matrix Exported from CompleteEase (Woollam)
    dMM = Cloude_Decomposition(MM)
    resultsLH.append(dMM)


# Polar plots of Mueller Matrices for LH Gammadions

plt.rcParams["font.family"] = "Arial"
fig, axs = plt.subplots(4, 4,figsize=(24,24), dpi = 300,subplot_kw=dict(projection='polar'))

y = []
for i in range(len(MM[0])):
    y.append(MM[0][i])
    
x = []
for i in range(len(azimuthLH)):
    x.append(azimuthLH[i])

X,Y = np.meshgrid(np.radians(x),y)

n=len(azimuthLH)
m=len(MM[0])

z1 = np.zeros((m,n))

z2 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z2[i][j] = resultsLH[j][1][i]

z3 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z3[i][j] = resultsLH[j][2][i]

z4 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z4[i][j] = resultsLH[j][3][i]

z5 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z5[i][j] = resultsLH[j][4][i]
z6 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z6[i][j] = resultsLH[j][5][i]

z7 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z7[i][j] = resultsLH[j][6][i]

z8 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z8[i][j] = resultsLH[j][7][i]

        
z9 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z9[i][j] = resultsLH[j][8][i]

z10 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z10[i][j] = resultsLH[j][9][i]

z11 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z11[i][j] = resultsLH[j][10][i]        
z12 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z12[i][j] = resultsLH[j][11][i]   
        
z13 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z13[i][j] = resultsLH[j][12][i]           
        
z14 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z14[i][j] = resultsLH[j][13][i]          
        
z15 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z15[i][j] = resultsLH[j][14][i]
        
z16 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z16[i][j] = resultsLH[j][15][i]  
z17 = np.zeros((m,n))

for j in range(len(azimuthLH)):
    for i in range(len(MM[0])):
        z17[i][j] = resultsLH[j][16][i]  
        

cm = plt.cm.seismic

axs[0,0].pcolormesh(X,Y,z1,vmin=-1, vmax=1, shading='gouraud', antialiased=True, cmap = cm)
fig.text(0.301, 0.832, 'nm',size = 16)

axs[0,0].set_rticks([200,250,300,350,400])

axs[0,0].set_rlim(200, 400)
axs[0,0].set_xlim(0, 2*np.pi)
axs[0,0].grid(color = 'gray', linestyle = '--', linewidth = 1)
axs[0,0].grid(True)
axs[0,0].tick_params(labelsize=16)

# M12
ax1 = axs[0,1].pcolormesh(X,Y,z2,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm) # was +/- 0.01 for AOI20
axs[0,1].set_rlim(200, 400)
axs[0,1].set_xlim(0, 2*np.pi)
axs[0,1].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[0,1].grid(False)
axs[0,1].set_rticks([])
axs[0,1].set_xticks([])
axs[0,1].set_title('m12',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax1, ax = axs[0,1])
cb.ax.tick_params(labelsize=18)

#M13
ax2 = axs[0,2].pcolormesh(X,Y,z3,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm)
axs[0,2].set_rlim(200, 400)
axs[0,2].set_xlim(0, 2*np.pi)
axs[0,2].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[0,2].grid(False)
axs[0,2].set_rticks([])
axs[0,2].set_xticks([])
axs[0,2].set_title('m13',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax2, ax = axs[0,2])
cb.ax.tick_params(labelsize=18)

#M14
ax3 = axs[0,3].pcolormesh(X,Y,z4,vmin=-0.05, vmax=0.05, shading='gouraud', antialiased=True, cmap = cm)
axs[0,3].set_rlim(200, 400)
axs[0,3].set_xlim(0, 2*np.pi)
axs[0,3].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[0,3].grid(False)
axs[0,3].set_rticks([])
axs[0,3].set_xticks([])
axs[0,3].set_title('m14',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax3, ax = axs[0,3])
cb.ax.tick_params(labelsize=18)

# M21
ax4 = axs[1,0].pcolormesh(X,Y,z5,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm) # was +/- 0.01 for AOI20
axs[1,0].set_rlim(200, 400)
axs[1,0].set_xlim(0, 2*np.pi)
axs[1,0].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[1,0].grid(False)
axs[1,0].set_rticks([])
axs[1,0].set_xticks([])
axs[1,0].set_title('m21',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax4, ax = axs[1,0])
cb.ax.tick_params(labelsize=18)

# M22
ax5 = axs[1,1].pcolormesh(X,Y,z6,vmin=-1.0, vmax=1.0, shading='gouraud', antialiased=True, cmap = cm)
axs[1,1].set_rlim(200, 400)
axs[1,1].set_xlim(0, 2*np.pi)
axs[1,1].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[1,1].grid(False)
axs[1,1].set_rticks([])
axs[1,1].set_xticks([])
axs[1,1].set_title('m22',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax5, ax = axs[1,1])
cb.ax.tick_params(labelsize=18)

# M23
ax6 = axs[1,2].pcolormesh(X,Y,z7,vmin=-0.05, vmax=0.05, shading='gouraud', antialiased=True, cmap = cm)
axs[1,2].set_rlim(200, 400)
axs[1,2].set_xlim(0, 2*np.pi)
axs[1,2].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[1,2].grid(False)
axs[1,2].set_rticks([])
axs[1,2].set_xticks([])
axs[1,2].set_title('m23',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax6, ax = axs[1,2])
cb.ax.tick_params(labelsize=18)

# M24
ax7 = axs[1,3].pcolormesh(X,Y,z8,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm)
axs[1,3].set_rlim(200, 400)
axs[1,3].set_xlim(0, 2*np.pi)
axs[1,3].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[1,3].grid(False)
axs[1,3].set_rticks([])
axs[1,3].set_xticks([])
axs[1,3].set_title('m24',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax7, ax = axs[1,3])
cb.ax.tick_params(labelsize=18)

# M31
ax8 = axs[2,0].pcolormesh(X,Y,z9,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm)
axs[2,0].set_rlim(200, 400)
axs[2,0].set_xlim(0, 2*np.pi)
axs[2,0].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[2,0].grid(False)
axs[2,0].set_rticks([])
axs[2,0].set_xticks([])
axs[2,0].set_title('m31',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax8, ax = axs[2,0])
cb.ax.tick_params(labelsize=18)

# M32
ax9 = axs[2,1].pcolormesh(X,Y,z10,vmin=-0.05, vmax=0.05, shading='gouraud', antialiased=True, cmap = cm)
axs[2,1].set_rlim(200, 400)
axs[2,1].set_xlim(0, 2*np.pi)
axs[2,1].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[2,1].grid(False)
axs[2,1].set_rticks([])
axs[2,1].set_xticks([])
axs[2,1].set_title('m32',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax9, ax = axs[2,1])
cb.ax.tick_params(labelsize=18)

# M33
ax10 = axs[2,2].pcolormesh(X,Y,z11,vmin=-1.0, vmax=1.0, shading='gouraud', antialiased=True, cmap = cm)
axs[2,2].set_rlim(200, 400)
axs[2,2].set_xlim(0, 2*np.pi)
axs[2,2].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[2,2].grid(False)
axs[2,2].set_rticks([])
axs[2,2].set_xticks([])
axs[2,2].set_title('m33',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax10, ax = axs[2,2])
cb.ax.tick_params(labelsize=18)

# M34
ax11 = axs[2,3].pcolormesh(X,Y,z12,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm) # was +/- 0.01 for AOI20
axs[2,3].set_rlim(200, 400)
axs[2,3].set_xlim(0, 2*np.pi)
axs[2,3].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[2,3].grid(False)
axs[2,3].set_rticks([])
axs[2,3].set_xticks([])
axs[2,3].set_title('m34',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax11, ax = axs[2,3])
cb.ax.tick_params(labelsize=18)

# M41
ax12 = axs[3,0].pcolormesh(X,Y,z13,vmin=-0.05, vmax=0.05, shading='gouraud', antialiased=True, cmap = cm)
axs[3,0].set_rlim(200, 400)
axs[3,0].set_xlim(0, 2*np.pi)
axs[3,0].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[3,0].grid(False)
axs[3,0].set_rticks([])
axs[3,0].set_xticks([])
axs[3,0].set_title('m41',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax12, ax = axs[3,0])
cb.ax.tick_params(labelsize=18)

# M42
ax13 = axs[3,1].pcolormesh(X,Y,z14,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm)
axs[3,1].set_rlim(200, 400)
axs[3,1].set_xlim(0, 2*np.pi)
axs[3,1].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[3,1].grid(False)
axs[3,1].set_rticks([])
axs[3,1].set_xticks([])
axs[3,1].set_title('m42',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax13, ax = axs[3,1])
cb.ax.tick_params(labelsize=18)

# M43
ax14 = axs[3,2].pcolormesh(X,Y,z15,vmin=-0.1, vmax=0.1, shading='gouraud', antialiased=True, cmap = cm) # was +/- 0.01 for AOI20
axs[3,2].set_rlim(200, 400)
axs[3,2].set_xlim(0, 2*np.pi)
axs[3,2].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[3,2].grid(False)
axs[3,2].set_rticks([])
axs[3,2].set_xticks([])
axs[3,2].set_title('m43',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax14, ax = axs[3,2])
cb.ax.tick_params(labelsize=18)

# M44
ax15 = axs[3,3].pcolormesh(X,Y,z16,vmin=-1, vmax=1, shading='gouraud', antialiased=True, cmap = cm)
axs[3,3].set_rlim(200, 400)
axs[3,3].set_xlim(0, 2*np.pi)
axs[3,3].grid( color = 'gray', linestyle = '--', linewidth = 1)
axs[3,3].grid(False)
axs[3,3].set_rticks([])
axs[3,3].set_xticks([])
axs[3,3].set_title('m44',fontsize=18,fontweight='bold')
cb = plt.colorbar(ax15, ax = axs[3,3])
cb.ax.tick_params(labelsize=18)

os.chdir(r"C:\Users\kmcpe\Documents\Leite_etal_bianisotropic\figures")
os.getcwd()
fig.savefig('MM_Azimuth_simulated_RH_gammadion_AOI20_per_pitch840_arm120_t40.png', dpi=fig.dpi)