import sys
import os
import pandas as pd
import numpy as np
curdir = os.getcwd()
os.chdir("..")
os.chdir("..")
print(os.path.abspath(os.curdir))
#exit()
jcm_root = "C:\\Users\\Colton Bruni\\Documents\\JCMsuite"
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import jcmwave
os.chdir("Github/Multivariable_Data_Struct/")
import matplotlib.pyplot as plt
import scipy
from scipy.ndimage import median_filter
folder = "Near_Field_60AOI_batch_13_14/project_results_azi"
#folder = "LH_Gammadion_Simulation_Racemic_Trp_500/project_results/"

#remember ----
# - old azi/wvl
# - new /wvl
# - new /trans_wvl
def ChiralDensityPlot(folder, trans_wvl, heights, azis):
    for azi in azis:
        for wvl in trans_wvl:
            plug = folder+azi+'/'+wvl
            first = True 
            for h in heights:
                cfb = jcmwave.loadcartesianfields(plug+'_ElectricChiralityDensity_xy_z'+h+'.jcm')
                chi = cfb['field'][1]
                cfb = jcmwave.loadcartesianfields(plug+'_MagneticChiralityDensity_xy_z'+h+'.jcm')
                chim = cfb['field'][1]
                C0 = chi[:,:,0].real+chim[:,:,0].real
                if first:
                   C1 = C0
                   first = False
                fig, ax = plt.subplots(figsize=(8,4))
                plt.subplot(1,1,1)
            data = C1/C0
            fdata = scipy.ndimage.median_filter(data,size=(4,4))
            plot1 = plt.pcolormesh(cfb['X']*10**9, cfb['Y']*10**9, C1/C0, cmap=plt.cm.seismic, shading='gouraud')
            cb = plt.colorbar(plot1)     
            plt.savefig("Normalized Pillar Chiral Density/azi"+str(azi)+str(wvl)+"normalized.png")
azis = [str(17 + i) for i in range(4)]
lambdas = [225 + i for i in range(1)]
heights = ['104', '120']
s_wvl = []
for l in lambdas:
   s_wvl.append('wvl' + str(l))
ChiralDensityPlot(folder, s_wvl, heights, azis) 
exit()



folder = "LH_Gammadion_Simulation_Racemic_Trp_500/project_results/"
cfb = jcmwave.loadcartesianfields(folder+'_ElectricChiralityDensity_xy_z90.jcm')
chi = cfb['field'][1]
cfb = jcmwave.loadcartesianfields(folder+'_MagneticChiralityDensity_xy_z90.jcm')
chim = cfb['field'][1]
import matplotlib.pyplot as plt
C1 = chi[:,:,0].real+chim[:,:,0].real
fig, ax = plt.subplots(figsize=(8,4))

plt.subplot(1,1,1)
plot1 = plt.pcolormesh(cfb['X']*10**9, cfb['Y']*10**9, C1, cmap=plt.cm.seismic, shading='gouraud')
cb = plt.colorbar(plot1)
plt.show()
folder = "LH_Gammadion_Simulation_Racemic_Trp_500/project_results/"
cfb = jcmwave.loadcartesianfields(folder+'_ElectricChiralityDensity_xy_z120.jcm')
chi = cfb['field'][1]
cfb = jcmwave.loadcartesianfields(folder+'_MagneticChiralityDensity_xy_z120.jcm')
chim = cfb['field'][1]
import matplotlib.pyplot as plt
C1 = chi[:,:,0].real+chim[:,:,0].real
fig, ax = plt.subplots(figsize=(8,4))

plt.subplot(1,1,1)
plot1 = plt.pcolormesh(cfb['X']*10**9, cfb['Y']*10**9, C1, cmap=plt.cm.seismic, shading='gouraud')
cb = plt.colorbar(plot1)
plt.show()
exit()

cfb = jcmwave.loadcartesianfields('Near_Field_60AOI_batch_13_14/project_results_azi13/wvl212_ElectricChiralityDensity_xy_z104.jcm')
chi = cfb['field'][1]
cfb = jcmwave.loadcartesianfields('Near_Field_60AOI_batch_13_14/project_results_azi13/wvl212_MagneticChiralityDensity_xy_z104.jcm')
chim = cfb['field'][1]
import matplotlib.pyplot as plt
C1 = chi[:,:,0].real+chim[:,:,0].real
fig, ax = plt.subplots(figsize=(8,4))

plt.subplot(1,1,1)
plot1 = plt.pcolormesh(cfb['X']*10**9, cfb['Y']*10**9, C1, cmap=plt.cm.seismic, shading='gouraud')
cb = plt.colorbar(plot1)
plt.show()

cfb = jcmwave.loadcartesianfields('Near_Field_60AOI_batch_13_14/project_results_azi13/wvl212_ElectricChiralityDensity_xy_z120.jcm')
chi = cfb['field'][1]
cfb = jcmwave.loadcartesianfields('Near_Field_60AOI_batch_13_14/project_results_azi13/wvl212_MagneticChiralityDensity_xy_z120.jcm')
chim = cfb['field'][1]
import matplotlib.pyplot as plt
C1 = chi[:,:,0].real+chim[:,:,0].real
lam = 212*1e-9
c0 = 299792458
n = (1)**0.5
omega = (2*np.pi*c0)/(lam)
eps0 = 8.854*1e-12
norm_A_sqr = 2
C = np.absolute((eps0*omega)/(2*c0))*norm_A_sqr
fig, ax = plt.subplots(figsize=(8,4))

plt.subplot(1,1,1)
plot1 = plt.pcolormesh(cfb['X']*10**9, cfb['Y']*10**9, C1, cmap=plt.cm.seismic, shading='gouraud')
cb = plt.colorbar(plot1)
plt.show()

#zero = pd.DataFrame(eleclow['field'][0][:,:,0].real + maglow['field'][0][:,:,0].real)
#one = pd.DataFrame(eleclow['field'][1][:,:,0].real + maglow['field'][1][:,:,0].real)

