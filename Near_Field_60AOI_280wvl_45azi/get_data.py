import jcmwave
import pandas as pd




 
eleclow = pd.DataFrame(jcmwave.loadcartesianfields('ElectricChirality_xy_z10.jcm')['field'])

elechigh = pd.DataFrame(jcmwave.loadcartesianfields('ElectricChiralityDensity_xy_z90.jcm')['field'])

maglow = pd.DataFrame(jcmwave.loadcartesianfields('MagneticChirality_xy_z10.jcm')['field'])

maghigh = pd.DataFrame(jcmwave.loadcartesianfields('MagneticChiralityDensity_xy_z90.jcm')['field'])

densities = [['eleclow', eleclow], ['elechigh', elechigh], ['maglow', maglow], ['maghigh', maghigh]]

for density in densities:
    density[1].to_excel(density[0]+'.xlsx')
 
