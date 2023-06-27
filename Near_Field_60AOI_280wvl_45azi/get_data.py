import jcmwave
import pandas as pd


eleclow = jcmwave.loadcartesianfields('./project_results/ElectricChiralityDensity_xy_z10.jcm')
maglow = jcmwave.loadcartesianfields('./project_results/MagneticChiralityDensity_xy_z10.jcm')
#pd.DataFrame(eleclow['field'][0][:,:,0].real + maglow['field'][0][:,:,0].real).to_excel('C0low.xlsx')
pd.DataFrame(eleclow['Y']).to_csv('eleclow_Y0.csv')
#pd.DataFrame(eleclow['X']).to_excel('eleclow_X0.xlsx')
#pd.DataFrame(eleclow['Z']).to_excel('eleclow_Z0.xlsx')

exit()
elechigh = pd.DataFrame(jcmwave.loadcartesianfields('./project_results/ElectricChiralityDensity_xy_z90.jcm')['field'])
elechighZero = pd.DataFrame(elechigh[0])
elechighOne = pd.DataFrame(elechigh[1])

maglow = pd.DataFrame(jcmwave.loadcartesianfields('./project_results/MagneticChiralityDensity_xy_z10.jcm')['field'])
maglowZero = pd.DataFrame(maglow[0])
maglowOne = pd.DataFrame(maglow[1])

maghigh = pd.DataFrame(jcmwave.loadcartesianfields('./project_results/MagneticChiralityDensity_xy_z90.jcm')['field'])
maghighZero = pd.DataFrame(maghigh[0])
maghighOne = pd.DataFrame(maghigh[1])

densities = [['eleclowOne', eleclowOne], ['eleclowZero', eleclowZero], ['elechighZero', elechighZero], ['elechighOne', elechighOne], ['maglowZero', maglowZero],['maglowOne', maglowOne], ['maghighZero', maghighZero], ['maghighOne', maghighOne]]

for density in densities:
    density[1].to_excel(density[0]+'.xlsx')
 
