import numpy as np
import math
from dataframe_format import Structure, cloude_decomp
import format_module
import json
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser
AOI = [60]
AZI = [0+i for i in range(46)]
WVL = [3000+1000*i for i in range(7)]

output_file = ".json"

keys = {}  # Create empty dictionary for keys

set = Structure(columns=['azi','aoi','svl','radius','pitch','mm', 'dmm','abs pillar','abs substrate','reflected flux','E','H'])
for wvl in WVL:
        
    format_module.save_json(output_file,dict)
    for aoi in AOI:
        parser = OptionParser()
        parser.add_option("-t", "--threads",
        action="store",type="int", dest="threads",
        help="number of threads to use")
        (options, args) = parser.parse_args()

        jcmwave.set_num_threads(options.threads)
        
        keys = {} # Create empty dictionary for keys
        results = [] # Create empty array for results
        tic = time.time() # use time() not clock() on linux system 
        
        # Set simulation parameters 
        keys = {
            'AOI': aoi,
            'radius': 1000,
            'vacuum_wavelength': wvl*1e-9,
            'uol': 1e-9,
            'display_triangulation' : 'no',
            'boundary' : 'Periodic',
            'info_level' : -1,
            'fem_degree' : 4,
            'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
            'thickness' : 1000,
            'pitch' : 3500, # pitch of square lattice
            'wvl' : wvl
            }

        # material properties
        keys['n_1'] = 1.00 # index refraction of Air
        
        store  = format_module.json_file(output_file).data
        store[aoi] = {}
 
        for azi in AZI:
            keys['sm_filename'] = '"'+'project_results/sm.jcm"'
            keys['azimuth'] = azi
            jcmwave.jcmt2jcm('boundary_conditions.jcmt', keys)
            jcmwave.jcmt2jcm('materials.jcmt', keys)
            jcmwave.jcmt2jcm('project.jcmpt', keys)
            jcmwave.jcmt2jcm('sources.jcmt', keys)
            jcmwave.jcmt2jcm('layout.jcmt', keys)
            jcmwave.solve('project.jcmp')
		
            # extract absorption data:
            fileName = './project_results/absorption.jcm'
            table = jcmwave.loadtable(fileName)
            n1 = np.nonzero(table['DomainId'] == 1) # row number for Material Id 1 (Air)
            n2 = np.nonzero(table['DomainId'] == 2) # row number for Material Id 2 (Ge Pillar)
            n3 = np.nonzero(table['DomainId'] == 3) # row number for Material Id 3 (CaF2 Substrate)

            # Power from plane wave
            #p_in = (0.5)*(Eamp**2/Z1)*(((3*np.sqrt(3))/2)*((keys['pitch']*keys['uol'])**2)))*math.cos(math.radians(keys['AOI'])) # ((Eamp V/m)/ohms)*m^2 = Watts
            #hexagon area
            #p_in = 1.0 * (((3*np.sqrt(3))/2)*((keys['pitch']*keys['uol'])**2))*math.cos(math.radians(keys['AOI'])) # Eamp is not used because PowerScaling is being used in sources.jcmt
            #square
            #p_in = 1.0 * (((keys['pitch']*keys['uol'])**2))*math.cos(math.radians(keys['AOI'])) # Eamp is not used because PowerScaling is being used in sources.jcmt
            p_in = 1.0 * (((keys['pitch']*keys['uol'])**2)) # Eamp is not used because PowerScaling is being used in sources.jcmt

            absS_pillar = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][0][n2]))/p_in)
            absP_pillar = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][1][n2]))/p_in)
         
            absS_substrate = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][0][n3]))/p_in)
            absP_substrate = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][1][n3]))/p_in)            

            # Gather Reflected Fourier Modes (Z)
            filename_fourierModes_r = './project_results/fourier_modes_r.jcm';
            fourierModes_r = jcmwave.loadtable(filename_fourierModes_r,format='named')
            powerFlux_r = jcmwave.convert2powerflux(fourierModes_r)
            
            # Reflected flux in normal direction
            reflectedS_flux = np.real(powerFlux_r['PowerFluxDensity'][0][:,2])[0]
            reflectedP_flux = np.real(powerFlux_r['PowerFluxDensity'][1][:,2])[0]

            filename_MM = './project_results/sm.jcm'
            table = jcmwave.loadtable(filename_MM)
            mm = []
            for i in range(4):
                row = []
                for j in range(4):
                    row.append(float(table['Mueller_xy'+str(1+i)+str(1+j)][0][0]))
                mm.append(row)
            dmm = cloude_decomp(mm)


            fileName = './project_results/electricfield.jcm'
            table = jcmwave.loadtable(fileName)
            n1 = np.nonzero(table['DomainId'] == 1) # row number for Material Id 1 (Air)
            n2 = np.nonzero(table['DomainId'] == 2) # row number for Material Id 2 (Ge Pillar)
            n3 = np.nonzero(table['DomainId'] == 3) # row number for Material Id 3 (CaF2 Substrate)
            print(n1)

            exit()
            fileName = './project_results/magneticfield.jcm'
            table = jcmwave.loadtable(fileName)
            n1 = np.nonzero(table['DomainId'] == 1) # row number for Material Id 1 (Air)
            n2 = np.nonzero(table['DomainId'] == 2) # row number for Material Id 2 (Ge Pillar)
            n3 = np.nonzero(table['DomainId'] == 3) # row number for Material Id 3 (CaF2 Substrate)

            entry = {
                    'azi' : 3,
                    'aoi' : 1, 
                    'wvl' : wvl, 
                    'radius' : keys['radius'], 
                    'pitch' : keys['pitch'], 
                    'mm': mm, 
                    'dmm' : dmm,
                    'abs pillar' : [absS_pillar, absP_pillar], 
                    'abs substrate' : [absS_substrate, absP_substrate], 
                    'reflected flux' : [reflectedS_flux, reflectedP_flux]
                    'E air' : 
                    'H air' : 
                    'E Ge' : 
                    'H Ge' : 
                    'E CaF2' : 
                    'H CaF2' : 
            }
            set.append(entry)

            set.save_csv(output_file)

	    
        toc = time.time() # use time() not clock() on linux system  
        t = toc-tic
        print ("Total runtime for WVL: "+str(wvl)+" AOI: "+str(aoi)+" aziumth: "+str(azi)+" - %6.4f s" % t)

