import numpy as np
import math
from dataframe_format import Structure, cloude_decomp
import format_module
import json
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser
AOI = [60]
AZI = [35+i for i in range(11)]
#AZI = [35]
#WVL = [6200, 6300]
WVL = [6200]
#WVL = [5450]

keys = {}  # Create empty dictionary for keys

set = Structure(columns=['azi','aoi','wvl','radius','pitch','mm', 'dmm','abs pillar','abs substrate','reflected flux','E','H'])
chirality_list = [False]
for is_Chiral in chirality_list:
    for wvl in WVL:
            
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
                'radius': 500,
                'vacuum_wavelength': wvl*1e-9,
                'uol': 1e-9,
                'display_triangulation' : 'no',
                'boundary' : 'Periodic',
                'info_level' : -1,
                'fem_degree' : 3,
                'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
                'thickness' : 1000,
                'pitch' : 4000, # pitch of square lattice
                'wvl' : wvl
                }

            # material properties
            keys['n_3'] = 1.00 # index refraction of Air

            Au_nk = np.loadtxt('data/Au_Babar_micron.nk') # You would need to change this to the silver data file I will send you
            wl_Au_data = []; n_Au_real = []; n_Au_imag = [] ### ALUMINUM NOT GOLD!
            for data in Au_nk:
                wl_Au_data.append(data[0]*1e-6) # e-10 for [ang], e-9 for [nm], e-6 for [um]
                n_Au_real.append(data[1])
                n_Au_imag.append(data[2])
		
            for azi in AZI:
                print("Wavelength: " + str(wvl))
                keys['vacuum_wavelength'] = keys['uol']*wvl
                Au_Material = np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_imag)
                keys['n_2'] = Au_Material
                keys['n_4'] = Au_Material
                keys['sm_filename'] = '"'+'project_results/sm.jcm"'
                jcmwave.jcmt2jcm('./boundary_conditions.jcmt', keys)
                jcmwave.jcmt2jcm('./materials.jcmt', keys)
                jcmwave.jcmt2jcm('./project.jcmpt', keys)
                jcmwave.jcmt2jcm('./sources.jcmt', keys)
                jcmwave.jcmt2jcm('./layout.jcmt', keys)
                jcmwave.solve('./project.jcmp')
		
                ## Gather Reflected Fourier Modes (+Z)
                filename_fourierModes_r = './project_results/fourier_modes_r.jcm';
                fourierModes_r = jcmwave.loadtable(filename_fourierModes_r,format='named')
                powerFlux_r = jcmwave.convert2powerflux(fourierModes_r)

                ## Reflected flux in normal direction
                P_s_r = np.sum(powerFlux_r['PowerFluxDensity'][0][:, 2]);
                P_p_r = np.sum(powerFlux_r['PowerFluxDensity'][1][:, 2]); 

                filename_MM = './project_results/sm.jcm'
                table = jcmwave.loadtable(filename_MM)
                mm = []
                for i in range(4):
                    row = []
                    for j in range(4):
                        row.append(float(table['Mueller_xy'+str(1+i)+str(1+j)][0][0]))
                    mm.append(row)
                dmm = cloude_decomp(mm)

                entry = {
                        'azi' : azi,
                        'aoi' : aoi, 
                        'wvl' : wvl, 
                        'pitch' : keys['pitch'], 
                        'mm': [mm], 
                        'dmm' : [dmm],
                        'reflected flux' : [[P_s_r, P_p_r]]
                }
                set.append(entry)
                set.save_csv('pillar_lattice_copies_1.csv')

	    
            toc = time.time() # use time() not clock() on linux system  
            t = toc-tic
            print ("Total runtime for WVL: "+str(wvl)+" AOI: "+str(aoi)+" aziumth: "+str(azi)+" - %6.4f s" % t)

