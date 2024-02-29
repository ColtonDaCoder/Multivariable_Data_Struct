import numpy as np
import math
from dataframe_format import Structure, cloude_decomp
import format_module
import json
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser
AOI = [60]
AZI = [2]
WVL = [5600]

keys = {}  # Create empty dictionary for keys

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
        keys['kappa'] = 0.00038 + 1j*0.00038        
        keys['n_1'] = 1.00 # index refraction of Air
         
        for azi in AZI:
            keys['sm_filename'] = '"'+'project_results/sm.jcm"'
            keys['azimuth'] = azi
            jcmwave.jcmt2jcm('boundary_conditions.jcmt', keys)
            jcmwave.jcmt2jcm('materials.jcmt', keys)
            jcmwave.jcmt2jcm('project.jcmpt', keys)
            jcmwave.jcmt2jcm('sources.jcmt', keys)
            jcmwave.jcmt2jcm('layout.jcmt', keys)
            jcmwave.solve('project.jcmp')
		 
        toc = time.time() # use time() not clock() on linux system  
        t = toc-tic
        print ("Total runtime for WVL: "+str(wvl)+" AOI: "+str(aoi)+" aziumth: "+str(azi)+" - %6.4f s" % t)

