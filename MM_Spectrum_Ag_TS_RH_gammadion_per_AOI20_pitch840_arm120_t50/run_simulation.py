import numpy as np
from dataframe_format import Structure, cloude_decomp
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--threads",action="store",type="int", dest="threads",help="number of threads to use")
(options, args) = parser.parse_args()
jcmwave.set_num_threads(options.threads)

tic = time.time() # use time() not clock() on linux system
AOI = 60
AZI = 22
WVL = [i*50+4000 for i in range(80)]
aoi = AOI
azi = AZI

set = Structure(columns=['azi','aoi','wvl','radius','pitch','mm', 'dmm','reflected flux'])

    
keys = {} # Create empty dictionary for keys
results = [] # Create empty array for results

# Set simulation parameters
keys = {
    'dielectric': 10, # This was the variable "spacer" before, but it isn't being used in the layout file
    'AOI':60, # This is the angle of incidence
    'gap':378,
    'uol': 1e-9,
    'display_triangulation' : 'no',
    'boundary' : 'Periodic', # Periodic or Transparent
    'info_level' : -1,
    'pitch': 7560,
    'fem_degree' : 2,
    'arm_thickness' : 1080,
    'angle': AZI, # This rotates the layout, so we will use it to change the azimuth
    'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
    'thickness' : 450, #thickness of the gammadion
    'film': 900, # Thickness of the flat film
    'MSL':900, #MaximumSideLength for the gammadion mesh (not the z mesh)
    'z_radius' : 45, # radius of curvature of gammadion in the z plane
    'z_radius_MSL' : 9 # maximum side length of z radius
    }

keys['delta_k'] = (2*np.pi)/(20*keys['pitch']*1e-9) # For the discritization of the continuous Fourier spectrum into modes (only for isolated setup)

# material properties
keys['n_3'] = 1.00 # index refraction of Air
    
Au_nk = np.loadtxt('../data/Au_Babar_micron.nk') # You would need to change this to the silver data file I will send you
wl_Au_data = []; n_Au_real = []; n_Au_imag = [] ### ALUMINUM NOT GOLD!
for data in Au_nk:
    wl_Au_data.append(data[0]*1e-6) # e-10 for [ang], e-9 for [nm], e-6 for [um]
    n_Au_real.append(data[1])
    n_Au_imag.append(data[2])
    
for wvl in WVL:
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
            'radius' : keys['radius'], 
            'pitch' : keys['pitch'], 
            'mm': [mm], 
            'dmm' : [dmm],
            'reflected flux' : [[P_s_r, P_p_r]]
    }
    set.append(entry)
    set.save_csv('Au_on_Au_Gammadion.csv')

toc = time.time() # use time() not clock() on linux system  
t = toc-tic
print ("Total runtime for WVL: "+str(wvl)+" AOI: "+str(aoi)+" aziumth: "+str(azi)+" - %6.4f s" % t)
