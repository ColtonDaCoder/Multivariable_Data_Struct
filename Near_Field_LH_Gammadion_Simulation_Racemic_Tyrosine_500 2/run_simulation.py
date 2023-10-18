import numpy as np
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser


keys = {}  # Create empty dictionary for keys
arms = [135]
for armi in arms:

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
		'gap':10,
		'uol': 1e-9,
		'display_triangulation' : 'no',
		'boundary' : 'Periodic',
		'info_level' : -1,
		'pitch': 500,
		'fem_degree' : 3,
		'arm_thickness' : armi,
		'angle': 0,
		'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
		'thickness' : 100, #thickness of the structure
		'adhesion_thickness': 5,
		'z_radius' : 5, # radius of curvature of dimer in z plane
		'z_radius_MSL' : 1, # maximum side length of z radius
        'aminoacid_thickness':20
		}

	tag_ = '_1nm_step_LH_Transm_Al_Gamm_' + str(keys['pitch']) + '_nm_pitch_' + str(keys['boundary']) + '_boundary' + 'NA_0.09_'  + '_FE_' + str(keys['fem_degree']) + '_arm_' + str(keys['arm_thickness'])+ '_gap_' + str(keys['gap']) + '_40nm_mesh_' + str(keys['aminoacid_thickness'])+ '_nm_racemic_tyr_odd'

	
			
	# material properties
	keys['n_3'] = 1.00 # index refraction of Air

	Fused_Silica_nk = np.loadtxt('../data/fused_silica_malitson.nk')
	wl_Fused_Silica_data = []; n_Fused_Silica_real = [];
	for data in Fused_Silica_nk:
	  wl_Fused_Silica_data.append(data[0]*1e-6) # e-10 for [ang], e-9 for [nm], e-6 for [um]
	  n_Fused_Silica_real.append(data[1])
		
	Al_nk = np.loadtxt('../data/Al_30A-s_FS_substrate.nk')
	wl_Al_data = []; n_Al_real = []; n_Al_imag = [] ### Aluminum data
	for data in Al_nk:
		wl_Al_data.append(data[0]*1e-9) # e-10 for [ang], e-9 for [nm], e-6 for [um]
		n_Al_real.append(data[1])
		n_Al_imag.append(data[2])
        
	Ti_nk = np.loadtxt('../data/Ti_Johnson_Christy.nk')
	wl_Ti_data = []; n_Ti_real = []; n_Ti_imag = [] ### Titanium data
	for data in Ti_nk:
		wl_Ti_data.append(data[0]*1e-6) # e-10 for [ang], e-9 for [nm], e-6 for [um]
		n_Ti_real.append(data[1])
		n_Ti_imag.append(data[2])
    
	Aminoacid_nk = np.loadtxt('../data/L-tyr_nx_kx.txt')
	wl_Aminoacid_data = []; n_Aminoacid_real = []; n_Aminoacid_imag = []
	for data in Aminoacid_nk:
		wl_Aminoacid_data.append(data[0]*1e-9) # e-10 for [ang], e-9 for [nm], e-6 for [um]
		n_Aminoacid_real.append(data[1])
		n_Aminoacid_imag.append(data[2])

			
	lambdas = np.linspace(228, 228, 1)*1e-9

	for keys['vacuum_wavelength'] in lambdas:
		print('Wavelength : %3.2f nm' % (keys['vacuum_wavelength']*1e9))
		keys['n_1'] = np.interp(keys['vacuum_wavelength'], wl_Fused_Silica_data, n_Fused_Silica_real)
		keys['n_2'] = np.interp(keys['vacuum_wavelength'], wl_Ti_data, n_Ti_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Ti_data, n_Ti_imag)
		keys['n_4'] = np.interp(keys['vacuum_wavelength'], wl_Al_data, n_Al_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Al_data, n_Al_imag)
		keys['n_5'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_imag)
		keys['sm_filename'] = '"'+'project_results/sm_%dnm.jcm"' % int(keys['vacuum_wavelength']*1e9)	#for saving MM file
		jcmwave.jcmt2jcm('./boundary_conditions.jcmt', keys)
		jcmwave.jcmt2jcm('./materials.jcmt', keys)
		jcmwave.jcmt2jcm('./project.jcmpt', keys)
		jcmwave.jcmt2jcm('./sources.jcmt', keys)
		jcmwave.jcmt2jcm('./layout.jcmt', keys)
		jcmwave.solve('./project.jcmp')
		

		
	toc = time.time() # use time() not clock() on linux system  
	t = toc-tic
	print ("Total runtime for "+tag_+": %6.4f s" % t)
