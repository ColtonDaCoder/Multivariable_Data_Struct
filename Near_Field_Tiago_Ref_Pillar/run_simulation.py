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
            'AOI': 60,
            'radius': 50,
            'uol': 1e-9,
            'display_triangulation' : 'no',
            'boundary' : 'Periodic',
            'info_level' : -1,
            'fem_degree' : 3,
            'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
            'thickness' : 50,
            'pitch' : 175, # pitch of square lattice (gammadion)
            'z_radius' : 5, # radius of curvature of dimer in z plane
            'z_radius_MSL' : 1, # maximum side length of z radius
            'azimuth' : 14
		}	
			
	# material properties
	keys['n_1'] = 1.00 # index refraction of Air

	Al_nk = np.loadtxt('Al_OMEL_mfp.nk')   #Al not GOLD!!!!!!!!!!!!!!!!!!!!!!!
	wl_Al_data = []; n_Al_real = []; n_Al_imag = []
	for data in Al_nk:	     
		wl_Al_data.append(data[0]*1e-9) # e-10 for [ang], e-9 for [nm], e-6 for [um]
		n_Al_real.append(data[1])
		n_Al_imag.append(data[2])

	Aminoacid_nk = np.loadtxt('L-tyr_nx_kx.txt')
	wl_Aminoacid_data = []; n_Aminoacid_real = []; n_Aminoacid_imag = []
	for data in Aminoacid_nk:
		wl_Aminoacid_data.append(data[0]*1e-9) # e-10 for [ang], e-9 for [nm], e-6 for [um]
		n_Aminoacid_real.append(data[1])
		n_Aminoacid_imag.append(data[2])
        
	Aminoacid_kappa = np.loadtxt('L-tyr_kappa.txt')
	wl_Aminoacid_kappa = []; kappa_Aminoacid_real = []; kappa_Aminoacid_imag = []
	for data in Aminoacid_kappa:
		wl_Aminoacid_kappa.append(data[0]*1e-9) # e-10 for [ang], e-9 for [nm], e-6 for [um]
		kappa_Aminoacid_real.append(data[1])
		kappa_Aminoacid_imag.append(data[2])

		
        
			
	lambdas = np.linspace(220, 220, 1)*1e-9

	for keys['vacuum_wavelength'] in lambdas:
		print('Wavelength : %3.2f nm' % (keys['vacuum_wavelength']*1e9))
		keys['n_2'] = np.interp(keys['vacuum_wavelength'], wl_Al_data, n_Al_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Al_data, n_Al_imag)
		keys['n_4'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_imag)
		keys['kappa'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_kappa, kappa_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_kappa, kappa_Aminoacid_imag)
		keys['sm_filename'] = '"'+'project_results/sm_%dnm.jcm"' % int(keys['vacuum_wavelength']*1e9)	#for saving MM file
		jcmwave.jcmt2jcm('./boundary_conditions.jcmt', keys)
		jcmwave.jcmt2jcm('./materials.jcmt', keys)
		jcmwave.jcmt2jcm('./project.jcmpt', keys)
		jcmwave.jcmt2jcm('./sources.jcmt', keys)
		jcmwave.jcmt2jcm('./layout.jcmt', keys)
		jcmwave.solve('./project.jcmp')
		

		
	toc = time.time() # use time() not clock() on linux system  
	t = toc-tic
	print ("Total runtime for : %6.4f s" % t)
