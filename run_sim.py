import numpy as np
import format_module
import json
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser
AOI = [i*4+20 for i in range(11)]
AZI = [i for i in range(46)]

WVL = [220, 222]

#AOI = [20]
#AZI = [34]

#WVL = [230]

output_folder = "kappa"
kappa = True

keys = {}  # Create empty dictionary for keys

set = Structure()
for wvl in WVL:
    output_file = output_folder+"/fe3_220_222.json"
    Au_nk = np.loadtxt('Al_OMEL_mfp.nk')   #Al not GOLD!!!!!!!!!!!!!!!!!!!!!!!
    wl_Au_data = []; n_Au_real = []; n_Au_imag = []
    for data in Au_nk:
        wl_Au_data.append(data[0]*1e-9) # e-10 for [ang], e-9 for [nm], e-6 for [um]
        n_Au_real.append(data[1])
        n_Au_imag.append(data[2])

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

    dict = {}
    format_module.save_json(output_file,dict)
    for aoi in AOI:
        parser = OptionParser()
        parser.add_option("-t", "--threads",
        action="store",type="int", dest="threads",
        help="number of threads to use")
        (options, args) = parser.parse_args()

        jcmwave.set_num_threads(options.threads)
        
        keys = {} # Create empty dictionary for keys
        dict[aoi] = {}
        results = [] # Create empty array for results
        tic = time.time() # use time() not clock() on linux system 
        # Set simulation parameters
        
        keys = {
            'AOI': aoi,
            'radius': 50,
            'vacuum_wavelength': wvl*1e-9,
            'uol': 1e-9,
            'display_triangulation' : 'no',
            'boundary' : 'Periodic',
            'info_level' : -1,
            'fem_degree' : 3,
            'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
            'thickness' : 50,
            'pitch' : 175, # pitch of square lattice (gammadion)
            'z_radius' : 5, # radius of curvature of dimer in z plane
            'z_radius_MSL' : 1 # maximum side length of z radius
            }

        # material properties
        keys['n_1'] = 1.00 # index refraction of Air
        
        store  = format_module.json_file(output_file).data
        store[aoi] = {}
 
        for azi in AZI:
            #print('Wavelength : %3.2f nm' % (keys['vacuum_wavelength']*1e9))
            keys['n_2'] = np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_imag)
            keys['n_4'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_imag)
            keys['kappa'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_kappa, kappa_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_kappa, kappa_Aminoacid_imag)
            #keys['sm_filename'] = '"'+'project_results/sm.jcm"'
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
            n2 = np.nonzero(table['DomainId'] == 2) # row number for Material Id 2 (Al Pillar)
            n3 = np.nonzero(table['DomainId'] == 3) # row number for Material Id 3 (Al Film)
            n4 = np.nonzero(table['DomainId'] == 4) # row number for Material Id 4 (Amino Acid)

            # Power from plane wave
            #p_in = (0.5)*(Eamp**2/Z1)*(((3*np.sqrt(3))/2)*((keys['pitch']*keys['uol'])**2)))*math.cos(math.radians(keys['AOI'])) # ((Eamp V/m)/ohms)*m^2 = Watts
            p_in = 1.0 * (((3*np.sqrt(3))/2)*((keys['pitch']*keys['uol'])**2))*math.cos(math.radians(keys['AOI'])) # Eamp is not used because PowerScaling is being used in sources.jcmt

            absR_pillar = np.real(table['ElectromagneticFieldAbsorption'][0][n2])/p_in
            absL_pillar = np.real(table['ElectromagneticFieldAbsorption'][1][n2])/p_in
         
            absR_film = np.real(table['ElectromagneticFieldAbsorption'][0][n3])/p_in
            absL_film = np.real(table['ElectromagneticFieldAbsorption'][1][n3])/p_in
            
            absR_amino = np.real(table['ElectromagneticFieldAbsorption'][0][n4])/p_in
            absL_amino = np.real(table['ElectromagneticFieldAbsorption'][1][n4])/p_in
                

            # Gather Reflected Fourier Modes (Z)
            filename_fourierModes_r = './project_results/fourier_modes_r.jcm';
            fourierModes_r = jcmwave.loadtable(filename_fourierModes_r,format='named')
            powerFlux_r = jcmwave.convert2powerflux(fourierModes_r)

            # Reflected flux in normal direction
            P_s_r = np.sum(powerFlux_r['PowerFluxDensity'][0][:, 2]);
            P_p_r = np.sum(powerFlux_r['PowerFluxDensity'][1][:, 2]); 
		
            filename_MM = './project_results/sm.jcm'
            table = jcmwave.loadtable(filename_MM)
            row1 = [
	        float(table['Mueller_xy11'][0]), 
		float(table['Mueller_xy12'][0]),
		float(table['Mueller_xy13'][0]),
		float(table['Mueller_xy14'][0])]
            row2 = [
		float(table['Mueller_xy21'][0]),
		float(table['Mueller_xy22'][0]),
		float(table['Mueller_xy23'][0]),
		float(table['Mueller_xy24'][0])]
            row3 = [
		float(table['Mueller_xy31'][0]),
		float(table['Mueller_xy32'][0]),
		float(table['Mueller_xy33'][0]),
		float(table['Mueller_xy34'][0])]
            row4 = [
		float(table['Mueller_xy41'][0]),
		float(table['Mueller_xy42'][0]),
		float(table['Mueller_xy43'][0]),
		float(table['Mueller_xy44'][0])]
            mm = [row1, row2, row3, row4]

            id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
            set.append(Tag(id_names, [azi, aoi, wvl, keys['radius'], keys['pitch'], keys['height'], kappa]), ["abs pillar", "abs film", "abs amino", "reflected flux", "mm"], [[absR_pillar, absL_pillar], [absR_film, absL_film], [absR_amino, absL_amino], [P_s_r, P_p_r], mm])
            set.save_json(output_file)
	    
        toc = time.time() # use time() not clock() on linux system  
        t = toc-tic
        print ("Total runtime for WVL: "+str(wvl)+" AOI: "+str(aoi)+" aziumth: "+str(azi)+" - %6.4f s" % t)

