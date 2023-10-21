import numpy as np
import math
from data_format import Structure, Tag
import format_module
import json
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser
AOI = [60]
AZI = [38+i for i in range(5)]
WVL = [275+i for i in range(6)]
kappa = True 

output_file = "60AOI_near.json"

keys = {}  # Create empty dictionary for keys

set = Structure()
for wvl in WVL:
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
            'z_radius_MSL' : 1, # maximum side length of z radius
            'isKappa' : kappa,
            'wvl' : wvl
            }

        # material properties
        keys['n_1'] = 1.00 # index refraction of Air
        
        store  = format_module.json_file(output_file).data
        store[aoi] = {}
 
        for azi in AZI:
            #print('Wavelength : %3.2f nm' % (keys['vacuum_wavelength']*1e9))
            keys['n_2'] = np.interp(keys['vacuum_wavelength'], wl_Al_data, n_Al_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Al_data, n_Al_imag)
            keys['n_4'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_data, n_Aminoacid_imag)
            keys['kappa'] = np.interp(keys['vacuum_wavelength'], wl_Aminoacid_kappa, kappa_Aminoacid_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Aminoacid_kappa, kappa_Aminoacid_imag)
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
            n2 = np.nonzero(table['DomainId'] == 2) # row number for Material Id 2 (Al Pillar)
            n3 = np.nonzero(table['DomainId'] == 3) # row number for Material Id 3 (Al Film)
            n4 = np.nonzero(table['DomainId'] == 4) # row number for Material Id 4 (Amino Acid)

            # Power from plane wave
            #p_in = (0.5)*(Eamp**2/Z1)*(((3*np.sqrt(3))/2)*((keys['pitch']*keys['uol'])**2)))*math.cos(math.radians(keys['AOI'])) # ((Eamp V/m)/ohms)*m^2 = Watts
            #hexagon area
            #p_in = 1.0 * (((3*np.sqrt(3))/2)*((keys['pitch']*keys['uol'])**2))*math.cos(math.radians(keys['AOI'])) # Eamp is not used because PowerScaling is being used in sources.jcmt

            #square
            #p_in = 1.0 * (((keys['pitch']*keys['uol'])**2))*math.cos(math.radians(keys['AOI'])) # Eamp is not used because PowerScaling is being used in sources.jcmt
            p_in = 1.0 * (((keys['pitch']*keys['uol'])**2)) # Eamp is not used because PowerScaling is being used in sources.jcmt

            absS_pillar = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][0][n2]))/p_in)
            absP_pillar = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][1][n2]))/p_in)
         
            absS_film = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][0][n3]))/p_in)
            absP_film = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][1][n3]))/p_in)
            
            absS_amino = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][0][n4]))/p_in)
            absP_amino = (np.real(np.sum(table['ElectromagneticFieldAbsorption'][1][n4]))/p_in)

            # Gather Reflected Fourier Modes (Z)
            filename_fourierModes_r = './project_results/fourier_modes_r.jcm';
            fourierModes_r = jcmwave.loadtable(filename_fourierModes_r,format='named')
            powerFlux_r = jcmwave.convert2powerflux(fourierModes_r)
            
            reflection_list = []
            total_P_s_r = np.real(powerFlux_r['PowerFluxDensity'][0][:,2])
            total_P_p_r = np.real(powerFlux_r['PowerFluxDensity'][1][:,2])
            # Reflected flux in normal direction
            for order in range(len(powerFlux_r['PowerFluxDensity'][0])):
                orders_list = [str(fourierModes_r['N1'][order]), str(fourierModes_r['N2'][order])]
                P_s_r = total_P_s_r[order]
                P_p_r = total_P_p_r[order]
                reflection_list.append([orders_list, [P_s_r, P_p_r]])

            filename_MM = './project_results/sm.jcm'
            table = jcmwave.loadtable(filename_MM)
            mm_orders = []	
            for num in range(len(table['Mueller_xy11'][0])):
                mm = []
                for i in range(4):
                    row = []
                    for j in range(4):
                        row.append(float(table['Mueller_xy'+str(1+i)+str(1+j)][0][num]))
                    mm.append(row)
                Kout = table['KOut'][num]
                thetaOut = table['ThetaOut'][num]
                phiOut = table['PhiOut'][num]
  
                order_list = [[[float(Kout[0].real), float(Kout[0].imag)], [float(Kout[1].real), float(Kout[1].imag)],[float(Kout[2].real), float(Kout[2].imag)]], [float(thetaOut), float(thetaOut)], [float(phiOut), float(phiOut)]]
                mm_orders.append([order_list, mm])

            id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']

            set.append(Tag(id_names, [azi, aoi, wvl, keys['radius'], keys['pitch'], 50, kappa]), ["abs pillar", "abs film", "abs amino", "reflected_diff_orders", "reflected flux", "mm_diff_orders", "mm"], [[absS_pillar, absP_pillar], [absS_film, absP_film], [absS_amino, absP_amino], len(reflection_list), reflection_list, len(mm_orders), mm_orders])
            #set.append(Tag(id_names, [azi, aoi, wvl, keys['radius'], keys['pitch'], 50, kappa]), ["abs film", "abs air", "reflected flux", "mm"], [[absR_film, absL_film], [absR_air, absL_air], [P_s_r, P_p_r], mm])
            set.save_json(output_file)
	    
        toc = time.time() # use time() not clock() on linux system  
        t = toc-tic
        print ("Total runtime for WVL: "+str(wvl)+" AOI: "+str(aoi)+" aziumth: "+str(azi)+" - %6.4f s" % t)

