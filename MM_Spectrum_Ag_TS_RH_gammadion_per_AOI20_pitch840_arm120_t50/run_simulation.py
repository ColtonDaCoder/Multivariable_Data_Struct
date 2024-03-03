import numpy as np
import jcmwave,time,imp,shutil,os 
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--threads",action="store",type="int", dest="threads",help="number of threads to use")
(options, args) = parser.parse_args()
jcmwave.set_num_threads(options.threads)
tic = time.time() # use time() not clock() on linux system
azi = []
azi = np.arange(0,182,2) # The azimuth range you want to sweep over
azimuth = []
for i in range(len(azi)):
    angle = int(azi[i])
    results = [] # create empty array to hold results for each azi
    print('Azimuth: %3.0f deg' % (angle))
        
    keys = {} # Create empty dictionary for keys
    results = [] # Create empty array for results

    # Set simulation parameters
    keys = {
        'dielectric': 10, # This was the variable "spacer" before, but it isn't being used in the layout file
        'AOI':60, # This is the angle of incidence
        'gap':42,
        'uol': 1e-9,
        'display_triangulation' : 'no',
        'boundary' : 'Periodic', # Periodic or Transparent
        'info_level' : -1,
        'pitch': 840,
        'fem_degree' : 2,
        'arm_thickness' : 120,
        'angle': angle, # This rotates the layout, so we will use it to change the azimuth
        'n_refinement_steps' : 0, # Currently we get non-physical results if this is >0
        'thickness' : 50, #thickness of the gammadion
        'film': 100, # Thickness of the flat film
        'MSL':100, #MaximumSideLength for the gammadion mesh (not the z mesh)
        'z_radius' : 5, # radius of curvature of gammadion in the z plane
        'z_radius_MSL' : 1 # maximum side length of z radius
        }

    keys['delta_k'] = (2*np.pi)/(20*keys['pitch']*1e-9) # For the discritization of the continuous Fourier spectrum into modes (only for isolated setup)
    tag_ = 'FE2_MSL100'
    #tag_ = '1nm_step_updated_'+str(keys['AOI'])+'AOI_RH_5nm_z-rounded_Tiago_Al_New_Mesh_MM_SP_Spectr_TS_Al_Gamm_' + str(keys['pitch']) + '_nm_pitch_' + str(keys['boundary']) + '_boundary' + 'NA_0.09_FE_' + str(keys['fem_degree']) + '_arm_' + str(keys['arm_thickness'])+ '_gap_' + str(keys['gap']) + '_50nm_mesh_' + str(keys['film'])+'nm_Al_sub_150nm_Air'

    # material properties
    keys['n_3'] = 1.00 # index refraction of Air
        
    Au_nk = np.loadtxt('../data/Au_Babar_micron.nk') # You would need to change this to the silver data file I will send you
    wl_Au_data = []; n_Au_real = []; n_Au_imag = [] ### ALUMINUM NOT GOLD!
    for data in Au_nk:
        wl_Au_data.append(data[0]*1e-6) # e-10 for [ang], e-9 for [nm], e-6 for [um]
        n_Au_real.append(data[1])
        n_Au_imag.append(data[2])
    
    lambdas = np.linspace(210, 400, 20)*1e-9
    n = len(lambdas)
    # Create empty arrays for the Mueller matrix data, the length of the arrays is equal to number of wavelengths we are scanning
    M11=[]
    M12=[]
    M13=[]
    M14=[]
    M21=[]
    M22=[]
    M23=[]
    M24=[]
    M31=[]
    M32=[]
    M33=[]
    M34=[]
    M41=[]
    M42=[]
    M43=[]
    M44=[]
    
    for keys['vacuum_wavelength'] in lambdas:
        print('Wavelength : %3.2f nm' % (keys['vacuum_wavelength']*1e9))
        lambda_int = int(keys['vacuum_wavelength']*1e9)
        keys['n_1'] = keys['n_3'] # index refraction of air
        keys['n_2'] = np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_imag) #Al not Au!
        keys['n_4'] = np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_real) + 1j*np.interp(keys['vacuum_wavelength'], wl_Au_data, n_Au_imag) #Al not Au!
        keys['sm_filename'] = '"'+'project_results/sm_%dnm.jcm"' % int(keys['vacuum_wavelength']*1e9)   #for saving MM file
        jcmwave.jcmt2jcm('./boundary_conditions.jcmt', keys)
        jcmwave.jcmt2jcm('./materials.jcmt', keys)
        jcmwave.jcmt2jcm('./project.jcmpt', keys)
        jcmwave.jcmt2jcm('./sources.jcmt', keys)
        jcmwave.jcmt2jcm('./layout.jcmt', keys)
        jcmwave.solve('./project.jcmp')

        # Constants
        #c0 = 299792458
        #mu0 = 4*np.pi*1e-7
        #eps0 = 1/(mu0*c0**2)
        #omega = 2*np.pi*c0/keys['vacuum_wavelength']
        #Z1 = np.sqrt(mu0/(eps0*keys['n_3']**2))
        
        # Power from plane wave
        #p_in = 2*(1/np.sqrt(2))**2/Z1*(keys['pitch']*keys['uol'])**2 ##2 in front because of amplitude
        #tot_p_in = p_in*2 # sum x and y polarization

        #Extract absorption data: (ONLY FOR CPL)
        #fileName = './project_results/electric_field_energy.jcm'
        #table = jcmwave.loadtable(fileName)
        #n4 = np.nonzero(table['DomainId'] == 4) # row number for Material Id 4 (Al Nanostructure)
        #n2 = np.nonzero(table['DomainId'] == 2) # row number for Material Id 2 (Al film)
        
        #Calculate Absorption (ONLY FOR CPL)
        #absR_Struc = -2*omega*np.imag(table['ElectricFieldEnergy'][0][n4][0])
        #absR_Film = -2*omega*np.imag(table['ElectricFieldEnergy'][0][n2][0])
        #absL_Struc = -2*omega*np.imag(table['ElectricFieldEnergy'][1][n4][0])
        #absL_Film = -2*omega*np.imag(table['ElectricFieldEnergy'][1][n2][0])

        #fileName = './project_results/absorption.jcm'
        #table = jcmwave.loadtable(fileName)
        #absR_Struc_JCM = np.real(table['ElectromagneticFieldAbsorption'][0][2]) 
        #absL_Struc_JCM = np.real(table['ElectromagneticFieldAbsorption'][1][2])
        #absR_Film_JCM = np.real(table['ElectromagneticFieldAbsorption'][0][0])
        #absL_Film_JCM = np.real(table['ElectromagneticFieldAbsorption'][1][0])

        
        ## Gather Reflected Fourier Modes (+Z)
        filename_fourierModes_r = './project_results/fourier_modes_r.jcm';
        fourierModes_r = jcmwave.loadtable(filename_fourierModes_r,format='named')
        powerFlux_r = jcmwave.convert2powerflux(fourierModes_r)

        ## Reflected flux in normal direction
        P_s_r = np.sum(powerFlux_r['PowerFluxDensity'][0][:, 2]);
        P_p_r = np.sum(powerFlux_r['PowerFluxDensity'][1][:, 2]); 

        filename_MM = './project_results/sm_' + str(int(keys['vacuum_wavelength']*1e9)) + 'nm.jcm'

        table = jcmwave.loadtable(filename_MM)
        M11.append(sum(table['Mueller_xy11'][0])) # Need to sum the different K-vectors, in the isolated case since we have a continuous Fourier spectrum
        M12.append(sum(table['Mueller_xy12'][0]))
        M13.append(sum(table['Mueller_xy13'][0]))
        M14.append(sum(table['Mueller_xy14'][0]))
        
        M21.append(sum(table['Mueller_xy21'][0]))
        M22.append(sum(table['Mueller_xy22'][0]))
        M23.append(sum(table['Mueller_xy23'][0]))
        M24.append(sum(table['Mueller_xy24'][0]))
        
        M31.append(sum(table['Mueller_xy31'][0]))
        M32.append(sum(table['Mueller_xy32'][0]))
        M33.append(sum(table['Mueller_xy33'][0]))
        M34.append(sum(table['Mueller_xy34'][0]))
        
        M41.append(sum(table['Mueller_xy41'][0]))
        M42.append(sum(table['Mueller_xy42'][0]))
        M43.append(sum(table['Mueller_xy43'][0]))
        M44.append(sum(table['Mueller_xy44'][0]))
        # save data to file
        #results.append([keys['vacuum_wavelength'], keys['pitch'], P_s_r, P_p_r,m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44])
        #np.savetxt('./my_results_' + tag + '.txt', results, header='wvl[m], pitch, Reflec_Pol-1, Reflec_Pol-2,m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44')
        #results.append([keys['vacuum_wavelength'], keys['pitch'], P_s_r, P_p_r,m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44])
    
    results_filename = "RH_gammadion_Azi" + str(angle) + ".txt"
    lambda_txt = "\t".join(str(int(lmbda*1e9)) for lmbda in lambdas) # Create txt string of lambdas seperated by spaces
    # Create 2D array of data for the azi, that follows the Woollam format
    results = [['','AOI', lambda_txt],
      ["m11",str(keys['AOI']),"\t".join(map(str, M11))],
      ["m12",str(keys['AOI']),"\t".join(map(str, M12))],
      ["m13",str(keys['AOI']),"\t".join(map(str, M13))],
      ["m14",str(keys['AOI']),"\t".join(map(str, M14))],
      ["m21",str(keys['AOI']),"\t".join(map(str, M21))],
      ["m22",str(keys['AOI']),"\t".join(map(str, M22))],
      ["m23",str(keys['AOI']),"\t".join(map(str, M23))],
      ["m24",str(keys['AOI']),"\t".join(map(str, M24))],
      ["m31",str(keys['AOI']),"\t".join(map(str, M31))],
      ["m32",str(keys['AOI']),"\t".join(map(str, M32))],
      ["m33",str(keys['AOI']),"\t".join(map(str, M33))],
      ["m34",str(keys['AOI']),"\t".join(map(str, M34))],
      ["m41",str(keys['AOI']),"\t".join(map(str, M41))],
      ["m42",str(keys['AOI']),"\t".join(map(str, M42))],
      ["m43",str(keys['AOI']),"\t".join(map(str, M43))],
      ["m44",str(keys['AOI']),"\t".join(map(str, M44))]]
    #print(results);  
    np.savetxt('./'+results_filename, results, delimiter='\t', newline='\n', fmt='%s') # Set format to string
 
toc = time.time() # use time() not clock() on linux system  
t = toc-tic
print ("Total runtime for "+tag_+": %6.4f s" % t)
