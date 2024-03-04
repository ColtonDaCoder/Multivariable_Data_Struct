import pandas as pd
import numpy as np
import ast
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def complete_MM_heatmap_plot(X, Y, x_name, y_name, kZ, no_kZ=None):
    fig, ax = plt.subplots(4,4, figsize=(10,8))
    X[0] = [i*(180/np.pi) for i in X[1]]
    for j in range(4):
        for i in range(4):
            val = i*4+j
            mm = str(str(j+1)+str(i+1))
            if not no_kZ == None: 
                Z = [[abs(kZ[val][i][j] - no_kZ[val][i][j])/no_kZ[val][i][j]*100 for j, azi in enumerate(kZ[val][i])] for i, aoi in enumerate(kZ[val])]  
            else:
                Z = [[kZ[val][i][j] for j, azi in enumerate(kZ[val][i])] for i, aoi in enumerate(kZ[val])]  
            Z = np.absolute(Z)
            #Z = np.log10(Z)
            Z = np.array(Z)
            norm_list = Z.flatten() 
            norm_list = np.delete(norm_list, np.where(norm_list == 1234)) 
            std = np.std(norm_list)
            mean = np.mean(norm_list)
            z_max = np.amax(norm_list)
            z_min = np.amin(norm_list)
            if abs(z_max) > abs(z_min):
                z_min = -z_max
            else:
                z_max = -z_min
            z_max = 3* std
            z_min = 0
            Z[Z == 1234] = z_min
            c = ax[j,i].pcolormesh(X[0], Y[0], Z, cmap=plt.cm.Reds, vmin=z_min, vmax=z_max)

            fmt = ticker.FormatStrFormatter("%.0f%%")
            cbar = fig.colorbar(c, ax=ax[i,j], format=fmt) 
            
            #for new pandas df
            ax[j,i].set_xlabel(x_name, fontsize=10)
            ax[j,i].set_ylabel(y_name, fontsize=10)

    plt.tight_layout(h_pad=1,w_pad=0.5)
    plt.show()

def plot2d(X, Y_list):
    colors = ['b','g','r','c','m','y']
    for i, Y in enumerate(Y_list):
        plt.scatter(X, Y, color=colors[i])
    plt.show()

def plot3d(X, Y, Z):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    colors = ['b','g','r','c','m','y']
    ax.scatter(X, Y, Z, color='r')
    max = np.amax(Y)
    min = np.amin(Y)
    ticks = np.arange(max, min, (max - min)/5)
    ax.set_yticks(ticks)
    plt.show()



def MM_heatmap_plot(X, Y, x_name, y_name, kZ):
    fig, ax = plt.subplots()
    X[0] = [i*(180/np.pi) for i in X[1]]

    Z = [[kZ[i][j] for j, azi in enumerate(kZ[i])] for i, aoi in enumerate(kZ)]  
    #Z = np.absolute(Z)
    #Z = np.log10(Z)
    Z = np.array(Z)
    norm_list = Z.flatten() 
    norm_list = np.delete(norm_list, np.where(norm_list == 1234)) 
    z_max = np.amax(norm_list)
    z_min = np.amin(norm_list)

    #Z[Z == 1234] = z_min
    c = ax.pcolormesh(X[0], Y[0], Z, cmap=plt.cm.Reds, vmin=z_min, vmax=z_max)

    fmt = ticker.FormatStrFormatter("%.2f")
    cbar = fig.colorbar(c, ax=ax, format=fmt) 
    
    #for new pandas df
    ax.set_xlabel(x_name, fontsize=10)
    ax.set_ylabel(y_name, fontsize=10)

    plt.tight_layout(h_pad=1,w_pad=0.5)
    plt.show()

#cloude decomposition matrix
def cloude_decomp(mm):
    #pauli matrices
    p = [
            [
                [1,0],
                [0,1]
            ],
            [
                [1,0],
                [0,-1]
            ],
            [
                [0,1],
                [1,0]
            ],
            [
                [0,-1j],
                [1j,0]
            ]
        ]
    A = np.array([np.conj(i).flatten() for i in p])
    #coherency matrix = (1/4) A(pauli(sub i) kronecker pauli(sub j))A^-1
    sigma = 0
    for i in range(4):
        for j in range(4):
            sigma = sigma + mm[i][j]*(A@np.kron(p[i],np.conj(p[j]))@np.linalg.inv(A))
    C = (1/4) * sigma

    w,v = np.linalg.eigh(C)
    Ψ1=np.transpose(v[:, 3])
    j11=Ψ1[0] + Ψ1[1]
    j12=Ψ1[2] - 1j*Ψ1[3]
    j21=Ψ1[2] + 1j*Ψ1[3]
    j22=Ψ1[0]-Ψ1[1]
    J=np.array([[j11,j12],[j21,j22]])
    A = (1/np.sqrt(2))*np.array([[1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1j, -1j, 0]])
    Ainv = np.linalg.inv(A)
    Jconj=np.matrix.conjugate(J)
    M = A@(np.kron(J,Jconj))@Ainv
    M = M.real
    return M 

def get_dmm(df, x, y, x_name, y_name):
    # x and y are lists of azi and wvl respectively 
    base_x,base_y = np.meshgrid(np.radians(x),y)
    elements = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]
    Z=[[ [] for i in range(len(y))] for e in elements]
    #Z = is len of elements, len of y
    for yi, yv in enumerate(y):
        for xi, xv in enumerate(x):
            mm = ast.literal_eval(df.loc[(df[x_name] == xv) & (df[y_name] == yv)].iloc[0]['mm'])
            dmm = cloude_decomp(mm)
            for mmi, element in enumerate(elements):
                sub = str(element)
                dmm_element = dmm[int(sub[0])-1][int(sub[1])-1]
                Z[mmi][yi].append(dmm_element)
                #Z[mmi][yi].append(1234)
                #Z[mmi][yi].append(kappa_value)
    X = [base_x for e in elements]
    Y = [base_y for e in elements]
    return X, Y, Z

def get_E(df, x, y, x_name, y_name, E_name):
    # x and y are lists of azi and wvl respectively 
    base_x,base_y = np.meshgrid(np.radians(x),y)
    elements = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]
    Z=[ [] for i in range(len(y))]
    #Z = is len of elements, len of y
    for yi, yv in enumerate(y):
        for xi, xv in enumerate(x):
            E = ast.literal_eval(df.loc[(df[x_name] == xv) & (df[y_name] == yv)].iloc[0][E_name])[1]
            Z[yi].append(E)
            #for mmi, element in enumerate(elements):
                #Z[mmi][yi].append(E)
    X = [base_x for e in elements]
    Y = [base_y for e in elements]
    return X, Y, Z

def getXY(df, x_name, y_name):
    X = dict()
    Y = dict()
    for element in df[x_name].values:
        X[int(element)] = None
    for element in df[y_name].values:
        Y[int(element)] = None
    return list(X.keys()), list(Y.keys())


file = "MM_Spectrum_Ag_TS_RH_gammadion_per_AOI20_pitch840_arm120_t50/Au_on_Au_Gammadion_fe2_MSL500.csv"

df = pd.read_csv(file)
x, y = getXY(df, 'azi', 'wvl')
X = []
Y = []
Y2 = []
X, Y, Z  = get_dmm(df, x, y, 'azi', 'wvl')
plt.plot([j[0] for j in Y[3]],[i[0] for i in Z[3]], 'bo')
plt.plot([j[0] for j in Y[12]],[i[0] for i in Z[12]], 'ro')
plt.xlabel("wvl (nm)")
plt.ylabel("cloude decomp MM (blue: 14; red: 41")
plt.show()
exit()
for index, row in df.iterrows():
    X.append(row['wvl'])
    print(ast.literal_eval(row['dmm']))
    #print(ast.literal_eval(row['dmm'])[3][0])
    exit()
    
plt.plot(X,Y, 'ro')
plt.plot(X,Y2, 'bo')
plt.title("AOI: 60, AZI: 22, 450 thick Au Gammadion on Au")
plt.xlabel("wvl (nm)")
plt.ylabel("reflected flux w/o PowerFluxScaling")
plt.show()
exit()
print(df)
exit()
for val in y:
    df
print(x,y)

exit()
X, Y, Z = get_dmm(df, x, y, 'azi','wvl')



file = "Far_field_MIR_60AOI_CaF2/high_far_field_60AOI_chiral.csv"

df = pd.read_csv(file)
x, y = getXY(df, 'azi', 'wvl')
x = x[1:]
y=y[1:]
#y=y[18:-10]
X, Y, Z = get_dmm(df, x, y, 'azi','wvl')

#complete_MM_heatmap_plot(X, Y, 'azi', 'wvl', Z)


file = "Far_field_MIR_60AOI_CaF2/high_far_field_60AOI_racemic.csv"

df = pd.read_csv(file)
x, y = getXY(df, 'azi', 'wvl')
#y=y[18:-10]
x = x[1:]

X, Y, no_kZ = get_dmm(df, x, y, 'azi','wvl')
complete_MM_heatmap_plot(X, Y, 'azi', 'wvl', Z, no_kZ)
exit()

#X_, Y_, Z = get_dmm(df, x, y, 'azi','wvl')
#complete_MM_heatmap_plot(X, Y, 'azi', 'wvl', Z)

file = "Far_field_MIR_60AOI_CaF2_2/high_far_field_60AOI_racemic_2.csv"

df2 = pd.read_csv(file)

file = "Far_field_MIR_60AOI_CaF2/high_far_field_60AOI_racemic_1_5.csv"

df3 = pd.read_csv(file)

df = pd.concat([df, df2, df3])

x, y = getXY(df, 'azi', 'wvl')
y = y[18:-11]
x = x[1:]

X, Y, no_kZ = get_dmm(df, x, y, 'azi','wvl')
complete_MM_heatmap_plot(X, Y, 'azi', 'wvl', Z, no_kZ)

exit()
#df = df.loc[(df["azi"] == 5)]
#Y0 = [ast.literal_eval(i)[0] for i in df["abs pillar"].values]
#Y1 = [ast.literal_eval(i)[1] for i in df["abs pillar"].values]
#plot3d(df["wvl"].values, Y0, df["azi"].values)

#kappa = np.loadtxt('LH_Gammadion_Simulation_Racemic_Trp_500/data/L-tyr_kappa.txt')
kappa = np.loadtxt('LH_Gammadion_Simulation_Racemic_Trp_500/data/L-tyr_nx_kx.txt')

print(kappa[0])
print(kappa[1])
plt.plot([i[0] for i in kappa],[i[1] for i in kappa])
#plt.plot([i[0] for i in kappa],[i[2] for i in kappa])
#plt.plot(kappa[0],kappa[2])
plt.show()



#X,Y,e = get_E(df, x, y, 'azi','wvl', 'abs pillar')
#MM_heatmap_plot(X, Y, 'azi', 'wvl', e)
