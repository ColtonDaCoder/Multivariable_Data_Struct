import json
import matplotlib.pyplot as plt
import numpy as np
import copy

class json_file:

    def __init__(self, filename, is_new=False):
        self.filename = filename
        if is_new:
            self.data = {}
        else:
            self.data = self.get_json(self.filename) 

    #returns json file of wavelength file
    def get_json(self, filename):
        with open(filename, 'r') as read_file:
            raw = json.load(read_file)
        return raw
 
    #if array is not given return entire mueller matrix
    #else return element indexed by array parameter
    def getMM(self, aoi, azi, element=None):
        mm = self.get_data()[str(aoi)][str(azi)]["mm"]
        if element == None:
            return mm
        else:
            sub = str(element)
            row, column = int(sub[0]),int(sub[1])
            return mm[row-1][column-1]

    def getdMM(self, aoi, azi, element=None):
        dMM = self.get_data()[str(aoi)][str(azi)]["dMM"]
        if element == None:
            return dMM
        else:
            sub = str(element)
            row, column = int(sub[0]),int(sub[1])
            return dMM[row-1][column-1]



    #using tiago DI
    def add_DI(self):        
        new_data = self.get_data()
        for aoi in self.get_data():
            for azi in self.data[aoi]:
                mm = np.asarray(self.getMM(aoi,azi)).flatten()
                di = np.sqrt(abs(sum(i**2 for i in mm) - mm[0]**2 ))/(np.sqrt(3)*mm[0])
                new_data[aoi][azi]["DI"] = di
        self.save_json(self.filename,new_data)

    def get_data(self):
        return copy.deepcopy(self.data)

    def getMM_col(self, element):
        new_data = self.get_data()
        for aoi in new_data:
            for azi in new_data[aoi]:
                new_data[aoi][azi]["dMM"] = self.getdMM(aoi,azi,element)
        return new_data

    def add_DMM(self):
        new_data = self.get_data()
        for aoi in self.get_data():
            for azi in self.data[aoi]:
                mm = np.asarray(self.getMM(aoi,azi))
                new_data[aoi][azi]["dMM"] = self.cloude_decomp(mm).tolist()
        self.save_json(self.filename,new_data) 

    #cloude decomposition matrix
    def cloude_decomp(self, mm):
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

    def save_json(self, filename, data): 
        with open(filename, 'w') as write_file:
            write_file.write(json.dumps(data, indent=4))

def save_json(filename, data): 
    with open(filename, 'w') as write_file:
        write_file.write(json.dumps(data, indent=4))

def tiago_to_json(filename):
    store = {}
    for i in range(35-20):
        aoi = str(i+20) + "AOI" 
        file = filename + aoi + ".txt"
        raw = np.loadtxt(file)
        store[aoi] = {}
        for j in range(46): 
            azi = raw[j][3]
            store[aoi][azi] = {}
            for k in range(16):
                store[aoi][azi]["mm"] = [raw[j][6:10].tolist(), raw[j][10:14].tolist(), raw[j][14:18].tolist(), raw[j][18:22].tolist()]
    return store

def old_convert_to_plot(file):
    all_data = json_file(file).data
    #remove last unfinished element
    all_data.popitem()
    y = [int(aoi[0:2]) for aoi in all_data]
    x = [float(azi) for azi in all_data["20"]]
    #len of X and Y is dependent on len of y
    base_x,base_y = np.meshgrid(np.radians(x),y)
    X = [base_x for i in range(16)]
    Y = [base_y for i in range(16)] 
    elements = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]
    Z=[[ [] for aoi in all_data ] for e in elements]
    #Z = is len of elements, len of y
    for i, aoi in enumerate(all_data): 
        for azi in all_data[aoi]:
            dMM = all_data[aoi][azi]["dMM"]
            for index, element in enumerate(elements):
                sub = str(element)
                Z[index][i].append(dMM[int(sub[0])-1][int(sub[1])-1])
    return X,Y,Z

def single_polar_plot(X, Y, Z, DECOMP=False):
    rTicks = [219,222,] 
    xTicks = [0, np.pi/8,np.pi/6,np.pi/4, np.pi/3] 

    fig, axis = plt.subplots(figsize=(10,10),subplot_kw=dict(projection='polar'))
    cm = plt.cm.Reds
    i = 0
    j = 3
    val = i*4+j
    Zmin = np.amin(Z[val])
    Zmax = np.amax(Z[val])
    cb = axis.pcolormesh(X[val],Y[val],Z[val],vmin=Zmin, vmax=Zmax, antialiased=True, cmap = cm) 
    axis.set_rticks(rTicks, fontsize=10) #AOI
    axis.set_xticks(xTicks, fontsize=10)
    axis.set_rlim(rTicks[0], rTicks[-1])
    axis.set_xlim(xTicks[0], xTicks[-1]) #azimuth
    axis.grid( color = 'gray', linestyle = '--', linewidth = 1 )
    title = str(j+1) + str(i+1)
    #if DECOMP and (i == 0 and j == 0):
        #title = "DI"
    fig.colorbar(cb,ax=axis,pad=0.2)

    axis.set_title(title,fontsize=20)
    plt.tight_layout(h_pad=1,w_pad=3)
    #plt.savefig(str(wvl)+"polar.png")
    plt.show()

def polar_plot(X, Y, Z, DECOMP=False):
    #file = "fe3_nokappa/noK_fe3_"+raw+".json"
    #X,Y,Z = convert_to_plot(file)
    #X = azimuth (xticks)
    #Y = AOI (rticks)
    #Z = dMM
    rTicks = [219,222,] 
    xTicks = [0,np.pi/8,np.pi/6,np.pi/4] 

    fig, axis = plt.subplots(4,4,figsize=(10,10),subplot_kw=dict(projection='polar'))
    cm = plt.cm.Reds
    for j in range(4):
        for i in range(4):
            val = i*4+j
            Zmin = np.amin(Z[val])
            Zmax = np.amax(Z[val])
            print(Zmin)
            #if DECOMP and (i == 0 and j == 0):
                #Zmin = 0.9
                #Zmax = 1.1
            cb = axis[j,i].pcolormesh(X[val],Y[val],Z[val],vmin=Zmin, vmax=Zmax, antialiased=True, cmap = cm) 
            #cb = axis[j,i].pcolormesh(X[val],Y[val],Z[val],vmin=Zmin, vmax=Zmax, shading='gouraud', antialiased=True, cmap = cm) 
            axis[j,i].set_rticks(rTicks, fontsize=10) #AOI
            axis[j,i].set_xticks(xTicks, fontsize=10)
            axis[j,i].set_rlim(rTicks[0], rTicks[-1])
            axis[j,i].set_xlim(xTicks[0], xTicks[-1]) #azimuth
            axis[j,i].grid( color = 'gray', linestyle = '--', linewidth = 1 )
            title = str(j+1) + str(i+1)
            #if DECOMP and (i == 0 and j == 0):
                #title = "DI"
            fig.colorbar(cb,ax=axis[j,i],pad=0.2)

            axis[j,i].set_title(title,fontsize=20)
    plt.tight_layout(h_pad=1,w_pad=3)
    #plt.savefig(str(wvl)+"polar.png")
    plt.show()

def complete_MM_heatmap_plot(X, Y, no_kZ, raw, diff_order, kZ=0):
    fig, ax = plt.subplots(4,4, figsize=(10,8))
    X[1] = [i*(180/np.pi) for i in X[1]]
    for j in range(4):
        for i in range(4):
            val = i*4+j
            mm = str(str(j+1)+str(i+1))

            #toggle for difference or magnitude
            if not kZ == 0:
                Z = [[(kZ[val][i][j] - no_kZ[val][i][j])/kZ[val][i][j] for j, azi in enumerate(kZ[val][i])] for i, aoi in enumerate(kZ[val])]   
                Z = np.absolute(Z)
                #Z = np.log10(Z)
            else:
                Z = [[no_kZ[val][i][j] for j, azi in enumerate(no_kZ[val][i])] for i, aoi in enumerate(no_kZ[val])]  
                Z = np.absolute(Z)
            Z = np.array(Z)
            norm_list = Z.flatten() 
            norm_list = np.delete(norm_list, np.where(norm_list == 1234)) 
            z_max = np.amax(norm_list)
            #z_max = np.amax(Z)
            #z_max = 1
            z_min = np.amin(norm_list)
            Z[Z == 1234] = z_min
            #z_min = -0.02
            #z_max = 0.05
            c = ax[j,i].pcolormesh(X[1][0], Y[1][0], Z, cmap=plt.cm.Reds, vmin=z_min, vmax=z_max)

            cbar = fig.colorbar(c, ax=ax[i,j]) 
            #cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
            ax[j,i].set_xlabel(X[0], fontsize=10)
            ax[j,i].set_ylabel(Y[0], fontsize=10)

            #ax[j,i].set_title(str(raw[0]) + ": " + str(raw[1]) + " MM: " + str(mm), + "diff " + diff_order)
    plt.tight_layout(h_pad=1,w_pad=0.5)
    #plt.savefig("magnitude_heatmap_fe3_full_MM/"+raw+"/Complete_MM.png")
    #plt.savefig("log_heatmap_fe3_"+raw+".png")
    plt.show()


def contour_plot(x, y, no_kZ, wvl, mm, kZ=None):
    sub = str(mm)
    val = (int(sub[0])-1)*4 + (int(sub[1])-1)
    X = [i*(180/np.pi) for i in x[val]]
    Y = y[val]

    #percent difference
    #Z = [[(kZ[val][i][j] - no_kZ[val][i][j])/kZ[val][i][j] for j, azi in enumerate(kZ[val][i])] for i, aoi in enumerate(kZ[val])]  
    Z = [[no_kZ[val][i][j] for j, aoi in enumerate(no_kZ[val][i])] for i, aoi in enumerate(no_kZ[val])]  
    print(no_kZ[0][0][0])
    Z = np.absolute(Z)
    #Z = np.log10(Z)
    fig, ax = plt.subplots()
    z_max = np.amax(Z)
    #z_max = 1
    z_min = np.amin(Z)

    c = ax.pcolormesh(X, Y, Z, cmap=plt.cm.Reds, vmin=z_min, vmax=z_max)
    #c = ax.imshow(Z, cmap=plt.cm.Reds)
    ax.set_title("wvl: " + str(wvl) + " MM: " + str(mm))

    # set the limits of the plot to the limits of the data
    ax.axis([30, np.amax(X), np.amin(Y), np.amax(Y)])

    fig.colorbar(c, ax=ax)

    ax.set_xlabel('azimuth')
    ax.set_ylabel('AOI')
    plt.show()
    #plt.text(80,62, raw+" "+str(mm)+"log(|(chiral - achiral/chiral)|)")
    #plt.savefig("magnitude_heatmap_fe3_full_MM/"+raw+"/"+raw+"_"+str(mm)+".png")


