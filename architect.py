import json_format as json
import numpy as np
from data_format import *
from format_module import *
import matplotlib as plt


def getXY(set, x_name, y_name):
    X = dict()
    Y = dict()
    for element in set.get_elements():
        X[int(element.get(x_name))] = None
        Y[int(element.get(y_name))] = None

    return list(X.keys()), list(Y.keys())

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



def azi_X_aoi_Y(set, wvl, kappa):
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
    radius = 50
    pitch = 175
    height = 50
    #remove last unfinished element
    #len of X and Y is dependent on len of y
    x, y = getXY(set, "azimuth", "AOI")
    
    base_x,base_y = np.meshgrid(x,y)
    elements = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]
    X = [base_x for e in elements]
    Y = [base_y for e in elements]
    Z=[[ [] for i in range(len(y))] for e in elements]
    #Z = is len of elements, len of y
    for yi, yv in enumerate(y):
        for xv in x:
            for mmi, element in enumerate(elements):
                try:
                    data_tag = Tag(id_names, [xv, yv, wvl, radius, pitch, height, kappa])
                    sub = str(element)
                    value = set.get(data_tag.hashable).get("dMM")[int(sub[0])-1][int(sub[1])-1]
                    Z[mmi][yi].append(value) 
                except:
                    print(str(xv) + " " + str(yv) + " " + str(element))
                    exit()
                    Z[mmi][yi].append(None)
                #Z[mmi][yi].append(kappa_value)
    return X, Y, Z


def azi_X_wvl_Y(set):
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
    radius = 50
    pitch = 175
    height = 50
    kappa = True
    #remove last unfinished element
    #len of X and Y is dependent on len of y
    x, y = getXY(set, "azimuth", "wvl")
    
    base_x,base_y = np.meshgrid(np.radians(x),y)
    elements = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]
    X = [base_x for e in elements]
    Y = [base_y for e in elements]
    Z=[[ [] for i in range(len(y))] for e in elements]
    #Z = is len of elements, len of y
    for yi, yv in enumerate(y):
        for xv in x:
            for mmi, element in enumerate(elements):
                try:
                    data_tag = Tag(id_names, [xv, 48, yv, radius, pitch, height, kappa])
                    #noKdata_tag = Tag(id_names, [xv, 56, yv, radius, pitch, height, False])
                    sub = str(element)
                    kappa_value = set.get(data_tag.hashable).get("dMM")[int(sub[0])-1][int(sub[1])-1]
                    #nokappa_value = set.get(noKdata_tag.hashable).get("dMM")[int(sub[0])-1][int(sub[1])-1]
                    #dif = (kappa_value-nokappa_value)/kappa_value 
                    #Z[mmi][yi].append(np.log10(np.absolute(dif)))
                    Z[mmi][yi].append(kappa_value)
                except:
                    print(str(xv) + " " + str(yv) + " " + str(element))
                    exit()
                    print(data_tag.hashable)
                    Z[mmi][yi].append(None)
                #Z[mmi][yi].append(kappa_value)
    return X, Y, Z

def get_sum(set):
    i = 0
    index = []
    store = []
    for hashable in set.data.keys():
        if set.get(hashable).get('wvl') == 222:
            break;
        pillar  = set.get(hashable).get("abs pillar")[0]
        film  = set.get(hashable).get("abs film")[0]
        amino  = set.get(hashable).get("abs amino")[0]
        amino  = set.get(hashable).get("abs amino")[0]
        reflect = set.get(hashable).get("reflected flux")[0]
        store.append(pillar+film+amino+reflect)
        #index.append(i)
        index.append(set.get(hashable).get("AOI"))
        i=i+1
    plt.pyplot.plot(index, store, 'bo')
    plt.pyplot.show()


def get_sep(set):
    i = 0
    index = []
    store = []
    inner = []
    total = []
    for hashable in set.data.keys():
        if set.get(hashable).get('wvl') == 222:
            break;
        pillar  = set.get(hashable).get("abs pillar")[0]
        film  = set.get(hashable).get("abs film")[0]
        amino  = set.get(hashable).get("abs amino")[0]
        reflect = set.get(hashable).get("reflected flux")[0]
        inner.append(pillar+film+amino+reflect)
        total.append(pillar+film+amino+reflect)
        #inner.append(film+reflect)
        #total.append(film+reflect)
        index.append(i)
        i=i+1
        if i == 46:
            i = 0
            store.append(inner)
            inner = []
    i=0
    for b in store:
        plt.pyplot.plot([a for a in range(46)], b, 'bo')
        plt.pyplot.title(i*4)
        plt.pyplot.yticks([0.9, 0.92,0.94, 0.96, 0.98, 1.0])
        plt.pyplot.show()
        i= i+1







id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
wvl = 220
#file = 'hashed.json'
file = "flat_test.json"
file = "fe3_NA_0_9.json"
set = Structure.from_json(file)
#X, Y, kZ = azi_X_aoi_Y(set, wvl, True)
get_sep(set)
#complete_MM_heatmap_plot(X, Y, kZ, wvl)
#file = "fe3_220_222_NA_2_5.json"
#set = Structure.from_json(file)
#get_sum(set)
#X, Y, kZ = azi_X_aoi_Y(set, "220")
#X, Y, no_kZ = azi_X_aoi_Y(set, 222, False)
#single_polar_plot(X, Y, kZ)

#for hashable in set.keys():
    #set.get(hashable).get('mm')
#X, Y, Z = azi_X_aoi_Y(set)
#contour_plot(X, Y, Z, 41)

#X, Y, Z = convert_to_plot(set)
#complete_MM_heatmap_plot(X, Y, Z, '220')

#wvl = 230
#set = Structure.from_json("FCC_220_230.json")
#old_data = json.get_json("FCC_no_kappa/fe3_"+str(wvl)+".json")
#kappa = False
#pitch = 175
#radius = 50
#height = 50
#print(old_data.get("20").get("30").get("DI"))
#for aoi in old_data:
    #for azi in old_data.get(aoi):
        #values = old_data.get(aoi).get(azi)
        #set.append(Tag(id_names, [azi, aoi, wvl, radius, pitch, height, kappa]), ["mm", "DI", "dMM"], [values.get("mm"), values.get("DI"), values.get("dMM")])
#set.save_json("FCC_220_230.json")

#data_tag = Tag(id_names, [30,20,220,radius, pitch, height, False])
#print(data_tag.hashable)
#print(new_set.data.get(data_tag.hashable))
#print(new_set.get_elements()[0].get("AOI"))
#print(len(Z[0][0]))
#print(len(Y))
#print(len(X))
#polar_plot(X, Y, Z)


#for entry in set.get_elements():
    #id_values = [] 
    #value_names = ['mm','DI','dMM']
    #values = []
    #for key in entry.keys():
        #if not key == 'dMM' and not key == 'DI' and not key == 'mm': 
            #id_values.append(entry.get(key)) 
        #else:
            #values.append(entry.get(key))
    #data_tag = Tag(id_names, id_values)
    #new_set.append(data_tag, value_names, values) 
#new_set.save_json("FCC_220_230.json")


    
    #data_tag = Tag(id_names, entry.
    #new_set.append()
    






#constraints for displaying data
#EX: two dynamic and the rest static to plot x, and y
#set.add_constraint('azimuth', [30,90], False)
#set.add_constraint('AOI', [20,60], False)
#set.add_constraint('kappa', False, True)
#file = "fe3_220_230.json"


