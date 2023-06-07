import json_format as json
import numpy as np
from data_format import *
from format_module import *
import matplotlib.pyplot as plt


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



def azi_X_aoi_Y(set, wvl, kappa, diff_order):
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
                    print(data_tag.hashable)
                    print(set.data)
                    sub = str(element)
                    print(sub)
                    value = set.get(data_tag.hashable).get("mm_information").get("dMM")[diff_order][int(sub[0])-1][int(sub[1])-1]
                    Z[mmi][yi].append(value) 
                except Exception as e:
                    print(e)
                    print(str(xv) + " " + str(yv) + " " + str(element))
                    exit()
                    Z[mmi][yi].append(None)
                #Z[mmi][yi].append(kappa_value)
    return X, Y, Z

def azi_X_wvl_Y(set, aoi, kappa, diff_order):
    diff_order_map = {('0', '0'): 0, ('1', '0'): 1, ('1', '-1'): 2, ('0', '-1'): 3, ('0', '1'): 4, ('1', '1'): 5}
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
    radius = 50
    pitch = 175
    height = 50
    #remove last unfinished element
    #len of X and Y is dependent on len of y
    x, y = getXY(set, "azimuth", "wvl")
    
    base_x,base_y = np.meshgrid(np.radians(x),y)
    elements = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]
    Z=[[ [] for i in range(len(y))] for e in elements]
    #Z = is len of elements, len of y
    order_map = dict()
    for yi, yv in enumerate(y):
        for xi, xv in enumerate(x):
            for mmi, element in enumerate(elements):
                try:
                    data_tag = Tag(id_names, [xv, aoi, yv, radius, pitch, height, kappa])
                    #noKdata_tag = Tag(id_names, [xv, 56, yv, radius, pitch, height, False])
                    sub = str(element)
                    data_point = set.get(data_tag.hashable)
                    for index, order in enumerate(data_point.get("reflected flux")):
                        order_tuple = (order[0][0], order[0][1])
                        order_map[order_tuple] = None 
                        if diff_order_map.get(order_tuple) == diff_order: 
                            dmm_element = data_point.get("mm_information")[index].get("dmm")[int(sub[0])-1][int(sub[1])-1]
                            break
                        else:
                            dmm_element = 1234

                    #nokappa_value = set.get(noKdata_tag.hashable).get("dMM")[int(sub[0])-1][int(sub[1])-1]
                    #dif = (kappa_value-nokappa_value)/kappa_value 
                    #Z[mmi][yi].append(np.log10(np.absolute(dif)))
                    Z[mmi][yi].append(dmm_element)
                except:
                    #print(str(xv) + " " + str(yv) + " " + str(element))
                    Z[mmi][yi].append(1234)
                #Z[mmi][yi].append(kappa_value)
    X = [base_x for e in elements]
    Y = [base_y for e in elements]
    print(order_map)
    return X, Y, Z

def azi_X_wvl_Y_reflect_sep_diff_orders(set, aoi, kappa, diff_order):
    diff_order_map = {('0', '0'): 0, ('1', '0'): 1, ('1', '-1'): 2, ('0', '-1'): 3, ('0', '1'): 4, ('1', '1'): 5}
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
    radius = 50
    pitch = 175
    height = 50
    #remove last unfinished element
    #len of X and Y is dependent on len of y
    x, y = getXY(set, "azimuth", "wvl")
    
    base_x,base_y = np.meshgrid(np.radians(x),y)
    Z=[ [] for i in range(len(y))]
    #Z = is len of elements, len of y
    order_map = dict()
    for yi, yv in enumerate(y):
        for xi, xv in enumerate(x):
            data_tag = Tag(id_names, [xv, aoi, yv, radius, pitch, height, kappa])
            data_point = set.get(data_tag.hashable)
            for index, order in enumerate(data_point.get("reflected flux")):
                order_tuple = (order[0][0], order[0][1])
                order_map[order_tuple] = None 
                if diff_order_map.get(order_tuple) == diff_order: 
                    DI = data_point.get("reflected flux")[index][1][0]
                    break
                else:
                    DI = 1234
            Z[yi].append(DI)
    X = base_x
    Y = base_y
    #print(order_map)
    return X, Y, Z, order_map

def azi_X_wvl_Y_reflect(set, aoi, kappa):
    diff_order_map = {('0', '0'): 0, ('1', '0'): 1, ('1', '-1'): 2, ('0', '-1'): 3, ('0', '1'): 4, ('1', '1'): 5}
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
    radius = 50
    pitch = 175
    height = 50
    #remove last unfinished element
    #len of X and Y is dependent on len of y
    x, y = getXY(set, "azimuth", "wvl")
    
    base_x,base_y = np.meshgrid(np.radians(x),y)
    Z=[ [] for i in range(len(y))]
    #Z = is len of elements, len of y
    order_map = dict()
    for yi, yv in enumerate(y):
        for xi, xv in enumerate(x):
            data_tag = Tag(id_names, [xv, aoi, yv, radius, pitch, height, kappa])
            data_point = set.get(data_tag.hashable)
            total_reflect = 0
            for index, order in enumerate(data_point.get("reflected flux")):
                order_tuple = (order[0][0], order[0][1])
                order_map[order_tuple] = None 
                total_reflect = total_reflect + data_point.get("reflected flux")[index][1][0]
            Z[yi].append(total_reflect)
    X = base_x
    Y = base_y
    #print(order_map)
    return X, Y, Z

def azi_X_wvl_Y_DI(set, aoi, kappa, diff_order):
    diff_order_map = {('0', '0'): 0, ('1', '0'): 1, ('1', '-1'): 2, ('0', '-1'): 3, ('0', '1'): 4, ('1', '1'): 5}
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']
    radius = 50
    pitch = 175
    height = 50
    #remove last unfinished element
    #len of X and Y is dependent on len of y
    x, y = getXY(set, "azimuth", "wvl")
    
    base_x,base_y = np.meshgrid(np.radians(x),y)
    Z=[ [] for i in range(len(y))]
    #Z = is len of elements, len of y
    order_map = dict()
    for yi, yv in enumerate(y):
        for xi, xv in enumerate(x):
            data_tag = Tag(id_names, [xv, aoi, yv, radius, pitch, height, kappa])
            data_point = set.get(data_tag.hashable)
            for index, order in enumerate(data_point.get("reflected flux")):
                order_tuple = (order[0][0], order[0][1])
                order_map[order_tuple] = None 
                if diff_order_map.get(order_tuple) == diff_order: 
                    DI = data_point.get("mm_information")[index].get("DI")
                    break
                else:
                    DI = 1234
            Z[yi].append(DI)
    X = base_x
    Y = base_y
    print(order_map)
    return X, Y, Z

def get_sum(set):
    i = 0
    x = []
    y = []
    z = []
    sum = []
    orders_list = []
    azi_list = []
    aoi_list = []
    wvl_list = []
    for hashable in set.data.keys():
        s_p = 0
        pillar  = set.get(hashable).get("abs pillar")[s_p]
        film  = set.get(hashable).get("abs film")[s_p]
        amino  = set.get(hashable).get("abs amino")[s_p]
        reflect = np.sum([i[1][s_p] for i in set.get(hashable).get("reflected flux")])
        sum.append(pillar+film+amino+reflect)
        orders_list.append(set.get(hashable).get("reflected_diff_orders"))
        azi_list.append(set.get(hashable).get("azimuth"))
        aoi_list.append(set.get(hashable).get("AOI"))
        wvl_list.append(set.get(hashable).get("wvl"))
        i=i+1
    #plt.scatter(azi_list, wvl_list, c=orders_list)
    plt.scatter(wvl_list, sum, c=azi_list)
    plt.title("1")
    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    #ax.scatter(azi_list, wvl_list, aoi_list, c=sum)
    #fig.colorbar(ax.collections[0])
    plt.show()
    #plt.scatter(azi_list, wvl_list, c=sum)
    #plt.show()

def get_reflect(set, aoi, azi):
    diff_order_map = {('0', '0'): 0, ('1', '0'): 1, ('1', '-1'): 2, ('0', '-1'): 3, ('0', '1'): 4, ('1', '1'): 5}
    i = 0
    x = []
    y = []
    z = []
    sum = []
    orders_list = []
    azi_list = []
    aoi_list = []
    wvl_list = []
    s_p = 0
    aoi = aoi
    azi = azi
    radius = 50
    pitch = 175
    height = 50
    kappa = False
    id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']

    complete_orders = [[] for i in range(6)]
    complete_wvl = [[] for i in range(6)]

    len_list = dict()

    for element in set.get_elements():
        len_list[int(element.get("wvl"))] = None
    len_list.popitem()
    diff_map = dict()
    for wvl in len_list.keys():
        data_tag = Tag(id_names, [azi, aoi, wvl, radius, pitch, height, kappa])
        data_point = set.get(data_tag.hashable).get("reflected flux")
        for order, value in enumerate(data_point):
            diff_tuple = (value[0][0], value[0][1])
            diff_map[diff_tuple] = None
            index = diff_order_map.get(diff_tuple)
            complete_orders[index].append(value[1][s_p])
            complete_wvl[index].append(wvl)
    print(diff_map)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'] 
    sum = []
    orders_list = []
    azi_list = []
    aoi_list = []
    wvl_list = []
    for hashable in set.data.keys():
        s_p = 0
        pillar  = set.get(hashable).get("abs pillar")[s_p]
        film  = set.get(hashable).get("abs film")[s_p]
        amino  = set.get(hashable).get("abs amino")[s_p]
        reflect = np.sum([i[1][s_p] for i in set.get(hashable).get("reflected flux")])
        sum.append(pillar+film+amino+reflect)
        orders_list.append(set.get(hashable).get("reflected_diff_orders"))
        azi_list.append(set.get(hashable).get("azimuth"))
        aoi_list.append(set.get(hashable).get("AOI"))
        wvl_list.append(set.get(hashable).get("wvl"))
        i=i+1
    for i in range(len(complete_wvl)):
        print(len(complete_orders[2]))
        plt.plot(complete_wvl[i], complete_orders[i], colors[i]+'o')
    plt.legend(diff_order_map.keys())
    plt.title(str(aoi) + " AOI; " + str(azi) + " azimuth")
    print(len(complete_wvl[1]))
    print(len(complete_orders[1]))

    #plt.scatter(complete_wvl[1], complete_orders[1], )
    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d') 
    #ax.scatter(azi_list, wvl_list, aoi_list, c=sum)
    #fig.colorbar(ax.collections[0])
    plt.show()



def get_sep(set):
    i = 0
    index = []
    store = []
    inner = []
    total = []
    for hashable in set.data.keys():
        pillar  = set.get(hashable).get("abs pillar")[0]
        film  = set.get(hashable).get("abs film")[0]
        amino  = set.get(hashable).get("abs amino")[0]
        reflect = np.sum(i[1][0] for i in set.get(hashable).get("reflected flux"))
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




def list_toDict(set):
    for key in set.data.keys():
        full_order_list = []
        for order in set.data.get(key).get("mm"):
            dimen = order[0]
            mm = order[1] 
            dmm = cloude_decomp(mm)
            map = dict()
            map["OutX"] = dimen[0][0][0]
            map["OutY"] = dimen[0][1][0]
            map["OutZ"] = dimen[0][2][0]
            map["OutTheta"] = dimen[1][0]
            map["OutPhi"] = dimen[2][0]
            order_map = dict()
            order_map["dimensions"] = map
            order_map["mm"] = mm
            order_map["dmm"] = dmm.tolist()
            full_order_list.append(order_map)
        set.data[key]["mm_information"] = full_order_list
        set.data[key].pop("mm")
    set.save_json("proper_60AOI.json")



id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']

file = 'proper_60AOI_DI.json'
set = Structure.from_json(file)
aoi = 60
wvl = 230

X, Y, Z = azi_X_wvl_Y(set, aoi, False, 0)
complete_MM_heatmap_plot(X,Y,Z, '', 0)

polar_plot(X,Y,Z)
exit()

#exit()
#get_reflect(set, aoi, 40)
#for i in range(24):
    #get_reflect(set, aoi,i*2)
#exit()
X, Y, kZ = azi_X_wvl_Y_reflect(set, aoi, False)
DI_heatmap_plot(['azi',X],['wvl',Y],kZ,"AOI: " + str(aoi) + ", Combined Reflect")

for i in range(3):
    #X, Y, kZ = azi_X_wvl_Y(set, aoi, False,i)
    #complete_MM_heatmap_plot(['azi', X], ['wvl', Y], kZ, ['aoi', aoi],i)
    #X, Y, kZ = azi_X_wvl_Y(set, aoi, False,i)
    X, Y, kZ, order = azi_X_wvl_Y_reflect_sep_diff_orders(set, aoi, False,i)
    DI_heatmap_plot(['azi',X],['wvl',Y],kZ,"AOI: " + str(aoi) + ", Diff Order: " + str(i))



#aoi = 40
#get_sum(set)
#for i in range(8):
    #X, Y, kZ = azi_X_wvl_Y(set, aoi, False, i)
    #complete_MM_heatmap_plot(['azi', X], ['wvl', Y], kZ, ['aoi', aoi], i)



    
#set.add_set(Structure.from_json(second))
#print(set.get("_0_20_224_50_175_50_False"))
#X, Y, kZ = azi_X_aoi_Y(set, wvl, True)
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


