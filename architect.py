import json_format as json
import numpy as np
from data_format import *
from format_module import *


def getXY(set, x_name, y_name):
    X = dict()
    Y = dict()
    for element in set.get_elements():
        X[int(element.get(x_name))] = None
        Y[int(element.get(y_name))] = None
    X.popitem()
    return list(X.keys()), list(Y.keys())

def convert_to_plot(set):
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
                    data_tag = Tag(id_names, [xv, 56, yv, radius, pitch, height, kappa])
                    noKdata_tag = Tag(id_names, [xv, 56, yv, radius, pitch, height, False])
                    sub = str(element)
                    kappa_value = set.get(data_tag.hashable).get("dMM")[int(sub[0])-1][int(sub[1])-1]
                    nokappa_value = set.get(noKdata_tag.hashable).get("dMM")[int(sub[0])-1][int(sub[1])-1]
                    dif = (kappa_value-nokappa_value)/kappa_value 
                    #Z[mmi][yi].append(np.log10(np.absolute(dif)))
                    Z[mmi][yi].append(nokappa_value)
                except:
                    Z[mmi][yi].append(None)
                #Z[mmi][yi].append(kappa_value)
    return X, Y, Z


id_names = ['azimuth','AOI','wvl','radius','pitch','height', 'kappa']

set = Structure.from_json("hashed.json")
X, Y, Z = convert_to_plot(set)
complete_MM_heatmap_plot("FCC","230")

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


