import json_format as json
from data_format import *


id_names = ['azimuth','AOI','wvl','radius','pitch','height']
data = Structure()

#append to data structure
mm = []
data.append(Tag(id_names, [30,45, 210, 50, 175, 100]), mm)



#constraints for displaying data
#EX: two dynamic and the rest static to plot x, and y
data.add_constraint('azimuth', [0])
data.add_constraint('AOI', [0])
data.add_constraint('wvl', [0])
data.add_constraint('radius', [0])
data.add_constraint('pitch', [0])
data.add_constraint('height', [0])








