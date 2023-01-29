import json_format as json
import variables_format as vars


#create large data
def create_large_data():
    pass

def create_vars():
    pass


def create_data(dynamic: list, static: list):
    interval = (0,2)
    data = {}
    possible = True
    for dlist in dynamic:
        for step in interval(dlist):
            indexes = dlist.get(step)
            for i in indexes:
                if not data.get(i):
                    for slist in static
                        #create large_data
                        #create variables (static.value) in variables format
                        if not large_data.get(i).get(slist) == s.value:
                            possible = False
                            break
                    if possible:
                        data.add(i)


