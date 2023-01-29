#------ ADD JSON DOWNLOADING -------
import json_format as json

#dictionary for identifiers of a tag
class Tag():
    def __init__(self, id_names: list, id_values: list):
        self.tag = dict()
        #makes dictionary of identifiers to search name for value
        for i in range(len(id_names)):
            self.tag[id_names[i]] = id_values[i]
        return tag 

#static constraint has range with length 1
#dynamic constraint has range with length 2
class Constraint():
    def __init__(self, name, range):
        self.name = name
        self.range = range 

    def check(self, element):
        if range[0] <= element.get(self.name) <= range[-1]:
            return True
        else:
            return False

class Structure():
    #give list of tag names
    def __init__(self,):
        self.data = [] 
        self.constraints = []

    def add_constraint(self, name, range):
        self.constraints.append(Constraint(name, range))

    def append(self, tag: Tag, value):
        self.data.append([tag, value])
    
    #returns elements that fit the dynamic and static constraints
    def get_elements(self,):
        element_list = []
        for element in data:
            for constraint in constraints:
                if constraint.check(element):
                    element_list.append(element)
        return element_list
