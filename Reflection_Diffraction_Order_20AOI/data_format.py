#------ ADD JSON DOWNLOADING -------
import json_format as json

class Tag():
#dictionary for identifiers of a tag
    def __init__(self, id_names: list, id_values: list):
        self.tag = dict()
        self.hashable = ""
        #makes dictionary of identifiers to search name for value
        for i in range(len(id_names)):
            self.tag[id_names[i]] = id_values[i]
            self.hashable = self.hashable + "_" + str(id_values[i])

#static constraint has range with length 1
#dynamic constraint has range with length 2
class Constraint():
    def __init__(self, name, range, is_static):
        self.name = name
        self.range = range
        self.is_static = is_static

    def check(self, element: dict()):
        if self.is_static:
            if int(self.range) == int(element.get(self.name)):
                return True
        else:
            if float(self.range[0]) <= float(element.get(self.name)) <= float(self.range[-1]):
                return True 
        return False

class Structure():
    #give list of tag names
    def __init__(self,):
        self.data = dict()
        self.constraints = []

    def add_constraint(self, name, range, is_static):
        self.constraints.append(Constraint(name, range, is_static))

    def append(self, data_tag: Tag, value_names: list, values: list):
        for i, name in enumerate(value_names):
            data_tag.tag[name] = values[i]
        self.data[data_tag.hashable] = data_tag.tag
    
    def save_json(self, filename):
        json.save_json(filename, self.data)

    def get(self, hashable):
        return self.data.get(hashable)
    
    #returns elements that fit the dynamic and static constraints
    def get_elements(self,):
        element_list = []
        for element in self.data:
            possible = True
            for constraint in self.constraints:
                if not constraint.check(element):
                    possible = False
                    break
            if(possible):
                element_list.append(self.data.get(element))
        return element_list

    def from_json(filename):
        new_structure = Structure()
        new_structure.data = json.get_json(filename)
        return new_structure
