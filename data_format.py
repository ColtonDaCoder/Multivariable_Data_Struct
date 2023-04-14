#------ ADD JSON DOWNLOADING -------
import json_format as json
import numpy as np

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

    #give list of tag names
    def __init__(self,):
        self.data = dict()
        self.constraints = []

    def full_decomp(self):
        for hashable in self.data.keys():
            self.data.get(hashable)['dMM'] = self.cloude_decomp(self.data.get(hashable).get('mm')).tolist()

    def full_decomp_diff_orders(self,):
        for hashable in self.data.keys():
            dMM = []
            for matrix in self.data.get(hashable).get('mm'):
                dMM.append(self.cloude_decomp(matrix).tolist())
            self.data.get(hashable)['dMM'] = dMM

    def add_set(self, set):
        for hashable in set.data.keys():
            self.data[hashable] = set.get(hashable)

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
