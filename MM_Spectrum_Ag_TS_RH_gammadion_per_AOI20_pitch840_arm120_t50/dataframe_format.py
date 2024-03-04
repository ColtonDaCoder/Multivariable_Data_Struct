import pandas as pd
import numpy as np
 
class Structure():
    def __init__(self, columns=[]):
        self.columns = columns 
        self.df = pd.DataFrame(columns=self.columns)

    def append(self,dict):
        self.df = pd.concat([self.df, pd.DataFrame(dict, index=[0])], ignore_index=True)

    def save_csv(self, filename=None):
        if filename is not None:
            self.filename = filename
        self.df.to_csv(self.filename, index=False)
    
    def get(self, search):
        # search format: [[name, [range of values],[name, [range of values]]]
        holder = Structure(columns=self.df.columns)
        return_struc = self
        for input in search:
            for sim in return_struc.df.to_dict(orient='records'):
                for value in input[1]:
                    if sim[input[0]] == value:
                        holder.append(sim)
            return_struc = holder 
            holder = Structure(columns=self.df.columns)
        return return_struc

    def from_csv(filename):
        df = pd.read_csv(filename)
        holder = Structure(df.columns)
        holder.df = df
        return holder
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
