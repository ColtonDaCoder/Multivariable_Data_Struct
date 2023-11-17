import pandas as pd
 
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

struct = Structure(columns=['azi','aoi','results'])

#adding results
results = Structure(columns=['mm', 'dmm', 'E'])
mm = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]
results.append({'mm':[mm], 'dmm' : [mm], 'E' : [5]})

#adding each spec + results
dict = {'azi' : 3,'aoi' : 1, 'results': results}
struct.append(dict)

struct.save_csv('test.csv')

print(pd.read_csv('test.csv'))

#print(struct.get([['azi', [1]]]).df)
