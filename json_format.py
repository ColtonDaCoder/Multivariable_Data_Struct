import json
import copy

def get_json(self, filename):
    with open(filename, 'r') as read_file:
        raw = json.load(read_file)
    return raw

def save_json(self, filename, data): 
    with open(filename, 'w') as write_file:
        write_file.write(json.dumps(data, indent=4))
    self.reload()

def get_data(self):
    return copy.deepcopy(self.data)
