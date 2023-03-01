import json
import copy

def get_json(filename):
    with open(filename, 'r') as read_file:
        raw = json.load(read_file)
    return raw

def save_json(filename, data): 
    with open(filename, 'w') as write_file:
        write_file.write(json.dumps(data, indent=4))

