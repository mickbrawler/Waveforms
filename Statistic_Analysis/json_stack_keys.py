import glob
import json
import argparse

# works like a charm, put in your path and
# make sure to include the path of the produced json in their merge name

# This function will be run on the jsons of output, RHO_MOD, and MAX_OS
# It updates all trial dictionaries, and since each dictionary has a different key 
# (which indicates the trial used), then will be added, till 1 dictionary for all
# trials is made

def json_stack_keys(jsons_path, merge_path_name):
    print("in stack_keys.py")
    # Include last / at end of path!
    files = glob.glob("{}*.json".format(jsons_path))

    count = 0
    for file in files:
        if count == 0:
            count += 1
            with open(file, "r") as f:
                C_dictionary = json.load(f)
        else:
            with open(file, "r") as f:
                C_dictionary.update(json.load(f))

    with open('{}.json'.format(merge_path_name), 'w') as f:
        json.dump(C_dictionary, f, indent=2)
    
    print("leaving stack_key")
#def json_stack_keys(jsons_path, merge_path_name):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--jsons_path', type=str)
    parser.add_argument('--merge_path_name', type=str)
    args = parser.parse_args()

    json_stack_keys(args.jsons_path,args.merge_path_name)
