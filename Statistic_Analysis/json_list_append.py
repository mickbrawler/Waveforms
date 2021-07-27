import glob
import json
import argparse

# This function will be run on the jsons of pot_threshold
# pot_thresholds is a 2d list, where each list is a stat's onsource peaks for a trial
# The lists will have each trials values added to them until we havethe onsource peaks for all trials 

def json_list_append(jsons_path, merge_path_name):
    print("entering list_append")
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
                C_dictionary_new = json.load(f)
        
            for i in C_dictionary:
                C_dictionary[str(i)] += C_dictionary_new[str(i)]

    with open('{}.json'.format(merge_path_name), 'w') as f:
        json.dump(C_dictionary, f, indent=2)
    print("leaving list_append")

# def json_list_append(jsons_path, merge_path_name):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--jsons_path', type=str)
    parser.add_argument('--merge_path_name', type=str)
    args = parser.parse_args()

    json_list_append(args.jsons_path,args.merge_path_name)
