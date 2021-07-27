import glob
import json
import argparse

# This function will be run on the jsons of MAX_BG_TEMP
# MAX_BG_TEMP always has the same keys, as they are the templates used per trial
# As such we update the max offsource value for each stat of each template of each trial
# We therefore have the crosses across all trials

def json_update_components(jsons_path, merge_path_name):
    print("inside components")    
    with open("Merged_jsons/Merged_Peaks.json", "r") as f:
        RHO_MOD = json.load(f)
    
    # way of getting the first key of a dictionary, used to get tempn and statn quick and dirty
    get_first_key = []
    get_first_key += RHO_MOD.keys()
    dict_key = get_first_key[0]
    
    # way of getting tempn and statn from RHO_MOD
    tempn = len(RHO_MOD[dict_key][0][0])
    statn = len(RHO_MOD[dict_key][0])
    
    # Include last / at end of path!
    files = glob.glob("{}*.json".format(jsons_path))

    count = 0
    for file in files:
#         trialn = "" # method of obtaining this file's trialnumber through its name
#         for character in file:
#             if character.isdigit():
#                 trialn += str(character)
#         trialn = int(trialn)

        if count == 0:
            count += 1
            with open(file, "r") as f:
                C_dictionary = json.load(f)
        else:
            with open(file, "r") as f:
                C_dictionary_new = json.load(f)
                
            for stat in range(statn):
                for j in range(tempn):
                
                    if (C_dictionary_new[str(j)][stat] > C_dictionary[str(j)][stat]): # every value in MAX_BG_TEMP dictionary changes from 0 to that templates max
                        C_dictionary[str(j)][stat] = C_dictionary_new[str(j)][stat]

    with open('{}.json'.format(merge_path_name), 'w') as f:
        json.dump(C_dictionary, f, indent=2)
    print("leaving components")

# def json_update_components(jsons_path, merge_path_name):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--jsons_path', type=str)
    parser.add_argument('--merge_path_name', type=str)
    args = parser.parse_args()

    json_update_components(args.jsons_path,args.merge_path_name)
