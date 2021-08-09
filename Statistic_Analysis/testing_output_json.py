import json


with open("all_dir/Merged_jsons/Merged_output.json", "r") as f:
    output = json.load(f)

#print(output["0"][0])
#print(output["0"][1])
#print(output["0"][2])
print(output["0"][3])
