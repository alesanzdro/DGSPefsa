
import sys
import json
import os.path
path_checkm = str(sys.argv[1])

#path_checkm = '/ALMEIDA/PROJECTS/BACTERIAS/DGSP/analysis/230210_LSPV002/pipeline/22_LM_03313/checkm/storage/bin_stats.analyze.tsv'

# Check if file exist
if os.path.isfile(path_checkm):
    # deserializes into dict 
    # and returns dict.
    with open(path_checkm) as f:
        data = json.load(f) 
    f = open(path_checkm, "r")
    data = f.read()

else:
    print ("ERROR - CHECKM FILE bin_stats.analyze.tsv DOES NOT EXIST")

# 

# JSON string
a = '{"name": "Bob", "languages": "English"}'
  


y = json.loads(a)
  
print("JSON string = ", y)
print()
  
  
  
# JSON file
f = open ('data.json', "r")
  
# Reading from file
data = json.loads(f.read())
  
# Iterating through the json
# list
for i in data['emp_details']:
    print(i)
  
# Closing file
f.close()