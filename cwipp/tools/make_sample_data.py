#November 2020
#nmoisseeva@eoas.ubc.ca
#this script makes some sample data to use for testing the module
import pickle
import numpy as np
import config

#load pickled plume data
with open('plumes.pkl','rb') as f:
     all_plumes = pickle.load(f)

#create a storage dictionary (to be replaced with json in pipeline)
fires = {}

#make a new dictionary containing input data only
for plume in all_plumes[:15]:
    fires[plume.name] = {}
    fires[plume.name]['sounding'] = np.load(config.wrfdir + 'profiles/profT0' + plume.name + '.npy')
    fires[plume.name]['I'] = plume.I
    fires[plume.name]['truth'] = plume.profile

np.save('sample_fires.npy',fires, allow_pickle=True)
