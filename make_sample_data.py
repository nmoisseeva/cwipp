#November 2020
#nmoisseeva@eoas.ubc.ca
#this script makes some sample data to use for testing the module
import pickle
import numpy as np

#load pickled plume data
with open('plumes.pkl','rb') as f:
     all_plumes = pickle.load(f)

#create a storage dictionary (to be replaced with json in pipeline)
fires = {}

#to TEST!!!!!!
#make a sounding dataset [z,T,U]
T0 = np.load(config.wrfdir + 'profiles/profT0' + name + '.npy')    #load initial temperature profile
U0 = np.load(config.wrfdir + 'profiles/profU0' + name + '.npy')    #load initial temperature profile
metlvls = np.arange(0,len(T0)*config.dz,config.dz)
sounding = [metlvls,T0,U0]

#make a new dictionary containing input data only
for plume in all_plumes[:15]:
    fires[plume.name] = {}
    fires[plume.name]['sounding'] = sounding
    fires[plume.name]['I'] = plume.I
    fires[plume.name]['truth'] = plume.profile

np.save('fires.npy',fires, allow_pickle=True)
