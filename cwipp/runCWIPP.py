# October 2020
#nmoisseeva@eoas.ubc.ca
# This code is a sample implementation of the model with test data

import numpy as np
import matplotlib.pyplot as plt
import imp

#import all common project variables
import config
import utils
import cwipp
import graphics

print('==========CWIPP PLUME-RISE PARAMETERIZATION SCHEME========')
#load inputs
inputs = np.load(config.input_data,allow_pickle=True).item()
print('Loading data: %s' %config.input_data)

#loop through all fires
for fire in inputs:
    print('Processing fire: %s' %fire)
    case = cwipp.MODplume(fire)
    case.get_sounding(inputs[case.name]['sounding'])
    case.I = inputs[case.name]['I']
    case.iterate(config.biasFit)
    case.classify()

    #FUTURE ADDITION: incroporate vertical profile code here to reconstruct full profile here
    
    #update data
    inputs[case.name]['zCL'] = case.zCL
'''
    #temporary plot for sanity check (use with sample_fires.npy)
    plt.figure()
    plt.title('TRUE PROFILE cs MODEL GUESS')
    plt.plot(inputs[case.name]['truth'], config.interpZ,color='grey', label='target smoke profile')
    plt.axhline(y = case.zCL, color='C1', label='modelled injection height')
    plt.gca().set(xlabel='smoke concentration', ylabel='height [m]')
    plt.legend()
    plt.show()
'''
#write back to fires dictionary
np.save('output.npy',inputs, allow_pickle=True)

print('...COMPLETED')
