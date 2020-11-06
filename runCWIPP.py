# October 2020
#nmoisseeva@eoas.ubc.ca
# This code is a sample implementation of the model with test data

import numpy as np
import matplotlib.pyplot as plt
import imp

#import all common project variables
import config
imp.reload(config)
import utils
imp.reload(utils)
import cwipp
imp.reload(cwipp)
import graphics
imp.reload(graphics)


inputs = np.load(config.input_data,allow_pickle=True).item()

for fire in inputs:
	print('Processing fire: %s' %fire)
	case = cwipp.MODplume(fire)
	case.get_sounding(inputs[case.name]['sounding'])
	case.I = inputs[case.name]['I']
	case.iterate(config.biasFit)
	case.classify()
	case.get_profile()

	#temporary plot for sanity check
	plt.figure()
	plt.plot(fires[case.name]['truth'], config.interpZ)
	plt.plot(case.profile,config.interpZ)
	plt.show()

	#write back to fires dictionary
	inputs[case.name]['profile'] = case.profile

# np.save('output.npy',inputs, allow_pickle=True)
