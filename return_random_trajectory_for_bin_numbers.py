def return_random_trajectory_for_bin_numbers(xbin,ybin):
	import numpy as np
	import pickle
	import matplotlib.pyplot as plt
	import random

	binwidth = 10
	[eccSOS2, conSOS2] = pickle.load(open('SumOfSquares2.pkl','rb'))
	_,xedges,yedges,_ = plt.hist2d(eccSOS2,conSOS2,bins = [range(0,260,binwidth),range(0,150,binwidth)])

	xbin_center = (xbin-1/2)*binwidth
	ybin_center = (ybin-1/2)*binwidth
	possible_values = []
	for i in range(len(eccSOS2)):
		if np.abs(eccSOS2[i]-xbin_center)<=binwidth/2 \
			and np.abs(conSOS2[i]-ybin_center)<=binwidth/2:
				possible_values.append(i)

	random_choice = random.choice(possible_values)
	return(random_choice)


