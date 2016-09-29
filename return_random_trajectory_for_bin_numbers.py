def return_random_trajectory_for_bin_numbers(xbin,ybin):
	import numpy as np
	import pickle
	import matplotlib.pyplot as plt
	import random

	xmax = 125
	ymax = 125
	numberofbins = 25
	xbinwidth = xmax/numberofbins
	ybinwidth = ymax/numberofbins
	[eccSOS2, conSOS2] = pickle.load(open('ALLSumOfSquares2.pkl','rb'))
	_,xedges,yedges,_ = plt.hist2d(eccSOS2,conSOS2,bins = [np.arange(0,xmax,xbinwidth),np.arange(0,ymax,ybinwidth)])

	xbin_center = (xbin-1/2)*xbinwidth
	ybin_center = (ybin-1/2)*ybinwidth
	possible_values = []
	for i in range(len(eccSOS2)):
		if np.abs(eccSOS2[i]-xbin_center)<=xbinwidth/2 \
			and np.abs(conSOS2[i]-ybin_center)<=ybinwidth/2:
				possible_values.append(i+1)
	random.seed()
	random_choice = random.choice(possible_values)
	return(random_choice)