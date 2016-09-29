def plot_2d_hist(**kwargs):
	import numpy as np 
	import pickle
	import matplotlib.pyplot as plt 
	import matplotlib

	matplotlib.rcParams['xtick.direction'] = 'out'
	matplotlib.rcParams['ytick.direction'] = 'out'
	matplotlib.rcParams['pdf.fonttype'] = 42

	PlotHistograms = kwargs.get('PlotHistograms',False)
	assert type(PlotHistograms) == bool, "PlotHistograms must be a boolean"

	[EccSumOfSquares, ConcSumOfSquares] = pickle.load(open('ALLSumOfSquares2.pkl','rb'))

	

	plt.figure()

	NumberOfBins = 25
	EccMax = 125
	ConcMax = 125
	if PlotHistograms == True: plt.subplot(2,2,3)
	plt.hist2d(EccSumOfSquares,ConcSumOfSquares,bins = [np.arange(0,EccMax,EccMax/NumberOfBins),np.arange(0,ConcMax,ConcMax/NumberOfBins)],\
					norm=matplotlib.colors.LogNorm())
	ax = plt.gca()
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	x0,x1 = plt.xlim(0,EccMax)
	y0,y1 = plt.ylim(0,ConcMax)
	ax.set_aspect(abs(x1-x0)/abs(y1-y0))
	plt.xlabel('Eccentric Velocity S.O.S.')
	plt.ylabel('Concentric Velocity S.O.S.')
	plt.set_cmap('bone_r')
	if PlotHistograms == False: plt.colorbar()
	
	if PlotHistograms == True:
		plt.subplot(2,2,4)
		plt.hist(ConcSumOfSquares,bins = np.arange(0,ConcMax,ConcMax/NumberOfBins), normed=True, color = '0.75',orientation = 'horizontal')
		ax = plt.gca()
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		X0,X1 = ax.get_xlim()
		plt.ylim(y0,y1)
		Y0,Y1 = ax.get_ylim()
		ax.set_aspect(abs(X1-X0)/abs(Y1-Y0))
		plt.xlabel('Frequency')

		plt.subplot(2,2,1)
		plt.hist(EccSumOfSquares,bins = np.arange(0,EccMax,EccMax/NumberOfBins), normed=True, color = '0.75')
		ax = plt.gca()
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		plt.xlim(x0,x1)
		X0,X1 = ax.get_xlim()
		Y0,Y1 = ax.get_ylim()
		ax.set_aspect(abs(X1-X0)/abs(Y1-Y0))
		plt.ylabel('Frequency')
	
	plt.show()
