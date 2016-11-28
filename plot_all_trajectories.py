class Spline:
	"""
	Initiate a class variable spline that has one break at x = x_break starting at x_initial and has
	the equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3.

	pp_func(X)
	~~~~~~~~~~~~~~~~~~~

	Takes in X array and outputs the piecewise polynomial associated with this spline.

	pp_deriv()
	~~~~~~~~~~~~~~~~~~~

	Takes in X array and outputs the piecewise polynomial associated with the spline's derivative.

	find_max_and_min()
	~~~~~~~~~~~~~~~~~~~

	Takes in the min and max values for both x and y and will find the maximum y values of the piecewise
	polynomial. To do this, first we find the extrema point (find_extrema) by inputing the x values that 
	set the derivate of the piecewise polynomial equal to zero (quadratic formula). Next we ensure that
	the zero values are in fact real (is_real). We then filter out the zeros that are not in the 
	appropriate domains (is_in_appropriate_domain). To see if these values are maximum or minimum, we 
	plug them back into the second derivative of the appropriate piecewise polynomial (second_deriv_is_neg()
	and second_deriv_is_pos(), respectively). Finally we determine the y value of these extrema by using 
	the class function self.pp_func().

	is_initial_slope_positive()
	~~~~~~~~~~~~~~~~~~~

	This takes in X and will check to make sure that for the first 2500 entries in X, that the derivative 
	of the piecewise polynomial (pp_deriv()) will be positive. Make sure that X is at least 2500 in length.

	is_within_bounds()
	~~~~~~~~~~~~~~~~~~~

	This checks to see if the maximum maximum value and the minimum mininum value calculated above will 
	fall between y_min and y_max. This makes use of self.find_max_and_min()

	"""

	def __init__(self,a,b,c,d,x_initial,x_break):
		import numpy as np
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.x_initial = x_initial
		self.x_break = x_break
		#self.all_values = {'A': a, 'B' : b, 'C' : c, 'D' : d, 'init' : x_initial, 'break' : x_break}
	def pp_func(self,X):
		import numpy as np
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
									[lambda X: self.a[0] + self.b[0]*(X-self.x_initial) + self.c[0]*(X-self.x_initial)**2 + self.d[0]*(X-self.x_initial)**3, \
									lambda X: self.a[1] + self.b[1]*(X-self.x_break) + self.c[1]*(X-self.x_break)**2 + self.d[1]*(X-self.x_break)**3])
		return(result)
	def pp_deriv(self,X):
		import numpy as np
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
									[lambda X: self.b[0] + 2*self.c[0]*(X-self.x_initial) + 3*self.d[0]*(X-self.x_initial)**2, \
									lambda X: self.b[1] + 2*self.c[1]*(X-self.x_break) + 3*self.d[1]*(X-self.x_break)**2])
		return(result)
	def find_max_and_min(self,x_min,x_max,y_min,y_max):
		def find_extrema():
			import numpy as np
			extrema_1 = np.float(self.x_initial + (- 2*self.c[0] + (4*self.c[0]**2 - 12*self.b[0]*self.d[0])**.5)/(6*self.d[0]))
			extrema_2 = np.float(self.x_initial + (- 2*self.c[0] - (4*self.c[0]**2 - 12*self.b[0]*self.d[0])**.5)/(6*self.d[0]))
			extrema_3 = np.float(self.x_break + (- 2*self.c[1] + (4*self.c[1]**2 - 12*self.b[1]*self.d[1])**.5)/(6*self.d[1]))
			extrema_4 = np.float(self.x_break + (- 2*self.c[1] - (4*self.c[1]**2 - 12*self.b[1]*self.d[1])**.5)/(6*self.d[1]))
			return(extrema_1,extrema_2,extrema_3,extrema_4)
		def is_real(x_value):
			result = not isinstance(x_value,complex)
			return(result)
		def is_in_appropriate_domain(x_value,x_min,x_max,segment_number):
			if segment_number == 1:
				result = x_value >= x_min and x_value <= self.x_break
			elif segment_number == 2:
				result = x_value >= self.x_break and x_value <= x_max
			return(result)
		def second_deriv_is_neg(x_value,segment_number):
			if segment_number == 1:
				x_not = self.x_initial
			elif segment_number == 2:
				x_not = self.x_break
			second_deriv =  2*self.c[segment_number-1] + 6*self.d[segment_number-1]*(x_value-x_not)
			result = second_deriv<0
			return(result)
		def second_deriv_is_pos(x_value,segment_number):
			if segment_number == 1:
				x_not = self.x_initial
			elif segment_number == 2:
				x_not = self.x_break
			second_deriv =  2*self.c[segment_number-1] + 6*self.d[segment_number-1]*(x_value-x_not)
			result = second_deriv>0
			return(result)
		def determine_if_max_or_min(extrema_1,extrema_2,extrema_3,extrema_4,x_min,x_max):
			import numpy as np
			maxima = []
			minima = []
			if is_real(extrema_1) and is_in_appropriate_domain(extrema_1,x_min,x_max,1):
				if second_deriv_is_neg(extrema_1,1):
					maxima.append(np.float(self.pp_func(extrema_1)))
				elif second_deriv_is_pos(extrema_1,1):
					minima.append(np.float(self.pp_func(extrema_1)))
			if is_real(extrema_2) and is_in_appropriate_domain(extrema_2,x_min,x_max,1):
				if second_deriv_is_neg(extrema_2,1):
					maxima.append(np.float(self.pp_func(extrema_2)))
				elif second_deriv_is_pos(extrema_2,1):
					minima.append(np.float(self.pp_func(extrema_2)))
			if is_real(extrema_3) and is_in_appropriate_domain(extrema_3,x_min,x_max,2):
				if second_deriv_is_neg(extrema_3,2):
					maxima.append(np.float(self.pp_func(extrema_3)))
				elif second_deriv_is_pos(extrema_3,2):
					minima.append(np.float(self.pp_func(extrema_3)))
			if is_real(extrema_4) and is_in_appropriate_domain(extrema_4,x_min,x_max,2):
				if second_deriv_is_neg(extrema_4,2):
					maxima.append(np.float(self.pp_func(extrema_4)))
				elif second_deriv_is_pos(extrema_4,2):
					minima.append(np.float(self.pp_func(extrema_4)))
			return(maxima,minima)
		extrema_1,extrema_2,extrema_3,extrema_4 = find_extrema()
		maxima, minima = determine_if_max_or_min(extrema_1,extrema_2,extrema_3,extrema_4,x_min,x_max)
		return(maxima,minima)
	def is_initial_slope_positive(self,X):
		result = min(self.pp_deriv(X[:2501]))>=0
		return(result)
	def is_within_bounds(self,x_min,x_max,y_min,y_max):
		import numpy as np
		maxima,minima = self.find_max_and_min(x_min,x_max,y_min,y_max)
		if len(maxima) == 0:
			maxima = y_max
		if len(minima) == 0:
			minima = y_min
		result = np.max(maxima) <= y_max and np.min(minima) >= y_min
		return(result)

def plot_all_trajectories(**kwargs):
	"""
	**kwargs
	=======================================

	PlotAngles is a boolean that will plot the 3 joint angles if True.
	Default is set to false.

	"""

	import numpy as np 
	import matplotlib.pyplot as plt 
	import pickle
	import matplotlib as mpl 
	from danpy.stuff import statusbar
	import time
	import ipdb
	from mpl_toolkits.mplot3d import Axes3D

	mpl.rcParams['pdf.fonttype'] = 42
	mpl.rcParams['xtick.direction'] = 'out'
	mpl.rcParams['ytick.direction'] = 'out'

	Range = kwargs.get('Range',None)
	if Range != None:
		assert type(Range) == list or type(Range)==np.ndarray, "Range must be a list/array"
		assert len(Range) == 2, "Range must be a 1x2 list/array"
	PlotConfigurationSpace = kwargs.get('PlotConfigurationSpace',False)
	assert type(PlotConfigurationSpace)==bool, "PlotConfigurationSpace must be either True/False."

	dt = 0.0001
	Time = np.arange(0,0.55,dt)
	HeightInInches = 71
	Height = HeightInInches*2.54
	ShoulderToElbowLength = 0.186*Height
	ForearmLength = 0.146*Height
	HandLength = 0.108*Height
	x = lambda a1,a2,a3: np.cumsum([0, ShoulderToElbowLength*np.sin(a1), ForearmLength*np.sin(a1+a2), HandLength*np.sin(a1+a2-a3)])
	y = lambda a1,a2,a3: np.cumsum([0, -ShoulderToElbowLength*np.cos(a1), -ForearmLength*np.cos(a1+a2), -HandLength*np.cos(a1+a2-a3)])

	if Range == None:
		Splines = pickle.load(open('ALLAngleSplines.pkl','rb'))
		AppropriateIndex = np.random.choice(range(0,len(Splines[0])),20)
		Title = 'Plot All Trajectories'
	else:
		Splines = pickle.load(open('ALLAngleSplines.pkl','rb'))
		[EccSumOfSquares, ConcSumOfSquares] = pickle.load(open('ALLSumOfSquares2.pkl','rb'))
		AppropriateIndex = np.where(np.logical_and(EccSumOfSquares>=Range[0],EccSumOfSquares<=Range[1]))[0]
		AppropriateIndex = np.random.choice(AppropriateIndex,20)
		Title = 'Plot within range ('+str(Range[0])+', '+str(Range[1])+')'

	fig1 = plt.figure(1)
	ax = plt.gca()
	if PlotConfigurationSpace == True:
		fig2 = plt.figure(2,figsize = (8,4))
		ax_config = fig2.gca(projection='3d')
		plt.figure(1)

	StartTime = time.time()
	
	for i in range(len(AppropriateIndex)):
		Angle1Spline = Splines[0][AppropriateIndex[i]]
		Angle2Spline = Splines[1][AppropriateIndex[i]]
		Angle3Spline = Splines[2][AppropriateIndex[i]]
		Angle1 = Angle1Spline.pp_func(Time)
		Angle2 = Angle2Spline.pp_func(Time)
		Angle3 = Angle3Spline.pp_func(Time)
		X = [ShoulderToElbowLength*np.sin(Angle1[j])+ForearmLength*np.sin(Angle1[j]+Angle2[j])+HandLength*np.sin(Angle1[j]+Angle2[j]-Angle3[j]) for j in range(len(Time))]
		Y = [-ShoulderToElbowLength*np.cos(Angle1[j])-ForearmLength*np.cos(Angle1[j]+Angle2[j])-HandLength*np.cos(Angle1[j]+Angle2[j]-Angle3[j]) for j in range(len(Time))]
		plt.plot(X,Y,'0.75')
		if PlotConfigurationSpace==True: 
			plt.figure(2,figsize = (8,4))
			plt.plot((180/np.pi)*Angle1,(180/np.pi)*Angle2,(180/np.pi)*Angle3,color = '0.75',lw = 1)
			plt.plot([(180/np.pi)*Angle1[-1]],[(180/np.pi)*Angle2[-1]],[(180/np.pi)*Angle3[-1]],marker = 'o', color = 'k')
			plt.plot([(180/np.pi)*Angle1[0]],[(180/np.pi)*Angle2[0]],[(180/np.pi)*Angle3[0]],marker = 'o', color = 'k')
			ax_config.set_xlabel('PSR\n(in degrees)')
			ax_config.set_ylabel('EFE\n(in degrees)')
			ax_config.set_zlabel('WFE\n(in degrees)')
			ax_config.set_xlim3d(-15, 145)
			ax_config.set_ylim3d(75, 140)
			ax_config.set_zlim3d(-95, 5)
			ax_config.view_init(24, -110)
			plt.figure(1)
		statusbar(i,len(AppropriateIndex),Title="Plot All Trajectories",StartTime=StartTime)
	plt.plot(x(Angle1[0],Angle2[0],Angle3[0]),y(Angle1[0],Angle2[0],Angle3[0]),'k',linestyle = '-',marker = 'o') # initial configuration
	plt.plot(x(Angle1[-1],Angle2[-1],Angle3[-1]),y(Angle1[-1],Angle2[-1],Angle3[-1]),'k',linestyle = '-',marker = 'o') # final configuruationax.spines['left'].set_color('none')
	plt.plot(6*np.sin(np.arange(0,2*np.pi,0.001)),6*np.cos(np.arange(0,2*np.pi,0.001)),'k')
	ax.spines['right'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['top'].set_color('none')
	plt.yticks([])
	plt.tick_params(
	    axis='both',          # changes apply to both axes
	    which='both',      # both major and minor ticks are affected
	    bottom='off',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    right = 'off',		# ticks along the right edge are off
	    left = 'off',		# ticks along the left edge are off
	    labelbottom='off') # labels along the bottom edge are off
	ax.set_aspect('equal', 'datalim')

	plt.show()
			
