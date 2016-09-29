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
def run_analysis_for_trial(shotnumber,xbin,ybin,N):
	import matplotlib.pyplot as plt 
	import numpy

	def statusbar(i,N,**kwargs):
		"""
		i is the current iteration (must be an int) and N is the length of 
		the range (must be an int). i must also be in [0,N). 
		
		~~~~~~~~~~~~~~
		**kwargs
		~~~~~~~~~~~~~~
		
		StartTime should equal time.time() and should be defined before your
		loop to ensure that you get an accurate representation of elapsed time.

		Title should be a str that will be displayed before the statusbar. Title
		should be no longer than 25 characters.

		~~~~~~~~~~~~~~

		NOTE: you should place a print('\n') after the loop to ensure you
		begin printing on the next line.

		"""
		import time
		StartTime = kwargs.get("StartTime",False)
		Title = kwargs.get("Title",'')

		assert type(i)==int, "i must be an int"
		assert type(N)==int, "N must be an int"
		assert N>i, "N must be greater than i"
		assert N>0, "N must be a positive integer"
		assert i>=0, "i must not be negative (can be zero)"
		assert type(Title) == str, "Title should be a string"
		assert len(Title) <= 25, "Title should be less than 25 characters"
		if Title != '': Title = ' '*(25-len(Title)) + Title + ': '
		statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
		if StartTime != False:
			print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) + 'sec        \r', end='')
		else:
			print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')

	def find_similar_shots(xbin,ybin,shotnumber):
		import numpy as np 
		import pickle
		import random
		import matplotlib.pyplot as plt 
		import time

		shotnumber -= 1
		dt = 0.0001
		Time = np.arange(0,0.55,dt)
		Splines = pickle.load(open('AllAngleSplines.pkl','rb'))
		Angle1Spline = Splines[0][shotnumber]
		Angle2Spline = Splines[1][shotnumber]
		Angle3Spline = Splines[2][shotnumber]
		def XY(A1,A2,A3,Time):
			HeightInInches = 71
			Height = HeightInInches*2.54
			ShoulderToElbowLength = 0.186*Height
			ForearmLength = 0.146*Height
			HandLength = 0.108*Height
			a1,a2,a3 = A1.pp_func(Time),A2.pp_func(Time),A3.pp_func(Time)
			X = [ShoulderToElbowLength*np.sin(a1[i])+ForearmLength*np.sin(a1[i]+a2[i])+HandLength*np.sin(a1[i]+a2[i]-a3[i]) for i in range(len(Time))]
			Y = [-ShoulderToElbowLength*np.cos(a1[i])-ForearmLength*np.cos(a1[i]+a2[i])-HandLength*np.cos(a1[i]+a2[i]-a3[i]) for i in range(len(Time))]
			return(np.array(X),np.array(Y))

		gx,gy = XY(Angle1Spline,Angle2Spline,Angle3Spline,Time)

		xmax = 125
		ymax = 125
		numberofbins = 25
		xbinwidth = xmax/numberofbins
		ybinwidth = ymax/numberofbins
		[eccSOS2, conSOS2] = pickle.load(open('ALLSumOfSquares2.pkl','rb'))
		_,xedges,yedges = np.histogram2d(eccSOS2,conSOS2,bins = [np.arange(0,xmax,xbinwidth),np.arange(0,ymax,ybinwidth)])

		xbin_center = (xbin-1/2)*xbinwidth
		ybin_center = (ybin-1/2)*ybinwidth
		possible_index_values = []
		for i in range(len(eccSOS2)):
			if np.abs(eccSOS2[i]-xbin_center)<=xbinwidth/2 \
				and np.abs(conSOS2[i]-ybin_center)<=ybinwidth/2:
					possible_index_values.append(i)
		def error(gx,gy,x,y):
			error_value = np.sum(((gx-x)**2+(gy-y)**2)**0.5)
			return(error_value)

		if len(possible_index_values) != 0:
			Error = []
			StartTime = time.time()
			for i in range(len(possible_index_values)):
				A1,A2,A3 = Splines[0][possible_index_values[i]],Splines[1][possible_index_values[i]],Splines[2][possible_index_values[i]]
				x,y = XY(A1,A2,A3,Time)
				Error.append(error(gx,gy,x,y))
				statusbar(i,len(possible_index_values),StartTime=StartTime,Title = "find_similar_shots")
			index = [x for (y,x) in sorted(zip(Error, possible_index_values))]
			#error = [Error[possible_index_values.index(y)] for y in index]
		else:
			index = 'None'		
		return(index)

	def plot_many_trajectories(reference_trial_number, additional_trials_index):
		"""
		additional_trials_index must be a list

		"""

		import numpy as np 
		import matplotlib.pyplot as plt 
		import pickle
		import matplotlib as mpl 
		import time
		import random

		mpl.rcParams['pdf.fonttype'] = 42
		mpl.rcParams['xtick.direction'] = 'out'
		mpl.rcParams['ytick.direction'] = 'out'

		if type(additional_trials_index) == int: additional_trials_index = [additional_trials_index]
		assert type(additional_trials_index) == list, "additional_trials_index must be a list"

		reference_trial_index = reference_trial_number - 1
		dt = 0.0001
		Time = np.arange(0,0.55,dt)
		Splines = pickle.load(open('AllAngleSplines.pkl','rb'))
		Angle1Spline = Splines[0][reference_trial_index]
		Angle2Spline = Splines[1][reference_trial_index]
		Angle3Spline = Splines[2][reference_trial_index]

		def XY(A1,A2,A3,Time):
			HeightInInches = 71
			Height = HeightInInches*2.54
			ShoulderToElbowLength = 0.186*Height
			ForearmLength = 0.146*Height
			HandLength = 0.108*Height
			a1,a2,a3 = A1.pp_func(Time),A2.pp_func(Time),A3.pp_func(Time)
			X = [ShoulderToElbowLength*np.sin(a1[i])+ForearmLength*np.sin(a1[i]+a2[i])+HandLength*np.sin(a1[i]+a2[i]-a3[i]) for i in range(len(Time))]
			Y = [-ShoulderToElbowLength*np.cos(a1[i])-ForearmLength*np.cos(a1[i]+a2[i])-HandLength*np.cos(a1[i]+a2[i]-a3[i]) for i in range(len(Time))]
			return(np.array(X),np.array(Y))

		Xref,Yref = XY(Angle1Spline,Angle2Spline,Angle3Spline,Time)
		N = len(additional_trials_index)
		if N<=5:
			row,col = 1,N
		else:
			col = 5
			if N%5 == 0:
				row = int(N/5)
			else:
				row = int((N-N%5)/5 + 1)
		plt.figure()
		ax = plt.gca()
		plt.suptitle('Trial Number ' + str(reference_trial_number))
		def remove_axes(ax):
			ax.spines['left'].set_color('none')
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
		StartTime = time.time()
		for i in range(N):
			plt.subplot(row,col,i+1)
			plt.plot(Xref,Yref,color = '#71C177',lw = 2)
			A1,A2,A3 = Splines[0][additional_trials_index[i]],Splines[1][additional_trials_index[i]],Splines[2][additional_trials_index[i]]
			X,Y = XY(A1,A2,A3,Time)
			plt.plot(X,Y,color = '#5D4EA1')
			ax = plt.gca()
			ax.set_title('Trial Number ' + str(additional_trials_index[i]+1))
			remove_axes(ax)
			statusbar(i,N,StartTime = StartTime,Title = "plot_many_trajectories")
		print('\n')

	index = find_similar_shots(xbin,ybin,shotnumber)
	if index == "None": 
		print("No shots in these bins")
	else:
		if len(index) < N: N = len(index)
		print('\n')
		plot_many_trajectories(shotnumber,index[:N])
		plt.show()
