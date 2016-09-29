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

def plot_three_trajectory_numbers(good,fair,poor):
	import matplotlib.pyplot as plt
	import numpy as np 
	import pickle

	dt = 0.0001
	Time = np.arange(0,0.55,dt)
	fig1 = plt.figure()
	ax = fig1.gca()
	
	HeightInInches = 71
	Height = HeightInInches*2.54
	ShoulderToElbowLength = 0.186*Height
	ForearmLength = 0.146*Height
	HandLength = 0.108*Height
	
	GoodAngle1 = pickle.load(open('AllAngleSplines.pkl','rb'))[0][good-1].pp_func(Time)
	GoodAngle2 = pickle.load(open('AllAngleSplines.pkl','rb'))[1][good-1].pp_func(Time)
	GoodAngle3 = pickle.load(open('AllAngleSplines.pkl','rb'))[2][good-1].pp_func(Time)
	GoodX = [ShoulderToElbowLength*np.sin(GoodAngle1[i])+ForearmLength*np.sin(GoodAngle1[i]+GoodAngle2[i])+HandLength*np.sin(GoodAngle1[i]+GoodAngle2[i]-GoodAngle3[i]) for i in range(len(Time))]
	GoodY = [-ShoulderToElbowLength*np.cos(GoodAngle1[i])-ForearmLength*np.cos(GoodAngle1[i]+GoodAngle2[i])-HandLength*np.cos(GoodAngle1[i]+GoodAngle2[i]-GoodAngle3[i]) for i in range(len(Time))]

	plt.plot(GoodX,GoodY,color = '#71C177',lw = 2)

	FairAngle1 = pickle.load(open('AllAngleSplines.pkl','rb'))[0][fair-1].pp_func(Time)
	FairAngle2 = pickle.load(open('AllAngleSplines.pkl','rb'))[1][fair-1].pp_func(Time)
	FairAngle3 = pickle.load(open('AllAngleSplines.pkl','rb'))[2][fair-1].pp_func(Time)
	FairX = [ShoulderToElbowLength*np.sin(FairAngle1[i])+ForearmLength*np.sin(FairAngle1[i]+FairAngle2[i])+HandLength*np.sin(FairAngle1[i]+FairAngle2[i]-FairAngle3[i]) for i in range(len(Time))]
	FairY = [-ShoulderToElbowLength*np.cos(FairAngle1[i])-ForearmLength*np.cos(FairAngle1[i]+FairAngle2[i])-HandLength*np.cos(FairAngle1[i]+FairAngle2[i]-FairAngle3[i]) for i in range(len(Time))]

	plt.plot(FairX,FairY,color = '#F37722',lw = 2)

	PoorAngle1 = pickle.load(open('AllAngleSplines.pkl','rb'))[0][poor-1].pp_func(Time)
	PoorAngle2 = pickle.load(open('AllAngleSplines.pkl','rb'))[1][poor-1].pp_func(Time)
	PoorAngle3 = pickle.load(open('AllAngleSplines.pkl','rb'))[2][poor-1].pp_func(Time)
	PoorX = [ShoulderToElbowLength*np.sin(PoorAngle1[i])+ForearmLength*np.sin(PoorAngle1[i]+PoorAngle2[i])+HandLength*np.sin(PoorAngle1[i]+PoorAngle2[i]-PoorAngle3[i]) for i in range(len(Time))]
	PoorY = [-ShoulderToElbowLength*np.cos(PoorAngle1[i])-ForearmLength*np.cos(PoorAngle1[i]+PoorAngle2[i])-HandLength*np.cos(PoorAngle1[i]+PoorAngle2[i]-PoorAngle3[i]) for i in range(len(Time))]

	plt.plot(PoorX,PoorY,color = '#5D4EA1',lw = 2)

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
	ax.set_title('Good Trial Number: ' + str(good) +\
					'\nFair Trial Number: ' + str(fair) +\
					'\nPoor Trial Number: ' + str(poor))
	ax.set_aspect('equal', 'datalim')
	
	plt.show()
