import ipdb
import numpy as np 
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
import pickle



EndTime = 0.55
ChangeInTime = 0.0001
Time = np.arange(0,EndTime+ChangeInTime,ChangeInTime,float)

PositionInX = []
PositionInY = []

HeightInInches = 71
DegreeToRadianFactor = np.pi/180

Height = HeightInInches*2.54
ShoulderToElbowLength = 0.186*Height
ForearmLength = 0.146*Height
HandLength = 0.108*Height

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
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.x_initial = x_initial
		self.x_break = x_break
		#self.all_values = {'A': a, 'B' : b, 'C' : c, 'D' : d, 'init' : x_initial, 'break' : x_break}
	def pp_func(self,X):
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
									[lambda X: self.a[0] + self.b[0]*(X-self.x_initial) + self.c[0]*(X-self.x_initial)**2 + self.d[0]*(X-self.x_initial)**3, \
									lambda X: self.a[1] + self.b[1]*(X-self.x_break) + self.c[1]*(X-self.x_break)**2 + self.d[1]*(X-self.x_break)**3])
		return(result)
	def pp_deriv(self,X):
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
									[lambda X: self.b[0] + 2*self.c[0]*(X-self.x_initial) + 3*self.d[0]*(X-self.x_initial)**2, \
									lambda X: self.b[1] + 2*self.c[1]*(X-self.x_break) + 3*self.d[1]*(X-self.x_break)**2])
		return(result)
	def find_max_and_min(self,x_min,x_max,y_min,y_max):
		def find_extrema():
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
		maxima,minima = self.find_max_and_min(x_min,x_max,y_min,y_max)
		if len(maxima) == 0:
			maxima = y_max
		if len(minima) == 0:
			minima = y_min
		result = np.max(maxima) <= y_max and np.min(minima) >= y_min
		return(result)

Angle1Splines = pickle.load(open('SplineClassObjects.pkl','rb'))[0]
Angle2Splines = pickle.load(open('SplineClassObjects.pkl','rb'))[1]
Angle3Splines = pickle.load(open('SplineClassObjects.pkl','rb'))[2]



def find_trajectories(xbin,ybin):
	"""
	xbin and ybin should be integers between 1 and 25. Will return all trial numbers that fall within that bin.
	"""
	Ecc_Sum_of_Squares = np.array(pickle.load(open('SumOfSquares.pkl','rb'))[0])
	Conc_Sum_of_Squares = np.array(pickle.load(open('SumOfSquares.pkl','rb'))[1])
	x_lower, x_upper = (xbin-1)*25, (xbin)*25
	y_lower, y_upper = (ybin-1)*25, (ybin)*25
	results = np.nonzero([Ecc_Sum_of_Squares[i] >= x_lower and \
						Ecc_Sum_of_Squares[i] <= x_upper and \
						Conc_Sum_of_Squares[i] >= y_lower and \
						Conc_Sum_of_Squares[i] <= y_upper for i in range(100000)])[0]
	return(results)

def position_in_x(Angle1, Angle2, Angle3):
	"""
	Takes in the values of Angle1, Angle2, and Angle3 and generates the x position of the endpoint.
	"""
	PositionInX = ShoulderToElbowLength*sin(Angle1) \
					+ ForearmLength*sin(Angle1+Angle2) \
					+ HandLength*sin(Angle1+Angle2-Angle3)
	return(PositionInX)
def position_in_y(Angle1,Angle2,Angle3):
	"""
	Takes in the values of Angle1, Angle2, and Angle3 and generates the y position of the endpoint.
	"""
	PositionInY = -ShoulderToElbowLength*cos(Angle1) \
					- ForearmLength*cos(Angle1+Angle2) \
					- HandLength*cos(Angle1+Angle2-Angle3)
	return(PositionInY)
def generate_posture_array(angle1,angle2,angle3):
	xresult = np.array([	0, \
							ShoulderToElbowLength*sin(angle1), \
							ShoulderToElbowLength*sin(angle1) + ForearmLength*sin(angle1+angle2), \
							ShoulderToElbowLength*sin(angle1) + ForearmLength*sin(angle1+angle2) + HandLength*sin(angle1+angle2-angle3)])
	yresult = np.array([	0, \
							-ShoulderToElbowLength*cos(angle1), \
							-ShoulderToElbowLength*cos(angle1) - ForearmLength*cos(angle1+angle2), \
							-ShoulderToElbowLength*cos(angle1) - ForearmLength*cos(angle1+angle2) - HandLength*cos(angle1+angle2-angle3) ])
	
	return(xresult,yresult)
def plot_random_trajectory(xbin,ybin,Angle1Splines,Angle2Splines,Angle3Splines,Time):
	TrialNumbers = find_trajectories(xbin,ybin)
	np.random.seed()
	RandomTrialNumber = np.random.choice(TrialNumbers)
	angle1 = Angle1Splines[RandomTrialNumber].pp_func(Time)
	angle2 = Angle2Splines[RandomTrialNumber].pp_func(Time)
	angle3 = Angle3Splines[RandomTrialNumber].pp_func(Time)
	X,Y = position_in_x(angle1,angle2,angle3), position_in_y(angle1,angle2,angle3)
	xcirc, ycirc = 5*sin(np.arange(0,2*np.pi,0.01)), 5*cos(np.arange(0,2*np.pi,0.01))
	plt.figure()
	plt.plot(X,Y, lw = 3)
	plt.plot(xcirc,ycirc,'k')
	for i in range(0,len(Time),np.int(len(Time)/10)):
		x_arm,y_arm = generate_posture_array(angle1[i],angle2[i],angle3[i])
		plt.plot(x_arm,y_arm, 'k')
	plt.xticks([])
	plt.yticks([])	
	plt.axis('equal')
	plt.title('Trial '+ str(RandomTrialNumber))
	plt.show()




