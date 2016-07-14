from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pickle


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

def animate_trajectory_number(TrialNumber,Angle1Splines,Angle2Splines,Angle3Splines):
	#------------------------------------------------------------
	# set up initial state and global variables

	dt = 0.0001
	Time = np.arange(0,0.55,dt)
	#shot_attempt = ThreeLinkSystem(22375,Angle1Splines,Angle2Splines,Angle3Splines,Time)

	Angle1 = Angle1Splines[TrialNumber].pp_func(Time)
	Angle2 = Angle2Splines[TrialNumber].pp_func(Time)
	Angle3 = Angle3Splines[TrialNumber].pp_func(Time)

	#------------------------------------------------------------
	# set up figure and animation
	fig = plt.figure()
	ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
	                     xlim=(-40, 80), ylim=(-60, 60))

	line, = ax.plot([], [], 'o-', lw=2)
	time_template = 'time = %.3fs'
	time_text = ax.text(0,-50, '', fontsize = 16)

	def init():
	    """initialize animation"""
	    line.set_data([], [])
	    time_text.set_text('')
	    return(line, time_text)

	def animate(i):
		"""perform animation step"""
		HeightInInches = 71
		Height = HeightInInches*2.54
		ShoulderToElbowLength = 0.186*Height
		ForearmLength = 0.146*Height
		HandLength = 0.108*Height
		x = np.cumsum([0, ShoulderToElbowLength*sin(Angle1[100*i]), ForearmLength*sin(Angle1[100*i]+Angle2[100*i]), HandLength*sin(Angle1[100*i]+Angle2[i]-Angle3[100*i])])
		y = np.cumsum([0, -ShoulderToElbowLength*cos(Angle1[100*i]), -ForearmLength*cos(Angle1[100*i]+Angle2[100*i]), -HandLength*cos(Angle1[100*i]+Angle2[i]-Angle3[100*i])])
		line.set_data(x,y)
		time_text.set_text(time_template % (i*dt*100))
		return(line, time_text) 

	# choose the interval based on dt and the time to animate one step
	from time import time
	t0 = time()
	animate(0)
	t1 = time()
	interval = 100

	ani = animation.FuncAnimation(fig, animate, frames = 55,
	                              interval=interval, repeat = True, blit=False, init_func=init)

	# save the animation as an mp4.  This requires ffmpeg or mencoder to be
	# installed.  The extra_args ensure that the x264 codec is used, so that
	# the video can be embedded in html5.  You may need to adjust this for
	# your system: for more information, see
	# http://matplotlib.sourceforge.net/api/animation_api.html
	#ani.save('TrialNumber'+np.str(TrialNumber)+'.mp4', fps=10)

	plt.show()