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

def plot_trajectory_number(trialnumber, **kwargs):
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

	mpl.rcParams['pdf.fonttype'] = 42

	PlotAngles = kwargs.get('PlotAngles',False)
	PlotVm = kwargs.get('PlotVm',False)
	assert type(PlotAngles) == bool, "PlotAngles must be a boolean"
	assert type(PlotVm) == bool, "PlotVm must be a boolean"

	dt = 0.0001
	Time = np.arange(0,0.55,dt)
	Splines = pickle.load(open('SplineClassObjects2.pkl','rb'))
	Angle1Spline = Splines[0][trialnumber]
	Angle2Spline = Splines[1][trialnumber]
	Angle3Spline = Splines[2][trialnumber]

	Angle1 = Angle1Spline.pp_func(Time)
	Angle2 = Angle2Spline.pp_func(Time)
	Angle3 = Angle3Spline.pp_func(Time)

	HeightInInches = 71
	Height = HeightInInches*2.54
	ShoulderToElbowLength = 0.186*Height
	ForearmLength = 0.146*Height
	HandLength = 0.108*Height
	x = lambda a1,a2,a3: np.cumsum([0, ShoulderToElbowLength*np.sin(a1), ForearmLength*np.sin(a1+a2), HandLength*np.sin(a1+a2-a3)])
	y = lambda a1,a2,a3: np.cumsum([0, -ShoulderToElbowLength*np.cos(a1), -ForearmLength*np.cos(a1+a2), -HandLength*np.cos(a1+a2-a3)])

	X = [ShoulderToElbowLength*np.sin(Angle1[i])+ForearmLength*np.sin(Angle1[i]+Angle2[i])+HandLength*np.sin(Angle1[i]+Angle2[i]-Angle3[i]) for i in range(len(Time))]
	Y = [-ShoulderToElbowLength*np.cos(Angle1[i])-ForearmLength*np.cos(Angle1[i]+Angle2[i])-HandLength*np.cos(Angle1[i]+Angle2[i]-Angle3[i]) for i in range(len(Time))]
	
	LengthOfIntervals = 500
	plt.cla()

	plt.figure
	ax = plt.gca()
	[plt.plot(x(Angle1[i],Angle2[i],Angle3[i]),y(Angle1[i],Angle2[i],Angle3[i]),'0.75') \
		for i in range(0,len(Time),LengthOfIntervals)]
	plt.plot(x(Angle1[-1],Angle2[-1],Angle3[-1]),y(Angle1[-1],Angle2[-1],Angle3[-1]),'0.75')
	plt.plot(X,Y,'r')
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
	ax.set_title('Trial Number: ' + str(trialnumber))
	ax.set_aspect('equal', 'datalim')

	if PlotAngles == True:
		plt.figure()
		plt.plot(Time,Angle1)
		plt.ylim((-15*np.pi/180, 145*np.pi/180))
		plt.xticks([0, Time[-1]],['Start', 'Release'])
		plt.yticks([-10*np.pi/180, 140*np.pi/180],['-10$^\circ$','140$^\circ$'])
		ax1 = plt.gca()
		ax1.spines['left'].set_position('zero')
		ax1.spines['right'].set_color('none')
		ax1.spines['top'].set_color('none')
		ax1.xaxis.set_ticks_position('bottom')
		ax1.yaxis.set_ticks_position('left')
		ax1.set_title('Trial Number: ' + str(trialnumber))
		plt.ylabel('Shoulder Angle')

		plt.figure()
		plt.plot(Time,Angle2)
		plt.ylim((75*np.pi/180, 140*np.pi/180))
		plt.xticks([0, Time[-1]],['Start', 'Release'])
		plt.yticks([80*np.pi/180, 135*np.pi/180],['80$^\circ$','135$^\circ$'])
		ax2 = plt.gca()
		ax2.spines['left'].set_position('zero')
		ax2.spines['right'].set_color('none')
		ax2.spines['top'].set_color('none')
		ax2.xaxis.set_ticks_position('bottom')
		ax2.yaxis.set_ticks_position('left')
		ax2.set_title('Trial Number: ' + str(trialnumber))
		plt.ylabel('Elbow Angle')

		plt.figure()
		plt.plot(Time,Angle3)
		plt.ylim((-95*np.pi/180, 5*np.pi/180))
		plt.xticks([0, Time[-1]],['Start', 'Release'])
		plt.yticks([-90*np.pi/180, 0*np.pi/180],['-90$^\circ$','0$^\circ$'])
		ax3 = plt.gca()
		ax3.spines['left'].set_position('zero')
		ax3.spines['right'].set_color('none')
		ax3.spines['top'].set_color('none')
		ax3.xaxis.set_ticks_position('bottom')
		ax3.yaxis.set_ticks_position('left')
		ax3.set_title('Trial Number: ' + str(trialnumber))
		plt.ylabel('Wrist Angle')

	def normalized_muscle_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		#R_tranpose Column 1
		r1DELTa = 19
		r1CB = 20
		r1DELTp = -8
		r1BIC = 15
		r1TRI = -15
		r1BRA = 0
		r1BRD = 0
		r1PRO = 0
		r1FCR = 0
		r1ECRB = 0
		r1ECRL = 0
		r1FCU = 0
		r1FDS = 0
		r1PL = 0
		r1ECU = 0
		r1EDM = 0
		r1EDC = 0
		r1APL = 0

		#R_tranpose Column 2
		r2DELTa = 0
		r2CB = 0
		r2DELTp = 0
		r2BIC = lambda Angle2: (24777*np.pi)/10000 + (1288265228720957*Angle2)/35184372088832 \
								- (2429*np.pi*Angle2)/125 + (68251*np.pi*(Angle2**2))/5000 \
								- (10427*np.pi*(Angle2**3))/5000 + (20571*(Angle2**2))/10000 \
								- (14043*(Angle2**3))/2500 + 84533/10000
		r2TRI = lambda Angle2: - (8759*(Angle2**3))/5000 + (93509*(Angle2**2))/10000 \
								- (88691*Angle2)/10000 - 863614486669217/35184372088832
		r2BRA = lambda Angle2: - (12667*(Angle2**3))/2000 + (30689*(Angle2**2))/1250 \
								- (4544779416463265*Angle2)/281474976710656 + 1139910323808397/70368744177664
		r2BRD = lambda Angle2: (28129*np.pi)/10000 - (23671*Angle2)/2000 - (57781*np.pi*Angle2)/10000 \
								+ (3629*np.pi*(Angle2**2))/1250 - (197*np.pi*(Angle2**3))/500 \
								+ (24636921970321*(Angle2**2))/549755813888 - (33739*(Angle2**3))/2500 + 38141/2500
		r2PRO = lambda Angle2: (3933*np.pi)/10000 - (10079*Angle2)/10000 - (13103*np.pi*Angle2)/1250 \
								+ (2597831076304493*np.pi*(Angle2**2))/70368744177664 + (2202*np.pi**2*Angle2)/625 \
								- (93111*np.pi*(Angle2**3))/2500 + (72987*np.pi*(Angle2**4))/5000 \
								- (20089*np.pi*(Angle2**5))/10000 - (4369*np.pi**2)/10000 \
								- (6847666938421497*(Angle2**2))/562949953421312 + (53151*(Angle2**3))/2500 \
								- (5503*(Angle2**4))/500 + (8763*(Angle2**5))/5000 \
								- (1466808324885735*np.pi**2*(Angle2**2))/140737488355328 + (51333*np.pi**2*(Angle2**3))/5000 \
								- (39919*np.pi**2*(Angle2**4))/10000 + (273*np.pi**2*(Angle2**5))/500 + 22081/2000
		r2FCR = 14
		r2ECRB = lambda Angle2: (8199*np.pi)/5000 + (44637*Angle2)/2500 - (5073*np.pi*Angle2)/10000 \
								- (471*np.pi*(Angle2**2))/5000 - (28827*(Angle2**2))/10000 - 1407/125
		r2ECRL = lambda Angle2: (74361*np.pi)/10000 + (72089699777459*Angle2)/4398046511104 \
								- (8783*np.pi*Angle2)/5000 + (371*np.pi**2*Angle2)/5000 - (1667*np.pi**2)/1250 \
								- 38517/5000
		r2FCU = 19
		r2FDS = 20
		r2PL = 25
		r2ECU = -23
		r2EDM = -10
		r2EDC = -20
		r2APL = 0

		#R_tranpose Column 3
		r3DELTa = 0
		r3CB = 0
		r3DELTp = 0
		r3BIC = 0
		r3TRI = 0
		r3BRA = 0
		r3BRD = 0
		r3PRO = 0
		r3FCR = lambda Angle3: (3199*Angle3)/2000 + 3301/250
		r3ECRB = lambda Angle3: (21411*Angle3)/10000 - 7562500789275879/562949953421312
		r3ECRL = lambda Angle3: (457*Angle3)/200 - 58583/5000
		r3FCU = lambda Angle3: (13307*(Angle3**2))/10000 + (1869*Angle3)/400 + 1578328710658497/140737488355328
		r3FDS = lambda Angle3: (2099*(Angle3**2))/2000 + (10641*Angle3)/10000 + 5824674283064289/562949953421312
		r3PL = lambda Angle3: (5011*(Angle3**2))/10000 + (13821*Angle3)/10000 + 3749/400
		r3ECU = lambda Angle3: (3883*Angle3)/1250 - 21289/2500
		r3EDM = lambda Angle3: (7603*Angle3)/2500 - 7791/625
		r3EDC = lambda Angle3: (693*Angle3)/400 - 35319/2500
		r3APL = lambda Angle3: 1171/2000 - (171*(Angle3**2))/2000 - (73*Angle3)/2500

		OptimalMuscleLength = np.array([(9.8),     (9.3),      (13.7),      (11.6),      (13.4),    (8.6),   \
	                      				(17.3),    (4.9),       (6.3),       (5.9),       (8.1),    (5.1),   \
	                       				(8.4),     (6.4),       (6.2),       (6.8),        (7.),    (7.1)] )
		def muscle_1_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1DELTa  									\
						- Angle2Spline.pp_deriv(Time)*r2DELTa  								\
						-Angle3Spline.pp_deriv(Time)*r3DELTa)  								\
	         	        /(10*OptimalMuscleLength[0])
			return(velocity)
		def muscle_2_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1CB  								\
	                    -Angle2Spline.pp_deriv(Time)*r2CB  									\
	                    -Angle3Spline.pp_deriv(Time)*r3CB)  								\
	                    /(10*OptimalMuscleLength[1])
			return(velocity)
		def muscle_3_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1DELTp  								\
	                    -Angle2Spline.pp_deriv(Time)*r2DELTp  								\
	                    -Angle3Spline.pp_deriv(Time)*r3DELTp)  								\
	                    /(10*OptimalMuscleLength[2])
			return(velocity)
		def muscle_4_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1BIC  								\
	                    -Angle2Spline.pp_deriv(Time)*r2BIC(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3BIC)  								\
	                    /(10*OptimalMuscleLength[3])
			return(velocity)
		def muscle_5_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1TRI  								\
	                    -Angle2Spline.pp_deriv(Time)*r2TRI(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3TRI)  								\
	                    /(10*OptimalMuscleLength[4])
			return(velocity)
		def muscle_6_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1BRA  								\
	                    -Angle2Spline.pp_deriv(Time)*r2BRA(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3BRA)  								\
	                    /(10*OptimalMuscleLength[5])
			return(velocity)
		def muscle_7_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1BRD  								\
	                    -Angle2Spline.pp_deriv(Time)*r2BRD(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3BRD)  								\
	                    /(10*OptimalMuscleLength[6])
			return(velocity)
		def muscle_8_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1PRO  								\
	                    -Angle2Spline.pp_deriv(Time)*r2PRO(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3PRO)  								\
	                    /(10*OptimalMuscleLength[7])
			return(velocity)
		def muscle_9_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1FCR  								\
	                    -Angle2Spline.pp_deriv(Time)*r2FCR  								\
	                    -Angle3Spline.pp_deriv(Time)*r3FCR(Angle3Spline.pp_func(Time)))  	\
	                    /(10*OptimalMuscleLength[8])
			return(velocity)
		def muscle_10_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1ECRB  								\
	                    -Angle2Spline.pp_deriv(Time)*r2ECRB(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3ECRB(Angle3Spline.pp_func(Time)))  	\
	                    /(10*OptimalMuscleLength[9])
			return(velocity)
		def muscle_11_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1ECRL  								\
	                    -Angle2Spline.pp_deriv(Time)*r2ECRL(Angle2Spline.pp_func(Time))  	\
	                    -Angle3Spline.pp_deriv(Time)*r3ECRL(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[10])
			return(velocity)
		def muscle_12_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1FCU  								\
	                    -Angle2Spline.pp_deriv(Time)*r2FCU  								\
	                    -Angle3Spline.pp_deriv(Time)*r3FCU(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[11])
			return(velocity)
		def muscle_13_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1FDS  								\
	                    -Angle2Spline.pp_deriv(Time)*r2FDS  								\
	                    -Angle3Spline.pp_deriv(Time)*r3FDS(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[12])
			return(velocity)
		def muscle_14_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1PL  									\
	                    -Angle2Spline.pp_deriv(Time)*r2PL  									\
	                    -Angle3Spline.pp_deriv(Time)*r3PL(Angle3Spline.pp_func(Time))) 		\
	                    /(10*OptimalMuscleLength[13])
			return(velocity)
		def muscle_15_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1ECU  								\
	                    -Angle2Spline.pp_deriv(Time)*r2ECU  								\
	                    -Angle3Spline.pp_deriv(Time)*r3ECU(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[14])
			return(velocity)
		def muscle_16_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1EDM  								\
	                    -Angle2Spline.pp_deriv(Time)*r2EDM  								\
	                    -Angle3Spline.pp_deriv(Time)*r3EDM(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[15])
			return(velocity)
		def muscle_17_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1EDC  								\
	                    -Angle2Spline.pp_deriv(Time)*r2EDC  								\
	                    -Angle3Spline.pp_deriv(Time)*r3EDC(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[16])
			return(velocity)
		def muscle_18_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
			velocity = (-Angle1Spline.pp_deriv(Time)*r1APL  								\
	                    -Angle2Spline.pp_deriv(Time)*r2APL  								\
	                    -Angle3Spline.pp_deriv(Time)*r3APL(Angle3Spline.pp_func(Time))) 	\
	                    /(10*OptimalMuscleLength[17])
			return(velocity)	
		NormalizedMuscleVelocity = np.array([ 	muscle_1_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_2_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_3_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_4_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_5_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_6_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_7_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_8_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_9_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_10_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_11_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_12_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_13_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_14_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_15_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_16_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_17_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
												muscle_18_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time) 	],	\
				                                float)
		return(NormalizedMuscleVelocity)
	def plot_normalized_muscle_velocity(Angle1Splines,Angle2Splines,Angle3Splines,Time):
		NormalizedMuscleVelocity = normalized_muscle_velocity(Angle1Splines,Angle2Splines,Angle3Splines,Time)
		plt.figure()
		plt.plot(Time,NormalizedMuscleVelocity.T,color = 'b')
		ax4 = plt.gca()
		ax4.spines['left'].set_position('zero')
		ax4.spines['right'].set_color('none')
		ax4.spines['top'].set_color('none')
		ax4.xaxis.set_ticks_position('bottom')
		ax4.yaxis.set_ticks_position('left')
		plt.xticks([0, Time[-1]],['Start', 'Release'])
		plt.ylabel('$\hat{v}_{m}}$')

		plt.figure()
		plt.plot(Time,NormalizedMuscleVelocity.T,color = 'b')
		ax4 = plt.gca()
		ax4.spines['left'].set_position('zero')
		ax4.spines['right'].set_color('none')
		ax4.spines['top'].set_color('none')
		ax4.xaxis.set_ticks_position('bottom')
		ax4.yaxis.set_ticks_position('left')
		plt.xticks([0, Time[-1]],['Start', 'Release'])
		plt.ylabel('$\hat{v}_{m}}$')
		plt.ylim(-5,5)

		
	if PlotVm == True: plot_normalized_muscle_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time) 
	
	plt.show()
			
