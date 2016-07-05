
import ipdb
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from numpy import sin,cos,tan,sqrt
from numpy.linalg import inv


EndTime = 0.55
ChangeInTime = 0.0001
Time = np.arange(0,EndTime+ChangeInTime,ChangeInTime,float)
#NumberOfAdditionalLoops = 99
#NumberOfTrials = 1000

PositionInX = []
PositionInY = []

HeightInInches = 71
DegreeToRadianFactor = np.pi/180

#global ShoulderToElbowLength,ForearmLength,HandLength,Height
Height = HeightInInches*2.54
ShoulderToElbowLength = 0.186*Height
ForearmLength = 0.146*Height
HandLength = 0.108*Height
ReleaseAngle = DegreeToRadianFactor*30
ShotDepth = 457.2

Angle1Initial = DegreeToRadianFactor*-10
Angle1Final = DegreeToRadianFactor*120

Angle2Initial = DegreeToRadianFactor*85
Angle2Final = DegreeToRadianFactor*84

Angle3Initial = DegreeToRadianFactor*0
Angle3Final = DegreeToRadianFactor*-35

def position_in_x(Angle1, Angle2, Angle3):
	PositionInX = ShoulderToElbowLength*sin(Angle1) \
					+ ForearmLength*sin(Angle1+Angle2) \
					+ HandLength*sin(Angle1+Angle2-Angle3)
	return(PositionInX)

def position_in_y(Angle1,Angle2,Angle3):
	PositionInY = -ShoulderToElbowLength*cos(Angle1) \
					- ForearmLength*cos(Angle1+Angle2) \
					- HandLength*cos(Angle1+Angle2-Angle3)
	return(PositionInY)

FinalPositionInX = position_in_x(Angle1Final,Angle2Final,Angle3Final)
FinalPositionInY = position_in_y(Angle1Final,Angle2Final,Angle3Final)

def displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth):
	FootLength = 0.152*Height
	BallRadius = 11.9
	BallDisplacementInX = BallRadius*cos(ReleaseAngle)
	HoopRadius = 22.9
	ChangeInX = ShotDepth +  FootLength - FinalPositionInX - BallDisplacementInX - HoopRadius
	return(ChangeInX/100.)

def displacement_in_y(FinalPositionInY,ReleaseAngle):
	BallRadius = 11.9
	BallDisplacementInY = BallRadius*sin(ReleaseAngle)
	HoopHeight = 304.8
	ChangeInY = HoopHeight - 0.87*Height - FinalPositionInY - BallDisplacementInY
	return(ChangeInY/100.)

def initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle, ShotDepth):
	# VelocityInX = Velocity*cos(ReleaseAngle)
	# Time = displacement_in_x(...)/VelocityInX = displacement_in_x(...)/(Velocity*cos(ReleaseAngle))
	# displacement_in_y(...) = VelocityInY*Time - (9.8/2)*Time**2
	# VelocityInY/VelocityInX = [Velocity*sin(ReleaseAngle)]/[Velocity*cos(ReleaseAngle)] = tan(ReleaseAngle)
	# VelocityInY*Time = VelocityInY*displacement_in_x(...)/VelocityInX = tan(ReleaseAngle)*displacement_in_x(...)
	# Time**2 = displacement_in_x(...)**2/VelocityInX**2 = (1/Velocity**2)*displacement_in_x(...)/(cos(ReleaseAngle)**2)
	# Velocity**2 = (9.8/2)*(displacement_in_x(...)/cos(ReleaseAngle)**2)/(displacement_in_x(...)*tan(ReleaseAngle)-displacement_in_y(...))
	Velocity = sqrt((9.8/2.)*(displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth)/cos(ReleaseAngle))**2 \
				/ (displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth)*tan(ReleaseAngle) \
					- displacement_in_y(FinalPositionInY,ReleaseAngle)))
	return(Velocity)

InitialProjectileVelocity = initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle,ShotDepth)
EndpointVelocity = [cos(ReleaseAngle)*InitialProjectileVelocity*100, sin(ReleaseAngle)*InitialProjectileVelocity*100, 0.0]

def jacobian(Angle1Final,Angle2Final,Angle3Final):
	Angle1,Angle2,Angle3 = sp.symbols('Angle1,Angle2,Angle3',real=True)
	x = ShoulderToElbowLength*sp.sin(Angle1) \
		+ ForearmLength*sp.sin(Angle1+Angle2) \
		+ HandLength*sp.sin(Angle1+Angle2-Angle3)
	y = -ShoulderToElbowLength*sp.cos(Angle1) \
		- ForearmLength*sp.cos(Angle1+Angle2) \
		- HandLength*sp.cos(Angle1+Angle2-Angle3)
	alpha = Angle1 + Angle2 - Angle3

	GeometricModel = sp.Matrix([x,y,alpha])
	SymbolicJacobianMatrix = GeometricModel.jacobian([Angle1,Angle2,Angle3])
	JacobianMatrix = SymbolicJacobianMatrix.subs([(Angle1,Angle1Final), (Angle2,Angle2Final), (Angle3,Angle3Final)]).evalf()
	return(np.array(JacobianMatrix).astype(float))

JacobianMatrix = jacobian(Angle1Final,Angle2Final,Angle3Final)
AngularVelocities = np.dot(inv(jacobian(Angle1Final,Angle2Final,Angle3Final)),EndpointVelocity)

Angle1Bounds = [Angle1Initial, DegreeToRadianFactor*140]
AngularVelocity1Initial = 0
AngularVelocity1Final = AngularVelocities[0]

Angle2Bounds = [0, DegreeToRadianFactor*135]
AngularVelocity2Initial = 0
AngularVelocity2Final = AngularVelocities[1]

Angle3Bounds = [DegreeToRadianFactor*-90, 0]
AngularVelocity3Initial = 0
AngularVelocity3Final = AngularVelocities[2]

def c_matrix(x1,x2,x3):
	"""
	Takes in the values of x1, x2, and x3 to create the C matrix needed to find the coefficients of a clamped
	cubic spline with only one break (i.e. Cx = y, where x is an array of c coefficients for the
	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns a matrix.
	"""
	C = np.matrix([	[	2*(x2-x1), 		(x2-x1), 			0			],   \
					[	(x2-x1), 		2*(x3-x1), 		(x3-x2)		],   \
					[	0,				(x3-x2),		2*(x3-x2)	] 	], \
					float)
	return(C)
def y_vector(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
	"""
	Takes in the values of (x1,y1), (x2,y2), and (x3,y3) to create the y array necessary for the clamped cubic
	spline matrix manipulation for one break only (i.e. Cx = y, where x is an array of c coefficients for the
	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns an array.

	"""
	y = np.array([	[	3*(y2-y1)/(x2-x1) - 3*initial_slope 				],  \
					[	3*(y3-y2)/(x3-x2) - 3*(y2-y1)/(x2-x1) 	],  \
					[	3*final_slope - 3*(y3-y2)/(x3-x2) 				] 	],  \
					float)
	return(y)
def c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
	"""
	Using matrix manipulations the equation Cx = y necessary for the c coefficients for a clamped cubic spline 
	with only one break (i.e. Cx = y, where x is an array of c coefficients for the piecewise polynomial 
	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) can be rearranged such that x = C.T*y. The values
	(x1,y1), (x2,y2), and (x3,y3) are the three points needed to the spline and initial_slope and final_slope
	are the endpoint conditions. Returns an array.
	"""
	C = c_matrix(x1,x2,x3)
	y = y_vector(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
	CCoefficients = np.dot(inv(C),np.array(y).astype(float))
	return(np.array(CCoefficients).astype(float))
def d_coefficients(x1,x2,x3,CCoefficients):
	"""
	Uses the c coefficients and the values of x1, x2, and x3 to find the d coefficients for the	piecewise 
	polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an array with
	three elements. Returns an array.
	"""
	DCoefficients = np.array([	[	(CCoefficients[1]-CCoefficients[0])/(3*(x2-x1))	],  \
								[	(CCoefficients[2]-CCoefficients[1])/(3*(x3-x2))	] ],  \
								float)
	return(DCoefficients)
def b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients):
	"""
	Uses the c and d coefficients and the values of (x1,y1), (x2,y2), and (x3,y3) to find the b coefficients for 
	the	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an 
	array with two or more elements and DCoefficients must be an array with two elements. Returns an array.
	"""
	BCoefficients = np.array([	((y2-y1)/(x2-x1)-CCoefficients[0]*(x2-x1) - DCoefficients[0]*((x2-x1)**2)),  \
								((y3-y2)/(x3-x2)-CCoefficients[1]*(x3-x2) - DCoefficients[1]*((x3-x2)**2)) 	]).astype(float)
	return(BCoefficients)
def test_b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients,expected_slope):
	"""
	Tests to make sure that the generated b coefficients match the expected slope. Uses the c and d coefficients 
	and the values of (x1,y1), (x2,y2), and (x3,y3) to find the b coefficients for the	piecewise polynomial 
	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an array with two or more 
	elements and DCoefficients must be an array with two elements. Returns TRUE if expected_slope equals b.
	"""
	B = b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients)
	result = B[0] == expected_slope
	return(result)
	assert B[0]==expected_slope, "First b coefficient (%f) does not equal initial slope (%f)." (B[0],expected_slope)
def a_coefficients(y1,y2):
	"""
	Uses the y values of (x1,y1) and (x2,y2) to find the a coefficients for the	piecewise polynomial equation 
	y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns an array.
	"""
	ACoefficients = np.array([	y1,    \
								y2  ]).astype(float)
	return(ACoefficients)
def spline_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
	"""
	Uses the values of (x1,y1), (x2,y2), and (x3,y3) to find the coefficients for the piecewise polynomial 
	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) for a clamped cubic spline with one break only. 
	Returns coefficient arrays A, B, C,and D.
	"""
	C = c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
	D = d_coefficients(x1,x2,x3,C)
	B = b_coefficients(x1,x2,x3,y1,y2,y3,C,D)
	A = a_coefficients(y1,y2)
	return(A,B,C,D)
def generate_random_point(xmin,xmax,ymin,ymax):
	"""
	Generates a random point in Cartesian space such that x is in [xmin, xmax] and y is in [ymin,ymax]. 
	Number is chosen at random from a uniform distribution. Returns two floats.
	"""
	x_rand = np.random.uniform(xmin,xmax)
	y_rand = np.random.uniform(ymin,ymax)
	return(x_rand,y_rand)
def test_endpoint_slope(b,c,d,x_n_minus_1,x_n,expected_slope):
	"""
	Takes in the cubic spline coefficients for the derivative of y = a + b*(x-x_n_minus_1) + c*(x-x_n_minus_1)**2 + d*(x-x_n_minus_1)**3 
	(y' = b + 2*c*(x-x_n_minus_1) + 3*d*(x-x_n_minus_1)**2)	for the last piecewise polynomial and tests to see if the expected slope at 
	the endpoint is equal to the actual	endpoint slope. The variable x_n_minus_1 is the initial value of the final piecewise polynomial 
	and x_n is the final data point. Returns TRUE if they are equal.

	"""
	actual_slope = b + 2*c*(x_n-x_n_minus_1) + 3*d*(x_n-x_n_minus_1)**2
	result = abs(actual_slope-expected_slope)<0.001
	return(result)
def test_for_discontinuity(a_n,b_n,c_n,d_n,x_n,x_n_plus_1,y_n_plus_1):
	"""
	Takes in the coefficients for a cubic spline polynomial y = a_n + b_n*(x-x_n) + c_n*(x-x_n)**2 + d_n*(x-x_n)**3
	and tests to see if the final y value for this piecewise polynomial is equal to the initial y value of the next 
	piecewise polynomial (i.e. when x = x_n_plus_1). The variable x_n is the initial x value of the preceding 
	polynomial, and x_n_plus_1 is the transition value from one polynomial to the next. y_n_plus_1 is the initial y
	value for the next piecewise polynomial. 
	"""
	y_n_final = a_n + b_n*(x_n_plus_1-x_n) + c_n*(x_n_plus_1-x_n)**2 + d_n*(x_n_plus_1-x_n)**3
	result = abs(y_n_final-y_n_plus_1)<0.001
	return(result)
def piecewise_for_2_intervals(x_initial,x_transition,A,B,C,D,X):
	"""
	This is to generate a piecewise polynomial array for a cubic spline that has one break ONLY. A, B, C, and D 
	must each have two elements corresponding to the coefficients of each polynomial. x_initial is the initial 
	x value for the first polynomial and x_transition is the final x value for the first polynomial AND the 
	initial x value for the second polynomial (thus being the transition x value). X is a 1D array of the
	independent variable. Returns an array with len(result)=len(X).
	"""
	result = np.piecewise(X,[X <= x_transition, X > x_transition], \
										[lambda X: A[0] + B[0]*(X-x_initial) + C[0]*(X-x_initial)**2 + D[0]*(X-x_initial)**3, \
										lambda X: A[1] + B[1]*(X-x_transition) + C[1]*(X-x_transition)**2 + D[1]*(X-x_transition)**3])
	return(result)
def is_within_bounds(Spline,ymin,ymax):
	"""
	This takes in a 1D Spline array and tests to see if the values are within the allowed bounds [ymin, ymax].
	Returns TRUE if all values are within bounds.

	"""
	result = max(Spline)<=ymax and min(Spline)>=ymin
	return(result)
def clamped_cubic_spline(x_initial,x_final,y_initial,y_final,initial_slope,final_slope,ymin,ymax,X):
	"""
	This will take in the initial and final values for both x and y, as well as the desired initial and final
	slopes and generate 1000 clamped cubic spline that produce y values that are within the bounds [ymin, ymax].
	Returns a list of arrays, each of len(X).
	"""
	i = 0
	while i < 1000:
		x_rand,y_rand = generate_random_point(x_initial,x_final,ymin,ymax)
		A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initial_slope,final_slope)
		assert test_b_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,C,D,initial_slope)
		assert test_endpoint_slope(B[1],C[1],D[1],x_rand,x_final,final_slope),"Problem with Endpoint Slope"
		assert test_for_discontinuity(A[0],B[0],C[0],D[0],x_initial,x_rand,A[1]), "Jump Discontinuity at t = %f!" %x_rand
		if i == 0:
			ClampedSpline = piecewise_for_2_intervals(x_initial,x_rand,A,B,C,D,X)
			if is_within_bounds(ClampedSpline, ymin, ymax):
				i+=1
		elif i == 1:
			NextClampedSpline = piecewise_for_2_intervals(x_initial,x_rand,A,B,C,D,X)
			if is_within_bounds(NextClampedSpline, ymin, ymax):
				ClampedSpline = np.concatenate(([ClampedSpline], [NextClampedSpline]), axis = 0)
				i+=1
		else:
			NextClampedSpline = piecewise_for_2_intervals(x_initial,x_rand,A,B,C,D,X)
			if is_within_bounds(NextClampedSpline, ymin, ymax):
				ClampedSpline = np.concatenate((ClampedSpline, [NextClampedSpline]), axis = 0)
				i+=1
	return(ClampedSpline)

Angle1Spline = clamped_cubic_spline(0,EndTime,Angle1Initial,Angle1Final,AngularVelocity1Initial, \
										AngularVelocity1Final,Angle1Bounds[0],Angle1Bounds[1],Time)
Angle2Spline = clamped_cubic_spline(0,EndTime,Angle2Initial,Angle2Final,AngularVelocity2Initial, \
										AngularVelocity2Final,Angle2Bounds[0],Angle2Bounds[1],Time)
Angle3Spline = clamped_cubic_spline(0,EndTime,Angle3Initial,Angle3Final,AngularVelocity3Initial, \
										AngularVelocity3Final,Angle3Bounds[0],Angle3Bounds[1],Time)

def plot_spline_results(Time,AngleSplines):
	"""
	For a 1D vector Time and an N-D list of arrays AngleSplines, this will plot all calculated splines on the
	same figure. Must close final figure in order to run next line of code.
	"""
	plt.figure()
	for i in range(0,AngleSplines.shape[0]):
		plt.plot(Time,AngleSplines[i])
	plt.show()

#plot_spline_results(Time,Angle1Spline)
#plot_spline_results(Time,Angle2Spline)
#plot_spline_results(Time,Angle3Spline)














