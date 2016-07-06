
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

FinalPositionInX = position_in_x(Angle1Final,Angle2Final,Angle3Final)
FinalPositionInY = position_in_y(Angle1Final,Angle2Final,Angle3Final)

def displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth):
	""" 
	Utilizes the final endpoint position in x, the release angle, the shot depth and known values of the players
	height, ball radius, and hoop radius to calculate the total distance the projectile will be from the basket
	in the x direction.
	"""
	FootLength = 0.152*Height
	BallRadius = 11.9
	BallDisplacementInX = BallRadius*cos(ReleaseAngle)
	HoopRadius = 22.9
	ChangeInX = ShotDepth +  FootLength - FinalPositionInX - BallDisplacementInX - HoopRadius
	return(ChangeInX/100.)
def displacement_in_y(FinalPositionInY,ReleaseAngle):
	""" 
	Utilizes the final endpoint position in y and the release angle as well as known values of the players
	height, ball radius, and hoop height to calculate the total distance the projectile will be from the basket
	in the y direction.
	"""
	BallRadius = 11.9
	BallDisplacementInY = BallRadius*sin(ReleaseAngle)
	HoopHeight = 304.8
	ChangeInY = HoopHeight - 0.87*Height - FinalPositionInY - BallDisplacementInY
	return(ChangeInY/100.)
def initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle, ShotDepth):
	"""
	Takes the final endpoint positions in x and y (initial projectile positions) as well as release angle and 
	shot depth to find the initial velocity of the projectile to make it into the basket. The following logic 
	helps explain the equation used:

	VelocityInX = Velocity*cos(ReleaseAngle)
	Time = displacement_in_x(...)/VelocityInX = displacement_in_x(...)/(Velocity*cos(ReleaseAngle))
	displacement_in_y(...) = VelocityInY*Time - (9.8/2)*Time**2
	VelocityInY/VelocityInX = [Velocity*sin(ReleaseAngle)]/[Velocity*cos(ReleaseAngle)] = tan(ReleaseAngle)
	VelocityInY*Time = VelocityInY*displacement_in_x(...)/VelocityInX = tan(ReleaseAngle)*displacement_in_x(...)
	Time**2 = displacement_in_x(...)**2/VelocityInX**2 = (1/Velocity**2)*displacement_in_x(...)/(cos(ReleaseAngle)**2)
	Velocity**2 = (9.8/2)*(displacement_in_x(...)/cos(ReleaseAngle)**2)/(displacement_in_x(...)*tan(ReleaseAngle)-displacement_in_y(...))
	"""
	Velocity = sqrt((9.8/2.)*(displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth)/cos(ReleaseAngle))**2 \
				/ (displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth)*tan(ReleaseAngle) \
					- displacement_in_y(FinalPositionInY,ReleaseAngle)))
	return(Velocity)

InitialProjectileVelocity = initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle,ShotDepth)
EndpointVelocity = [cos(ReleaseAngle)*InitialProjectileVelocity*100, sin(ReleaseAngle)*InitialProjectileVelocity*100, 0.0]

def jacobian(Angle1Final,Angle2Final,Angle3Final):
	"""
	Generates the Jacobian matrix for the geometric model of a three link planar system (i.e. 3 dof's).
	Final joint angles are then substituted to create the Jacobian matrix of the endpoint in the final release 
	posture.
	"""
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
	result = abs(B[0]-expected_slope)< 0.001
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
def is_slope_positive(x_initial,x_transition,B,C,D,X):
	"""
	Takes in the initial x values for the two piecewise polynomials made by the clamped cubic splines with one
	break and evaluates if the slope is always positive. Returns TRUE if slope is positive for the first 250 ms
	of a 550 ms shot.
	"""
	derivative = np.piecewise(X,[X <= x_transition, X > x_transition], \
										[lambda X: B[0] + 2*C[0]*(X-x_initial) + 3*D[0]*(X-x_initial)**2, \
										lambda X: B[1] + 2*C[1]*(X-x_transition) + 3*D[1]*(X-x_transition)**2])
	result = min(derivative[:2501])>=0
	return(result)
def clamped_cubic_spline(x_initial,x_final,y_initial,y_final,initial_slope,final_slope,ymin,ymax,X,**options):
	"""
	This will take in the initial and final values for both x and y, as well as the desired initial and final
	slopes and generate 1000 clamped cubic spline that produce y values that are within the bounds [ymin, ymax].
	Returns a list of arrays, each of len(X). Options allows for slope limitations on shoulder rotations such
	that the derivative of the spline is always positive to match observations (slope = "Shoulder").
	"""
	i = 0
	while i < 1000:
		x_rand,y_rand = generate_random_point(x_initial,x_final,ymin,ymax)
		A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initial_slope,final_slope)
		assert test_b_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,C,D,initial_slope), "Initial slope does not match the expected value"
		assert test_endpoint_slope(B[1],C[1],D[1],x_rand,x_final,final_slope),"Problem with Endpoint Slope"
		assert test_for_discontinuity(A[0],B[0],C[0],D[0],x_initial,x_rand,A[1]), "Jump Discontinuity at t = %f!" %x_rand
		if options.get("angle") == "Shoulder":
			options_slope_condition = is_slope_positive(x_initial,x_final,B,C,D,X)
		else:
			options_slope_condition = True

		if i == 0:
			ClampedSpline = piecewise_for_2_intervals(x_initial,x_rand,A,B,C,D,X)
			if is_within_bounds(ClampedSpline, ymin, ymax) and options_slope_condition:
				i+=1
		elif i == 1:
			NextClampedSpline = piecewise_for_2_intervals(x_initial,x_rand,A,B,C,D,X)
			if is_within_bounds(NextClampedSpline, ymin, ymax) and options_slope_condition:
				ClampedSpline = np.concatenate(([ClampedSpline], [NextClampedSpline]), axis = 0)
				i+=1
		else:
			NextClampedSpline = piecewise_for_2_intervals(x_initial,x_rand,A,B,C,D,X)
			if is_within_bounds(NextClampedSpline, ymin, ymax) and options_slope_condition:
				ClampedSpline = np.concatenate((ClampedSpline, [NextClampedSpline]), axis = 0)
				i+=1
	return(ClampedSpline)

Angle1Spline = clamped_cubic_spline(0,EndTime,Angle1Initial,Angle1Final,AngularVelocity1Initial, \
										AngularVelocity1Final,Angle1Bounds[0],Angle1Bounds[1],Time,\
										angle = "Shoulder")
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
def moment_arm_matrix(Angle1,Angle2,Angle3):
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
	r2BIC = lambda Angle2: (24777*np.pi)/10000 + (1288265228720957*Angle2)/35184372088832 - (2429*np.pi*Angle2)/125 + (68251*np.pi*(Angle2**2))/5000 - (10427*np.pi*(Angle2**3))/5000 + (20571*(Angle2**2))/10000 - (14043*(Angle2**3))/2500 + 84533/10000
	r2TRI = lambda Angle2: - (8759*(Angle2**3))/5000 + (93509*(Angle2**2))/10000 - (88691*Angle2)/10000 - 863614486669217/35184372088832
	r2BRA = lambda Angle2: - (12667*(Angle2**3))/2000 + (30689*(Angle2**2))/1250 - (4544779416463265*Angle2)/281474976710656 + 1139910323808397/70368744177664
	r2BRD = lambda Angle2: (28129*np.pi)/10000 - (23671*Angle2)/2000 - (57781*np.pi*Angle2)/10000 + (3629*np.pi*(Angle2**2))/1250 - (197*np.pi*(Angle2**3))/500 + (24636921970321*(Angle2**2))/549755813888 - (33739*(Angle2**3))/2500 + 38141/2500
	r2PRO = lambda Angle2: (3933*np.pi)/10000 - (10079*Angle2)/10000 - (13103*np.pi*Angle2)/1250 + (2597831076304493*np.pi*(Angle2**2))/70368744177664 + (2202*np.pi**2*Angle2)/625 - (93111*np.pi*(Angle2**3))/2500 + (72987*np.pi*(Angle2**4))/5000 - (20089*np.pi*(Angle2**5))/10000 - (4369*np.pi**2)/10000 - (6847666938421497*(Angle2**2))/562949953421312 + (53151*(Angle2**3))/2500 - (5503*(Angle2**4))/500 + (8763*(Angle2**5))/5000 - (1466808324885735*np.pi**2*(Angle2**2))/140737488355328 + (51333*np.pi**2*(Angle2**3))/5000 - (39919*np.pi**2*(Angle2**4))/10000 + (273*np.pi**2*(Angle2**5))/500 + 22081/2000
	r2FCR = 14
	r2ECRB = lambda Angle2: (8199*np.pi)/5000 + (44637*Angle2)/2500 - (5073*np.pi*Angle2)/10000 - (471*np.pi*(Angle2**2))/5000 - (28827*(Angle2**2))/10000 - 1407/125
	r2ECRL = lambda Angle2: (74361*np.pi)/10000 + (72089699777459*Angle2)/4398046511104 - (8783*np.pi*Angle2)/5000 + (371*np.pi**2*Angle2)/5000 - (1667*np.pi**2)/1250 - 38517/5000
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

	R_Matrix = np.matrix([	[r1DELTa, 		r1CB, 			r1DELTp, 		r1BIC, 			r1TRI, 			r1BRA, 	\
							r1BRD, 			r1PRO, 			r1FCR, 			r1ECRB, 		r1ECRL, 		r1FCU, 	\
							r1FDS, 			r1PL,   		r1ECU, 			r1EDM, 			r1EDC, 			r1APL], \
							[r2DELTa, 		r2CB, 			r2DELTp, 		r2BIC(Angle2), 	r2TRI(Angle2), 	r2BRA(Angle2),  \
							r2BRD(Angle2), 	r2PRO(Angle2), 	r2FCR,			r2ECRB(Angle2), r2ECRL(Angle2), r2FCU,  \
							r2FDS, 			r2PL, 			r2ECU, 			r2EDM, 			r2EDC, 			r2APL], \
							[r3DELTa, 		r3CB, 			r3DELTp, 		r3BIC, 			r3TRI, 			r3BRA, \
							r3BRD, 			r3PRO, 			r3FCR(Angle3), 	r3ECRB(Angle3), r3ECRL(Angle3), r3FCU(Angle3), \
							r3FDS(Angle3), 	r3PL(Angle3), 	r3ECU(Angle3), 	r3EDM(Angle3), 	r3EDC(Angle3), 	r3APL(Angle3)]	], \
							float)
	return(R_Matrix)

	print(moment_arm_matrix(0,0,0))


#plot_spline_results(Time,Angle1Spline)
#plot_spline_results(Time,Angle2Spline)
#plot_spline_results(Time,Angle3Spline)














