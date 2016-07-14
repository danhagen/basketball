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
Angle1Conditions = [AngularVelocity1Initial, Angle1Initial, Angle1Final, AngularVelocity1Final]

Angle2Bounds = [0, DegreeToRadianFactor*135]
AngularVelocity2Initial = 0
AngularVelocity2Final = AngularVelocities[1]
Angle2Conditions = [AngularVelocity1Initial, Angle1Initial, Angle1Final, AngularVelocity1Final]

Angle3Bounds = [DegreeToRadianFactor*-90, 0]
AngularVelocity3Initial = 0
AngularVelocity3Final = AngularVelocities[2]
Angle3Conditions = [AngularVelocity3Initial, Angle3Initial, Angle3Final, AngularVelocity3Final]

def c_coefficients(RandomMiddleAngle,RandomMiddleTime,AngleConditions):
	C = np.matrix([	[2*(RandomMiddleTime-0), 	(RandomMiddleTime-0), 		0							],   \
					[(RandomMiddleTime-0), 		2*(EndTime-0), 				(EndTime-RandomMiddleTime)	],   \
					[0,							(EndTime-RandomMiddleTime),	2*(EndTime-RandomMiddleTime)] ], \
					float)
	y = np.matrix([	[3*(RandomMiddleAngle-AngleConditions[1])/(RandomMiddleTime-0) - 3*AngleConditions[0]],  \
					[3*(AngleConditions[2]-RandomMiddleAngle)/(EndTime-RandomMiddleTime) - 3*(RandomMiddleAngle-AngleConditions[1])/(RandomMiddleTime-0)],  \
					[3*AngleConditions[3] - 3*(AngleConditions[2]-RandomMiddleAngle)/(EndTime-RandomMiddleTime)] ],  \
					float)
	CCoefficients = np.dot(inv(C),y)
	c_coefficients_float = np.array(CCoefficients).astype(float)
	ipdb.set_trace()
	return(c_coefficients_float)

def d_coefficients(RandomMiddleTime,CCoefficients):
	DCoefficients = np.array([	[(CCoefficients[1]-CCoefficients[0])/(3*(RandomMiddleTime-0))],  \
								[(CCoefficients[2]-CCoefficients[1])/(3*(EndTime-RandomMiddleTime))] ],  \
								float)
	return(DCoefficients)

def b_coefficients(RandomMiddleTime,RandomMiddleAngle,AngleConditions,CCoefficients,DCoefficients):
	BCoefficients = np.array([((RandomMiddleAngle-AngleConditions[1])/(RandomMiddleTime-0)-CCoefficients[0]*(RandomMiddleTime-0) - DCoefficients[0]*((RandomMiddleTime-0)**2)),  \
							((AngleConditions[2]-RandomMiddleAngle)/(EndTime-RandomMiddleTime)-CCoefficients[1]*(EndTime-RandomMiddleTime) - DCoefficients[1]*((EndTime-RandomMiddleTime)**2)) ]).astype(float)
	return(BCoefficients)

def a_coefficients(RandomMiddleAngle, AngleConditions):
	ACoefficients = np.array([	np.float(AngleConditions[1]),    \
								RandomMiddleAngle  ])
	return(ACoefficients)

def cubic_spline(AngleBounds, AngleConditions, Time):
	RandomMiddleTime = np.random.uniform(0,EndTime)
	RandomMiddleAngle = np.random.uniform(AngleBounds[0],AngleBounds[1])
	C = c_coefficients(RandomMiddleTime,RandomMiddleAngle,AngleConditions)
	D = d_coefficients(RandomMiddleTime,C)
	B = b_coefficients(RandomMiddleTime,RandomMiddleAngle,AngleConditions,C,D)
	A = a_coefficients(RandomMiddleAngle, AngleConditions)
	AngleSpline = np.piecewise(Time,[Time < RandomMiddleTime, Time >= RandomMiddleTime], \
					[lambda Time: A[0] + B[0]*(Time-0) + C[0]*(Time-0)**2 + D[0]*(Time-0)**3, \
					lambda Time: A[1] + B[1]*(Time-RandomMiddleTime) + C[1]*(Time-RandomMiddleTime)**2 + D[1]*(Time-RandomMiddleTime)**3])
	return(AngleSpline,B)

Angle1Spline,B = cubic_spline(Angle1Bounds,Angle1Conditions,Time)
plt.plot(Time,Angle1Spline)
plt.show()











