
import numpy as np
import sympy as sp
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

Angle1Initial = DegreeToRadianFactor*-10
Angle1Final = DegreeToRadianFactor*120

Angle2Initial = DegreeToRadianFactor*85
Angle2Final = DegreeToRadianFactor*84

Angle3Initial = DegreeToRadianFactor*0
Angle3Final = DegreeToRadianFactor*-35

#global ShoulderToElbowLength,ForearmLength,HandLength,Height
Height = HeightInInches*2.54
ShoulderToElbowLength = 0.186*Height
ForearmLength = 0.146*Height
HandLength = 0.108*Height
ReleaseAngle = DegreeToRadianFactor*30
ShotDepth = 457.2


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
#print(AngularVelocities)













