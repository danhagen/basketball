import numpy as np
from numpy import sin,cos,tan,sqrt

EndTime = 0.55
ChangeInTime = 0.0001
Time = np.arange(0,EndTime+ChangeInTime,ChangeInTime,float)
#NumberOfAdditionalLoops = 99
#NumberOfTrials = 1000

PositionInX = []
PositionInY = []

HeightInInches = 71
DegreeToRadianFactor = 3.14/180

Angle1Initial = DegreeToRadianFactor*-10
Angle1Final = DegreeToRadianFactor*120

Angle2Initial = DegreeToRadianFactor*85
Angle2Final = DegreeToRadianFactor*84

Angle3Initial = DegreeToRadianFactor*0
Angle3Final = DegreeToRadianFactor*-35

Height = HeightInInches*2.54
ShoulderToElbowLength = 0.186*Height
ForearmLength = 0.146*Height
HandLength = 0.108*Height

ReleaseAngle = DegreeToRadianFactor*30
ShotDepth = 457.2

FinalPositionInX = ShoulderToElbowLength*sin(Angle1Final) \
					+ ForearmLength*sin(Angle1Final+Angle2Final) \
					+ HandLength*sin(Angle1Final+Angle2Final-Angle3Final)
FinalPositionInY = -ShoulderToElbowLength*cos(Angle1Final) \
					- ForearmLength*cos(Angle1Final+Angle2Final) \
					- HandLength*cos(Angle1Final+Angle2Final-Angle3Final)
def displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth,Height):
	FootLength = 0.152*Height
	BallRadius = 11.9
	BallDisplacementInX = BallRadius*cos(ReleaseAngle)
	HoopRadius = 22.9
	ChangeInX = ShotDepth +  FootLength - FinalPositionInX - BallDisplacementInX - HoopRadius
	return(ChangeInX/100.)
def displacement_in_y(FinalPositionInY,ReleaseAngle,Height):
	BallRadius = 11.9
	BallDisplacementInY = BallRadius*sin(ReleaseAngle)
	HoopHeight = 304.8
	ChangeInY = HoopHeight - 0.87*Height - FinalPositionInY - BallDisplacementInY
	return(ChangeInY/100.)
def initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle, Height, ShotDepth):
	# VelocityInX = Velocity*cos(ReleaseAngle)
	# Time = displacement_in_x(...)/VelocityInX = displacement_in_x(...)/(Velocity*cos(ReleaseAngle))
	# displacement_in_y(...) = VelocityInY*Time - (9.8/2)*Time**2
	# VelocityInY/VelocityInX = [Velocity*sin(ReleaseAngle)]/[Velocity*cos(ReleaseAngle)] = tan(ReleaseAngle)
	# VelocityInY*Time = VelocityInY*displacement_in_x(...)/VelocityInX = tan(ReleaseAngle)*displacement_in_x(...)
	# Time**2 = displacement_in_x(...)**2/VelocityInX**2 = (1/Velocity**2)*displacement_in_x(...)/(cos(ReleaseAngle)**2)
	# Velocity**2 = (9.8/2)*(displacement_in_x(...)/cos(ReleaseAngle)**2)/(displacement_in_x(...)*tan(ReleaseAngle)-displacement_in_y(...))
	Velocity = sqrt((9.8/2.)*(displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth,Height)/cos(ReleaseAngle))**2 \
				/ (displacement_in_x(FinalPositionInX,ReleaseAngle,ShotDepth,Height)*tan(ReleaseAngle) \
					- displacement_in_y(FinalPositionInY,ReleaseAngle,Height)))
	return(Velocity)

InitialProjectileVelocity = initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle, Height,ShotDepth)
print InitialProjectileVelocity