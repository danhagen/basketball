import numpy as np
from numpy import cos,sin,sqrt

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

FinalPositionInX = ShoulderToElbowLength*sin(Angle1Final) \
					+ ForearmLength*sin(Angle1Final+Angle2Final) \
					+ HandLength*sin(Angle1Final+Angle2Final-Angle3Final)
FinalPositionInY = -ShoulderToElbowLength*cos(Angle1Final) \
					- ForearmLength*cos(Angle1Final+Angle2Final) \
					- HandLength*cos(Angle1Final+Angle2Final-Angle3Final)

def initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle, Height):
	velocity = sqrt(-490.*((434.3+0.152*Height-FinalPositionInX-11.9*cos(ReleaseAngle))**2)/  \
                                ((((cos(ReleaseAngle))**2.0)*(304.8 - 0.87*Height-FinalPositionInY))  \
                                -(sin(ReleaseAngle)*cos(ReleaseAngle)*(434.3+0.152*Height-FinalPositionInX))))
	return(velocity)

InitialProjectileVelocity = initial_projectile_velocity(FinalPositionInX, FinalPositionInY, ReleaseAngle, Height)
print InitialProjectileVelocity