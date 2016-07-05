function [ AngularVelocities ] = InverseJacobian( EndpointVelocity,FinalAngles,ShoulderToElbowLength,ForearmLength,HandLength )
%   Inputs EndpointVelocity vector, FinalAngles vector,
%   ShoulderToElbowLength, ForearmLength, and HandLength and returns the
%   AngularVelocities associated with that EndpointVelocity by means of 
%   the InverseJacobianMatrix operation. Angles must be in radians and all
%   lengths must be in meters. EndpointVelocity Vector should be 3x1 with 0
%   as the third entry.
%   Modified 02/04/16 by Dan Hagen

syms Angle1 Angle2 Angle3 PositionInX PositionInY AngleAlpha  
syms SymbolicGeometricModel SymbolicJacobianMatrix SymbolicInverseJacobianMatrix

PositionInY = -ShoulderToElbowLength.*cos(Angle1)...
                -ForearmLength.*cos(Angle1 + Angle2)...
                -HandLength.*cos(Angle1 + Angle2 - Angle3);
PositionInX = ShoulderToElbowLength.*sin(Angle1)...
                +ForearmLength.*sin(Angle1 + Angle2)...
                +HandLength.*sin(Angle1 + Angle2 - Angle3); 
AngleAlpha = Angle1 + Angle2 - Angle3;

SymbolicGeometricModel = [PositionInX;PositionInY;AngleAlpha];
SymbolicJacobianMatrix = jacobian(SymbolicGeometricModel,[Angle1,Angle2,Angle3]);
SymbolicInverseJacobianMatrix = inv(SymbolicJacobianMatrix);
RealNumberInverseJacobianMatrix = subs(SymbolicInverseJacobianMatrix,[Angle1,Angle2,Angle3],[FinalAngles(1),FinalAngles(2),FinalAngles(3)]);

InverseJacobianMatrix = double(RealNumberInverseJacobianMatrix);
AngularVelocities =InverseJacobianMatrix*EndpointVelocity;

end

