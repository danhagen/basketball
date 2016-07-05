function [ AngleSplines, AngleSplineStructures ] = CubicSpline( AngleConditions, MaximumAngle, MinimumAngle, Time, NumberOfTrials )
%   Calculates the Physiologically Feasible Cubic Spline Interpolations for
%   For NumberOfTrials given AngleConditions (Initial Angular Velocity,   
%   Initial Angle, Final Angle, Final Angular Velocity) and Time. Angles 
%   must be in radians. If end angular velocity conditions are not
%   specified, function will assume zero velocities at both endpoints.
%   Updated 02/03/16 5:00 PM

if length(AngleConditions) == 2
    AngleConditions = [0 AngleConditions(1) AngleConditions(2) 0];
else
end

EndTime = Time(end);
InitialAngularVelocity = AngleConditions(1);
InitialAngle = AngleConditions(2);
FinalAngle = AngleConditions(3);
FinalAngularVelocity = AngleConditions(4);

AngleSplines = zeros(NumberOfTrials, length(Time));
AngleSplineStructures = struct('form','pp',...
                                'breaks',num2cell(zeros(NumberOfTrials,1)),...
                                'coefs',num2cell(zeros(NumberOfTrials,1)),...
                                'pieces',num2cell(ones(NumberOfTrials,1)),...
                                'order',num2cell(ones(NumberOfTrials,1)),...
                                'dim',num2cell(zeros(NumberOfTrials,1)));
                            
I = 1:NumberOfTrials;

rng('shuffle');
RandomMiddleTime = EndTime.*rand(NumberOfTrials,1);
RandomMiddleAngle = MinimumAngle +(MaximumAngle-MinimumAngle).*rand(NumberOfTrials,1);

for i = 1:NumberOfTrials 
    AngleSplineStructures(i) = csape([0 RandomMiddleTime(i) EndTime], ...
        [InitialAngularVelocity InitialAngle RandomMiddleAngle(i) FinalAngle FinalAngularVelocity]);
    AngleSplines(i,:) = ppval(AngleSplineStructures(i), Time);
    while (max(AngleSplines(i,:))>MaximumAngle) || (min(AngleSplines(i,:))<MinimumAngle)
        rng('shuffle');
        RandomMiddleTime(i) = EndTime.*rand();
        RandomMiddleAngle(i) = MinimumAngle +(MaximumAngle-MinimumAngle).*rand();
        
        AngleSplineStructures(i) = csape([0 RandomMiddleTime(i) EndTime], ...
            [InitialAngularVelocity InitialAngle RandomMiddleAngle(i) FinalAngle FinalAngularVelocity]);
        AngleSplines(i,:) = ppval(AngleSplineStructures(i), Time);
    end 
end

% % Filtering out Splines greater than MaximumAngle
% I = I(max(AngleSplines,[],2)<=MaximumAngle);
% AngleSplines = AngleSplines(max(AngleSplines,[],2)<=MaximumAngle,:);
% 
% % Filtering out Splines less than MinimumAngle
% I = I(min(AngleSplines,[],2)>=MinimumAngle);
% AngleSplines = AngleSplines(min(AngleSplines,[],2)>=MinimumAngle,:);
% 
% AngleSplineStructures = AngleSplineStructures(I);
end

