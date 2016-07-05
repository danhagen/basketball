%% Modified 03/11/16 4:21 PM

clear all
close all
clc

EndTime = .55;
ChangeInTime = .0001;
Time = 0:ChangeInTime:EndTime;
NumberOfAdditionalLoops = 99;
NumberOfTrials = 1000;

PositionInX = [];
PositionInY = [];

HeightInches = 71; % Height in in
DegreesToRadianFactor = pi/180;
Angle1Initial = DegreesToRadianFactor*-10;
Angle1Final = DegreesToRadianFactor*120;

Angle2Initial = DegreesToRadianFactor*85;
Angle2Final = DegreesToRadianFactor*84;

Angle3Initial = DegreesToRadianFactor*0;
Angle3Final = DegreesToRadianFactor*-35;

Height = HeightInches *2.54; % Height in cm
ShoulderToElbowLength = .186*Height;% Length of upper arm (shoulder to elbow) cm
ForearmLength = .146*Height; % Length of forearm cm
HandLength = .108*Height; % Length of hand

InitialAngles = [Angle1Initial Angle2Initial Angle3Initial]';
FinalAngles = [Angle1Final Angle2Final Angle3Final]';

ReleaseAngle = 30; % in degrees
ReleaseAngle = DegreesToRadianFactor*ReleaseAngle; % in radians

for LoopNumber = 1:NumberOfAdditionalLoops+1    
    
    AngularVelocities = zeros(1,3);
    
    FinalPositionInX = ShoulderToElbowLength*sin(FinalAngles(1))...
                        +ForearmLength*sin(FinalAngles(1)+FinalAngles(2))...
                        +HandLength*sin(FinalAngles(1)+FinalAngles(2)-FinalAngles(3)); % in cm
    FinalPositionInY = -ShoulderToElbowLength*cos(FinalAngles(1))...
                        -ForearmLength*cos(FinalAngles(1)+FinalAngles(2))...
                        -HandLength*cos(FinalAngles(1)+FinalAngles(2)-FinalAngles(3)); % in cm
    InitialProjectileVelocity = sqrt(-490.*((434.3+0.152*Height-FinalPositionInX-11.9*cos(ReleaseAngle))^2)/...
                                ((((cos(ReleaseAngle))^2)*(304.8 - 0.87*Height-FinalPositionInY))...
                                -(sin(ReleaseAngle)*cos(ReleaseAngle)*(434.3+0.152*Height-FinalPositionInX)))); % in cm/s
    EndpointVelocity = [cos(ReleaseAngle)*InitialProjectileVelocity sin(ReleaseAngle)*InitialProjectileVelocity 0]'; % cm/s
    AngularVelocities(1,:) = InverseJacobian(EndpointVelocity,FinalAngles,...
                            ShoulderToElbowLength,ForearmLength,HandLength);
    
    Angle1Maximum = DegreesToRadianFactor*140;
    Angle1Minimum = Angle1Initial;
    AngularVelocity1Final = AngularVelocities(1);
    AngularVelocity1Initial = 0;
    Angle1Conditions = [AngularVelocity1Initial Angle1Initial Angle1Final AngularVelocity1Final];
    
    Angle2Maximum = DegreesToRadianFactor*135;
    Angle2Minimum = DegreesToRadianFactor*80;
    AngularVelocity2Final = AngularVelocities(2);
    AngularVelocity2Initial = 0;
    Angle2Conditions = [AngularVelocity2Initial Angle2Initial Angle2Final AngularVelocity2Final];

    Angle3Maximum = DegreesToRadianFactor*0;
    Angle3Minimum = DegreesToRadianFactor*-90;
    AngularVelocity3Final = AngularVelocities(3);
    AngularVelocity3Initial = 0;
    Angle3Conditions = [AngularVelocity3Initial Angle3Initial Angle3Final AngularVelocity3Final];

    [~, Angle1SplineStructures] = CubicSpline(Angle1Conditions, Angle1Maximum, Angle1Minimum, Time, NumberOfTrials);
    Angle1SplineStructures = LimitingSlope(Angle1SplineStructures, Time, 0);
    [~, Angle2SplineStructures] = CubicSpline(Angle2Conditions, Angle2Maximum, Angle2Minimum, Time, NumberOfTrials);
    [~, Angle3SplineStructures] = CubicSpline(Angle3Conditions, Angle3Maximum, Angle3Minimum, Time, NumberOfTrials);

    rng('shuffle');
    TrialOrder = randi([1 length(Angle1SplineStructures)],1,NumberOfTrials);
 
    Angle1SplineStructures = Angle1SplineStructures(TrialOrder(1,:));

    NormalizedMuscleVelocity = NormalizedMomentArmMatrix(Angle1SplineStructures,Angle2SplineStructures,Angle3SplineStructures,Time);
    
    [~, GoodEccentricContractionTrials, ~, BadEccentricContractionTrials] = FilteringEccentricTrialsBy5MuscleLengthsPerSecondCriterion(NormalizedMuscleVelocity,520);
    [~, GoodConcentricContractionTrials, ~, BadConcentricContractionTrials] = FilteringConcentricTrialsBy5MuscleLengthsPerSecondCriterion(NormalizedMuscleVelocity,520);
    
    [MaxEccContractionVelocity,MaxConcContractionVelocity] = ...
                    MaximumContractionVelocities(NormalizedMuscleVelocity);
    [SumOfMaxEccContractionVelocitySquared,SumOfMaxConcContractionVelocitySquared] = ...
                    SumOfMaximumContractionVelocitiesSquared(NormalizedMuscleVelocity);

    GoodMaxEccContractionVelocity = MaxEccContractionVelocity(GoodEccentricContractionTrials);
    BadMaxEccContractionVelocity = MaxEccContractionVelocity(BadEccentricContractionTrials);
    
    GoodSumOfMaxEccContractionVelocitySquared = SumOfMaxEccContractionVelocitySquared(GoodEccentricContractionTrials);
    BadSumOfMaxEccContractionVelocitySquared = SumOfMaxEccContractionVelocitySquared(BadEccentricContractionTrials);
    
    GoodMaxConcContractionVelocity = MaxConcContractionVelocity(GoodConcentricContractionTrials);
    BadMaxConcContractionVelocity = MaxConcContractionVelocity(BadConcentricContractionTrials);
    
    GoodSumOfMaxConcContractionVelocitySquared = SumOfMaxConcContractionVelocitySquared(GoodConcentricContractionTrials);
    BadSumOfMaxConcContractionVelocitySquared = SumOfMaxConcContractionVelocitySquared(BadConcentricContractionTrials);
 
    BestEccentricShot = find(SumOfMaxEccContractionVelocitySquared == min(SumOfMaxEccContractionVelocitySquared)); 
    BestConcentricShot = find(SumOfMaxConcContractionVelocitySquared == min(SumOfMaxConcContractionVelocitySquared));
    BestShot = find((SumOfMaxEccContractionVelocitySquared+SumOfMaxConcContractionVelocitySquared)...
                    == min(SumOfMaxEccContractionVelocitySquared+SumOfMaxConcContractionVelocitySquared));
   
    filename = ['LoopNumber' num2str(LoopNumber) '.mat'];
    save(filename,'Angle1SplineStructures','Angle2SplineStructures','Angle3SplineStructures',...
                    'Time','BestEccentricShot','BestConcentricShot','BestShot');

    clear('FinalPositionInX','FinalPositionInY','InitialProjectileVelocity',...
        'EndpointVelocity','AngularVelocities','TrialOrder',...
        'Angle1SplineStructures','Angle2SplineStructures',...
        'Angle3SplineStructures','NormalizedMuscleVelocity',...
        'GoodEccentricContractionTrials','BadEccentricContractionTrials',...
        'MaxEccContractionVelocity','SumOfMaxEccContractionVelocitySquared',...
        'GoodMaxEccContractionVelocity','BadMaxEccContractionVelocity',...
        'GoodSumOfMaxEccContractionVelocitySquared','BadSumOfMaxEccContractionVelocitySquared',...
        'GoodConcentricContractionTrials','BadConcentricContractionTrials',...
        'MaxConcContractionVelocity','SumOfMaxConcContractionVelocitySquared',...
        'GoodMaxConcContractionVelocity','BadMaxConcContractionVelocity',...
        'GoodSumOfMaxConcContractionVelocitySquared','BadSumOfMaxConcContractionVelocitySquared',...
        'BestEccentricShot','BestConcentricShot','BestShot');
end
clear all
    
    
        
    



