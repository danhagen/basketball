function [ EccentricMuscleVelocitiesLessThan5MuscleLengthsPerSecond, LessThan5MuscleLengthsPerSecondTrialIndex, ...
            EccentricMuscleVelocitiesMoreThan5MuscleLengthsPerSecond, MoreThan5MuscleLengthsPerSecondTrialIndex ] ...
            = FilteringEccentricTrialsBy5MuscleLengthsPerSecondCriterion( NormalizedMuscleVelocity, CutOffTime )
%	Inputs a NumberOfMuscle x TimeLength x NumberOfTrials Matrix and the 
%   CutOffTime and returns a matrix of Trials that meet the 5 muscle
%   lengths per second criterion and a matrix of Trials that fail to meet
%   the criterion as well as the index value. CutOff time is used to eliminate late stage eccentric
%   contractions evident in all trajectories (in ms).

    CutOffTime=CutOffTime*10;
    Index = 1:size(NormalizedMuscleVelocity,2);
    EccentricMuscleVelocitiesLessThan5MuscleLengthsPerSecond = NormalizedMuscleVelocity(:,:,sum(max(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)>5)==0);
    LessThan5MuscleLengthsPerSecondTrialIndex = Index(sum(max(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)>5)==0);
    EccentricMuscleVelocitiesMoreThan5MuscleLengthsPerSecond = NormalizedMuscleVelocity(:,:,sum(max(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)>5)~=0);
    MoreThan5MuscleLengthsPerSecondTrialIndex = Index(sum(max(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)>5)~=0);
end

