function [ ConcentricMuscleVelocitiesMoreThanNeg5MuscleLengthsPerSecond, MoreThanNeg5MuscleLengthsPerSecondTrialIndex, ...
            ConcentricMuscleVelocitiesLessThanNeg5MuscleLengthsPerSecond, LessThanNeg5MuscleLengthsPerSecondTrialIndex ] ...
            = FilteringConcentricTrialsBy5MuscleLengthsPerSecondCriterion( NormalizedMuscleVelocity, CutOffTime )
%	Inputs a NumberOfMuscle x TimeLength x NumberOfTrials Matrix and the 
%   CutOffTime and returns a matrix of Trials that meet the -5 muscle
%   lengths per second criterion and a matrix of Trials that fail to meet
%   the criterion as well as the index value. CutOff time is used to
%   eliminate late stage concentric contractions evident in all
%   trajectories (in ms).

    CutOffTime=CutOffTime*10;
    Index = 1:size(NormalizedMuscleVelocity,2);
    ConcentricMuscleVelocitiesMoreThanNeg5MuscleLengthsPerSecond = NormalizedMuscleVelocity(:,:,sum(min(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)<-5)==0);
    MoreThanNeg5MuscleLengthsPerSecondTrialIndex = Index(sum(min(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)<-5)==0);
    ConcentricMuscleVelocitiesLessThanNeg5MuscleLengthsPerSecond = NormalizedMuscleVelocity(:,:,sum(min(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)<-5)~=0);
    LessThanNeg5MuscleLengthsPerSecondTrialIndex = Index(sum(min(NormalizedMuscleVelocity(:,1:CutOffTime,:),[],2)<-5)~=0);
end


