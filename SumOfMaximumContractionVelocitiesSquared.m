function [ SumOfMaxEccContractionVelocitySquared,SumOfMaxConcContractionVelocitySquared ] = SumOfMaximumContractionVelocitiesSquared( NormalizedMuscleVelocity )
%Returns the sum of maximum contraction velocities squared for both
%eccentrically and concentrically contracting muscles.
%NormalizedMuscleVelocity must be a NxMxP matrix where N denotes muscle, M
%denotes time, and P denotes trial number.
%Created 2/29/16. Modified 2/29/16.
if nargin > 1
    error('myfuns:SumOfMaximumContractionVelocitiesSquared:TooManyInputs', ...
        'requires 1 input');
end
if nargin == 0
    error('myfuns:SumOfMaximumContractionVelocitiesSquared:NotEnoughInputs', ...
        'requires 1 input');
end

SumOfMaxEccContractionVelocitySquared = permute(sum(max((NormalizedMuscleVelocity(:,1:5200,:)...
                                                .*double(NormalizedMuscleVelocity(:,1:5200,:)>=0))...
                                                ,[],2).^2),[3,2,1]);
SumOfMaxConcContractionVelocitySquared = permute(sum(min((NormalizedMuscleVelocity(:,1:5200,:)...
                                                .*double(NormalizedMuscleVelocity(:,1:5200,:)<=0))...
                                                ,[],2).^2),[3,2,1]);
end

