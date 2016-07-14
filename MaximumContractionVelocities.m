function [ MaxEccContractionVelocity,MaxConcContractionVelocity ] = MaximumContractionVelocities( NormalizedMuscleVelocity )
%Returns maximum values for both eccentrically and concentrically
%contracting muscles. NormalizedMuscleVelocity must be a NxMxP matrix where
%N denotes muscle, M denotes time, and P denotes trial number.
%Created 2/29/16. Modified 2/29/16.
if nargin > 1
    error('myfuns:MaximumContractionVelocities:TooManyInputs', ...
        'requires 1 input');
end
if nargin == 0
    error('myfuns:MaximumContractionVelocities:NotEnoughInputs', ...
        'requires 1 input');
end
MaxEccContractionVelocity = permute(max(max((NormalizedMuscleVelocity(:,1:5200,:)...
                                    .*double(NormalizedMuscleVelocity(:,1:5200,:)>=0))...
                                    ,[],2)),[3,2,1]);

MaxConcContractionVelocity = permute(min(min((NormalizedMuscleVelocity(:,1:5200,:)...
                                    .*double(NormalizedMuscleVelocity(:,1:5200,:)<=0))...
                                    ,[],2)),[3,2,1]);
end

