function [NewSplineStructures] = LimitingSlope( SplineStructures, Time, SlopeLimitation )
%   Takes a SplineStructures, Time, and SlopeLimitation and returns only  
%   those splines that meet the slope limitations.
I = 1:length(SplineStructures);
SplineStructuresDerivative = struct('form','pp',...
                                    'breaks',num2cell(zeros(length(SplineStructures),1)),...
                                    'coefs',num2cell(zeros(length(SplineStructures),1)),...
                                    'pieces',num2cell(ones(length(SplineStructures),1)),...
                                    'order',num2cell(ones(length(SplineStructures),1)),...
                                    'dim',num2cell(zeros(length(SplineStructures),1)));
SplineDerivative = zeros(length(SplineStructures),length(Time));
for i = 1:length(SplineStructures)
    SplineStructuresDerivative(i) = ppdiff(SplineStructures(i));
    SplineDerivative(i,:) = ppval(SplineStructuresDerivative(i),Time);     
end

% SlopeLimitation is a lower bound. For upperbound switch to <= 
I = I(sum(SplineDerivative(:,1:2500)>=SlopeLimitation,2)==size(SplineDerivative(:,1:2500),2));
NewSplineStructures = SplineStructures(I,:);

end


