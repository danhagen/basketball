function [Count] = HistogramCount(NumberOfBins)
%Generates 2D histogram count for AllSumOfSquareHistogramData.m. Input the
%number of bins desired (default is 20).

    % Check number of inputs.
    if nargin > 1
        error('myfuns:PlotXY:TooManyInputs', ...
            'requires at most 1 input');
    end

    % Fill in unset optional values.
    switch nargin
        case 0
            NumberOfBins = 20;
    end
    
    load('AllSumOfSquareHistogramData.mat');

    XUpperBounds = 0:500/NumberOfBins:500;
    YUpperBounds = 0:500/NumberOfBins:500;

    BinNumberX = zeros(100000,1);
    BinNumberY = zeros(100000,1); 
    for i = 1:length(BinNumberX) 
        BinNumberX(i) = sum(AllSumOfMaxEccContractionVelocitySquared(i)>XUpperBounds); 
        BinNumberY(i) = sum(AllSumOfMaxConcContractionVelocitySquared(i)>YUpperBounds); 
    end
    Count = zeros(length(XUpperBounds),length(YUpperBounds));
    for i = 1:length(XUpperBounds)
        for j = 1:length(YUpperBounds)
            for k = 1:100000
                if BinNumberX(k)==i && BinNumberY(k)==j
                    Count(j,i) = Count(j,i) + 1;
                end
            end
        end
    end
end


