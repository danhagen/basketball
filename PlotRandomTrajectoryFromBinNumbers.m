function [] = PlotRandomTrajectoryFromBinNumbers(BinNumberX,BinNumberY)

    load('AllSumOfSquareHistogramData.mat')
    XBinCenter = 6.25 + (BinNumberX-1)*12.5;
    YBinCenter = 6.25 + (BinNumberY-1)*12.5;
    I = 1:100000;
    X = zeros(1,100000); 
    Y = zeros(1,100000);
    
    for i = 1:100000 
        if abs(AllSumOfMaxEccContractionVelocitySquared(i)-XBinCenter)<=12.5 
            X(i) = 1; 
        end
        if abs(AllSumOfMaxConcContractionVelocitySquared(i)-YBinCenter)<=12.5 
            Y(i) = 1; 
        end
    end
    
    I = I(X+Y==2);
    rng('shuffle');
    RandomTrialNumber = I(randi(length(I)));
    PlotXY(RandomTrialNumber);
end

