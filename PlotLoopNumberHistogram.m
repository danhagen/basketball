clear all
clc

NumberOfLoops = 100;
NumberOfTrials = 1000;
AllSumOfMaxEccContractionVelocitySquared = zeros(NumberOfLoops*NumberOfTrials,1);
AllSumOfMaxConcContractionVelocitySquared = zeros(NumberOfLoops*NumberOfTrials,1);

for i = 1:NumberOfLoops
    FileName = ['.\LoopNumberTrials\LoopNumber' num2str(i) '.mat'];
    load(FileName,'Angle1SplineStructures','Angle2SplineStructures',...
        'Angle3SplineStructures','Time');
    NormalizedMuscleVelocity = NormalizedMomentArmMatrix(Angle1SplineStructures,Angle2SplineStructures,Angle3SplineStructures,Time);
    [SumOfMaxEccContractionVelocitySquared,SumOfMaxConcContractionVelocitySquared] = ...
                SumOfMaximumContractionVelocitiesSquared(NormalizedMuscleVelocity);
    AllSumOfMaxEccContractionVelocitySquared((i-1)*NumberOfTrials+1:i*NumberOfTrials,1)=SumOfMaxEccContractionVelocitySquared;
    AllSumOfMaxConcContractionVelocitySquared((i-1)*NumberOfTrials+1:i*NumberOfTrials,1)=SumOfMaxConcContractionVelocitySquared;
    clear('Angle1SplineStructures','Angle2SplineStructures',...
        'Angle3SplineStructures','Time');
end

figure; 
    subplot(2,1,1);
        bar(hist(AllSumOfMaxEccContractionVelocitySquared,30),'b'); 
        title('Histogram of ALL Eccentric Trials Sum of Squares');
        ylabel('Frequency');
        xlabel('Sum of Squares Value');
    subplot(2,1,2);
        bar(hist(AllSumOfMaxConcContractionVelocitySquared,30),'r'); 
        title('Histogram of ALL Concentric Trials Sum of Squares');
        ylabel('Frequency');
        xlabel('Sum of Squares Value');
        
figure;
    hist3([AllSumOfMaxEccContractionVelocitySquared AllSumOfMaxConcContractionVelocitySquared],[100 100]);
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    xlabel('Sum of Maximum Eccentric Contraction Velocities Squared');
    ylabel('Sum of Maximum Concentric Contraction Velocities Squared');
    zlabel('Frequency');
    axis off
    
save('AllSumOfSquareHistogramData.mat','AllSumOfMaxEccContractionVelocitySquared',...
    'AllSumOfMaxConcContractionVelocitySquared');
    

        

