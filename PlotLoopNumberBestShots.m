clear all
clc

load('.\LoopNumberTrials\LoopNumber1.mat','Time');
NumberOfLoops = 100;
NumberOfTrials = 1000;
AllBestEccentricShots = zeros(2,length(Time),NumberOfLoops);
AllBestConcentricShots = zeros(2,length(Time),NumberOfLoops );
AllBestShots = zeros(2,length(Time),NumberOfLoops);

for i = 1:NumberOfLoops
    FileName = ['.\LoopNumberTrials\LoopNumber' num2str(i) '.mat'];
    load(FileName,'Angle1SplineStructures','Angle2SplineStructures','Angle3SplineStructures',...
                    'BestEccentricShot','BestConcentricShot','BestShot');
    [BestEccPositionInX, BestEccPositionInY] = PlotXY(BestEccentricShot+(i-1)*1000);
    [BestConcPositionInX, BestConcPositionInY] = PlotXY(BestConcentricShot+(i-1)*1000);
    [BestPositionInX, BestPositionInY] = PlotXY(BestEccentricShot+(i-1)*1000);
    AllBestEccentricShots(:,:,i) = [BestEccPositionInX;BestEccPositionInY];
    AllBestConcentricShots(:,:,i) = [BestConcPositionInX;BestConcPositionInY];
    AllBestShots(:,:,i) = [BestPositionInX;BestPositionInY];
    clear('Angle1SplineStructures','Angle2SplineStructures','Angle3SplineStructures',...
            'BestEccentricShot','BestConcentricShot','BestShot');
end

AllBestEccentricShotsInX = permute(AllBestEccentricShots(1,:,:),[2,3,1]);
AllBestEccentricShotsInY = permute(AllBestEccentricShots(2,:,:),[2,3,1]);

figure('Name','Best Eccentric Shots');
    plot(AllBestEccentricShotsInX,AllBestEccentricShotsInY);
    xlabel('x (cm)','FontName','AvantGarde');
    ylabel('y (cm)','FontName','AvantGarde');
    set(gca,'TickDir',      'out',...
            'TickLength',   [0.02 0.02],...
            'XTick',        -20:20:60,...
            'YTick',        -60:20:60,...
            'XColor',       [.3 .3 .3],...
            'YColor',       [.3 .3 .3],...
            'LineWidth',        1,...
            'box',          'off');
    axis equal;
    set(gcf,'PaperPositionMode','auto');
    
AllBestConcentricShotsInX = permute(AllBestConcentricShots(1,:,:),[2,3,1]);
AllBestConcentricShotsInY = permute(AllBestConcentricShots(2,:,:),[2,3,1]);
    
figure('Name','Best Concentric Shots');
    plot(AllBestConcentricShotsInX,AllBestConcentricShotsInY);
    xlabel('x (cm)','FontName','AvantGarde');
    ylabel('y (cm)','FontName','AvantGarde');
    set(gca,'TickDir',      'out',...
            'TickLength',   [0.02 0.02],...
            'XTick',        -20:20:60,...
            'YTick',        -60:20:60,...
            'XColor',       [.3 .3 .3],...
            'YColor',       [.3 .3 .3],...
            'LineWidth',        1,...
            'box',          'off');
    axis equal;
    set(gcf,'PaperPositionMode','auto');

AllBestShotsInX = permute(AllBestShots(1,:,:),[2,3,1]);
AllBestShotsInY = permute(AllBestShots(2,:,:),[2,3,1]);    
    
figure('Name','Best Overall Shots');
    plot(AllBestShotsInX,AllBestShotsInY);
    xlabel('x (cm)','FontName','AvantGarde');
    ylabel('y (cm)','FontName','AvantGarde');
    set(gca,'TickDir',      'out',...
            'TickLength',   [0.02 0.02],...
            'XTick',        -20:20:60,...
            'YTick',        -60:20:60,...
            'XColor',       [.3 .3 .3],...
            'YColor',       [.3 .3 .3],...
            'LineWidth',        1,...
            'box',          'off');
    axis equal;
    set(gcf,'PaperPositionMode','auto');

save('AllBestShots.mat','AllBestEccentricShots','AllBestConcentricShots',...
        'AllBestShots');

        

