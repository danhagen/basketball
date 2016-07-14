function [PositionInX,PositionInY] = PlotXY(TrialNumber,ShotType)
%PlotXY will plot TrialNumber in XY plane. Can return the values for
%PositionInX and PositionInY.

    % Check number of inputs.
    if nargin > 2
        error('myfuns:PlotXY:TooManyInputs', ...
            'requires at most 2 inputs');
    end

    % Fill in unset optional values.
    switch nargin
        case 1
            if isnumeric(TrialNumber) == 0
                error('myfuns:PlotXY:NotEnoughInputs', ...
                        'Must input TrialNumber');
            end
            ShotType = '';
        case 2
            if isnumeric(ShotType) == 1
                error('myfuns:PlotXY:IncorrectFormat', ...
                        'ShotType must be a string');
            end
    end
    
    HeightInches = 71;
    Height = HeightInches*2.54;
    ShoulderToElbowLength = .186*Height;
    ForearmLength = .146*Height;
    HandLength = .108*Height;
    LoopNumber = sum(TrialNumber-[0:1000:100000]>0);
    TrialNumber = TrialNumber - (LoopNumber-1)*1000;
    NumberOfIntervals = 12;
    FileName = ['.\LoopNumberTrials\LoopNumber' num2str(LoopNumber) '.mat'];
    load(FileName,'Angle1SplineStructures','Angle2SplineStructures',...
                'Angle3SplineStructures','Time');
    Angle1 = ppval(Angle1SplineStructures(TrialNumber),Time);
    Angle2 = ppval(Angle2SplineStructures(TrialNumber),Time);
    Angle3 = ppval(Angle3SplineStructures(TrialNumber),Time);
    LengthOfIntervals = int16(length(Angle1)/NumberOfIntervals);
    XFunction = @(Angle1,Angle2,Angle3) ShoulderToElbowLength*sin(Angle1)+ForearmLength*sin(Angle1 + Angle2)+HandLength*sin(Angle1 + Angle2 - Angle3); % in cm
    YFunction = @(Angle1,Angle2,Angle3) -ShoulderToElbowLength*cos(Angle1)-ForearmLength*cos(Angle1 + Angle2)-HandLength*cos(Angle1 + Angle2 - Angle3); % in cm
    TrialX = XFunction(Angle1,Angle2,Angle3);
    TrialY = YFunction(Angle1,Angle2,Angle3);
    if nargout == 0
        figure('Name',ShotType);
            plot(TrialX,TrialY,'r'); hold on;
            Angle=0:0.01:2*pi; 
            XCircle=5*cos(Angle);
            YCircle=5*sin(Angle);
            plot(XCircle,YCircle,'k','LineWidth',2);
            title(['Trial Number ' num2str(TrialNumber + (LoopNumber-1)*1000)],'FontName','AvantGarde','FontSize',14);
            scatter(0,0,'ko','LineWidth',2);
            line(   [0 ...
                    ShoulderToElbowLength*sin(Angle1(1))...
                    ShoulderToElbowLength*sin(Angle1(1))+ForearmLength*sin(Angle1(1)+Angle2(1))...
                    ShoulderToElbowLength*sin(Angle1(1))+ForearmLength*sin(Angle1(1)+Angle2(1))+HandLength*sin(Angle1(1)+Angle2(1)-Angle3(1))],...
                    [0 ...
                    -ShoulderToElbowLength*cos(Angle1(1))...
                    -ShoulderToElbowLength*cos(Angle1(1))-ForearmLength*cos(Angle1(1)+Angle2(1))...
                    -ShoulderToElbowLength*cos(Angle1(1))-ForearmLength*cos(Angle1(1)+Angle2(1))-HandLength*cos(Angle1(1)+Angle2(1)-Angle3(1))],...
                    'Color','k','LineWidth',2);%,'Marker','o','MarkerFaceColor','w'); 
            for j = 1:NumberOfIntervals
                line(   [0 ...
                        ShoulderToElbowLength*sin(Angle1(j*LengthOfIntervals))...
                        ShoulderToElbowLength*sin(Angle1(j*LengthOfIntervals))+ForearmLength*sin(Angle1(j*LengthOfIntervals)+Angle2(j*LengthOfIntervals))...
                        ShoulderToElbowLength*sin(Angle1(j*LengthOfIntervals))+ForearmLength*sin(Angle1(j*LengthOfIntervals)+Angle2(j*LengthOfIntervals))+HandLength*sin(Angle1(j*LengthOfIntervals)+Angle2(j*LengthOfIntervals)-Angle3(j*LengthOfIntervals))],...
                        [0 ...
                        -ShoulderToElbowLength*cos(Angle1(j*LengthOfIntervals))...
                        -ShoulderToElbowLength*cos(Angle1(j*LengthOfIntervals))-ForearmLength*cos(Angle1(j*LengthOfIntervals)+Angle2(j*LengthOfIntervals))...
                        -ShoulderToElbowLength*cos(Angle1(j*LengthOfIntervals))-ForearmLength*cos(Angle1(j*LengthOfIntervals)+Angle2(j*LengthOfIntervals))-HandLength*cos(Angle1(j*LengthOfIntervals)+Angle2(j*LengthOfIntervals)-Angle3(j*LengthOfIntervals))],...
                        'Color','k','LineWidth',2);%,'Marker','o','MarkerFaceColor','w'); 
            end 
            scatter(0,0,'ko','LineWidth',2,'MarkerFaceColor','w');
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

    else if nargout == 2
            PositionInX = TrialX;
            PositionInY = TrialY;
        else if nargout == 1
                error('myfuns:PlotXY:NotEnoughOutputs', ...
                    'Outputs must be PositionInX and PositionInY or empty');
            else if nargout > 2
                    error('myfuns:PlotXY:TooManyOutputs', ...
                        'Outputs must be PositionInX and PositionInY or empty');
                end      
            end
        end
    end
end

