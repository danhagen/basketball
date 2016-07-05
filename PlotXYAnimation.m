function [] = PlotXYAnimation(TrialNumber)
%PlotXYAnimation will plot TrialNumber in XY plane and generate an .avi
% file.

    % Check number of inputs.
    if nargin > 1
        error('myfuns:PlotXYAnimation:TooManyInputs', ...
            'requires 1 input');
    else if nargin == 0
            error('myfuns:PlotXYAnimation:NotEnoughInputs', ...
                'requires 1 input');
        end
    end
    
    % Check number of outputs.
    if nargout > 1
        error('myfuns:PlotXYAnimation:TooManyOutputs', ...
            'no outputs required');
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
    XFunction = @(Angle1,Angle2,Angle3) ShoulderToElbowLength*sin(Angle1)+ForearmLength*sin(Angle1 + Angle2)+HandLength*sin(Angle1 + Angle2 - Angle3); % in cm
    YFunction = @(Angle1,Angle2,Angle3) -ShoulderToElbowLength*cos(Angle1)-ForearmLength*cos(Angle1 + Angle2)-HandLength*cos(Angle1 + Angle2 - Angle3); % in cm
    TrialX = XFunction(Angle1,Angle2,Angle3);
    TrialY = YFunction(Angle1,Angle2,Angle3);
    figure;
        Animation = VideoWriter(['TrialNumber' num2str(TrialNumber + (LoopNumber-1)*1000)  '.avi']);
        open(Animation);
        set(gca,'nextplot','replacechildren');

        for i=1:50:length(Time)
            plot(TrialX,TrialY,'r:');
            hold on
            Angle=0:0.01:2*pi;
            XCircle=5*cos(Angle);
            YCircle=5*sin(Angle);
            plot(XCircle,YCircle,'k','LineWidth',2);
            title(['Trial Number ' num2str(TrialNumber + (LoopNumber-1)*1000)],'FontName','AvantGarde','FontSize',14);
            line(   [0 ...
                    ShoulderToElbowLength*sin(Angle1(i))...
                    ShoulderToElbowLength*sin(Angle1(i))+ForearmLength*sin(Angle1(i)+Angle2(i))...
                    ShoulderToElbowLength*sin(Angle1(i))+ForearmLength*sin(Angle1(i)+Angle2(i))+HandLength*sin(Angle1(i)+Angle2(i)-Angle3(i))],...
                    [0 ...
                    -ShoulderToElbowLength*cos(Angle1(i))...
                    -ShoulderToElbowLength*cos(Angle1(i))-ForearmLength*cos(Angle1(i)+Angle2(i))...
                    -ShoulderToElbowLength*cos(Angle1(i))-ForearmLength*cos(Angle1(i)+Angle2(i))-HandLength*cos(Angle1(i)+Angle2(i)-Angle3(i))],...
                    'Color','k','LineWidth',2,'Marker','o','MarkerFaceColor','w'); 
            scatter(0,0,'ko','LineWidth',2);
            hold off
            xlim([-20 80]);
            ylim([-60 60]);
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
            Frame = getframe;
            writeVideo(Animation,Frame);
        end
        set(gcf,'PaperPositionMode','auto');
        close(Animation); 
end

