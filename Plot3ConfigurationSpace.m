function [ ] = Plot3ConfigurationSpace( GoodTrialNumber, FairTrialNumber, PoorTrialNumber )
%Plots the Configuration Space for 3 trajectories (Good, Fair, and Poor).
%Must input 3 trajectories in order to execute.

% Check number of inputs.
    if nargin > 3
        error('myfuns:Plot3ConfigurationSpace:TooManyInputs', ...
            'requires 3 inputs');
    end
    if nargin == 2
        error('myfuns:Plot3ConfigurationSpace:NotEnoughInputs', ...
            'requires 3 inputs or default 0 inputs');
    end
    if nargin == 1
        error('myfuns:Plot3ConfigurationSpace:NotEnoughInputs', ...
            'requires 3 inputs or default 0 inputs');
    end
    if nargin == 0
        GoodTrialNumber = 90773;
    end
        
    TrialNumber = [GoodTrialNumber FairTrialNumber PoorTrialNumber];
    Color = ['b' 'g' 'r'];
    figure;
        for i = 1:3
            LoopNumber = sum(TrialNumber(i)-[0:1000:100000]>0);
            TrialNumber(i) = TrialNumber(i) - (LoopNumber-1)*1000;
            FileName = ['.\LoopNumberTrials\LoopNumber' num2str(LoopNumber) '.mat'];
            load(FileName, 'Angle1SplineStructures','Angle2SplineStructures',...
                'Angle3SplineStructures','Time');
            Angle1 = ppval(Angle1SplineStructures(TrialNumber(i)),Time);
            Angle2 = ppval(Angle2SplineStructures(TrialNumber(i)),Time);
            Angle3 = ppval(Angle3SplineStructures(TrialNumber(i)),Time);
            plot3(Angle1,Angle2,Angle3,...
                'Color',        Color(i),...
                'LineWidth',    2); 
                hold on;
            scatter3([Angle1(1) Angle1(end)],...
                [Angle2(1) Angle2(end)],...
                [Angle3(1) Angle3(end)],...
                'MarkerEdgeColor',  'k',...
                'MarkerFaceColor',  'k',...        
                'LineWidth',       2);
            xlabel('PSR','FontName','Times');
            ylabel('EFE','FontName','Times');
            zlabel('WFE','FontName','Times');
            set(gca,'TickDir',      'out',...
                    'TickLength',   [0.02 0.02],...
                    'XColor',       [.3 .3 .3],...
                    'YColor',       [.3 .3 .3],...
                    'ZColor',       [.3 .3 .3],...
                    'LineWidth',    1,...
                    'box',          'off');
            grid on;
            hold on;
        end
        axis equal
        xlim([0 pi]);
        ylim([pi/2 pi]);
        zlim([-pi/2 0]);
        set(gcf,'PaperPositionMode','auto');
        hold off;
end

