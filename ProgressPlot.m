classdef ProgressPlot < handle
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        curIdx = 1;
        hFig
        hAxis1
        hAxis2
        hLine1
        hLine2
        hLine3
        hLine4
        hLine5
        hLeg
    end
    
    methods
        % Constructor
        function obj = ProgressPlot()
            obj.hFig = figure();
            set(obj.hFig, 'color', 'w')
            set(obj.hFig, 'units', 'inches')
            set(obj.hFig, 'position', [1 2.5 9 4])
            
            
            subplot(1,2,1);
            obj.hAxis1 = get(gca);
            obj.hLine1 = plot(nan(1e4,1), nan(1e4,1), '.b-');
            hold on
            obj.hLine2 = plot(nan(1e4,1), nan(1e4,1), '.r-');
            hold off
            obj.hAxis1.XLim = [0 5];
            obj.hAxis1.YLim = [0 1];
            xlabel('Iteration')
            ylabel('f(x)')
            grid on
            obj.hLeg = legend('Train','Validate','location','northeast');
            
            subplot(1,2,2);
            
            obj.hAxis2 = get(gca);
            obj.hLine3 = plot(nan(1e4,1), nan(1e4,1), '.b-');
            hold on;
            obj.hLine4 = plot(nan(1e4,1), nan(1e4,1), '.g-');
            obj.hLine5 = plot(nan(1e4,1), nan(1e4,1), '.r-');
            hold off
            
            obj.hAxis2.XLim = [0 5];
            xlabel('Iteration')
            ylabel('R^2')
            legend('Surge', 'Heave', 'Pitch', 'location', 'south')
            grid on
            ylim([0 1])
        end
        
        function update(obj, iter, f1, f2, rsq1, rsq2, rsq3)
            % add data to plot
            obj.hLine1.XData(obj.curIdx) = iter;
            obj.hLine1.YData(obj.curIdx) = f1;
            obj.hLine2.XData(obj.curIdx) = iter;
            obj.hLine2.YData(obj.curIdx) = f2;
            
            obj.hLine3.XData(obj.curIdx) = iter;
            obj.hLine4.XData(obj.curIdx) = iter;
            obj.hLine5.XData(obj.curIdx) = iter;
            obj.hLine3.YData(obj.curIdx) = rsq1;
            obj.hLine4.YData(obj.curIdx) = rsq2;
            obj.hLine5.YData(obj.curIdx) = rsq3;
            
            % if its the first time, set the y axis properly
            if obj.curIdx == 1
                obj.hAxis1.YLim(2) = (0.01*ceil(110*max([f1 f2])));
            end
            % if for some reason f grows, expand the y axis
            if max([f1 f2]) > obj.hAxis1.YLim(2)
                obj.hAxis1.YLim(2) = 1.1*ceil(f);
            end
            
            % expand x range if needed
            if iter > obj.hAxis1.XLim(2), 
                obj.hAxis1.XLim(2) = ceil(1.2*obj.hAxis1.XLim(2));
                obj.hAxis2.XLim(2) = ceil(1.2*obj.hAxis1.XLim(2));
            end

            drawnow()
            obj.curIdx = obj.curIdx + 1;
        end
    end
    
end

