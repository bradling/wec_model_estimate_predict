classdef DisturbanceEstimator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        nStates
        type % linear or EKF
        parameters % if EKF, what parameters are being estimated
        m
        k
        Q
        R
        b      % Only needed if not being estimated
        feFreq %   "     "    "  "    "       " 
        
    end
    
    properties (SetAccess = public)
        dt = 0.5;
    end
    
    methods
        function obj = DisturbanceEstimator(WecSystemModel, Q, R, parameters, parmValues)
            % second input is a cell array of parameters to be estimated.
            % Possible values are {'feFreq', 'dampCoeff'}. If it is empty
            % or not provided then a linear kalman filter is used.
            
            % TODO - better input error checking
				
				if ~isa(WecSystemModel, 'WecSystemModel')
					error('input object is not a WecSystemModel')
				end
            
            % input handling
            if isstruct(parmValues)
                fields = fieldnames(parmValues);
                for ii = 1:length(fields)
                    switch fields{ii}
                        case 'b'
                            obj.b = parmValues.(fields{ii});
                        case 'feFreq'
                            obj.feFreq = parmValues.(fields{ii});
                        otherwise
                            error('bad parameter value')
                    end
                end
            end %if

            % get parameters that are to be estimated
            errorCase = and(~isempty(parameters), ...
                any(~ismember(parameters, {'feFreq', 'dampCoeff'} )));
            if errorCase == true
                error('bad estimated parameter specification');
            end
            clear errorCase
            obj.parameters = parameters;
            obj.nStates = 4 + length(obj.parameters);
            if obj.nStates == 4
                obj.type = 'linear';
                %obj.b = 1000;
                %obj.feFreq = 2*pi/10;
            else
                obj.type = 'ekf';
            end
            
            obj.m = WecSystemModel.mass + WecSystemModel.Ainf;
            obj.k = WecSystemModel.kHyd;
            
            % need to check that these sizes are right
            obj.Q = Q;
            obj.R = R;
            
            
        end %DisturbanceEstimator
        
        function [estimates, varargout] = calc_estimation(obj, measurements, initEst, P0)
            % The core function. Makes estimations given the model and
            % everything else.
            if strcmp(obj.type, 'linear') == true
                [xHat, xp] = kalman_filter(obj, measurements, initEst, P0);
                estimates.zHat    = xHat(1,:)';
                estimates.zDotHat = xHat(2,:)';
                estimates.feHat   = xHat(3,:)';
            else % then it must be ekf
                [xHat, xp] = ekf(obj, measurements, initEst, P0);
                estimates.zHat    = xHat(1,:)';
                estimates.zDotHat = xHat(2,:)';
                estimates.feHat   = xHat(3,:)';
                
                myCase = sprintf('%i%i', ismember({'dampCoeff', 'feFreq'}, obj.parameters));
                switch myCase
                    case '10'
                        estimates.dampCoeff = xHat(5,:)';
                    case '01'
                        estimates.feFreq    = xHat(5,:)';
                    case '11'
                        estimates.dampCoeff = xHat(6,:)';
                        estimates.feFreq    = xHat(5,:)';
                end%switch
                
                
            end % if
            if nargout == 2
                varargout{1} = xp;
            end
        end % calc_estimation
        
    end % public methods
    
    methods (Access = private)
        
        % Note: Algorithm comes from http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
        function [xHat, xp] = kalman_filter(obj, y, initEst, P0)
            if size(y,2) < size(y,1)
                y = y';
            end
            A = calc_a(obj, zeros(4,1));
            C = calc_c(obj);
            I = eye(obj.nStates);
            
            xHat = nan(obj.nStates, size(y,2));
            xp   = nan(obj.nStates, size(y,2));
            xHat(:,1) = initEst;
            xp(:,1)   = initEst;
            P = P0;
            for ii = 1:size(y,2)-1
                xp(:,ii+1) = A*xHat(:,ii);  
                P = A*P*A' + obj.Q;
                K = P*C' / (C*P*C' + obj.R);
                xHat(:, ii+1) = xp(:,ii+1) + K*(y(:,ii+1) - C*xp(:,ii+1));  
                P = (I - K*C)*P;
            end
            
        end %kalman_filter
        
        % Note: Algorithm comes from http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
        function [xHat, xp] = ekf(obj, y, initEst, P0)
            if size(y,2) < size(y,1)
                y = y';
            end
            
            C = calc_c(obj);
            I = eye(obj.nStates);            
            
            xHat = nan(obj.nStates, size(y,2));
            xp   = nan(obj.nStates, size(y,2));
            xHat(:,1) = initEst;
            xp(:,1)   = initEst;
            P = P0;
            for ii = 1:size(y,2)-1
                [A, F] = calc_a(obj, xHat(:,ii));
                xp(:,ii+1) = F * xHat(:,ii);
                P = A*P*A'+ obj.Q;
                K = P*C' / (C*P*C' + obj.R);
                xHat(:,ii+1) = xp(:,ii+1) + K*(y(:,ii+1) - C*xp(:,ii+1));
                P = (I - K*C)*P;
            end %loop
            
            
        end %ekf
        
        
        function C = calc_c(obj)
            C = [eye(2) zeros(2,obj.nStates - 2)];
        end
        
        
        [A, varargout] = calc_a(obj, estState)
        
    end % private methods
    
end

