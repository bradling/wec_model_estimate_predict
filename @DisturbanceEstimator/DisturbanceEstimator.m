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
                            error('bad parameter')
                    end
                end
            end %if

            
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
        
        function estimates = calc_estimation(obj, measurements, initEst, P0)
            % The core function. Makes estimations given the model and
            % everything else.
            if strcmp(obj.type, 'linear') == true
                xHat = kalman_filter(obj, measurements, initEst, P0);
                estimates.zHat    = xHat(1,:)';
                estimates.zDotHat = xHat(2,:)';
                estimates.feHat   = xHat(3,:)';
            else % then it must be ekf
                xHat = ekf(obj, measurements, initEst, P0);
                estimates.zHat    = xHat(1,:)';
                estimates.zDotHat = xHat(2,:)';
                estimates.feHat   = xHat(3,:)';
                switch ismember({'dampCoeff', 'feFreq'}, obj.parameters)
                    case [1 0]
                        estimates.dampCoeff = xHat(5,:)';
                    case [0 1]
                        estimates.feFreq    = xHat(5,:)';
                    case [1 1]
                        estimates.dampCoeff = xHat(6,:)';
                        estimates.feFreq    = xHat(5,:)';
                end%switch
            end % if
        end % calc_estimation
        
    end % public methods
    
    methods (Access = private)
        
        % Note: Algorithm comes from http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
        function xHat = kalman_filter(obj, y, initEst, P0)
            if size(y,2) < size(y,1)
                y = y';
            end
            A = calc_a(obj, zeros(4,1));
            C = calc_c(obj);
            I = eye(obj.nStates);
            
            xHat = nan(obj.nStates, size(y,2));
            xHat(:,1) = initEst;
            P = P0;
            for ii = 1:size(y,2)-1
                xt = A*xHat(:,ii);
                P = A*P*A' + obj.Q;
                K = P*C' / (C*P*C' + obj.R);
                xHat(:, ii+1) = xt + K*(y(:,ii) - C*xt);
                P = (I - K*C)*P;
            end
            
        end %kalman_filter
        
        % Note: Algorithm comes from http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
        function xHat = ekf(obj, y, initEst, P0)
            if size(y,2) < size(y,1)
                y = y';
            end
            
            C = calc_c(obj);
            I = eye(obj.nStates);            
            
            xHat = nan(obj.nStates, size(y,2));
            xHat(:,1) = initEst;
            P = P0;
            for ii = 1:size(y,2)-1
                [A, F] = calc_a(obj, xHat(:,ii));
                
                xt = F * xHat(:,ii);
                P = A*P*A'+ obj.Q;
                K = Pp*C' / (C*P*C' + obj.R);
                xHat(:,ii+1) = xt + K*(y(:,ii) - C*xt);
                P = (I - K*C)*P;
            end %loop
            
        end %ekf
        
        
        function C = calc_c(obj)
            C = [eye(2) zeros(2,obj.nStates - 2)];
        end
        
        
        [A, varargout] = calc_a(obj, estState)
        
    end % private methods
    
end

