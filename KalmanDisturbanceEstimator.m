classdef KalmanDisturbanceEstimator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        origSys
        augSys
        nDisturbances
        feFreq
        type = 'undefined'
        Q
        R
    end
    
    methods
        function obj = KalmanDisturbanceEstimator(varargin)
            % obj = KalmanDisturbanceEstimator(sysModel)
            %  sysModel is a structure defining the continuous state space 
            %  model used to make the predictions in the kalman filter.
            %  the state space model is inthe form:
            %    x(k+1) = a x(k) + b u(k) + d w(k)
            %    y(k)   = c x(k)
            %
            % sysModel has fields:
            % a, b, c, d
            %  NOTE: d maps the disturbance to the state.
            %        d DOES NOT map the input to the output. See form above
            if nargin > 0
                obj.origSys = varargin{1};
                obj.nDisturbances = size(obj.origSys.d, 2);
            end
        end
        
        function harmonic_model(obj, feFreq, dt)
            if and(numel(feFreq) ~= 1, numel(feFreq) ~= obj.nDisturbances)
                error('Incorrect number of assumed frequencies')
            end
            if and(numel(feFreq) == 1, obj.nDisturbances > 1) 
                obj.feFreq = repmat(feFreq, 1, obj.nDisturbances);
            else
                obj.feFreq = feFreq(:)';
            end
            obj.type = 'harmonic';
            
            n = size(obj.origSys.a, 1);
            m = obj.nDisturbances;
            a = [obj.origSys.a  obj.origSys.d        zeros(n,m) ; 
                 zeros(m,n)     zeros(m,m)           eye(m,m)   ; 
                 zeros(m,n)    -diag(obj.feFreq.^2)  zeros(m,m) ];
            b = [obj.origSys.b ; zeros(2*m, size(obj.origSys.b, 2)) ];
            c = [obj.origSys.c zeros(size(obj.origSys.c, 1), 2*m) ];
            d = zeros(size(c,1), size(b, 2));
            
            % convert to a discrete system for use with the kalman filter.
            sys = c2d(ss(a,b,c,d), dt, 'foh');
            obj.augSys.a = sys.a;
            obj.augSys.b = sys.b;
            obj.augSys.c = sys.c;
            obj.augSys.dt = dt;
        end
        
        function define_covariances(obj, Q, R)
            % Q is process noise covariance, R is measurement noise covariance
            if strcmp(obj.type, 'undefined'), 
                error('Need to define disturbance model first'); 
            end
            
            szQ = size(Q);
            szR = size(R);
            nStates = size(obj.augSys.a,1);
            nMeasure = size(obj.augSys.c,1);
            if or(szQ(1) ~= nStates, szQ(2) ~= nStates)
                error('Q is incorrect size');
            elseif or(szR(1) ~= nMeasure, szR(2) ~= nMeasure)
                error('R is incorrect size');
            end
            obj.Q = Q;
            obj.R = R;
        end
        
        function [xHat, xp, innov] = calc_estimates(obj, y, p0, x0)
            % each column of row is associated with a single time point
            n = size(obj.augSys.a, 1);
            A = obj.augSys.a;
            C = obj.augSys.c;
            I = eye(n);
            
            xHat = nan(n, size(y,2));
            xp   = nan(n, size(y,2));
            innov = nan(size(y,1), size(y,2));
            xHat(:,1) = x0;
            xp(:,1)   = x0;
            innov(:,1) = y(:,1) - C*x0;
            P = p0;
            for ii = 1:size(y,2)-1
                xp(:,ii+1) = A*xHat(:,ii);  
                P = A*P*A' + obj.Q;
                K = P*C' / (C*P*C' + obj.R);
                innov(:,ii+1) = y(:,ii+1) - C*xp(:,ii+1);
                xHat(:, ii+1) = xp(:,ii+1) + K*(y(:,ii+1) - C*xp(:,ii+1));  
                P = (I - K*C)*P;
            end
        end
    end
    
end

