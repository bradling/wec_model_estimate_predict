classdef AdaptiveKalmanDisturbanceEstimator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        origSys
        augSys
        linearSys
        nDisturbances
        Q
        R
        idxDist
        idxDistDot
        idxOmega
    end
    
    methods
        function obj = AdaptiveKalmanDisturbanceEstimator(varargin)
            % obj = AdaptiveKalmanDisturbanceEstimator(sysModel)
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
        
        function [Fa, linA] = get_a(obj, state)
            n = size(obj.origSys.a, 1);
            m = obj.nDisturbances;
            deltaDot = state(obj.idxDistDot);
            omegaT = state(obj.idxOmega);
            % calculate the nonlinear system A matrix
            f = [obj.origSys.a  obj.origSys.d   zeros(n,m)  zeros(n,m); 
                  zeros(m,n)     zeros(m,m)      eye(m,m)    zeros(m,m); 
                  zeros(m,n)    -diag(omegaT)    zeros(m,m)  zeros(m,m);
                  zeros(m,n)     zeros(m,m)      zeros(m,m)  zeros(m,m)]; 
            Fa = expm(f .* obj.augSys.dt);
            % calculate the linearized system Jacobian to use as A
            J = [obj.origSys.a  obj.origSys.d   zeros(n,m)  zeros(n,m); 
                 zeros(m,n)     zeros(m,m)      eye(m,m)    zeros(m,m); 
                 zeros(m,n)    -diag(omegaT)    zeros(m,m)  diag(deltaDot);
                 zeros(m,n)     zeros(m,m)      zeros(m,m)  zeros(m,m)];
            linA = expm(J .* obj.augSys.dt);
        end
        
        function harmonic_model(obj, dt)
            
            n = size(obj.origSys.a, 1);
            m = obj.nDisturbances;
            a = [obj.origSys.a  obj.origSys.d  zeros(n,m)  zeros(n,m); 
                 zeros(m,n)     zeros(m,m)     eye(m,m)    zeros(m,m); 
                 zeros(m,n)    -diag([1 1 1])  zeros(m,m)  zeros(m,m);
                 zeros(m,n)     zeros(m,m)     zeros(m,m)  zeros(m,m)];
            b = [obj.origSys.b ; zeros(3*m, size(obj.origSys.b, 2)) ];
            c = [obj.origSys.c zeros(size(obj.origSys.c, 1), 3*m) ];
            d = zeros(size(c,1), size(b, 2));
            
            % convert to a discrete system for use with the kalman filter.
            %sys = c2d(ss(a,b,c,d), dt, 'foh');
            obj.augSys.a = a;
            obj.augSys.b = b;
            obj.augSys.c = c;
            obj.augSys.dt = dt;
            
            obj.idxDist = (1:obj.nDisturbances) + n;
            obj.idxDistDot = obj.idxDist + obj.nDisturbances;
            obj.idxOmega = obj.idxDistDot + obj.nDisturbances;
            
            % now calc linearized system
        end

        
        function define_covariances(obj, Q, R)
            % Q is process noise covariance, R is measurement noise covariance
            
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
        
        function [results, xHat, xp, innov] = calc_estimates(obj, results, p0, x0)
            % each column of row is associated with a single time point
            y = [results.sDotNoisy; results.sNoisy];
            
            results.estimator = obj;
            
            n = size(obj.augSys.a, 1);
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
                [F, J] = obj.get_a(xHat(:,ii));
                xp(:,ii+1) = F*xHat(:,ii);  
                P = J*P*J' + obj.Q;
                K = P*C' / (C*P*C' + obj.R);
                innov(:,ii+1) = y(:,ii+1) - C*xp(:,ii+1);
                xHat(:, ii+1) = xp(:,ii+1) + K*(y(:,ii+1) - C*xp(:,ii+1));  
                P = (I - K*C)*P;
            end
            results.feHat = xHat(obj.idxDist,:);
        end
    end
    
end

