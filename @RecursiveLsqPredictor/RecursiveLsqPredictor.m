classdef RecursiveLsqPredictor
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nIn
        nFwd
        alpha
        beta
    end
    
    methods
        function obj = RecursiveLsqPredictor(varargin)
            % predictor = RecursiveLsqPredictor(nIn, nFwd, alpha)
            if nargin > 0
                obj.nIn = varargin{1};
                obj.nFwd = varargin{2};
                obj.alpha = varargin{3};
            end
        end
        
        function results = make_predictions(obj, results, P0, beta0, varargin)
            % results = make_predictions(RecursiveLsqPredictor, results,
            %                            P0, beta0, flag)
            %
            % P0 can be a scalar or a square matrix. In the case where it
            % is a scalar we will assume that P0 is a diagonal matrix with
            % diagonal values = P0.
            % if it is a matrix it must be size (nIn x nIn), and must be a
            % positive definate symmetric matrix.
            %
            % beta0 is the initial prediction vector.
            %
            % results must be a ResultsStruct class.
            %
            % flag is an optional input. if it is set to: 'fe', true values
            % of excitation force will be used to make predictions.
            %
            % note: It is currenty set up that P0 and beta0 are the same
            % for each mode of motion
            % TODO:  allow these to be cell arrays with each entry
            % corresponding to a differnet mode of motion.
            %
 
            
            % do some input checking/handling.
            if ~isa(results, 'ResultsStruct'), error('Input not of class ResultsStruct'); end
            if numel(P0) == 1, P0 = P0.*eye(obj.nIn); end
            
            
            % initialize yHat.
            yHat = nan(size(results.s,1), size(results.s,2), obj.nFwd);
            
            beta = nan(obj.nIn, size(results.s, 1), size(results.s,2));
            
            % loop through each mode of motion
            for jj = 1:size(results.s, 1)
                P = P0;
                beta(:,jj, obj.nIn+1) = beta0;
                
                y = results.feHat(jj,:);
                if nargin == 5, 
                    if strcmp(varargin{1}, 'fe'),
                        y = results.fe(jj,:);
                    end
                end
                x = y(1:obj.nIn);
                % perform first prediction
                yHat(jj, obj.nIn+1,1) = x * beta(:,jj, obj.nIn+1);
                if obj.nFwd > 1
                    % propagate predictions forward
                    xTemp = x;
                    for kk = 2:obj.nFwd
                        xTemp = [xTemp(2:end) yHat(jj, obj.nIn+1, kk-1)];
                        yHat(jj, obj.nIn+1, kk) = xTemp * beta(:,jj, obj.nIn+1);
%                        yHat(jj, obj.nIn+1,kk) = [x(kk:obj.nIn) shiftdim(yHat(jj,obj.nIn+1,1:kk-1), 1)] * beta;
                    end
                end
                % do all the rest of the predictions.
                for ii = (obj.nIn+2):(length(y) - obj.nFwd)
                    beta(:,jj, ii) = beta(:,jj, ii-1) + (P*x') / (obj.alpha + x*P*x') * (y(ii-1) - x*beta(:,jj, ii-1));
                    % Update covariance
                    P = (1/obj.alpha) .*  ( P - (P*(x'*x)*P) ./ (obj.alpha + x*P*x') );
                    x = y(ii-obj.nIn:ii-1);
                    yHat(jj,ii,1) = x * beta(:,jj, ii);
                    if obj.nFwd > 1
                        xTemp = x;
                        for kk = 2:obj.nFwd
                            % TODO
                            % ********************************************
                            % Need a case for when nFwd > nIn
                            % ********************************************
                            %if nFwd <=
                            %yHat(jj,ii,kk) = [x(kk:obj.nIn) shiftdim(yHat(jj,ii,1:kk-1), 1)] * beta;
                            xTemp = [xTemp(2:end) yHat(jj, ii, kk-1)];
                            yHat(jj, ii, kk) = xTemp * beta(:,jj, ii);
                        end
                    end
% This is not currently implemented. It can easily be added in if I decide
% I need it
%                    % Reset P if need be
%                     if pResetFlag == true
%                         if mod(i,resetFreq(1)) == 0
%                             P = resetFreq(2) .* P;
%                         end
%                     end
                end
                
 %               obj.beta = beta;
                results.feStar = yHat;
                results.predictor = obj;
                results.userData.beta = beta;
            end
        end % make_predictions()
    end
    
end

