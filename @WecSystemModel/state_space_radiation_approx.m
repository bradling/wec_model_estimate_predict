function [varargout] = state_space_radiation_approx(obj, modelOrder, varargin)
% Constructs a state space approximation to the radiation impulse response
% function for the given WecSystemModel instance. Minimizes the error using
% FMINCON.
%
% [...] = state_space_radiation_approx(WecSystemModel, modelOrder)
%   will run optimixation with my chosen default values. modelOrder is the
%   number of states the approximation will be.
%
% [...] = state_space_radiation_approx(WecSystemModel, modelOrder, options)
%   options is a structure containing any possible options to pass to
%   fmincon. Additional possible fields are: "plotFlag" true/false, "X0",
%   starting point for optimization, must be a vector with length
%   2*modelOrder, and "Weight", which weighs each point in the IRF in the
%   cost function to be optimized. Must be the same length as the IRF.
%
% [ssRad] = state_space_radiation_approx(...) outputs the resulting state
%   space model
%
% [ssRad, outMessage] = state_space_radiation_approx(...) also outputs the
%   exit message from FMINCON

% TODO -
%   Maybe I can calculate the gradient analytically???

if modelOrder < 0 || mod(modelOrder,1) ~=0
    error('modelOrder must be a positive integer!')
end

% get the radiation IRF
nVals = length(obj.radIrf.t);

% Get inputs
if nargin > 2
    [options, plotFlag, X0, weight] = get_inputs(varargin{1});
else
    [options, plotFlag, X0, weight] = get_inputs();
end

% configure approximation model
C = [ zeros(1,modelOrder-1) 1];

% Run optimixation
[xstar, ~, exitFlag, outMessage] = fminunc(@cost, X0, options);

if exitFlag <= 0
    warning('Optimization may not have converged. Consider running again')
end


% Calculate resulting state space model
ssRad.A = [ [zeros(1,modelOrder-1) ; eye(modelOrder-1)] -xstar(1:modelOrder)'];
ssRad.B = xstar(modelOrder+1:end)';
ssRad.C = [ zeros(1,modelOrder-1) 1];

% set obj property that holds my state space model.
obj.ssRad = ssRad;

% Set my outputs
if nargout > 1
    varargout{1} = ssRad;
end
if nargout == 2
    varargout{2} = outMessage;
end

% If we want to plot the comparison - do it here
if plotFlag == true
    radApprox = nan(size(obj.radIrf.irf));
    for jj = 1:nVals
        radApprox(jj) = ssRad.C*expm(obj.radIrf.t(jj).*ssRad.A)*ssRad.B;
    end

    h = obj.disp_irf('rad');
    set(h, 'name', 'Radiation IRF SS Approximation')
    hold on
    plot(obj.radIrf.t, radApprox, 'r*')
    hold off
    legend('Actual IRF', 'SS Approximation', 'location', 'northeast')
    title({'Radiation Force State Space Approximation', ...
        sprintf('Model Order: %i', modelOrder)});
end


% my cost function
function Q = cost(x)
    A = [ [zeros(1,modelOrder-1) ; eye(modelOrder-1)] -x(1:modelOrder)'];
    B = x(modelOrder+1:end)';
    q = nan(nVals,1);
    for ii = 1:nVals
        q(ii) = weight(ii) .* ...
            (obj.radIrf.irf(ii) - C * expm(obj.radIrf.t(ii).*A) * B ) ^2;
    end
    Q = sum(q);
end

% This function figures out all the inputs. Sorry its ugly but it works?
    function [options, plotFlag, X0, weight] = get_inputs(varargin)
        % set default values
        options.MaxFunEvals = 1200;
        options.LargeScale  = 'off';
        options.Display     = 'off';
        
        X0 = [3.*ones(1,modelOrder) repmat(max(obj.radIrf.irf)/(3.*modelOrder),1,modelOrder)];
        weight = linspace(1,1, nVals);
        plotFlag = true;
        
        if nargin > 0
            if isfield(varargin{1}, 'plotFlag')
                plotFlag = varargin{1}.plotFlag;
                varargin{1} = rmfield(varargin{1}, 'plotFlag');
            end
            if isfield(varargin{1}, 'X0')
                if length(varargin{1}.X0) == 2*modelOrder
                    X0 = varargin{1}.X0;
                    varargin{1} = rmfield(varargin{1}, 'X0');
                else
                    error('X0 is the wrong size')
                end
            end
            if isfield(varargin{1}, 'Weight')
                if length(varargin{1}.Weight) == nVals
                    weight = varargin{1}.Weight;
                    varargin{1} = rmfield(varargin{1}, 'Weight');
                else
                    error('Weight is the wrong size')
                end
            end
            fields = fieldnames(options);
            for ii = 1:length(fields)
                if ~isfield(varargin{1}, fields{ii})
                    varargin{1}.(fields{ii}) = options.(fields{ii});
                end
            end
            clear options
            options = varargin{1};
        end
    end

end