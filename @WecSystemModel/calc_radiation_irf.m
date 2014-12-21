function [rad] = calc_radiation_irf(obj, varargin)
% Given a .hyd file, calculates the causal radiation impulse response
%
%
% References:   * Ruehl, 2011, Eqn, 2.7
%               * Falnes, 2002, Eqn 2.172-2.174
%
% Inputs:
%   mag   - magnitude of radiation damping response
%   omega  - corresponding frequency vector. must be in RAD/S
%   tspan - t vector used for output impulse response function
% Outputs
%   rad - Structure with two fields:
%     > .irf - contains irf values
%     > .t   - corresponding time vector. Same as tspan input if used
%
% [rad] = calc_radiation_irf(mag, omega, tspan(opt) )

% Version 1.2
%   8/25/2014 - BL
%   1. Squashed the bug. switched coefficent calc from time to omega as the
%      equation specifys.
%   2. Added a zero point for magnitude of radiation response
% Version 1.1
%   8/18/2014 - BL
%   Updated to use simpson's rule for better accuracy. If number of
%   subdivisions is odd it will fall back on trapeziodal rule and issue a
%   warning.
% Version 1.0
%   Created, 8/15/2014, Author: B.L.

% Equations:
%
% $$ f_{r}(t)_{irf} = \frac{2}{\pi} \int_{0}^{\infty} f_{r}(\omega)cos(\omega t) d \omega$$
%




% -------------------------------------------------------------------------
% Algorithm

% input handling
if nargin > 1
    rad.t = varargin{1};
else
    rad.t = 0:0.05:20;
end

omega = obj.radFreq(:,1);
mag   = obj.radFreq(:,2);



% Add zero point
omega = [0 ; omega];
mag   = [0 ; mag];

% idx = find(mag == max(mag));
% if length(idx) == 1;
%     idxFlip = idx;
%     mag = [mag(1:idxFlip) ; flipud(mag(1:idxFlip-1))];
% elseif mod(length(idx),2) == 1
%     idxFlip = idx(floor(length(idx)/2));
%     mag = [mag(1:idxFlip) ; flipud(mag(1:idxFlip-1))];
% else
%     idxFlip = idx((length(idx)/2));
%     mag = [mag(1:idxFlip) ; flipud(mag(1:idxFlip))];
% end
% domega = floor( (omega(2) - omega(1)) .* 1000 ) ./ 1000;
% omega = domega .* [ 0:length(mag)-1 ]';

% Implement integration - Trapezoidal rule
rad.irf = nan(size(rad.t));


if mod(length(omega)-1,2) ~= 0 % Use trapezoidal rule
%    warning('Odd Number of subdivisions. Implementing via Trapeziodal Rule');
    coeff = ( omega(end)-omega(1) ) / ( 2 * (length(omega)-1) );
    for ii= 1:length(rad.t)
        f =  ( mag .* cos( omega.*rad.t(ii) ) );
        rad.irf(ii) = (2/pi) .* coeff .* (sum(f) + sum(f(2:end-1)));
    end
else % Use Simpson's composite rule
    coeff = ( omega(end)-omega(1) ) / ( 3 * (length(omega)-1) );
    for ii = 1:length(rad.t)
        f = ( mag .* cos( omega.*rad.t(ii) ) );
        rad.irf(ii) = (2/pi) .* coeff .* (sum(f) + sum(f(2:end-1)) + 2.*sum(f(2:2:end-1)));
    end
end


