function feIrf = calc_fe_irf(obj)
% Calculate the Fe IRF Function
%
% REFERENCE: falnes, ocean waves and oscillating systems
%            pp. 142, 33
%
% Inputs:
%   phase - vector of phase shifts of excitation force
%   mag   - magnitude of excitation force response
%   freq  - corresponding frequency vector. must be in RAD/S
%   tspan - t vector used for output impulse response function
%   zero_value - value used for magnitude at zero. default is 1.1 * mag(0)+
% Outputs
%   irf - Structure with two fields:
%     > .irf - contains irf values
%     > .t   - corresponding time vector. Same as tspan input if used
%
% [irf] = calc_irf(phase, mag, freq, tSpan(opt), zero_value(opt) )

% Version 1.2
%   8/25/2014 - BL
%   Bug Fixed. Results now match expected freq response
% Version 1.1
%   8/18/2014 - BL
%   Updated to use simpson's rule for better accuracy. If number of
%   subdivisions is odd it will fall back on trapeziodal rule and issue a
%   warning.
% Version 1.0
%   Create 8/12/2014. Author: B.L.

% Equations
% 
% $$F_{e}(t) = \frac{1}{2\pi}\int_{-\infty}^{\infty}fe(\omega)e^{i \omega t} d \omega$$
%
%
% $$F_{e}(t) = \frac{1}{2\pi} 
% \int_{-\infty}^{\infty}R(\omega)cos(\omega t) +
% I(\omega)sin(\omega t) d\omega + 
% \frac{i}{2\pi}\int_{-\infty}^{\infty}R(\omega)sin(\omega t) +
% I(\omega)cos(\omega t) d\omega$$




% -------------------------------------------------------------------------
% INput case handling
mag   = obj.feFreq(:,2);
phase = obj.feFreq(:,3);
freq  = obj.feFreq(:,1);
t = obj.feIrf.t;

zero_value = 1.1*mag(1);



% Convert to two sided responses
% Assumes Re part is even, Im part is odd
omega = [-flipud(freq); 0;  freq];

realPart = mag .* cos(phase);
imagPart = mag .* sin(phase);

realPart = [ flipud(realPart) ; zero_value ; realPart];
imagPart = [-flipud(imagPart) ; 0   ; imagPart];

%fe_freq = complex(realPart,imagPart);


% Performs integration via trapeziodal rule
fe_irf = zeros(size(t));
if mod(length(omega)-1,2) ~= 0 % Use trapezoidal rule
    warning('Odd Number of subdivisions. Implementing via Trapeziodal Rule');
    coeff = ( omega(end)-omega(1) ) / ( 2 * (length(omega)-1) );
    for p = 1:length(t)
        fx = (1/(2*pi))  .* ( realPart.*cos(omega.*t(p)) - imagPart.*sin(omega.*t(p)) ) + ...
            (1i/(2*pi)) .* ( realPart.*sin(omega.*t(p)) - imagPart.*cos(omega.*t(p)) );
        fe_irf(p) = coeff .* (sum(fx) + sum(fx(2:end-1)));
    end
else % Use Simpson's composite rule
    coeff = ( omega(end)-omega(1) ) / ( 3 * (length(omega)-1) ); % Simpson's rule
    for p = 1:length(t)
        fx = (1/(2*pi))  .* ( realPart.*cos(omega.*t(p)) - imagPart.*sin(omega.*t(p)) ) + ...
            (1i/(2*pi)) .* ( realPart.*sin(omega.*t(p)) - imagPart.*cos(omega.*t(p)) );
        fe_irf(p) = coeff .* (sum(fx) + sum(fx(2:end-1)) + 2.*sum(fx(2:2:end-1)));
    end
end

% To fix numerical inacuracies in complex conjugate math
% A better fix would be to 
fe_irf = real(fe_irf);

% Output into a structure
feIrf.irf = fe_irf;
feIrf.t   = t;
