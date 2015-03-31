function [t, eta, varargout] = jonswap_timeseries(Hs, Ts, gamma, tspan, dt)
% calculates a jonswap water surface elevation time series
%
% Source: 
% Random Seas and Design of Maritime Structures, 3rd., Yoshimi Goda, pp. 35



% Frequency values
% bin width = 0.05, range = [0 3] (Hz)
df = 0.001;
f = df/2:df:(0.75-df/2);

% constants in JONSWAP Spectrum Eqn
Bj = 0.0624*(1.094 - 0.01915*log(gamma)) / ...
    (0.23 + 0.0336*gamma - 0.185 *(1.9 + gamma)^(-1) );
Tp = Ts / (1 - 0.132*(gamma + 0.2)^(-0.559));
sigma = 0.09 .* ones(size(f));
sigma(f <= (1/Tp)) = 0.07;

% calculate Spectrum
S = (Bj * Hs^2 * Tp^(-4)) .* f.^(-5) .* exp(-1.25 .* (Tp.*f).^(-4)) ...
    .* gamma .^ ( exp(-(Tp.*f - 1).^2 ./ (2*sigma.^2)) );
a = sqrt(2 .* S .* df);

% Generate time series elevation data
t = tspan(1):dt:tspan(2);
eta = zeros(size(t));
for ii = 1:length(S)
    eta = eta + a(ii) .* sin(2.*pi.*f(ii).*t + 2*pi*rand(1));
end

% output the spectrum as well
if nargout > 2
    varargout{1}.f = f;
    varargout{2}.S = S;
end

