function [phi] = calc_phase_shift(x1, x2, dt, omega)
% CALC_PHASE_SHIFT calculates phase shift between two time domain sinusoids
% performs calculation in time domain by looking at delta t between peaks
%
% inputs: x1, x2:  time domain signals. must have same freq and sample rate
%         dt:      the sample rate of the signals
%         omega:   the angular frequency of the inputs in [rad/s]
%
% output: phi: phase in radians of x2 wrt x1. 
%                 i.e. phi(x1) = 0; phi(x2) = phi
%
% [phi] = calc_phase_shift(x1, x2, dt, omega)

% created 2/11/2015, author: B. Ling

if length(x1) ~= length(x2); error('x1 and x1 not same size'); end


[~, locs1] = findpeaks(x1);%, 'MinPeakHeight', 0.95*max(x1));
[~, locs2] = findpeaks(x2);%, 'MinPeakHeight', 0.95*max(x2));


nPeaks = min([length(locs1) length(locs2)]) - 2;

locs1 = locs1(2:nPeaks);
locs2 = locs2(2:nPeaks);


phi = (mean(locs1 - locs2) .* dt) .* omega;

% [s1, f1, phi1] = spectral_analysis(x1, dt);
% [s2, f2, phi2] = spectral_analysis(x2, dt);
% 
% phi = phi2(s2 == max(s2)) - phi1(s1 == max(s1));

if phi > pi; phi = phi - 2*pi; end
if phi < -pi; phi = phi + 2*pi; end


function [ S, F, varargout] = spectral_analysis( x, ts )

% Check input sizes
if size(x,2) == 1 && size(x,1) > 1, x = x'; end
if size(ts,2) == 1 && size(ts,1) > 1, ts = ts';end
if size(ts,1) == 1 && size(x,1) > 1, ts = linspace(ts,ts,size(x,1))';end

% Run FFT
N    = size(x,2);
X    = fft(x,[],2);
Xmag = abs(X./N);

% Convert FFT results to spectral results
% If N is ODD
if mod(N,2) == 1
    mid = ceil(N/2);
    S = [Xmag(:,1)  2.*Xmag(:,2:mid)./sqrt(2)] .^2;
% if N is even - extra center term at the fold
else
    mid = N/2;
    S = [Xmag(:,1)  2.*Xmag(:,2:mid)./sqrt(2)  Xmag(:,mid+1)./sqrt(2)] .^2;
end

% calculate freq vector
F = (1./ts) * ( (0:(size(S,2))-1) ./ (N) );

% calculate phase angles
if nargout == 3
    varargout{1} = angle(X(1:length(S)));
end
