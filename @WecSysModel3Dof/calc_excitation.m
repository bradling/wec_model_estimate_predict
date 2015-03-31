function [fe, feTime, varargout] = calc_excitation(obj, eta, dt)
% Calculates the time domain excitation force via convolution.
% If the sampling rate of the irf is greater than eta, then eta values will
% be upsampled to match by using linear interpolation. If they have the
% same sampling rate, it will be left be.
%
% NOTES: dt/dtIrf MUST BE AN INTEGER. OTHERWISE THIS WILL BREAK.
%
% Inputs:
%   eta: vector containing water surface elevation
%   dt: scalar - sample interval of eta in seconds
%
% [Fe, feTime] = calc_excitation(WecModel, eta, dt)



% Created 7/28/2014
% Bradley Ling
%
% Updated 8/15: need to scale irf to the sampling time. Should implement
% trapeziodal integration for the future.
%
% Updated 8/18: Implemented trapeziodal integration for convolution.
% Tested results with many uniform waves and compared to frequency domain
% results.
%
% Updated 8/22: changed interpolation to cubic

% Updated 8/25: Squashed some bugs with trimming output. Also optimized
% this process to use less memory and calcs.
%  
% References:
% Ruehl Thesis 2011; Bosma et. al. 2013



% -------------------------------------------------------------------------
% Lets interpolate eta so it has the same sampling rate as the irf
%  CUBIC INTERPOLATION CURRENTLY IMPLEMENTED
%
% NOTES: ** dt/dtIrf MUST BE AN INTEGER. OTHERWISE THIS WILL BREAK.
%        ** the time range of the irf must be symmetrical. 


% I want eta to be a row vector
szEta = size(eta);
if szEta(1) > szEta(2), eta = eta';end

irf = obj.feIrf;
dtIrf = irf.t(2) - irf.t(1);

ratio = round(dt/dtIrf);

% Check to make sure ratio is acceptable
if abs(ratio - dt/dtIrf) > 1e-5
    error('sampling rates are not compatable for interpolation')
elseif ratio ~= 1
    etaInterpTime = dtIrf .* ( 0:(length(eta)*ratio - (ratio - 1)) );
    etaInterp     = interp1(dt .* (0:length(eta)-1), eta, etaInterpTime, 'pchip');
else
    % don't do interpolation
    %warning('No Interpolation performed')
    etaInterp = eta;
    etaInterpTime = dtIrf .* [0:length(eta)-1];
end



% Implement Convolution via trapeziodal rule
fe = nan(size(irf.irf, 1), size(etaInterp, 2));
halfWidth = (size(irf.irf, 2) - 1 ) / 2;
startIdx = halfWidth+1;
endIdx   = length(etaInterp)-halfWidth;
idxIrf = -halfWidth:halfWidth;


coeff = (irf.t(end) - irf.t(1)) / (4 .* halfWidth); 
for it = startIdx:endIdx;
    for jj = 1:size(irf.irf, 1);
        fx = irf.irf(jj, :) .* etaInterp(idxIrf + it);
        fe(jj, it) = coeff * sum( [fx fx(2:end-1)] );
    end
end


% Trim the output
trimmedTime = etaInterpTime(~isnan(fe(1,:)));
idxKeep = ~isnan(fe(1,:));
fe = fe(:, idxKeep);

feTime = floor( (min(trimmedTime):dt:max(trimmedTime)) ./ dt) .* dt;

idx = dsearchn(trimmedTime', feTime');
fe = fe(:, idx);
if nargout == 3, 
    temp = etaInterp(idxKeep);
    varargout{1} = temp(idx);
end



