function verify_state_space_model(obj, varargin)
% runs the state space simulation at a sweep of single sine waves.
% compares results to frequency-response results from Aqwa.
%
% verify_state_space_model(wecSystemModel)
%
% verify_state_space_model(wecSystemModel, plotFlag) plots results if
%   plotFlag is true. Default is true.

if isempty(obj.ssRad)
    error('State Space Radiation Approximation not calculated.');
end

% set bGen to zero
bGen = obj.bGen;
obj.bGen = 0;

plotFlag = true;
if nargin == 2
    plotFlag = varargin{1};
end

freq = 0.1:0.1:4;
raox = nan(size(freq));
frMag = nan(size(freq));
feMag = nan(size(freq));
fePhase = nan(size(freq));
zdotMax = nan(size(freq));

dt = 0.1;
t = 0:dt:1000;

dlt = '';
for ii = 1:length(freq)
    msg = sprintf('Processing %i of %i',ii,length(freq));
    fprintf([dlt msg])
    dlt = repmat('\b',1,length(msg));
    eta = sin(freq(ii).*t);
    simresults = obj.run_state_space_simulation(eta, dt);
    
    raox(ii) = max(simresults.z(6000:end));
    frMag(ii) = max(simresults.fr(6000:end));
    feMag(ii) = max(simresults.fe(6000:end));
    zdotMax(ii) = max(simresults.zdot(6000:end));
    
    % get phase shift too!
    [~, feLocs] = findpeaks(simresults.fe(4750:end));
    etaIdx = dsearchn(t' ,simresults.t(4750:end)');
    [~, etaLocs] = findpeaks(eta(etaIdx));
    fePhase(ii) = mean( freq(ii) .* ...
        (simresults.t(4749+feLocs(end-4:end)) - t(etaIdx(etaLocs(end-4:end)))));

    % Move phase value to between (pi, pi)
    if fePhase(ii) > pi
        fePhase(ii) = fePhase(ii) - 2*pi;
    elseif fePhase(ii) < -pi
        fePhase(ii) = fePhase(ii) + 2*pi;
    end
end
fprintf(dlt)

% save these summary results.
obj.verify_results.raoRef = abs(obj.feFreq(:,2) ./ ...
    (obj.kHyd  - (obj.mass + obj.addMass(:,2)) .* obj.addMass(:,1).^2 + ...
    1j .* obj.radFreq(:,2) .* obj.radFreq(:,1)) );
obj.verify_results.raoSim = raox;
obj.verify_results.feSim  = feMag;
obj.verify_results.frMag  = frMag;
obj.verify_results.raoZdotSim = zdotMax;
obj.verify_results.fePhase = fePhase;

% Check the validity of the frequency response of the radiation force.
frMagRef = nan(size(freq));
for ii = 1:length(freq)
    if freq(ii) < obj.radFreq(1,1) || freq(ii) > obj.radFreq(end,1)
        frMagRef(ii) = nan;
    else
        R = interp1(obj.radFreq(:,1), obj.radFreq(:,2), freq(ii));
        A = interp1(obj.addMass(:,1), obj.addMass(:,2), freq(ii));
        frMagRef(ii) = abs( (R + 1j*freq(ii)*(A - obj.Ainf)) * zdotMax(ii) );
    end
end


% Compare RAOs
if plotFlag == true
    figure
    set(gcf,'color','w')
    set(gcf,'name','Heave Position RAO Verification')
    plot(2.*pi ./ obj.radFreq(:,1), obj.verify_results.raoRef, 'g', 2.*pi ./ freq, raox, 'b')
    legend('Aqwa', 'SS Time Domain', 'location', 'northeast')
    xlabel('Period (sec)')
    title('Position RAO')
    xlim([0 40])
    grid on

    figure
    set(gcf,'color','w')
    set(gcf,'name','Radiation Freq. Response Verification')
    plot(2.*pi ./ freq, frMagRef, 'g', 2.*pi ./ freq, frMag, 'b')
    legend('Aqwa', 'SS Time Domain', 'location', 'northeast')
    xlabel('Period (sec)')
    title('Radiation Force')
    xlim([0 40])
    grid on

    figure
    set(gcf,'color','w')
    set(gcf,'name','Excitation Freq. Response Verification')
    
    subplot(2,1,1)
    plot(2.*pi ./ obj.feFreq(:,1), obj.feFreq(:,2), 'g', 2.*pi ./ freq, feMag, 'b')
    legend('Aqwa', 'SS Time Domain', 'location', 'northeast')
    xlabel('Period (sec)')
    title('Excitation Force')
    xlim([0 40])
    grid on

    subplot(2,1,2)
    plot(2.*pi ./ obj.feFreq(:,1), obj.feFreq(:,3), 'g', 2.*pi ./ freq, fePhase, 'b')
    xlim([0 40])
    xlabel('Period (sec)')
    ylabel('Phase Shift (rad)')
    grid on
end

% reset bGen
obj.bGen = bGen;
