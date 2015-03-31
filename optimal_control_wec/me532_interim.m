% me532_interim.m
clc; clear all; close all;
warning('off','all')
addpath('../')

% Load Wec Model

% Define Wec Model

hydFile= '../../2015_disturbance_estimation_prediction/hydrodynamic_model/hyd_files/cylinder_d10_h5_fb2_lesswidefreq.hyd';

sysOpt.mass = 140076.5;
sysOpt.kHyd = 1025 * 9.81 * pi * 25;
sysOpt.feIrf.t = -25:0.05:25;
sysOpt.radIrf.t = 0:0.05:20;
sysOpt.bGen = 1e5;

wecModel = WecSystemModel(hydFile, sysOpt);
options.plotFlag = true;
options.Weights = linspace(1.2, 1, length(wecModel.radIrf.t)) .^ 2;
wecModel.state_space_radiation_approx(3, options);



%%
% look at Average power produced at different frequencies and generator
% damping ratios

bGen = [0 1e4 1e5 3e5 1e6 3e6 1e7];

% wave periods
%T = 2:1:30;
%omega = 2.* pi ./ T;
omega = [0.05, 0.1:0.1:2];
T = 2.*pi ./ omega;
dt = 0.05;
t = 0:dt:1000;
idxEnd = floor(length(t)/4);

powerMag   = nan(length(T), length(bGen));
powerPhase = nan(length(T), length(bGen));
posMag     = nan(length(T), length(bGen));
posPhase   = nan(length(T), length(bGen));
velMag     = nan(length(T), length(bGen));

for ii = 1:length(omega)
    tic;
    eta = 0.5 .* sin(omega(ii) .* t);
    for jj = 1:length(bGen)
        wecModel.set_bgen(bGen(jj));
        simResults = wecModel.run_state_space_simulation(eta,dt);
        power = simResults.zDot .^2 .* bGen(jj);
        powerMag(ii, jj) = mean(power(idxEnd:end));
%        powerPhase(ii, jj) = calc_phase_shift(simResults.eta(idxEnd:end), ...
%                                              power(idxEnd:end), ...
%                                              dt, omega(ii));                                   
        posMag(ii, jj) = 2 .* max(simResults.z(idxEnd:end));
        posPhase(ii, jj) = calc_phase_shift(simResults.eta(idxEnd:end), ...
                                            simResults.z(idxEnd:end), ...
                                            dt, omega(ii));
        velMag(ii, jj) = 2 .* max(simResults.zDot(idxEnd:end));
    end
    toc
end


legEntry = cell(size(bGen));
legEntry{1} = sprintf('b_{g}=0 kNs/m');
for ii = 2:length(bGen)
    legEntry{ii} = sprintf('b_{g}=%.g kNs/m', bGen(ii)./1e3);
    legEntry{ii}(end-8:end-7) = '';
end

h = figure;
h.Color = 'w';
h.Units = 'inches';
h.Position = [1 1 6 4];

hsp1=subplot('Position', [0.1300 0.5838 0.6150 0.3412]);
plot(omega, posMag)
grid on
ylabel('Position RAO (m/m)')
set(hsp1,'xtick',0:0.2:2);

hsp2=subplot('Position', [0.1300 0.1100 0.6150 0.3412]);
plot(omega, (180/pi).*posPhase)
grid on
xlabel('Frequency (rad/s)')
ylabel('Phase w.r.t. \eta, (deg)')
set(hsp2,'xtick',0:0.2:2);
set(hsp2,'ytick',-180:45:180);
set(hsp2,'ylim',[-90 45])

hl = legend(legEntry);
set(hl, 'fontsize', 8);
hl.Position = [0.7650 0.3330 0.2205 0.3342];


h4 = figure;
h4.Color = 'w';
h4.Units = 'inches';
h4.Position = [1 1 6 4];


plot(omega, powerMag./1e3)
grid on
xlabel('Frequency (rad/s)')
ylabel('Average Power per H^2 (kW/m^2)','interpreter','tex')
hl = legend(legEntry);
set(hl, 'fontsize', 8);

% subplot(2,2,4)
% plot(omega, powerPhase)
% grid on
% xlabel('Frequency (rad/s)')
% ylabel('Phase Shift w.r.t. \eta (deg)')


%%
% Look at eigenvalues of system.

% uncontrolled - no generator force

wecModel.set_bgen(0);
SS = wecModel.construct_state_space_model();
lambda = eig(SS.A);

% With a generator
wecModel.set_bgen(1e5);
SS2 = wecModel.construct_state_space_model();
lambda2 = eig(SS2.A);
wecModel.set_bgen(1e6);
SS3 = wecModel.construct_state_space_model();
lambda3 = eig(SS3.A);


fprintf('Eigenvalues of system A matrix with no generator:\n')
disp(lambda);

fprintf('Eigenvalues with bGen = 1e5:\n')
disp(lambda2);

h2 = figure;
h2.Color = 'w';
h2.Units = 'inches';
h2.Position = [1.5 0.5 5 4.5];

p1 = plot(real(lambda), imag(lambda), 'o', real(lambda2), imag(lambda2), 'x', ...
     real(lambda3), imag(lambda3), '+', 'linewidth', 2);
grid on
xlabel('Real')
ylabel('Imag')
legend('b_g = 0 kNs/m', 'b_g = 100 kNs/m', 'b_g = 1000  kNs/m', 'location', 'northwest')
%title('Eigenvalues of A')

warning('on','all')


%% Plot system response in a random Pierson-Moskowitz Spectrum
clearvars t eta dt


Hs = 3.3;
Ts = 9;
gamma = 1.5;
rng('default') 
[t, eta] = jonswap_timeseries(Hs, Ts, gamma, [0 1000], 0.1);


b = [3e5 1e6];
wecModel.set_bgen(b(1));
simResults1 = wecModel.run_state_space_simulation(eta, t(2)-t(1));
power1 = (simResults1.zDot.^2 .* b(1)) ./ 1e3;

wecModel.set_bgen(b(2));
simResults2 = wecModel.run_state_space_simulation(eta, t(2)-t(1));
power2 = (simResults2.zDot.^2 .* b(2)) ./ 1e3;


h3 = figure;
h3.Color = 'w';
h3.Units = 'inches';
h3.Position = [2 0.7 6 5];

xrange = [800 920];
%meanPower = mean( power1(t >= xrange(1) & t <= xrange(2)));

subplot(2,1,1);
hp31 = plot(t, eta, simResults1.t, simResults1.z, 'r', ...
    simResults2.t, simResults2.z, 'linewidth', 1.3);
hp31(3).Color = [0.3 0.7 0.2];
grid on
xlim(xrange)
ylabel('Position (m)')
legend('Water Surface', 'WEC Position, b_g = 300 kNs/m', ...
    'WEC Position, b_g = 1000 kNs/m', 'location', 'northeast')

subplot(2,1,2)
%plot(simResults.t, power, 'b', xrange, [meanPower meanPower], 'r');
hp32 = plot(simResults1.t, power1, 'r', simResults2.t, power2, ...
    'linewidth', 1.3);
hp32(2).Color = [0.3 0.7 0.2];
grid on
xlim(xrange)
ylabel('Power Produced (kW)')
%legend('Instantaneous Power', 'Mean Power', 'location', 'northeast')
legend('b_g = 300 kNs/m', 'b_g = 1000 kNs/m')
xlabel('Time (sec)')