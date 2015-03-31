% initialize_model()


addpath('../');
clc; close all;

% Define wec model and calc its radiation approximation
hydFile= '../../2015_disturbance_estimation_prediction/hydrodynamic_model/hyd_files/cylinder_d10_h5_fb2_lesswidefreq.hyd';
sysOpt.mass = 140076.5;
sysOpt.kHyd = 1025 * 9.81 * pi * 25;
sysOpt.feIrf.t = -25:0.05:25;
sysOpt.radIrf.t = 0:0.05:20;
sysOpt.bGen = 0;

wecModel = WecSystemModel(hydFile, sysOpt);
options.plotFlag = true;
options.Weights = linspace(1.2, 1, length(wecModel.radIrf.t)) .^ 2;
wecModel.state_space_radiation_approx(3, options);


% build time series sea state

Hs = 3.3;
Ts = 9;
gamma = 1.5;
rng('default');
dt = 0.1;
[tEta, eta] = jonswap_timeseries(Hs, Ts, gamma, [0 1000], dt);

% calculate excitation force
[fe, feTime, feEta] = wecModel.calc_excitation(eta, dt);


% Run open loop simulation
% uncontrolled generator
bgen = 5e5;
wecModel.set_bgen(bgen);
simOutConstGen = wecModel.run_state_space_simulation(eta, dt);
powerConstGen = (simOutConstGen.zDot.^2 .* bgen) ./ 1e3;
wecModel.set_bgen(0);

avgPowerConstGen = mean(powerConstGen)./1000;

