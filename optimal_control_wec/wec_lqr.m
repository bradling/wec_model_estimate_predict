% Test LQR type power maximizing control
%function varargout = wec_lqr(varargin)

%% First Define WEC Model

% Need access to WecSystemModel class
addpath('../');
clc; close all;

if nargin == 0
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
else
    wecModel = varargin{1};
end



%% Define sea state

Hs = 3.3;
Ts = 9;
gamma = 1.5;
rng('default');
dt = 0.1;
[tEta, eta] = jonswap_timeseries(Hs, Ts, gamma, [0 1000], dt);


%% Run open loop simulation with no and some generator damping


% No generator - just a buoy
simOutUnc = wecModel.run_state_space_simulation(eta, dt);

% uncontrolled generator
bgen = 5e5;
wecModel.set_bgen(bgen);
simOutConstGen = wecModel.run_state_space_simulation(eta, dt);
powerConstGen = (simOutConstGen.zDot.^2 .* bgen) ./ 1e3;
wecModel.set_bgen(0);

fprintf('Mean Power (bgen=%.0f): %.3f kW\n',bgen./1e3, mean(powerConstGen));


%% Calculate Controller

% calculate the excitation force
[fe, feTime, feEta] = wecModel.calc_excitation(eta, dt);

%options = optimoptions('fmincon','GradObj','off', 'TolFun', 1, 'MaxFunEvals', 50);
%[xStar, fStar] = fmincon(@lqr_controller, [-9 1e-4 1], diag([1 -1 -1]), [0;0;0], [], [], [], [], [], options);
%    function avgPower = lqr_controller(x)


% Define my Q, R, and S matrices
s1 = 9;
r1 = .000008;
qp = 1;


% dependent term to ensure Qtilda is P.D.
qv = s1^2/r1;

% convert to matrices
Q = diag([0 0 0 qp qv]);
R = [r1];
S = [0 0 0 0 s1]';

% get state space model
ss = wecModel.construct_state_space_model();
ss.B = ss.B; %TODO - what is the D matrix - map disturbance to state trajectory
n = size(ss.A, 1);




% Change of variables to use LQR with cross terms
Rinv = inv(R);
Qt = Q - S*Rinv*S';
At = ss.A - ss.B*Rinv*S';

% solve infinate time lqr problem to find PI and K
[PI, ~, ~] = care(At, ss.B, Qt, R);
b = calc_feedfwd_term(At, ss.B, Rinv, PI, feTime);
simLqr = run_lqr_sim(ss, R, S, PI, b, fe, feTime, dt);


fprintf('Mean Power (bgen=%.0f): %.3f kW\n',bgen./1e3, mean(powerConstGen));
fprintf('Mean Power (LQR Control): %.3f kW\n', mean(simLqr.power)./1000);


scrsz = get( groot, 'Screensize' );

h1 = figure;
set(h1, 'color', 'w')
set(h1, 'name', 'LQR Control Plots')
set(h1, 'position', [0.1*scrsz(3) 0.05*scrsz(4) 0.45*scrsz(3) 0.8*scrsz(4)])

subplot(4,1,1)
plot(feTime, simLqr.z, feTime, simOutConstGen.z, t, eta)
grid on
xlim([800 860])
ylabel('Position (m)')
legend('LQR', 'b_g=500','\eta','location','best')

subplot(4,1,2)
plot(feTime, simLqr.zDot, feTime, simOutConstGen.zDot)
grid on
xlim([800 860])
ylabel('Velocity (m/s)')

subplot(4,1,3)
plot(feTime, simLqr.power./1e3, feTime, powerConstGen)
grid on
xlim([800 860])
ylabel('Power (kW)')

subplot(4,1,4)
plot(feTime, simLqr.uc./1e3)
grid on
xlim([800 960])
ylabel('Generator Force (kN)')
xlabel('Time (sec)')
%end


