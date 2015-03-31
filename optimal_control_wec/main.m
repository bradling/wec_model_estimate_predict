% test indirect optimal controller
%function meanPower = run_indirect_opt_control(varargin)



if exist('wecModel', 'var') == 0
    if exist('model_data.mat', 'file') == 2
        load('model_data.mat')
    else
        initialize_model;
    end
end


% time span to run optimal control on
tspan = [600 920];
idx = dsearchn(feTime', tspan(1)):dsearchn(feTime', tspan(end));

x0 = simOutConstGen.state(:, idx(1));

% Configure controller
Q = diag([0 0 0 1e5 0]);
R = [0.4e-6];
S = zeros(5,1);
S(end) = 1.2;

% for hard controller limits on Fgen
umin = -1e6;
umax = 1e6;


ss = wecModel.construct_state_space_model();


warning('off','all')
results = indirect_opt(ss.A, ss.B, ss.B, x0, Q, R, S, fe(idx), feTime(idx), dt);
resultsLimits = indirect_with_limits(ss.A, ss.B, ss.B, x0, Q, 0.5.*R, S, umin, umax, fe(idx), feTime(idx), dt);
warning('on','all')

meanPower = mean(results.power)./1e3;
meanPowerLimits = mean(resultsLimits.power)./1e3;


fprintf('\nResults:\n')
fprintf('Case     Max x   Max Fg (kN)  Mean Power (kW)\n')
fprintf('Unc.     %4.2f     %8.0f       %8.1f\n', ...
    max(abs(simOutConstGen.z(idx))), ...
    max(abs(-bgen .* simOutConstGen.zDot(idx)))./1e3, ...
    mean(powerConstGen(idx)));
fprintf('Contr.   %4.2f     %8.0f       %8.1f\n', max(abs(results.z)), ...
    max(abs(results.uc))./1e3, meanPower);
fprintf('Limits.  %4.2f     %8.0f       %8.1f\n', max(abs(resultsLimits.z)), ...
    max(abs(resultsLimits.uc))./1e3, meanPowerLimits);

%fprintf('Mean Power (unc)         : %.3f kW\n', mean(powerConstGen(idx)));
%fprintf('Mean Power (con)         : %.3f kW\n', meanPower);
%fprintf('Mean Power (hard limits) : %.3f kW\n', meanPowerLimits);
    
scrsz = get( groot, 'Screensize' );

h1 = figure;
set(h1, 'color', 'w')
set(h1, 'name', 'Power Optimizing Control Plots')
set(h1, 'position', [0.1*scrsz(3) 0.05*scrsz(4) 0.45*scrsz(3) 0.8*scrsz(4)])

subplot(3,1,1)
plot(feTime(idx), results.z, ...
     feTime(idx), resultsLimits.z, ...
     feTime(idx), simOutConstGen.z(idx), tEta, eta, 'g')
grid on
xlim([tspan(1)+120 tspan(1)+180])
ylabel('Position (m)')
%legend('Controlled', 'Hard F_g limits', 'const. b_g', 'eta', 'location', 'best')
legend('Case 2', 'Case 3', 'Case 1', 'Water Surf.', 'location', 'best');


subplot(3,1,2)
plot(feTime(idx), results.power./1e3, ...
     feTime(idx), resultsLimits.power./1e3, ...
     feTime(idx), powerConstGen(idx))
grid on
ylabel('Power (kW)')
xlim([tspan(1)+120 tspan(1)+180])


subplot(3,1,3)
plot(feTime(idx), results.uc ./1e3, ...
     feTime(idx), resultsLimits.uc ./1e3, ...
     feTime(idx), -bgen .* simOutConstGen.zDot(idx) ./ 1e3)
grid on
ylabel('F_g (kN)')
xlim([tspan(1)+120 tspan(1)+180])
xlabel('Time (sec)')

% plot(feTime(idx), results.bgen, feTime([idx([1 end])), [bgen bgen])
% grd on
% ylabel('bgen')


figure
set(gcf, 'color', 'w')
set(gcf, 'name', 'Fe and Velocities')
plot(feTime(idx), results.zDot, ...
     feTime(idx), resultsLimits.zDot, ...
     feTime(idx), simOutConstGen.zDot(idx), feTime(idx), fe(idx)./1e6, 'g')
xlabel('Time (sec)')
ylabel('Velocity // fe/1e6')
grid on
legend('Controlled', 'Hard F_g limits', 'Cont. b', 'Fe', 'location', 'best')


h2 = figure;
set(h2, 'color', 'w')
set(h2, 'units', 'inches')
set(h2, 'position', [1.6 0.3 5 4])
set(h2, 'name', 'Generator Damping Ratio')

plot(feTime(idx), results.bgen./1e3, feTime(idx), resultsLimits.bgen./1e3)
xlabel('Time (sec)')
ylabel('b_g (kNs/m)')
grid on
xlim([tspan(1)+120 tspan(1)+180])
legend('Case 2', 'Case 3', 'location', 'northeast')



