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
Q = diag([0 0 0 0 0]);
R = [0.4e-6];
S = zeros(5,1);
S(end) = 1.2;
umin = -0.5e6;
umax = 1e6;


ss = wecModel.construct_state_space_model();
nRuns = 7;
svals = logspace(2,5,nRuns);

nRuns = 1;
svals = 1;
meanPower = nan(size(svals));

for ii = 1:nRuns
    fprintf('Iter %i of %i, ',ii, nRuns)
    S(end) = svals(ii);
    %results = indirect_opt(ss.A, ss.B, ss.B, x0, Q, R, S, fe(idx), feTime(idx), dt); 
    results = indirect_with_limits(ss.A, ss.B, ss.B, x0, Q, R, S, umin, umax, fe(idx), feTime(idx), dt); 

    meanPower(ii) = mean(results.power)./1e3;
    fprintf('  mean power: %.3f kW\n', meanPower(ii));
end

if length(svals) > 1
    figure
    semilogx(svals, meanPower)
    xlabel('S')
    ylabel('meanPower (kW)')
    grid on
    
    S(end) = svals(meanPower == max(meanPower));
    results = indirect_opt(ss.A, ss.B, ss.B, x0, Q, R, S, fe(idx), feTime(idx), dt);
else
    fprintf('\nResults:\n')
    fprintf('Mean Power (unc): %.3f kW\n', mean(powerConstGen(idx)));
    fprintf('Mean Power (con): %.3f kW\n', meanPower);
end

scrsz = get( groot, 'Screensize' );

h1 = figure;
set(h1, 'color', 'w')
set(h1, 'name', 'Power Optimizing Control Plots')
set(h1, 'position', [0.1*scrsz(3) 0.05*scrsz(4) 0.45*scrsz(3) 0.8*scrsz(4)])

subplot(3,1,1)
plot(feTime(idx), results.z, feTime(idx), simOutConstGen.z(idx), tEta, eta, 'g')
grid on
xlim(tspan)
ylabel('Position (m)')
legend('Controlled', 'const. b_g', 'eta', 'location', 'best')

subplot(3,1,2)
plot(feTime(idx), results.power./1e3, feTime(idx), powerConstGen(idx))
grid on
xlabel('Time (sec)')
ylabel('Power (kW)')

subplot(3,1,3)
plot(feTime(idx), results.uc , feTime(idx), -bgen .* simOutConstGen.zDot(idx))
grid on
ylabel('Fgen')
% plot(feTime(idx), results.bgen, feTime([idx([1 end])), [bgen bgen])
% grd on
% ylabel('bgen')


figure
set(gcf, 'color', 'w')
set(gcf, 'name', 'Fe and Velocities')
plot(feTime(idx), results.zDot, feTime(idx), simOutConstGen.zDot(idx), feTime(idx), fe(idx)./1e6, 'g')
xlabel('Time (sec)')
ylabel('Velocity // fe/1e6')
grid on
legend('Controlled', 'Cont. b', 'Fe', 'location', 'best')



