function output = optimize_3dof_estimator()

% basic commands for constructing and using a 3DOF WEC Model


% initialize model
ansysFilename = './ansys_data/aqwa_results_disk_v4/ANALYSIS.LIS';
ansysFilename = 'hydDataModified_01_to_40.mat';
ansysFilename = './ansys_data/two_body_wec_fixed/ANALYSIS.LIS';
model = WecSysModel3Dof(ansysFilename);

% display frequency response plot
%hFreq = model.disp_freq;

% construct radiation approximation and plot results
model.construct_rad_approx(6, 'false');



% simulate under regular wave conditions.

dt = 0.1;
eta = sin(2*pi*0.1*(0:dt:1000));
results = model.run_sim(eta, dt);





%% Test in irregular wave conditions

Hs = 1;
Ts = 10;
gamma = 1.5;
rng('default');
dt = 0.1;
[tEta, eta] = jonswap_timeseries(Hs, Ts, gamma, [0 1000], dt);

noise = [0.02  .* randn(2,size(results.sDot,2)) ;
         0.002 .* randn(1,size(results.sDot,2)) ;
         0.02  .* randn(2,size(results.sDot,2)) ;
         0.002 .* randn(1,size(results.sDot,2)) ];

results = model.run_sim(eta, dt);
y = [results.sDot ; results.s] + noise;

% valdation model
Hs = 1.2;
Ts = 11;
gamma = 1.5;
dt = 0.1;
[~, etaValidate] = jonswap_timeseries(Hs, Ts, gamma, [0 300], dt);
valResults = model.run_sim(etaValidate, dt);


valResults = model.run_sim(etaValidate, dt);
valResults = results;

noise = [0.02  .* randn(2,size(valResults.sDot,2)) ;
         0.002 .* randn(1,size(valResults.sDot,2)) ;
         0.02  .* randn(2,size(valResults.sDot,2)) ;
         0.002 .* randn(1,size(valResults.sDot,2)) ];

yVal = [valResults.sDot ; valResults.s] + noise;




figure
plot(results.t, results.eta, results.t, results.s(2,:))
xlim([800 860])
xlabel('Time')
ylabel('Position')
title('Heave Position')
legend('Eta', 'x', 'location', 'best')


% run estimator

R = [11000 0      58000 ;
     0     480000 0 ; 
     58000 0      260000];
Ainf = [model.Ainf(1) model.Ainf(4) model.Ainf(5) ;
        model.Ainf(4) model.Ainf(2) model.Ainf(6) ;
        model.Ainf(5) model.Ainf(6) model.Ainf(3) ];
Jinv = inv(model.inertia + Ainf);
sysEst.b = zeros(6,1);
sysEst.d = [1e6.*Jinv ; zeros(3)];
sysEst.c = eye(6);
sysEst.a = [-Jinv*R -Jinv*model.kHyd ; eye(3) zeros(3)];

estimator = KalmanDisturbanceEstimator(sysEst);
estimator.harmonic_model(0.11*2*pi, dt);




%Q = diag([0.01 0.01 0.001 0.01 0.01 0.001 0.5 0.1 0.7 0 0 0]);
R = diag([0.02 0.02 0.002 0.02 0.02 0.002].^(0.5));

% find max excitation forces to scale them in my objective function.
maxFe = max(results.fe, [], 2)./1e6;
maxFeValidation = max(valResults.fe, [], 2)./1e6;

    function mse = run_estimator(x)
        %Q = diag([0.01 0.01 0.001 0.01 0.01 0.001 x 0 0 0]);
        Q = diag([x 0 0 0]);
        estimator.define_covariances(Q, R);
        [xHat, ~, ~] = estimator.calc_estimates(y, diag(repmat(1e4,12,1)), zeros(12,1));
        mse = mean([(results.fe(1,:)./1e6 - xHat(7,:)) ./ maxFe(1) ...
                    (results.fe(2,:)./1e6 - xHat(8,:)) ./ maxFe(2) ...
                    (results.fe(3,:)./1e6 - xHat(9,:)) ./ maxFe(3)] .^2);
%         rsq = [regression_analysis(results.fe(1,:)./1e6, xHat(7,:)) ; 
%                regression_analysis(results.fe(2,:)./1e6, xHat(8,:)) ; 
%                regression_analysis(results.fe(3,:)./1e6, xHat(9,:)) ];  
%         obj = 1-mean(rsq);
    end

    function stop = output_function(x, optimValues, state)
        stop = false;
        
        % perform validation check.
        Q = diag([x 0 0 0]);
        estimator.define_covariances(Q, R);
        [xHatVal, ~, ~] = estimator.calc_estimates(yVal, diag(repmat(1e4,12,1)), zeros(12,1));
        rsq = [regression_analysis(valResults.fe(1,:)./1e6, xHatVal(7,:)) ; 
               regression_analysis(valResults.fe(2,:)./1e6, xHatVal(8,:)) ; 
               regression_analysis(valResults.fe(3,:)./1e6, xHatVal(9,:)) ];  
%         obj = 1-mean(rsq);
        mse = mean([(valResults.fe(1,:)./1e6 - xHatVal(7,:)) ./ maxFeValidation(1) ...
                    (valResults.fe(2,:)./1e6 - xHatVal(8,:)) ./ maxFeValidation(2) ...
                    (valResults.fe(3,:)./1e6 - xHatVal(9,:)) ./ maxFeValidation(3)] .^2);
        
        
        progress_plot.update(optimValues.iteration, optimValues.fval, ...
                             mse, rsq(1), rsq(2), rsq(3));
    end
     
 
progress_plot = ProgressPlot;
options.MaxFunEvals = 500;
options.TolX = 0.00002;
options.TolFun = 1e-7;
options.Display = 'iter';
options.PlotFcns = @optimplotfval;
options.OutputFcn = @output_function;
%[xstar,fval,exitflag,out] = fminsearch(@run_estimator, [0.5 0.1 0.5], options);
%[xstar,fval,exitflag,out] = fminsearch(@run_estimator, [0.01 0.01 0.001 0.01 0.01 0.001 0.5 0.1 10], options);
Qinit = [0.03 0.01 0.001 0.03 0.01 0.001 0.5 0.1 10];
Qinit = [0.002 0.012 0.002 0.04 0.01 0.0002 0.001 0.1 25];
[xstar,fval,exitflag,out] = fminsearch(@run_estimator, Qinit, options);

%Qbest = diag([0.01 0.01 0.001 0.01 0.01 0.001 xstar 0 0 0]);
Qbest = diag([xstar 0 0 0]);

estimator.define_covariances(Qbest, R);
[xHat, xp, innov] = estimator.calc_estimates(y, diag(repmat(1e4,12,1)), zeros(12,1));

output.xHat = xHat;
output.xp = xp;
output.innov = innov;
output.Qbest = Qbest;
output.results = results;
output.exitflag = exitflag;
output.optimizeOut = out;
output.mseBest = fval;
output.model = model;

h=figure;
set(h,'color', 'w')
set(h,'units','inches')
set(h, 'position', [1 1 5 4.5])
set(h, 'name', 'RSQ Performance')

subplot(2,2,1)
plot(results.fe(3,:), xHat(9,:).*1e6, 'r.',...
     results.fe(2,:), xHat(8,:).*1e6, 'g.',...
     results.fe(1,:), xHat(7,:).*1e6, 'b.');
xlim(max(abs(results.fe(:))) .* [-1 1])
ylim(max(abs(results.fe(:))) .* [-1 1])
rsqAll = regression_analysis(...
            [results.fe(1,:)  results.fe(2,:) results.fe(3,:)], ...
            [xHat(7,:) xHat(8,:) xHat(9,:)] .*1e6);
title(sprintf('Combined. R^2 = %.3f',rsqAll))

subplot(2,2,2)
plot(results.fe(1,:), xHat(7,:) .*1e6, 'b.')
xlim(max(abs(results.fe(1,:))) .* [-1 1])
ylim(max(abs(results.fe(1,:))) .* [-1 1])
rsqSurge = regression_analysis(results.fe(1,:), xHat(7,:) .*1e6);
title(sprintf('surge. R^2 = %.3f',rsqSurge))


subplot(2,2,3)
plot(results.fe(2,:), xHat(8,:) .*1e6, 'g.')
xlim(max(abs(results.fe(2,:))) .* [-1 1])
ylim(max(abs(results.fe(2,:))) .* [-1 1])
rsqHeave = regression_analysis(results.fe(2,:), xHat(8,:) .*1e6);
title(sprintf('Heave. R^2 = %.3f',rsqHeave ))


subplot(2,2,4)
plot(results.fe(3,:), xHat(9,:) .*1e6, 'r.')
xlim(max(abs(results.fe(3,:))) .* [-1 1])
ylim(max(abs(results.fe(3,:))) .* [-1 1])
rsqPitch = regression_analysis(results.fe(3,:), xHat(9,:) .*1e6);
title(sprintf('Pitch. R^2 = %.3f',rsqPitch))


%error = [results.sDot ; results.s ; results.fe./1e6] - xp(1:9,:);
%Q = [cov(error') zeros(9,3) ; zeros(3,12)];
%[xHat, xp, innov] = estimator.calc_estimates(y, diag(repmat(1e4,12,1)), zeros(12,1));


figure;plot(1e6.*xHat(8,:)); hold on; plot(results.fe(2,:));title('heave Excitation');%xlim([6000 7000])
figure;plot(1e6.*xHat(7,:)); hold on; plot(results.fe(1,:));title('surge Excitation');%xlim([6000 7000])
figure;plot(1e6.*xHat(9,:)); hold on; plot(results.fe(3,:));title('pitch Excitation');%xlim([6000 7000])
figure;plot(xHat(4,:)); hold on; plot(results.s(1,:));title('surge displacement');%xlim([6000 7000])
figure;plot(xHat(5,:)); hold on; plot(results.s(2,:));title('heave displacement');%xlim([6000 7000])
figure;plot(xHat(6,:)); hold on; plot(results.s(3,:));title('pitch displacement');%xlim([6000 7000])





end