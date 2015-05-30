function output = optimize_3dof_nonlin_estimator_ss_rad(model, options, train, validate)

% TODO - modify this function to streamline with total workflow.

orderRadApprox = options.approxOrder;
% basic commands for constructing and using a 3DOF WEC Model




%% Test in irregular wave conditions

noiseVar = [0.02 0.02 0.002];

cnt = numel(train);
for ii = 1:cnt
    if isempty(train(ii).s)
        train(ii) = model.run_sim(train(ii));
    end
    fprintf('runsim train %.0f\n', ii);
    train(ii) = train(ii).add_noise(noiseVar, noiseVar);
end

cnt = numel(validate);
for ii = 1:cnt
    if isempty(validate(ii).s)
        validate(ii) = model.run_sim(validate(ii));
    end
    fprintf('runsim validate %.0f\n', ii);
    validate(ii) = validate(ii).add_noise(noiseVar, noiseVar);
end



%--------------------------------------------------------------------------
% run estimator

% first construct the simplified system
if isfield(options, 'R')
    Rapprox = options.R;
    Ainf = [model.Ainf(1) model.Ainf(4) model.Ainf(5) ;
        model.Ainf(4) model.Ainf(2) model.Ainf(6) ;
        model.Ainf(5) model.Ainf(6) model.Ainf(3) ];
    Jinv = inv(model.inertia + Ainf);
    ssApprox.b = zeros(6,1);
    ssApprox.d = [1e6.*Jinv ; zeros(3)];
    ssApprox.c = eye(6);
    ssApprox.a = [-Jinv*(Rapprox-model.userData.pto.get_b_matrix) Jinv*(model.userData.mooring.get_k_matrix - model.kHyd) ;
                  eye(3) zeros(3)];
else
    estModel = model.copy;
    estModel.construct_rad_approx(orderRadApprox);
    ssApprox = estModel.construct_state_space_model();

    Ainf = [estModel.Ainf(1) estModel.Ainf(4) estModel.Ainf(5) ;
            estModel.Ainf(4) estModel.Ainf(2) estModel.Ainf(6) ;
            estModel.Ainf(5) estModel.Ainf(6) estModel.Ainf(3) ];
    Jinv = inv(estModel.inertia + Ainf);
    % Scale the excitation forces for better numerical performance
    ssApprox.d = 1e6.*ssApprox.b;
    ssApprox.b = ssApprox.b;
    ssApprox.a(1:3,1:3) = ssApprox.a(1:3,1:3) + Jinv*model.userData.pto.get_b_matrix;
    ssApprox.a(1:3,4:6) = ssApprox.a(1:3,4:6) + Jinv*model.userData.mooring.get_k_matrix;
end



estimator = AdaptiveKalmanDisturbanceEstimator(ssApprox);
estimator.harmonic_model(0.5);


% specify measurement noise covariance
R = diag([noiseVar noiseVar].^(0.5));



    % This is my main objective function
    function obj = run_estimator(x)
        Q = diag([abs(x(1:end-3)') 0 0 0 abs(x(end-2:end)')]);
        estimator.define_covariances(Q, R);
        cnt = numel(train);
        for jj = cnt:-1:1
            [resultsTemp(jj), ~, ~, ~] = estimator.calc_estimates(train(jj), ...
                                  diag([repmat(1e4,size(estimator.augSys.a,1)-3, 1) ; 10 ; 10 ; 10]), ...
                                  [zeros(size(estimator.augSys.a,1)-3, 1); [1;1;1].*2.*pi./9 ]);
            resultsTemp(jj).feHat = resultsTemp(jj).feHat .* 1e6; % fix the scaling
        end
        [rsq, ~] = resultsTemp.calc_estimation_results;  
        obj = sum(1-rsq);
    end

    % This function is for the validation set.
    function stop = output_function(x, optimValues, state)
        stop = false;
        
        % perform validation check.
        Q = diag([abs(x(1:end-3)') 0 0 0 abs(x(end-2:end)')]);
        
        estimator.define_covariances(Q, R);
        cnt = numel(validate);
        for jj = cnt:-1:1
            [resultsVal(jj), ~, ~, ~] = estimator.calc_estimates(validate(jj), ...
                                  diag([repmat(1e4,size(estimator.augSys.a,1)-3, 1) ; 10 ; 10 ; 10]), ...
                                  [zeros(size(estimator.augSys.a,1)-3, 1); [1;1;1].*2.*pi./9 ]);
            resultsVal(jj).feHat = resultsVal(jj).feHat .* 1e6; % fix the scaling
        end
        [rsq, ~] = resultsVal.calc_estimation_results;  
        obj = sum(1-rsq);
         
        % Update plot
        progress_plot.update(optimValues.iteration, optimValues.fval, ...
                             obj, rsq(1), rsq(2), rsq(3), x);
        
        % stop optimixation if validation is getting worse
        if progress_plot.curIdx - progress_plot.minValIdx > 30
            stop = true;
        end
    end
     

% =========================================================================
% Run optimizer

if isfield(options, 'Qinit')
    Qinit = options.Qinit;
else
    % initial guess
    % Qinit = repmat(0.5, 1, 9+6*orderRadApprox);
    % Qinit = [0.5 0.5 0.01 0.5 0.5 0.01 repmat(0.1, 1, 3+6*orderRadApprox)];
    % Qinit = [0.5 0.5 0.02 0.2 0.4 0.1 0.05 0.05 0.07 0.1 0.3 0.1 0.03 0.1 0.1 0.1 0.2 0.05 0.1 0.05 0.25];
    Qinit = [0.0001 45 2.5 19 250 0.6 31 27 12 51 38 38 4 40 8 16 21 0.2 0.5 59 93]; % for 2nd order rad approx.
    % Qinit = [0.01 200 1 10 200 0.1 repmat(10, 1, 3+6*orderRadApprox)]; % first try for 4th order rad approx
    % Qinit = [0.01 50 1 10 50 0.1 repmat(10, 1, 3+6*orderRadApprox)];
    % Qinit = [0.01 10 0.005 7 45 0.1 repmat(12, 1, 6*orderRadApprox) 1 10 100];
    % Qinit = [0.01 8 0.005 6 35 0.004 repmat(12, 1, 6*orderRadApprox) 1 10 100];
    % configure options
end
progress_plot = ProgressPlot;
options.MaxFunEvals = 2000;
options.TolX = 0.00005;
options.TolFun = 1e-3;
options.Display = 'iter';
options.PlotFcns = @optimplotfval;
options.OutputFcn = @output_function;

% Run optimization
[xstar,fval,exitflag,out] = fminsearch(@run_estimator, Qinit, options);

% rerun the best result so I can plot results.
Qbest = diag([abs(xstar(1:end-3)') 0 0 0 abs(xstar(end-2:end)')]);
estimator.define_covariances(Qbest, R);
cnt = numel(validate);
for ii = cnt:-1:1
    [valResults(ii), xhat{ii}, xp{ii}, innov{ii}] = estimator.calc_estimates(validate(ii), ...
        diag([repmat(1e4,size(estimator.augSys.a,1)-3, 1) ; 10 ; 10 ; 10]), ...
        [zeros(size(estimator.augSys.a,1)-3, 1); [1;1;1].*2.*pi./9 ]);
    valResults(ii).feHat = valResults(ii).feHat .* 1e6;
end


% data to output
output.resultsValidation = valResults;
output.estimator = estimator;
output.xHat = xhat;
output.xp = xp;
output.innov = innov;
output.Qbest = Qbest;
output.exitflag = exitflag;
output.optimizeOut = out;
output.mseBest = fval;
output.model = model;


%==========================================================================
% Plot some results

% h=figure;
% set(h,'color', 'w')
% set(h,'units','inches')
% set(h, 'position', [1 1 5 4.5])
% set(h, 'name', 'RSQ Performance')
% 
% subplot(2,2,1)
% plot(valResults.fe(3,:), valResults.feHat(3,:), 'r.',...
%      valResults.fe(2,:), valResults.feHat(2,:), 'g.',...
%      valResults.fe(1,:), valResults.feHat(1,:), 'b.');
% xlim(max(abs(valResults.fe(:))) .* [-1 1])
% ylim(max(abs(valResults.fe(:))) .* [-1 1])
% rsqAll = regression_analysis(...
%             [valResults.fe(1,:)  valResults.fe(2,:) valResults.fe(3,:)], ...
%             [xhat(n+1,:) xhat(n+2,:) xhat(n+3,:)] .*1e6);
% title(sprintf('Combined. R^2 = %.3f',rsqAll))
% 
% subplot(2,2,2)
% plot(valResults.fe(1,:), valResults.feHat(1,:), 'b.')
% xlim(max(abs(valResults.fe(1,:))) .* [-1 1])
% ylim(max(abs(valResults.fe(1,:))) .* [-1 1])
% rsqSurge = regression_analysis(valResults.fe(1,:), valResults.feHat(1,:));
% title(sprintf('surge. R^2 = %.3f',rsqSurge))
% 
% 
% subplot(2,2,3)
% plot(valResults.fe(2,:), valResults.feHat(2,:), 'g.')
% xlim(max(abs(valResults.fe(2,:))) .* [-1 1])
% ylim(max(abs(valResults.fe(2,:))) .* [-1 1])
% rsqHeave = regression_analysis(valResults.fe(2,:), valResults.feHat(2,:));
% title(sprintf('Heave. R^2 = %.3f',rsqHeave ))
% 
% 
% subplot(2,2,4)
% plot(valResults.fe(3,:), valResults.feHat(3,:), 'r.')
% xlim(max(abs(valResults.fe(3,:))) .* [-1 1])
% ylim(max(abs(valResults.fe(3,:))) .* [-1 1])
% rsqPitch = regression_analysis(valResults.fe(3,:), valResults.feHat(3,:));
% title(sprintf('Pitch. R^2 = %.3f',rsqPitch))
% 
% 
% %error = [results.sDot ; results.s ; results.fe./1e6] - xp(1:9,:);
% %Q = [cov(error') zeros(9,3) ; zeros(3,12)];
% %[xHat, xp, innov] = estimator.calc_estimates(y, diag(repmat(1e4,12,1)), zeros(12,1));
% 
% 
% figure;plot(valResults.feHat(2,:)); hold on; plot(valResults.fe(2,:));title('heave Excitation');%xlim([6000 7000])
% figure;plot(valResults.feHat(1,:)); hold on; plot(valResults.fe(1,:));title('surge Excitation');%xlim([6000 7000])
% figure;plot(valResults.feHat(3,:)); hold on; plot(valResults.fe(3,:));title('pitch Excitation');%xlim([6000 7000])
% figure;plot(xhat(4,:)); hold on; plot(valResults.s(1,:));title('surge displacement');%xlim([6000 7000])
% figure;plot(xhat(5,:)); hold on; plot(valResults.s(2,:));title('heave displacement');%xlim([6000 7000])
% figure;plot(xhat(6,:)); hold on; plot(valResults.s(3,:));title('pitch displacement');%xlim([6000 7000])
% 
% 
% 


end