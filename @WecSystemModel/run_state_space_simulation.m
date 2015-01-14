% run simulation using state space model
function simResults = run_state_space_simulation(obj, eta, dt)

% Calculate excitation force
[simResults.fe, feTime] = calc_excitation(obj, eta, dt);

% ramp in fe
simResults.fe(1:200) = linspace(0,1,200) .* simResults.fe(1:200);

% Get state space model
SS = construct_state_space_model(obj);

% Run Simulation
X0 = zeros(size(SS.A,1), 1);
tvals = (dt .* (0:length(feTime)-1)) + feTime(1);
[simResults.t, yout] = ode45(@eom, tvals, X0);
simResults.t = simResults.t';

% calculate output equation
simResults.z = nan(1,length(simResults.t));
simResults.zDot = nan(1,length(simResults.t));
simResults.fr = nan(1,length(simResults.t));
simResults.fk = nan(1,length(simResults.t));
for ii = 1:length(simResults.t)
    simResults.z(:,ii) = SS.C(1,:)*yout(ii,:)';
    simResults.zDot(:,ii) = SS.C(2,:)*yout(ii,:)';
    simResults.fr(ii) = -yout(ii,3) * yout(ii,5);
    simResults.fk(ii) = -yout(ii,4) * obj.kHyd;
end

% This is my equations of motion
    function dx = eom(t,x)
        u = interp1(feTime, simResults.fe, t);
        dx = SS.A*x + SS.B * u;
    end
end