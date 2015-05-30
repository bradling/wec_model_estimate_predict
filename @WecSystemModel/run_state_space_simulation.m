% run simulation using state space model
function simResults = run_state_space_simulation(obj, eta, dt, varargin)
%
% simResults = run_state_space_simulation(obj, eta, dt). where dt is the
% sample time for eta.
% 
% simResults = run_state_space_simulation(obj, eta, dt, dtSim) : gives a
% different sample time for results than was input with eta.



% Calculate excitation force
[simResults.fe, feTime, simResults.eta] = calc_excitation(obj, eta, dt);

tvals = feTime;
if nargin == 4
    tvals = feTime(1):varargin{1}:feTime(end);
end

% ramp in fe
simResults.fe(1:200) = linspace(0,1,200) .* simResults.fe(1:200);

% Get state space model
SS = construct_state_space_model(obj);

% Run Simulation
X0 = zeros(size(SS.A,1), 1);
%tvals = (dt .* (0:length(feTime)-1)) + feTime(1);
[simResults.t, yout] = ode45(@eom, tvals, X0);
simResults.t = simResults.t';
simResults.state = yout';

% calculate output equation
simResults.z = nan(1,length(simResults.t));
simResults.zDot = nan(1,length(simResults.t));
simResults.fr = nan(1,length(simResults.t));
simResults.fk = nan(1,length(simResults.t));
for ii = 1:length(simResults.t)
    simResults.z(:,ii) = SS.C(1,:)*yout(ii,:)';
    simResults.zDot(:,ii) = SS.C(2,:)*yout(ii,:)';
    simResults.fr(ii) = obj.ssRad.C * yout(ii,1:3)' * yout(ii,5);
 %   simResults.fr(ii) = -yout(ii,3) * yout(ii,5);
    simResults.fk(ii) = -yout(ii,4) * obj.kHyd;
end

% This is my equations of motion
    function dx = eom(t,x)
        idx = find(floor(t/dt)*dt == feTime);
        if (t-feTime(idx)) > 1e-8
            u = (simResults.fe(idx+1) - simResults.fe(idx))*(t-feTime(idx))/(feTime(idx+1)-feTime(idx)) + simResults.fe(idx);
        else
            u = simResults.fe(idx);
        end
        dx = SS.A*x + SS.B * u;
    end
end