function simOut = run_lqr_sim(ss, R, S, PI, b, fe, feTime, dt)

% preliminary stuff
%dt = feTime(2) - feTime(1);
Rinv = inv(R);

% run simulation
x0 = zeros(size(ss.A,1), 1);
soln = ode45(@eom, [feTime(1) feTime(end)], x0);

% extract output of interest
state = deval(soln, feTime);
simOut.z = state(4,:);
simOut.zDot = state(5,:);

y = nan(3,length(feTime));
for ii = 1:length(feTime);
    y(:,ii) = eom(feTime(ii), state(:,ii), 'output');
end
simOut.uc = y(1,:);
simOut.power = y(2,:);
simOut.bgen = y(3,:);

% TODO calculate instantaneous power
%y = nan(size(ss.C,1), length(feTime));
%for ii = 1:length(feTime)
%    y(:,ii) = ss.C
%end

    function dx = eom(t, x, varargin)
        % calculate excitation force
        idx = find(floor(t/dt)*dt == feTime);
        if (t-feTime(idx)) > 1e-8
            ud = (fe(idx+1) - fe(idx))*(t-feTime(idx))/(feTime(idx+1)-feTime(idx)) + fe(idx);
            bval = (b(:,idx+1) - b(:,idx)).*(t-feTime(idx))/(feTime(idx+1)-feTime(idx)) + b(:,idx);
        else
            ud = fe(idx);
            bval = b(:,idx);
        end
        
        % calculate control force
        uc = -(Rinv*ss.B'*PI + Rinv*S')*x - Rinv*ss.B'*bval;
        
        % eom
        dx = ss.A*x + ss.B*uc + ss.B*ud;
        
        % if we call with a flag, output the control, power gen, and bgen
        % y = [uc, power, bgen]^T
        if nargin == 3
            dx = [uc, -uc * x(end), uc/x(end)]';
        end
    end
end