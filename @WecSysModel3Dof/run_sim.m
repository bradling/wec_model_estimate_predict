function results = run_sim(obj, eta, dt)
% run state space simulation
% dt must be the sampling rate of eta. In the future consider allowing a
% different output sampling rate


ss = obj.construct_state_space_model();
[fe, feTime, etaSim] = calc_excitation(obj, eta, dt);

% ramp up fe
%fe(1,:) = [linspace(0,1,5000*dt) ones(1, size(fe,2)-5000*dt)] .* fe(1,:);
%fe(2,:) = [linspace(0,1,5000*dt) ones(1, size(fe,2)-5000*dt)] .* fe(2,:);
%fe(3,:) = [linspace(0,1,5000*dt) ones(1, size(fe,2)-5000*dt)] .* fe(3,:);


simOut = ode45(@eom, feTime([1 end]), zeros(size(ss.a,1), 1));
zeta = deval(simOut, feTime);

results.t    = feTime - feTime(1);
results.eta  = etaSim;
results.s    = zeta(4:6,:);
results.sDot = zeta(1:3,:);
results.fe   = fe;

    function dy = eom(t, y)
        u = [ cust_interp(t, feTime, fe(1,:), dt) ;
              cust_interp(t, feTime, fe(2,:), dt) ;
              cust_interp(t, feTime, fe(3,:), dt) ];
        dy = ss.a * y + ss.b * u;
    end

end


function y = cust_interp(ti, t, y, dt)

idx = find(floor(ti/dt)*dt == t);
if (ti-t(idx)) > 1e-8
    y = (y(idx+1) - y(idx)) * (ti-t(idx)) / (t(idx+1) - t(idx)) + y(idx);
else
    y = y(idx);
end

end