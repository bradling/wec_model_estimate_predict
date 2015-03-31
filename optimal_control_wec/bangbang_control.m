function results = bangbang_control(A,B,D,S,x0,ulim,fe,feTime,dt)
% calculates the bang bang control for actuator limits and u = Fgen
%
% note ulim = [u_min u_max]


n = size(A,1);

% solve bvp
solinit = bvpinit(feTime, zeros(2*n,1));
soln = bvp4c(@bvp, @res, solinit);

y = deval(soln, feTime);
results.state = y(1:n,:);
results.z = y(n-1,:);
results.zDot = y(n,:);
results.uc = calc_control(y);
results.bgen = results.uc ./ results.zDot;
results.power = -results.zDot .* results.uc;

    
    % the function defining the ODE/BVP to be solved
    function dy = bvp(t,y)
        dy = nan(size(y));
        % calc fe
        feVal = cust_interp(t, feTime, fe, dt);
        
        % calc u - bang bang
        if (y(n) + y(end)) > 0, u = ulim(2);
        else u = ulim(1);
        end
        
        % differential eqn.
        dy(1:n) = A*y(1:n) + B*u + D*feVal;
        dy(n+1:end) = -A'*y(n+1:end) - S*u;
    end

    % residual function for boundary conditions
    function r = res(yi, yf)
        r = [(yi(1:n)' - x0'),  yf(n+1:end)'];
    end

    % for calculating control after solving bvp
    function get_u(state)
        u = ulim(2) .* ones(size(state,2));
        u( (state(n,:) + state(end,:)) > 0 ) = ulim(1);
    end

end
