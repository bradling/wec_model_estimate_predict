function results = indirect_with_limits( A, B, D, x0, Q, R, S, umin, umax, fe, feTime, dt )



% preliminaries
Rinv = inv(R);
n = size(A,1);

% solve bvp
initGuess = [zeros(n,1) 100*ones(n,1)];
solinit = bvpinit(feTime, zeros(2*n,1));
soln = bvp4c(@bvp, @res, solinit);

y = deval(soln, feTime);
results.state = y(1:n,:);
results.z = y(n-1,:);
results.zDot = y(n,:);
results.uc = calc_control(y);
results.ucSat = calc_control_sat_bgen(y);
results.bgen = -results.uc ./ results.zDot;
results.power = -results.zDot .* results.uc;


[~, y2] = ode45(@eom, feTime, x0);

results.power2 = -y2(:,n)' .* results.ucSat;
results.z2 = y2(:,n-1)';
results.zDot2 = y2(:,n)';
results.bgen2 = -results.ucSat ./ results.zDot2;

    function [dy] = bvp(t, y)
        % boundary value problem to solve resulting from FONCs
        % y = [x ; lambda]
        dy = nan(size(y));
        fe_i = cust_interp(t, feTime, fe, dt);
        
        u = calc_u(y(1:n), y(n+1:end));
        dy(1:n) = A*y(1:n) + B*u + D * fe_i;
        dy(n+1:end) = -2*Q*y(1:n) - 2*S*u - A'*y(n+1:end);
    end

    function dy = eom(t, x)
        fe_i = cust_interp(t, feTime, fe, dt);
        u = cust_interp(t, feTime, results.ucSat, dt);
        
        dy = A*x + B*u + D*fe_i;
    end

    function u = calc_u(x, lambda)
        % nominal value
        u = -Rinv*(x'*S + 0.5*B'*lambda);
        % impose limits
        if     u > umax, u = umax;
        elseif u < umin, u = umin; 
        end
    end

    function r = res(yo, yf)
        r = [(yo(1:n)' - x0'),  yf(n+1:end)'];
    end

    function u = calc_control(y)
        u = nan(1, size(y,2));
        for ii = 1:size(y,2)
            u(ii) = calc_u(y(1:n,ii), y(n+1:end,ii));
        end
    end

    function u = calc_control_sat_bgen(y)
        max_bgen = 4e7;
        u = nan(1, size(y,2));
        for ii = 1:size(y,2)
            u(ii) = calc_u(y(1:n,ii), y(n+1:end,ii));
            if abs(u(ii) ./ y(n, ii)) > max_bgen
                u(ii) = max_bgen * y(n, ii) * sign(u(ii));
            end
        end
    end
end

