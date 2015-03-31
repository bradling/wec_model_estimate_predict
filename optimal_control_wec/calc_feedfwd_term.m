function b = calc_feedfwd_term(At, B, Rinv, PI, tvals)

    b0 = zeros(size(At,1), 1);
    soln = ode45(@ffeq, [tvals(end) tvals(1)], b0);
    
    b = deval(soln, tvals);
    
    function db = ffeq(t,b)
        db = -(At - B*Rinv*B'*PI)' * b - PI*B;
    end
end
