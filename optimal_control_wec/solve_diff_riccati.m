function PI = solve_diff_riccati(A, B, R, Q, PI_tf, t)

% calc R^(-1) - it gets called often, so precalc speeds things up
rInv = inv(R);

% solve the diff riccati eqn bvp
n = size(A,1);
nn = numel(A);
solinit = bvpinit(linspace(t(1), t(end)), ones(1, nn));
soln = bvp4c( @riccati, @resid, solinit);

% evaluate solution at requested t values
%PI = deval(soln, t);
PI = reshape(deval(soln, t), [n n length(t)]);

    % Riccati differential equation
    function dPi = riccati(t, pm)
        p = reshape(pm, [n n]);
        dPiMat = p*A + A'*p - p*B*rInv*B'*p + Q;
        dPi = reshape(dPiMat', [nn, 1]);
    end
    

    % Residual function
    function res = resid(ya, yb)
        res = yb;
    end

end
