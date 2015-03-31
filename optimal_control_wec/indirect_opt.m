function [ results ] = indirect_opt( A, B, D, x0, Q, R, S, fe, feTime, dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
results.bgen = -results.uc ./ results.zDot;
results.power = -results.zDot .* results.uc;


    function [dy] = bvp(t, y)
        % boundary value problem to solve resulting from FONCs
        % y = [x ; lambda]
        fe_i = cust_interp(t, feTime, fe, dt);
        
        dy = [ (A - (0.5).*B*Rinv*S')   , -(0.5).*B*Rinv*B' ; 
               ((0.5).*S*Rinv*S' - 2*Q) , (0.5).*S*Rinv*B'-A' ] * y ...
             + [D ; zeros(n,1)] * fe_i;
    end

    function r = res(yo, yf)
        r = [(yo(1:n)' - x0'),  yf(n+1:end)'];
    end

    function u = calc_control(y)
        u = nan(1, size(y,2));
        for ii = 1:size(y,2)
            u(ii) = -0.5.*[ Rinv*S' , Rinv*B' ] * y(:,ii);
        end
    end
end

