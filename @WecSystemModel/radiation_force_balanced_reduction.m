function [varargout] = radiation_force_balanced_reduction(obj, tol, varargin)



% if modelOrder < 0 || mod(modelOrder,1) ~=0
%     error('modelOrder must be a positive integer!')
% end
plotFlag = true;
if nargin > 2
    plotFlag = varargin{1};
end

% get the radiation IRF
dt = obj.radIrf.t(2) - obj.radIrf.t(1);

[sys,~, svh] = imp2ss(obj.radIrf.irf, dt, 1, 1, tol);

varargout{1} = sys;

obj.ssRad.A = sys.A;
obj.ssRad.B = sys.B;
obj.ssRad.C = sys.C;
obj.ssRad.D = sys.D;


if nargout > 1
    varargout{2} = svh;
end



% If we want to plot the comparison - do it here
if plotFlag == true
    
    [radApprox, tApprox] = impulse(sys, obj.radIrf.t);
    
    h = obj.disp_irf('rad');
    set(h, 'name', 'Radiation IRF SS Approximation')
    hold on
    plot(tApprox, radApprox, 'r*')
    hold off
    legend('Actual IRF', 'SS Approximation', 'location', 'northeast')
    title({'Radiation Force State Space Approximation', ...
        sprintf('dt: %.4f, tol: %.2e, order: %.0n', dt, tol, size(sys.A,1))});
end


