function [rsq] = regression_analysis(results,targets,varargin)
% assumes output and target are 1-D arrays of the same dimension
% Calculates the correlation coefficient between results and targets.
% NOTE: function is designed for calculating regression coefficients as an
% error measure, therefore the regression model assumes f(x) = x + 0.
%
% set plotFlag = true to get scatter plot. If this input is omitted,
% default is to not plot results.
% 
% [rsq] = regression_analysis(results,targets,plotFlag)

% Change history:
% 10-7-14: Changed regression to assume f(x) = x


if nargin == 2
    plotflag = false;
else
    plotflag = varargin{1};
end


% Assume F(x) = x:
yresid = targets - results;
SSresid = sum(yresid.^2);
SStotal = (length(targets) - 1) * var(targets);
rsq = 1 - SSresid/SStotal;


if plotflag == true
    figure
    set(gcf,'color','w')
    
    x = linspace(min(min([results targets])), max(max([results targets])),10);
    plot(targets,results,'.', x,polyval([1 0],x),'r', 'linewidth',3);
    xlabel('Observed Values')
    ylabel('Predicted Values')
    legend('Predictions','Linear Regression','location','best')
    grid on
    title(sprintf('R^{2} = %.3f',rsq));

end
