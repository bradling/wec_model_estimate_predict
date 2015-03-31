function [ssRad, order] = rad_ss_approx(obj, order)
% calculates the approximate radiation state space realization

nTerms = size(obj.radFreq, 2);
ssRad = cell(1, nTerms);
dt = obj.radIrf.dt;

for ii = 1:nTerms
    ss = imp2ss(dt .* obj.radIrf.irf(ii, :), -dt, 1, 1);
    ssc = d2c(ss, 'tustin');
    ssRad{ii} = balred(ssc, order);
    clear ss
end