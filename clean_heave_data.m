function [idx] = clean_heave_data()
% this function loads the awacs data and gets rid of bad timeseries values



dataFile = '../2015_disturbance_estimation_prediction/wave_data/awac_timeSeriesData.mat';
data = load(dataFile);

% Define sea state
eta = data.awacs.heave(:,:);


diffEta = diff(eta,1,2);
maxDiff = max(diffEta,[],2);

idxGood = abs(maxDiff) < (1.5 .* data.awacs.Hm0);

idx = 1:size(eta,1);
idx = idx(idxGood);


