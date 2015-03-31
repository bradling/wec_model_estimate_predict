function [x,y] = format_input(heave,inLen,outLen)
% Converts a time series dataset into input-output pairs for prediction
% Works on a single time series as well as a matrix of time series'
% NOTE - each state is a ROW
%
% [x, y] = format_input(heave, inLen, outLen)


% input is 60 data points
% output is 1 data points
% inLen  = 60;
% outLen = 1;


% Num of different states
numStates = size(heave,1);

% Number of points available
numPts = size(heave,2);

% Number of patterns in data set
numPatt = numPts - inLen - outLen + 1;


x = nan(inLen,numPatt * numStates);
y = nan(outLen,numPatt * numStates);
w = 1;
k = 1;
for s = 1:numStates
    for j = 1:numPatt
        i = j*k - k + 1;
        x(:,w) = heave(s,i:i+(inLen-1));
        y(:,w) = heave(s,(i+inLen):(i+inLen+outLen-1));
        w = w + 1;
    end
end