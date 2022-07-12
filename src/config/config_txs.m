%% * System
% % number of transmit antennas
% nTxs = 4;
% number of tags
nTags = 2;
% number of available states at tags (i.e., modulation order)
nStates = 2;
% constellation diagram at tags
constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);
% carrier frequency
frequency = 2e8;
% amplitude scatter ratio at tags
scatterRatio = 0.5;
% backscatter/primary symbol duration ratio
symbolRatio = 20;
% average transmit power
transmitPower = db2pow(6);
% average noise power
noisePower = db2pow(-100);
% antenna gain
rxGain = db2pow(3);
% layout and distance
directDistance = 10;
coverage = 1;
[forwardDistance, backwardDistance] = layout(directDistance, nTags, coverage);
% path loss exponents
directExponent = 2.6;
forwardExponent = 2.4;
backwardExponent = 2;
% Ricean factors
directFactor = 5;
forwardFactor = 5;
backwardFactor = 10;

%% * Algorithm
% relative priority of primary link
weightSet = flip([1 : -0.1 : 0.5, 0.45 : -0.05 : 0.25, 0.2 : -0.025 : 0]);
% weightSet = ([1 : -0.1 : 0.5, 0.45 : -0.05 : 0.25, 0.2 : -0.025 : 0]);
% number of weights on rate region boundary
nWeights = length(weightSet);
% number of quantization bins
nBins = 2 ^ 6;

%% * Variable
% Variable = struct('nTxs', num2cell(2 .^ (1 : 3)));
Variable = struct('nTxs', num2cell(2 .^ (0 : 3)));
% Variable = struct('nTxs', num2cell(2 .^ (3 : -1 : 0)));
% Variable = struct('nTxs', {3, 4, 8});
% Variable = struct('nTxs', num2cell(8));
nVariables = length(Variable);

%% * PBS
% number of instances
nInstances = 3e2;
