%% * System
% number of transmit antennas
nTxs = 1;
% number of scatter antennas
nSxs = 1;
% number of receive antennas
nRxs = 1;
% number of tags
nTags = 4;
% number of available states at tags (i.e., modulation order)
nStates = 2;
% constellation diagram at tags
constellation = pskmod(transpose(0 : nStates - 1), nStates);
% amplitude scatter ratio at tags
scatterRatio = 0.5;
% backscatter/primary symbol duration ratio
symbolRatio = 20;
% average transmit power
transmitPower = db2pow(6);
% average noise power
noisePower = db2pow(-70);
% layout and distance
directDistance = 10;
coverage = 2;
[forwardDistance, backwardDistance] = layout(directDistance, nTags, coverage);
% path loss exponents
directExponent = 2.6;
forwardExponent = 2.4;
backwardExponent = 2;
% Ricean factors
directFactor = 5;
forwardFactor = 5;
backwardFactor = 5;

%% * Algorithm
% relative priority of primary link
weightSet = [0 : 0.01 : 0.05, 0.075, 0.1 : 0.05 : 0.25, 0.3 : 0.1 : 0.5, 0.75, 1];
% number of weights on rate region boundary
nWeights = length(weightSet);
% number of quantization bins
nBins = 2 ^ 9;

%% * Variable
Variable = struct('Distribution', {'exhaustion', 'kkt', 'cooperation', 'cooperation', 'cooperation', 'cooperation', 'uniform'}, 'Recovery', {[], [], [], 'marginalization', 'decomposition', 'randomization', []});
nVariables = length(Variable);

%% * PBS
% number of instances
nInstances = 1e3;
