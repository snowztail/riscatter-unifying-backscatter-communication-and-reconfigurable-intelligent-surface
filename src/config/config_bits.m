%% * System
% number of transmit antennas
nTxs = 4;
% number of scatter antennas
nSxs = 1;
% number of receive antennas
nRxs = 1;
% number of tags
nTags = 8;
% number of available states at tags (i.e., modulation order)
nStates = 2;
% constellation diagram at tags
constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);
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
weight = 0;
% number of weights on rate region boundary
nWeights = 1;
% % number of quantization bins
% nBins = 2 ^ 9;

%% * Variable
Variable = struct('nBits', num2cell(8 : 10));
nVariables = length(Variable);

%% * PBS
% number of instances
nInstances = 3e3;
