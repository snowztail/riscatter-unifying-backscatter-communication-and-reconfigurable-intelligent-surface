%% * System
% number of transmit antennas
nTxs = 4;
% number of scatter antennas
nSxs = 1;
% number of receive antennas
nRxs = 1;
% number of tags
nTags = 2;
% number of available states at tags (i.e., modulation order)
nStates = 4;
% constellation diagram at tags
constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);
% amplitude scatter ratio at tags
% scatterRatio = 0.5;
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
nBins = 2 ^ 8;

%% * Variable
Variable = struct('scatterRatio', num2cell(0.25 : 0.25 : 1));
nVariables = length(Variable);

%% * PBS
% number of instances
nInstances = 3e3;
