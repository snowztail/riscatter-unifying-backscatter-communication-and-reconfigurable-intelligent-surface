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
frequency = 9e8;
% amplitude scatter ratio at tags
scatterRatio = 0.5;
% backscatter/primary symbol duration ratio
symbolRatio = 10;
% average transmit power
transmitPower = db2pow(6);
% average noise power
noisePower = db2pow(-100);
% antenna gain
rxGain = db2pow(3);
% layout and distance
[directDistance, forwardDistance, backwardDistance] = layout(nTags);
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
weightSet = flip([0 : 1e-2 : 5e-2, 1e-1 : 5e-2 : 2.5e-1, 3e-1 : 1e-1 : 1]);
% number of weights on rate region boundary
nWeights = length(weightSet);


%% * Variable
Variable.nTxs = 2 .^ (0 : 3);
nVariables = length(Variable.nTxs);


%% * PBS
% number of instances
nInstances = 1e3;
