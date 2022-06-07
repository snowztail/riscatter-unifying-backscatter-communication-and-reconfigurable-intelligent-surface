%% * System
% number of transmit antennas
nTxs = 4;
% number of tags
nTags = 2;
% % number of available states at tags (i.e., modulation order)
% nStates = 2;
% % constellation diagram at tags
% constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);
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
weightSet = [0 : 1e-3 : 2e-2, 1];
% number of weights on rate region boundary
nWeights = length(weightSet);


%% * Variable
Variable.nStates = 2 .^ (1 : 4);
nVariables = length(Variable.nStates);


%% * PBS
nInstances = 3e2;
