%% * System
% number of transmit antennas
nTxs = 2;
% number of tags
nTags = 1;
% modulation order
nStates = 4;
% input and output alphabet size
[nInputs, nOutputs] = deal(nStates ^ nTags);
% harvest-backscatter ratio
reflectRatio = 1e-4;
% symbol period ratio
symbolRatio = 10;
% average transmit power
txPower = 10;
% average noise power
noisePower = db2pow(-70);
% backscatter symbol candidates
constellation = qammod(0 : nStates - 1, nStates) ./ max(abs(qammod(0 : nStates - 1, nStates)));

%% * Algorithm
% number of weights
nWeights = 20;
% primary-backscatter weight pairs
weightSet = [linspace(0, 0.1, nWeights - 1), 1];
% number of output discretization bins
nBins = 2 ^ 8;
% minimum rate gain per iteration
tolerance = 1e-6;
