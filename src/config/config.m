%% * System
% carrier frequency
frequency = 9e8;
% number of transmit antennas
nTxs = 3;
% number of tags
nTags = 2;
% number of available states at tags (i.e., modulation order)
nStates = 2;
% number of possible input and output tuples (i.e., alphabet size)
[nInputs, nOutputs] = deal(nStates ^ nTags);
% amplitude scatter ratio at tags
scatterRatio = 0.5;
% backscatter/primary symbol duration ratio
symbolRatio = 10;
% constellation diagram at tags
constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);
% constellation = transpose(exp(1i * 2 * pi * (0 : nStates - 1) / nStates));
% average transmit power
% transmitPower = db2pow(6);
transmitPower = 1e3;
% average noise power
noisePower = db2pow(-160);
% path loss exponents
directExponent = 2.6;
forwardExponent = 2.4;
backwardExponent = 2;
% angle of arrivals
directAoa = 2 * pi * rand;
forwardAoa = 2 * pi * rand(nTags, 1);
backwardAoa = 2 * pi * rand(nTags, 1);
% Ricean factors
directFactor = 5;
forwardFactor = 5;
backwardFactor = 5;

%% * Layout
% coordinates
ap = [5; 0];
user = [0; 0];
[tag(1, :), tag(2, :)] = pol2cart(2 * pi * rand(1, nTags), sqrt(1 * rand(1, nTags)));
% distances
directDistance = norm(ap - user);
forwardDistance = vecnorm(ap - tag);
backwardDistance = vecnorm(tag - user);

%% * Algorithm
% relative priority of primary link
weightSet = 0 : 0.05 : 1;
% number of weights on rate region boundary
nWeights = length(weightSet);
