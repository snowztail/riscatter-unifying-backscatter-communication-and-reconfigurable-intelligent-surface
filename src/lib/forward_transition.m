function [forwardTransition, sortedChannel, mapIndex] = forward_transition(symbolRatio, transmitPower, noisePower, equivalentChannel)
	% Function:
    %   - compute the channel forward transition probability for AmBC with energy detection
    %
    % Input:
	%	- symbolRatio: the (integer) ratio of backscatter symbol period to legacy symbol period
    %   - transmitPower: average transmit power at the AP
	%	- noisePower: average noise power at the user
	%	- equivalentChannel [nStates * 1]: the equivalent channel candidates corresponding to tag modulation states
    %
    % Output:
	%	- forwardTransition [nInputs * nOutputs]: the channel transition probability matrix
	%	- sortedChannel [nStates * 1]: sorted equivalent channel candidates with ascending receive energy
	%	- mapIndex: the index of mapping from constellation counterclockwise order to energy ascending order
    %
    % Comment:
    %   - for a specific channel realization, the forward transition probability matrix is obtained numerically
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 10

	% * Get data
	nStates = size(equivalentChannel, 1);
	[nInputs, nOutputs] = deal(nStates);

	% * Estimate the typical received signal energy levels
	energyLevel = energy_level(nOutputs, equivalentChannel, transmitPower, noisePower);

	% * Sort the energy level in ascending order
	[energyLevel, mapIndex] = sort(energyLevel);
	sortedChannel = equivalentChannel(mapIndex);

	% * Compute the energy detection threshold
	detectionThreshold = detection_threshold(symbolRatio, nOutputs, energyLevel);

	% * Calculate the forward transition probability
	forwardTransition = forward_transition_local(symbolRatio, nInputs, nOutputs, energyLevel, detectionThreshold);
end


function [energyLevel] = energy_level(nOutputs, equivalentChannel, transmitPower, noisePower)
	energyLevel = zeros(nOutputs, 1);
	for iOutput = 1 : nOutputs
		energyLevel(iOutput) = abs(equivalentChannel(iOutput)) ^ 2 * transmitPower + noisePower;
	end
end

function [detectionThreshold] = detection_threshold(symbolRatio, nOutputs, energyLevel)
	detectionThreshold = zeros(nOutputs + 1, 1);
	[detectionThreshold(1), detectionThreshold(end)] = deal(0, inf);
	for jOutput = 2 : nOutputs
		iOutput = jOutput - 1;
		detectionThreshold(jOutput) = symbolRatio * (energyLevel(iOutput) * energyLevel(jOutput)) / (energyLevel(iOutput) - energyLevel(jOutput)) * log(energyLevel(iOutput) / energyLevel(jOutput));
	end
end

function [forwardTransition] = forward_transition_local(symbolRatio, nInputs, nOutputs, energyLevel, detectionThreshold)
	forwardTransition = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergyDistribution = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ energyLevel(iInput))) ./ (energyLevel(iInput) .^ symbolRatio .* gamma(symbolRatio));
		% conditionalEnergyDistribution = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ energyLevel(iInput))) ./ (hpf(energyLevel(iInput) .^ hpf(symbolRatio) .* hpf(gamma(symbolRatio))));
		% conditionalEnergyDistribution = @(z) vpa(z .^ (symbolRatio - 1) .* exp(-z ./ energyLevel(iInput))) ./ (energyLevel(iInput) .^ symbolRatio .* gamma(symbolRatio));
		% conditionalEnergyDistribution = @(z) (z .^ vpa(symbolRatio - 1) .* exp(-z ./ vpa(energyLevel(iInput)))) ./ (vpa(energyLevel(iInput)) .^ vpa(symbolRatio) .* vpa(gamma(symbolRatio)));
		% conditionalEnergyDistribution = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ energyLevel(iInput))) ./ double(vpa(energyLevel(iInput)) .^ vpa(symbolRatio) .* vpa(gamma(symbolRatio)));
		% fplot(conditionalEnergyDistribution{iInput});
		for iOutput = 1 : nOutputs
			forwardTransition(iInput, iOutput) = integral(conditionalEnergyDistribution, detectionThreshold(iOutput), detectionThreshold(iOutput + 1));
		end
	end
end
