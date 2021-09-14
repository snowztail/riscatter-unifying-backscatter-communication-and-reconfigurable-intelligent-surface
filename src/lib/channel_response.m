function [channel] = channel_response(nTxs, nRxs, distance, frequency, propagationMode)
    % Function:
    %   - generate channel frequency response at specified distance and frequency
    %
    % Input:
    %   - nTxs: number of transmit antennas
    %   - nRxs: number of receive antennas
    %   - distance: distance between the transmitter and receiver
    %   - frequency: the carrier frequency
    %   - propagationMode: 'los' | 'nlos'
    %
    % Output:
    %   - channel [nTxs * nRxs]: one realization of channel frequency response
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 14

	% * Compute path loss by piecewise model
    pathLoss = path_loss(distance, frequency);

	% * Generate tap coefficients by TGn model
    [tapGain, tapDelay] = tap_tgn(nTxs, nRxs, propagationMode);

	% * Compute fading by tapped delay line model
	fading = fading_local(nTxs, nRxs, tapGain, tapDelay, frequency);

	% * Combine for channel response
    channel = fading / sqrt(pathLoss);
end


function [pathLoss] = path_loss(distance, frequency)
    % Function:
    %   - simulate the path loss at specified distance and frequency
    %
    % Input:
    %   - distance: distance between the transmitter and receiver
    %   - frequency: the carrier frequency
    %
    % Output:
    %   - pathLoss: large-scale signal attenuation
    %
    % Comment:
    %   - use piecewise path loss model of IEEE TGn channel model D
    %   - consists of a free-space model up to a breakpoint distance and typical urban loss model onwards
    %
    % Reference:
    %   - V. Erceg et al., "TGn channel models," in Version 4. IEEE 802.11–03/940r4, May 2004.
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 14

	% * Define the breakpoint distance and path loss exponents
	breakpointDistance = 10;
	freeSpaceExponent = 2;
	urbanExponent = 3.5;

	% * Compute path loss based on the piecewise model
    if distance <= breakpointDistance
        pathLoss = (4 * pi * distance * frequency / physconst('lightspeed')) ^ freeSpaceExponent;
    else
        pathLoss = (4 * pi * breakpointDistance * frequency / physconst('lightspeed')) ^ freeSpaceExponent * (distance / breakpointDistance) ^ urbanExponent;
    end
end

function [tapGain, tapDelay] = tap_tgn(nTxs, nRxs, propagationMode)
    % Function:
    %   - generate tapped delay line based on IEEE TGn channel model D
    %
    % Input:
    %   - nTxs: number of transmit antennas
    %   - nRxs: number of receive antennas
    %   - propagationMode: 'los' | 'nlos'
    %
    % Output:
    %   - tapGain {nTaps * 1}[nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps * 1]: tap delay
    %
    % Comment:
    %   - use power delay profile and path loss of IEEE TGn channel model D
    %   - for each tap, the LoS component is fixed while the NLoS component follows i.i.d. CSCG distribution with zero-mean and tap-power variance
    %   - LoS and NLoS have different Ricean factors
	%	- entries of LoS matrix are of unit-modulus
    %
    % Reference:
    %   - V. Erceg et al., "TGn channel models," in Version 4. IEEE 802.11–03/940r4, May 2004.
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 14

    % * Declare parameters
    nClusters = 3;
    nTaps = 18;

	% * Define tap delay and tap power
    tapDelay = 1e-9 * transpose([0 10 20 30 40 50 60 70 80 90 110 140 170 200 240 290 340 390]);
    tapPower = zeros(nTaps, nClusters);
    tapPower(:, 1) = db2pow(transpose([0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -9.0 -11.1 -13.7 -16.3 -19.3 -23.2 -inf -inf]));
    tapPower(:, 2) = db2pow(transpose([-inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -6.6 -9.5 -12.1 -14.7 -17.4 -21.9 -25.5 -inf]));
    tapPower(:, 3) = db2pow(transpose([-inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -18.8 -23.2 -25.2 -26.7]));

	% * Define Ricean factor
	switch propagationMode
	case 'los'
		riceanFactor = db2pow([3; -inf(nTaps - 1, 1)]);
	case 'nlos'
		riceanFactor = db2pow(-inf(nTaps, 1));
	end

    % * Generate LoS matrix
    losMatrix = exp(1i * 2 * pi * rand(nTxs, nRxs));

	% * Generate tap gains and sum over clusters
	tapGain = cell(nTaps, 1);
	for iTap = 1 : nTaps
		tapGain{iTap} = zeros(nTxs, nRxs);
		for iCluster = 1 : nClusters
			losGain = sqrt(riceanFactor(iTap) / (riceanFactor(iTap) + 1)) * losMatrix;
			nlosGain = sqrt(1 / (riceanFactor(iTap) + 1)) * sqrt(1 / 2) * (randn(nTxs, nRxs) + 1i * randn(nTxs, nRxs));
			tapGain{iTap} = tapGain{iTap} + sqrt(tapPower(iTap, iCluster)) * (losGain + nlosGain);
		end
	end
end

function [fading] = fading_local(nTxs, nRxs, tapGain, tapDelay, frequency)
    % Function:
	%	- compute small-scale fading at specified frequency
    %
    % Input:
    %   - nTxs: number of transmit antennas
    %   - nRxs: number of receive antennas
    %   - tapGain {nTaps * 1}[nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps * 1]: tap delay
    %   - frequency: the carrier frequency
    %
    % Output:
    %   - fading [nTxs * nRxs]: the small-scale fading
    %
    % Comment:
    %   - based on the tapped delay line model
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 14

	% * Declare variables
	nTaps = size(tapGain, 1);
	fading = zeros(nTxs, nRxs);

	% * Compute fading based on the tapped delay line model
	for iTx = 1 : nTxs
		for iRx = 1 : nRxs
			for iTap = 1 : nTaps
				fading(iTx, iRx) = fading(iTx, iRx) + tapGain{iTap}(iTx, iRx) * exp(1i * 2 * pi * frequency * tapDelay(iTap));
			end
		end
	end
end
