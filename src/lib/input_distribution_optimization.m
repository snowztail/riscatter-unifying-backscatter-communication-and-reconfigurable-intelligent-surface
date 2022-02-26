function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_optimization(nTags, dmtc, weight, symbolRatio, snr, tolerance)
	% Function:
	%	- optimize the tag input distribution to characterize the capacity region of user and tags
    %
    % Input:
	%	- nTags: number of tags
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight [(nTags + 1) * 1]: the relative priority of the user and tags
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (bps/Hz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (bpcu)
    %
    % Comment:
    %   - Obtain the joint input probability matrix by optimization (corresponding to the optimal input with full transmit cooperation)
	%	- Extract rank-1 solution from the joint input probability matrix by randomization or marginalization
	%	- the discrete memoryless MAC is given in joint (equivalent point-to-point) form P(y | x_1, ..., x_K), instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 23

	% * Declare default tolerance
	arguments
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
		tolerance = eps;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);
	nOutputs = size(dmtc, 2);

	% * Initialization
	indexCombination = combvec_nested(1 : nStates, nTags);
	powerSet = power_set(1 : nTags);
	nCases = length(powerSet);
	inputArraySet = cell(nCases, 1);
	rateBound = cell(nCases, 1);

	cvx_begin
		% * Create variables
		for iCase = 1 : nCases
			inputArraySet{iCase} = strcat('inputArray', sprintf('%d', powerSet{iCase}));;
			eval(['variable', ' ', inputArraySet{iCase}, '(', num2str(repmat(nStates,[1, length(powerSet{iCase})])), ')']);
		end
		variables rate(nTags, 1);

		% * Express upper bounds of sum rate on all subsets
		jointArray = eval(inputArraySet{iCase});
		for iCase = 1 : nCases
			if iCase < nCases
				complementArray = eval(inputArraySet{nCases - iCase});
			else
				complementArray = [];
			end
			rateBound{iCase} = rate_bound(dmtc, powerSet{iCase}, jointArray, complementArray);








			eval([rateBound{iCase} '= [2 3 4 5]']);
		end

		for iOutput = 1 : nOutputs
			pp12 = transpose(vec(transpose(P))) * dmtc(:, iOutput);
% 			f12t2(iOutput) = rel_entr(pp12, 1);
			f12t2(iOutput) = - entr(pp12);
			ft1(iOutput) = transpose(vec(transpose(P))) * entr(dmtc(:, iOutput));

			for iState = 1 : nStates
				iIndex = indexCombination(1, :) == iState;
				pp2 = P(iState, :) * dmtc(iIndex, iOutput);
				f2t2(iState, iOutput) = rel_entr(pp2, p1(iState));
			end

			for jState = 1 : nStates
				jIndex = indexCombination(2, :) == jState;
				pp1 = transpose(P(:, jState)) * dmtc(jIndex, iOutput);
				f1t2(jState, iOutput) = rel_entr(pp1, p2(jState));
			end
		end
		f1 = - sum(ft1) - sum(sum(f1t2));
		f2 = - sum(ft1) - sum(sum(f2t2));
		f12 = - sum(ft1) - sum(f12t2);

		maximize R1 + R2
		subject to
				0 <= R1 <= f1;
				0 <= R2 <= f2;
				R1 + R2 <= f12;
				P * ones(nStates, 1) == p1;
				P' * ones(nStates, 1) == p2;
				ones(1, nStates) * P * ones(nStates, 1) == 1;
				P == semidefinite(nStates);
	cvx_end
end
