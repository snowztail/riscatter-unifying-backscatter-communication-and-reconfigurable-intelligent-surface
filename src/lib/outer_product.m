function [outerProduct] = outer_product(matrix)
	% Function:
    %	- compute the outer product of all column vectors of the input matrix
    %
    % Input:
	%	- matrix: columns as the input vectors
    %
    % Output:
	%	- outerProduct: the outer product of input vectors
    %
    % Comment:
    %	- modified from https://stackoverflow.com/a/14356614/8891100
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 02

	nVectors = size(matrix, 2);
    outerProduct = matrix(:, 1);
	for iVector = 2 : nVectors
		vector = permute(matrix(:, iVector), circshift(1 : (ndims_modified(outerProduct) + 1), ndims_modified(outerProduct)));
		outerProduct = repmat(outerProduct, size(vector)) .* repmat(vector, size(outerProduct));
	end
end
