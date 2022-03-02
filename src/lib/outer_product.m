function [outerProduct] = outer_product(varargin)
	% Function:
    %   - compute the outer product of arbitrary number of vectors
    %
    % Input:
	%	- varargin: the cell containing input vectors
    %
    % Output:
	%	- outerProduct: the outer product of input vectors
    %
    % Comment:
    %   - Modified from https://stackoverflow.com/a/14356614/8891100
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 02

    outerProduct = varargin{1};
	for iVector = 2 : numel(varargin)
		vector = permute(varargin{iVector}, circshift(1 : (ndims_modified(outerProduct) + 1), ndims_modified(outerProduct)));
		outerProduct = repmat(outerProduct, size(vector)) .* repmat(vector, size(outerProduct));
% 		outerProduct = outerProduct .* permute(varargin{iVector}, circshift(1 : (ndims_modified(outerProduct) + 1), ndims_modified(outerProduct)));
% 		outerProduct = outerProduct .* permute(varargin{iVector}, circshift(1 : (ndims_modified(outerProduct) + ndims_modified(varargin{iVector})), [0, ndims_modified(outerProduct)]));
	end
end
