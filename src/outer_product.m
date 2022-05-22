function [outerProduct] = outer_product(matrix)
	nVectors = size(matrix, 2);
    outerProduct = matrix(:, 1);
	for iVector = 2 : nVectors
		vector = permute(matrix(:, iVector), circshift(1 : (dimension_number(outerProduct) + 1), dimension_number(outerProduct)));
		outerProduct = repmat(outerProduct, size(vector)) .* repmat(vector, size(outerProduct));
	end
end
