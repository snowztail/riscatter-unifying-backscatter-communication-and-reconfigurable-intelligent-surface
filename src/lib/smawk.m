function [p] = smawk(matrix, r, c)
	p = zeros(1, length(r));
	[c] = reduce(matrix, r, c);
	if length(r) == 1
		p = c;
	else
		[p(2 : 2 : end)] = smawk(matrix, r(2 : 2 : end), c);
		j = 1;
		for i = 1 : 2 : length(r)
			p(i) = c(j);
			if i < length(r)
				u = p(i + 1);
			else
				u = Inf;
			end
			while j <= length(r) && c(j) <= u
				if matrix(r(i), c(j)) < matrix(r(i), p(i))
					p(i) = c(j);
				end
				j = j + 1;
			end
			j = j - 1;
		end
	end
end


function [c] = reduce(matrix, r, c)
	i = 1;
	while length(r) < length(c)
		if matrix(r(i), c(i)) <= matrix(r(i), c(i + 1))
			if i < length(r)
				i = i + 1;
			elseif i == length(r)
				c(i + 1) = [];
			end
		else
			c(i) = [];
			if i > 1
				i = i - 1;
			end
		end
	end
end
