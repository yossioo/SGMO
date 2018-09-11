function output = unit_vector(vector_or_matrix, direction)
%UNIT_VECTOR Normalizes given vector or matrix
%   Direction is relevant for matrices, 1 normalizes columns; 2
%   normalizes rows
if nargin <2
    direction = 1;
end
if isvector(vector_or_matrix)
    % this is a vector
    output = vector_or_matrix./norm(vector_or_matrix);
else
    output = vector_or_matrix./vecnorm(vector_or_matrix,2,direction);
end
end

