function logical_indices = max_ind(input_array)
%MAX_IND yields logical indices of maximal array elements 
%   That's it pretty much...
max_value = max(input_array(:));
logical_indices = max_value == input_array;
end

