function [ind] = is_Polygon_in_array(P_mkII,Cell_Array)
%IS_POLYGON_IN_ARRAY Summary of this function goes here
%   Detailed explanation goes here
ind = [];
for i = 1:numel(Cell_Array)
    if P_mkII.Name == Cell_Array{i}.Name
        ind = i;
        break
    end
end

end

