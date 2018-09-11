function output_polyshape = get_unified_poly(T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
warning off MATLAB:polyshape:boundary3Points
output_polyshape = polyshape([0 0; 0 0; 0 0]);
for i = 1:numel(T)
    output_polyshape = union(output_polyshape,...
        polybuffer(T{i}.Shape,1e-5,'JointType','square'));
%     output_polyshape = polyshape([output_polyshape.Vertices; T{i}.Shape.Vertices]);
end
output_polyshape = simplify(polyshape(round(output_polyshape.Vertices,4)));
% [x,y] = boundary(output_polyshape);
% output_polyshape = polyshape(x,y);

% [xlim,ylim] = boundingbox(output_polyshape);
% CH = convhull(output_polyshape);
% diff = subtract(CH, output_polyshape);

% output_polyshape = polyshape([0 0; 0 0; 0 0]);
% for i = 1:numel(T)
%     output_polyshape = polyshape([output_polyshape.Vertices; T{i}.Shape.Vertices]);
% end
warning on MATLAB:polyshape:boundary3Points
end

