function stack_poly_to_vertex_w_edge(poly_parent,poly_child,vertex_parent_id,vertex_child_id, forward)
%STACK_POLY_TO_VERTEX_W_EDGE Summary of this function goes here 
%   Stacks a child object given parent object. 
%   Poly parent/child - P_mkII objects
%   IDs - index of vertices on each polygons to coincide
%   Side - if > 0, child will lean on following parent's edge, else
%   will lean on leading parents edge.

V1 = poly_parent.Edges(vertex_parent_id,:);
V2 = poly_child.Edges(vertex_child_id,:);

if forward > 0
    base_direction = poly_parent.find_vector_from_vertex(vertex_parent_id,1);
    align_direction = poly_child.find_vector_from_vertex(vertex_child_id,-1);
else
    base_direction = poly_parent.find_vector_from_vertex(vertex_parent_id,-1);
    align_direction = poly_child.find_vector_from_vertex(vertex_child_id,1);
end

% V1 = mean(PolyList{1}.Edges([e_ids(1),e_ids(1)+1],:));
% base_direction = PolyList{1}.find_vector_from_vertex(e_ids(1),1);
% align_direction = -PolyList{2}.find_vector_from_vertex(e_ids(2),1);
% V2 = PolyList{2}.Edges(e_ids(2),:);

poly_child.rotate(-angle_of_2_vec(base_direction, align_direction), V2);
poly_child.translate(V1-V2);

end

