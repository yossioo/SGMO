%% This script implements the search for new configuration of the objects
warning off MATLAB:polyshape:boundary3Points
%#ok<*SAGROW>
 %#ok<*UNRCH>
% clc;
clearvars -except P finger_d Filtered_Contacts_poly_ind Filtered_Contacts DEBUG

STACK_SHARPEST = true; % false
STACK_WIDEST = false;% true
STACK_CLOSE_TO_ORIGIN = false; % Implement this someday
FILL_CONCAVITY_1by1 = false;
FILL_CONCAVITY_THAT_FITS = true;
%% In case we have one object only
if numel(P) == 1
    PolyList{1} = P{1};
    return
end
%% First we wish to find the root object:
max_concave_angle = 0;
Root_Object = [];

for object = P
    p = object{1};
    conc_vert_angles = 360-p.V_angles(p.Concave_ind);
    if max(conc_vert_angles) > max_concave_angle
        Root_Object = p;
        PolyList{1} = p;
    end
end


if isempty(Root_Object) % no concave vertices found
    % Look for couple of inner angles higher than 90
    poly_list = {};
    vert_list = {};
    for object = P
        p = object{1};
        if find(p.V_angles > 90)
            poly_list{end+1} = p;
            diff_w_90 = p.V_angles - 90;
            diff_w_90(diff_w_90<0) = 180;
            ind = find(min(diff_w_90) == diff_w_90);
            vert_list{end+1} = ind(1);
        end
    end
    if numel(poly_list) > 1 % Found 2 such vertices
        % We wish to align the vertices together and rotate objects to
        % fit ------- What about case with more than 2 such objects ?
        edge_lengths_1 = repmat(poly_list{1}.E_lengths,3,1);
        edge_lengths_2 = repmat(poly_list{2}.E_lengths,3,1);
        adj_edge(1,1) = edge_lengths_1(vert_list{1}+poly_list{1}.N_e-1);
        adj_edge(1,2) = edge_lengths_1(vert_list{1}+poly_list{1}.N_e);
        adj_edge(2,1) = edge_lengths_2(vert_list{2}+poly_list{2}.N_e-1);
        adj_edge(2,2) = edge_lengths_2(vert_list{2}+poly_list{2}.N_e);
        
        [i,j] = find(adj_edge == max(adj_edge(:)));
        i = i(1);  % In case some lengths are equal
        j = j(1);
        Root_Object = poly_list{i(1)};
        PolyList(1) = poly_list(i);
        PolyList(2) = poly_list(3-i);
        
        V1 = PolyList{1}.Edges(vert_list{i},:);
        V2 = PolyList{2}.Edges(vert_list{3-i},:);
        
        if j == 2% Longest edge is the following edge on the root object
            edge_base_direction = diff(PolyList{1}.Edges([vert_list{i},vert_list{i}+1],:));
            
            % select leading edge on child object
            if vert_list{3-i} == 1
                edge_2_align_direction = diff(PolyList{2}.Edges([end,end-1],:));
            else
                edge_2_align_direction = diff(PolyList{2}.Edges([vert_list{3-i},vert_list{3-i}-1],:));
            end
        else % j==1 - on root object longest edge is the leading edge
            
            if vert_list{i} == 1
                edge_base_direction = diff(PolyList{1}.Edges([end,end-1],:));
            else
                edge_base_direction = diff(PolyList{1}.Edges([vert_list{i},vert_list{i}-1],:));
            end
            edge_2_align_direction = diff(PolyList{2}.Edges([vert_list{3-i},vert_list{3-i}+1],:));
        end
        % First we rotate about the contacting vertex
        PolyList{2}.rotate(-angle_of_2_vec(edge_base_direction, edge_2_align_direction), V2)
        
        % Then translate
        PolyList{2}.translate(V1-V2);
        
        
    else % No luck
        % We will look for 2 longest edges and stack them together
        % Find longest edges
        PolyList = {};
        e_lengths = [0 0];
        e_ids= [0 0];
        for object = P
            p = object{1};
            max_e_length = max(p.E_lengths);
            max_e_id = find(max_e_length == p.E_lengths);
            
            if isempty(PolyList)
                PolyList{1} = p;
                e_lengths(1) = max_e_length;
                e_ids(1) = max_e_id(1);
            else
                if max_e_length>e_lengths(1)
                    PolyList{2} = PolyList{1};
                    PolyList{1} = p;
                    e_lengths(2) = e_lengths(1) ;
                    e_ids(2) = e_ids(1);
                    e_lengths(1) = max_e_length;
                    e_ids(1) = max_e_id(1);
                elseif max_e_length>e_lengths(2)
                    PolyList{2} = p;
                    e_lengths(2) = max_e_length;
                    e_ids(2) = max_e_id(1);
                end
            end
        end
        
        % Now we have 2 objects with longest edges
        V1 = mean(PolyList{1}.Edges([e_ids(1),e_ids(1)+1],:));
        base_direction = PolyList{1}.find_vector_from_vertex(e_ids(1),1);
        align_direction = -PolyList{2}.find_vector_from_vertex(e_ids(2),1);
        V2 = PolyList{2}.Edges(e_ids(2),:);
        PolyList{2}.rotate(-angle_of_2_vec(base_direction, align_direction), V2);
        PolyList{2}.translate(V1-V2);
    end
    
end

%% Stacking process
clear p
% Remove items that are stacked from the initial cell array
for i = 1:numel(PolyList)
    ind = is_Polygon_in_array(PolyList{i},P);
    if ind
        P(ind) = [];
        disp(PolyList{i}.Name + " removed from P")
    end
end

P2={};
% disp('P contains:')
% for i = 1:numel(P)
%     fprintf("P{%d}: %s\n",i,P{i}.Name)
% end

while ~isempty(P)
    % Unify the shape
    uni = get_unified_poly(PolyList);
    UnifiedPolygon = Polygon_mkII(round(uni.Vertices,3));
    outer_angles = 360-UnifiedPolygon.V_angles;
%     vertex_distances_from_origin = vecnorm(uni.Vertices,2,2);
    outer_angles(~UnifiedPolygon.Concave_ind) = 0;
    conc_angles = 360-UnifiedPolygon.V_angles(UnifiedPolygon.Concave_ind);
    biggest_conc_ind = find(max(outer_angles) == outer_angles,1);
    outer_angles(~UnifiedPolygon.Concave_ind) = 1e3;
        
    if STACK_SHARPEST
        % Search for shape with sharpest angle
        min_angle = 180;
        p_i = 0;
        v_ind = 0;
        for i = 1:numel(P)
            if min(P{i}.V_angles) < min_angle
                p_i = i;
                min_angle = min(P{i}.V_angles);
                v_ind = find(P{i}.V_angles == min_angle,1);
            end
        end

    elseif STACK_WIDEST
        % Search for shape with widest angle
        max_angle = 0;
        p_i = 0;
        v_ind = 0;
        for i = 1:numel(P)
            if max(P{i}.V_angles) > max_angle
                p_i = i;
                v_ind = find(max(P{i}.V_angles) > max_angle,1);
                max_angle = max(P{i}.V_angles);
            end
        end
    end
    
    if FILL_CONCAVITY_1by1
        
    elseif FILL_CONCAVITY_THAT_FITS
        
        if STACK_SHARPEST
            if find(min_angle<=conc_angles)
                concavity_left_after_stacking = outer_angles-min_angle;
                concavity_left_after_stacking(concavity_left_after_stacking<0) = 1e3;
                conc_ind = find(concavity_left_after_stacking == min(concavity_left_after_stacking),1);
                stack_poly_to_vertex_w_edge(UnifiedPolygon, P{p_i}, conc_ind, v_ind,1)
            else
                % Push back the object, maybe later will be some concavity to
                % use.
                P2{end+1} = P{p_i};
                P(p_i) = [];
            end
        elseif STACK_WIDEST
            if find(max_angle<=conc_angles)
                % Can fit a poly somewhere
                concavity_left_after_stacking = outer_angles-max_angle;
                concavity_left_after_stacking(concavity_left_after_stacking<0) = 1e3;
                conc_ind = find(concavity_left_after_stacking == min(concavity_left_after_stacking),1);
                stack_poly_to_vertex_w_edge(UnifiedPolygon, P{p_i}, conc_ind, v_ind,1)
            else
                % Push back the object, maybe later will be some concavity to
                % use.
                P2{end+1} = P{p_i};
                P(p_i) = [];
            end
        end
        PolyList{end+1} = P{p_i};
        PolyList{end}.Name = strcat(num2str(numel(PolyList))," - ",PolyList{end}.Name );
        P(p_i) = [];
        
    end
    if numel(P) > 0
        disp("Not done yet")
        disp('P contains:')
        for i = 1:numel(P)
            fprintf("P{%d}: %s\n",i,P{i}.Name)
        end
    end
    if DEBUG
        f = figure(14);     clf;
        f.Name = "Stacking process";
        PolyList{1}.plot(PolyList{1}.Polygon_color,false); hold on;
        
        axis equal
        grid on
        pause(0.5)
        for i = 2:numel(PolyList)
            PolyList{i}.plot(PolyList{i}.Polygon_color,false)
            
            axis equal
            grid on
            pause(0.5)
        end
        axis equal
        grid on
        hold off
        %     drawnow
    end
end

while ~isempty(P2)
    % If P2 is not empty, there are some contacts that did not fit in
    % any concavity.
    % We try different way of stacking.
    break
end
% So far so good. Think over the logic of stacking inside of most
% fitting concavity, or the largest one or smth.


