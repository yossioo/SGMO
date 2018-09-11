% clc;
clearvars -except DEBUG Fingers finger_d PolyList TableVarNames
Move_from_vertex_ratio = 0.05; % set to positive if want to stay away from the vertices
%% This scripts evaluates grasps obtained in previous steps
% The script iterates over finger groups and checks for feasible
% variants. For each group evaluation is performed, and contacts are
% shifted if needed

%% First we derive the number of configuration groups
% If there are more than 1 contact groups we test all
N_p_i = numel(PolyList);
Groups_num_per_p = zeros(1,N_p_i);

for p_i = 1:N_p_i
    Groups_num_per_p(p_i) = max(Fingers.ContactGroup(Fingers.PolygonNum == p_i));
end

TotalGroups = prod(Groups_num_per_p);
GroupCombinations = ones(TotalGroups,N_p_i);
if TotalGroups > 1
    for p_i = 1:N_p_i
        %     Groups_num(p_i) = max(Fingers.ContactGroup(Fingers.PolygonNum == p_i));
        groups = Fingers.ContactGroup(Fingers.PolygonNum == p_i);
        GroupCombinations(:,p_i) = repmat(groups,TotalGroups/numel(groups),1);
    end
end
clear p_i groups


%% Test whether a group forms force closure contact and attempt to fix it

% BestFingersGroup = table([],[],[],[],[],[],[],...
%         'VariableNames', TableVarNames);
Best_Group_GQM = 0;
U = Polygon_mkII(get_unified_poly(PolyList(:)));

for combination = GroupCombinations'
    
    SelectedFingers = table([],[],[],[],[],[],[],...
        'VariableNames', TableVarNames);
    for p_i = 1:N_p_i
        f_i = Fingers.PolygonNum == p_i & Fingers.ContactGroup == combination(p_i);
        SelectedFingers = vertcat(SelectedFingers, Fingers(f_i,:));  %#ok<*AGROW>
    end
    clear p_i
    contact_set = SelectedFingers.ContactVector;
    % First check whether there 4 or more contacts
    N_cs = numel(contact_set);
    if N_cs >= 4
        % Get wrench convex hull
        [W_CH,W] = W_CH_from_Contacts(contact_set,U.Center,sqrt(U.Area));
        %         K = convhull([W_CH'; [0 0 0]]);
        TR = delaunayTriangulation(W_CH');
        figure(300), clf
        quiver3(0*W(1,:),0*W(1,:),0*W(1,:),W(1,:),W(2,:),W(3,:),'AutoScale','off');
        hold on, grid on;
        quiver3(0*W_CH(1,:),0*W_CH(1,:),0*W_CH(1,:),W_CH(1,:),W_CH(2,:),W_CH(3,:),'AutoScale','off');
        tetramesh(TR,'FaceAlpha',.2)
        axis equal
        grid on
        zlabel('\tau_z','FontSize',20)
        xlabel('f_x','FontSize',20)
        ylabel('f_y','FontSize',20)
        points = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]*1e-3;
        if find(isnan(TR.pointLocation(points))) %ismember(size(W_CH,2)+1,K)
            % the C-Hull did not contain the origin
            %
            % We wish to iterate over the finger contacts to see
            % whether some can make the CH to contain the origin
            %%% If no contact adjustment can solve it , skip the
            %%% configuration (or save it for double contact test?)
            c_indices = 1:N_cs;
            for c_i = c_indices
                p_i = SelectedFingers.PolygonNum(c_i);
                cts = contact_set(c_indices~=c_i);% Select other contacts
                W_CH_temp = reorder_W(W_CH_from_Contacts(cts,U.Center,sqrt(U.Area)));
                c = contact_set(c_i); % Select this contact
                
                % Extracting the allowed region for given contact
                edgeRange = SelectedFingers.EdgeRange(c_i,:);
                
                p1 = PolyList{p_i}.point_from_edgePosition(SelectedFingers.EdgeNum(c_i),0);
                p2 = PolyList{p_i}.point_from_edgePosition(SelectedFingers.EdgeNum(c_i),1);
                w = [c.direction_vector(:), c.direction_vector(:);...
                    cross2d(p1 - U.Center, c.direction_vector(:))/sqrt(U.Area),...
                    cross2d(p2 - U.Center, c.direction_vector(:))/sqrt(U.Area)];
                
                modified_EGW = slice_cone(-W_CH_temp,w);
                if ~isempty(modified_EGW)
                    % Next step is extracting the allowed finger position along the edge
                    L = diff(w(3,1:2));
                    inner_p1 = (modified_EGW(3,1) - w(3,1))/L;
                    inner_p2 = (modified_EGW(3,2) - w(3,1))/L;
                    
                    modified_EGW_range = [inner_p1,inner_p2];
                    modified_EGW_range(modified_EGW_range>(1-Move_from_vertex_ratio)) = 1-Move_from_vertex_ratio;
                    modified_EGW_range(modified_EGW_range<Move_from_vertex_ratio) = Move_from_vertex_ratio;
                    
                    % We find optimal location by checking the maximum of
                    % the minimal of CH face distances. Each face is
                    % checked: normal direction created, and then distance
                    % from the origin is checked.
                    pos_range = linspace(0,1);
                    torq_range = linspace(w(3,1), w(3,2));
                    d_Full = zeros(size(W_CH_temp,2),numel(pos_range));
                    points_along_the_line = [repmat(c.direction_vector(:),1,100);torq_range];
                    
                    V1 = W_CH_temp(:,end)-points_along_the_line;
                    V2 = W_CH_temp(:,1)-points_along_the_line;
                    d_Full(size(W_CH_temp,2),:) = dot(-points_along_the_line, cross(V1,V2)/norm(cross(V1,V2)));
                    
                    for face_i = 1:size(W_CH_temp,2)-1
                        V1 = W_CH_temp(:,face_i)-points_along_the_line;
                        V2 = W_CH_temp(:,face_i+1)-points_along_the_line;
                        d_Full(face_i,:) = dot(-points_along_the_line, cross(V1,V2)/norm(cross(V1,V2)));
                    end
                    d = round(min(d_Full,[],1),3); % Can round this, to eliminate numerical errors.
                    optimal_pos = pos_range(find(d == max(d),1));
                    optimal_pos(optimal_pos>max(modified_EGW_range)) = max(modified_EGW_range);
                    optimal_pos(optimal_pos<min(modified_EGW_range)) = min(modified_EGW_range);
%                     optimal_pos(optimal_pos>max(edgeRange)) = max(edgeRange);
%                     optimal_pos(optimal_pos<min(edgeRange)) = min(edgeRange);
                    
                    
                    if max(d) > Best_Group_GQM
                        contact_point = PolyList{p_i}.point_from_edgePosition(Fingers.EdgeNum(c_i),optimal_pos);
                        
                        BestFingersGroup = SelectedFingers;
                        BestFingersGroup.OptimalPosition(c_i) = optimal_pos;
                        BestFingersGroup.ContactVector(c_i).point_on_the_line = contact_point;
                        
                        Best_Group_GQM = max(d);
                    end
                end
                
                
            end
            
            clearvars n e1 e2 w modified_EGW
        else
            BestFingersGroup = SelectedFingers;% The origin is inside
            Best_Group_GQM = GQM_from_W(W_CH_from_Contacts(contact_set));

        end
    else
        % There are less then 4 contacts; we need to find some more 
        
    end
    % less then 4
    % Add contacts to the unified polygon?
    
end

% Display images
if DEBUG
    
    t = linspace(0,2*pi);
    x = finger_d/2*cos(t); y = finger_d/2*sin(t);
    figure(96); clf
    PolyList{1}.plot(); hold on; axis equal; grid on;
    for i = 2:numel(PolyList)
        PolyList{i}.plot()
        text(PolyList{i}.Center(1)-5,PolyList{i}.Center(2),num2str(i))
    end
    axis manual
    p = zeros(numel(contact_set),2);
    for c_i = 1:numel(BestFingersGroup.ContactVector)
        f = BestFingersGroup.ContactVector(c_i);
        f.plot_contact('b')
        p(c_i,:) = f.get_finger_center(finger_d);
        c = fill(p(c_i,1)+x,p(c_i,2)+y,'b','FaceAlpha',.2);
        %     s = scatter(p(1),p(2),100,'b', 'filled','MarkerFaceAlpha',0.2);
%         f.draw_inf_line('k')
    end
end


