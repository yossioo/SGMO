%%

% BACKUP
% THIS ITERATES OVER POLYGONS












% clc;
clearvars -except DEBUG Fingers finger_d PolyList TableVarNames
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

%% Test whether a group forms force closure contact

t = linspace(0,2*pi);
x = finger_d/2*cos(t); y = finger_d/2*sin(t);
U = Polygon_mkII(get_unified_poly(PolyList(:)));
for combination = GroupCombinations'
    contact_set = [];
    for p_i = 1:N_p_i
        f_i = Fingers.PolygonNum == p_i & Fingers.ContactGroup == combination(p_i);
        contact_set = vertcat(contact_set, Fingers.ContactVector(f_i));  %#ok<*AGROW>
    end
    % First check whether there 4 or more contacts
    N_cs = numel(contact_set);
    if N_cs >= 4
        % Get wrench convex hull
        [W_CH,W] = W_CH_from_Contacts(contact_set,U.Center,sqrt(U.Area));
        K = convhull([W_CH'; [0 0 0]]);
        figure(300), clf
        quiver3(0*W(1,:),0*W(1,:),0*W(1,:),W(1,:),W(2,:),W(3,:),'AutoScale','off');
        hold on, grid on; 
        quiver3(0*W_CH(1,:),0*W_CH(1,:),0*W_CH(1,:),W_CH(1,:),W_CH(2,:),W_CH(3,:),'AutoScale','off');
        
        if ismember(size(W_CH,2)+1,K)
            % the C-Hull did not contain the origin
            %
            % We wish to iterate over the finger contacts to see
            % whether some can make the CH to contain the origin
            %%% If no contact adjustment can solve it , skip the
            %%% configuration (or save it for double contact test?)
            p_indices = 1:N_p_i;
            for p_i = 1:N_p_i
                f_i = Fingers.PolygonNum == p_i & Fingers.ContactGroup == combination(p_i);
                cts = contact_set(p_indices~=p_i);% Select contacts for other polygons
                c = contact_set(p_indices==p_i); % Select contacts relevant to this polygon
                U = Polygon_mkII(get_unified_poly(PolyList(:)));
                W_CH_temp = W_CH_from_Contacts(cts,U.Center,sqrt(U.Area));
                ind_contact = Fingers.PolygonNum == p_i & Fingers.ContactGroup == combination(p_i);
                % Extracting the allowed region for given contact
                edgeRange = Fingers.EdgeRange(ind_contact,:);
                p1 = PolyList{p_i}.point_from_edgePosition(Fingers.EdgeNum(ind_contact),0);
                p2 = PolyList{p_i}.point_from_edgePosition(Fingers.EdgeNum(ind_contact),1);
                w = [c.direction_vector(:), c.direction_vector(:);...
                    cross2d(p1 - U.Center, c.direction_vector(:))/sqrt(U.Area),...
                    cross2d(p2 - U.Center, c.direction_vector(:))/sqrt(U.Area)];
                
                % Now we test whether given contact moved can complete
                % the grasp
                 modified_EGW = slice_cone(-W_CH_temp,w);
                 
                  if ~isempty(modified_EGW)
                        % Next step is extracting the allowed finger position along the edge
                        L = diff(w(3,1:2));
                        inner_p1 = (modified_EGW(3,1) - w(3,1))/L;
                        inner_p2 = (modified_EGW(3,2) - w(3,1))/L;

                        % We find optimal location by checking the maximum of
                        % the minimal of CH face distances. Each face is
                        % checked: normal direction created, and then distance
                        % from the origin is checked.
                        pos_range = linspace(inner_p1,inner_p2);
                        torq_range = linspace(modified_EGW(3,1), modified_EGW(3,2));
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
                        d = min(d_Full,[],1);
                        optimal_pos = pos_range(find(d == max(d),1));
                        optimal_pos(optimal_pos>max(edgeRange)) = max(edgeRange);
                        optimal_pos(optimal_pos<min(edgeRange)) = min(edgeRange);
                        if max(d)<= 0.001
                            continue
                        end
                        contact_point = PolyList{p_i}.point_from_edgePosition(Fingers.EdgeNum(p_i),optimal_pos);
                        Fingers.OptimalPosition(f_i) = optimal_pos;
%                         f =  ContactVector(contact_point,...
%                             PolyList{p_i}.find_normal_at_point(contact_point),finger_d,p_i);
                        Fingers.ContactVector(f_i).point_on_the_line = contact_point;
                        Fingers.Group_GQM(f_i) = max(d);
                  end
                  
                  clearvars n e1 e2 w modified_EGW
                  
            end
        else
            % The CH contains the origin and hence the grasp is
            % immobilizing.
        end
    else
        % less then 4
        % Add contacts to the unified polygon?
        
    end
    PolyList{1}.plot(); hold on; axis equal; grid on;
    for i = 2:numel(PolyList)
        PolyList{i}.plot()
        text(PolyList{i}.Center(1)-5,PolyList{i}.Center(2),num2str(i))
    end
    axis manual
    p = zeros(numel(contact_set),2);
    for c_i = 1:numel(Fingers.ContactVector)
        f = Fingers.ContactVector(c_i);
        f.plot_contact('b')
        p(c_i,:) = f.get_finger_center(finger_d);
        c = fill(p(c_i,1)+x,p(c_i,2)+y,'b','FaceAlpha',.2);
        %     s = scatter(p(1),p(2),100,'b', 'filled','MarkerFaceAlpha',0.2);
        f.draw_inf_line('k')
    end

end
