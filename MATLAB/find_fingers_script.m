% clc;
clearvars -except DEBUG PolyList Filtered_Contacts Filtered_Contacts_poly_ind finger_d
warning off MATLAB:polyshape:boolOperationFailed
%#ok<*AGROW>

Move_from_vertex_ratio = 0.1; % set to positive if want to stay away from the vertices
%%
TableVarNames = {'PolygonNum','ContactGroup','EdgeNum','EdgeRange','OptimalPosition','Group_GQM','ContactVector'};
Fingers = table([],[],[],[],[],[],[],...
    'VariableNames', TableVarNames);

% Serves to test whether the origin is inside of a polygon
dx = 1e-2;
Qx = dx*cosd([0;120;240]);
Qy = dx*sind([0;120;240]);
% Iterating over polygons
for p_i = 1:numel(PolyList)
    p = PolyList{p_i};
    
    
    % Creating wrench space vectors for relevant contacts
    I_contacts = Filtered_Contacts(Filtered_Contacts_poly_ind == p_i);
    N = numel(I_contacts);
    W = zeros(3,N);
    for c_i = 1:N
        W(:,c_i) = round([I_contacts{c_i}.direction_vector(:);
            I_contacts{c_i}.cross_around_point(p.Center) / p.Area^0.5],2);
    end
    try
        TR = delaunayTriangulation(W');
        [~,P] = freeBoundary(TR);
        FullTR = delaunayTriangulation(P);
        origin_is_already_in_the_CH  = pointLocation(FullTR,[0 0 0]);
        if origin_is_already_in_the_CH && GQM_from_W(W) > 1e-2
            continue
        end
    catch
    end
    if DEBUG
        figure(16);clf; %
        p.plot(); hold on; axis equal; grid on;
        for c_i = 1:N
            I_contacts{c_i}.plot_contact()
        end
    end
    
    % Now we wish to determine whether all vectors lie in the same
    % plane. If dot projections of 3rd and next vectors on the normal
    % of the plane of first 2 - all vectors in plane have size 0 then
    % all vectors are coplanar and  2 marginal can be selected.
    cp = zeros(1,size(W,2)); % cross product
    dp = zeros(1,size(W,2)); % dot  product
    for w_i = 3:size(W,2)
        dp(w_i) = round(dot(cross(W(:,1),W(:,2)),W(:,w_i)),2);
        cp(w_i) = round(dot(cross(W(:,1),W(:,2)),cross(W(:,1),W(:,w_i))),2);
    end
    
    % if there is some dot product that is non-zeros - means that the
    % vectors are not coplanar, and convex hull can be built.
    if find(abs(dp)>1e-2)
        % There are some vectors out of the plane of first 2 vectors
        % Find convex hull indices:
        TR = delaunayTriangulation([W';[0 0 0]]);
        [~,P] = freeBoundary(TR);
        [is,index] = ismember([0 0 0], P, 'rows');
        if is
            P(index,:) = [];
        end
        W_CH = P';
    else
        % Select 2 marginal
        min_cp = min(cp);
        max_cp = max(cp);
        
        if min_cp == 0
            W_min = W(:,1);
        else
            W_min = W(:,find(cp==min_cp,1));
        end
        
        if max_cp == 0
            W_max = W(:,2);
        else
            W_max = W(:,find(cp==max_cp,1));
        end
        
        % The convex hull is only 2 marginal vectors
        W_CH = [W_min,W_max];
    end
    W_CH = reorder_W(W_CH);
    if size(W_CH, 2) > 2  % more than 2 vectors - searching for 1 finger to complete
        
        % Iterating over edges
        for e_i = 1:p.N_e
            n = p.Inner_normals(e_i,:);
            e1 = (1-Move_from_vertex_ratio)*p.Edges(e_i,:) +...
                Move_from_vertex_ratio * p.Edges(e_i+1,:);
            e2 = Move_from_vertex_ratio*p.Edges(e_i,:) +...
                (1-Move_from_vertex_ratio) * p.Edges(e_i+1,:);
            % Suggested EGW
            w = [n, cross2d(e1-p.Center,n)/p.Area^0.5;
                n, cross2d(e2-p.Center,n)/p.Area^0.5]';
            % We wish to confine the EGW inside of the inverted cone.
            modified_EGW = slice_cone(-W_CH,w);
            
            % If there is such variant, we found that for some region,
            % fingers placed on the edge will complete the grasp
            if ~isempty(modified_EGW)
                % Next step is extracting the allowed finger position along the edge
                L = diff(w(3,1:2));
                p1_ = (modified_EGW(3,1) - w(3,1))/L;
                p2_ = (modified_EGW(3,2) - w(3,1))/L;
                % Remember to remap the position to be normalized
                % along the edge length and not only to the allowed
                % region lenth
                p1_ = Move_from_vertex_ratio+p1_*(1-2*Move_from_vertex_ratio);
                p2_ = Move_from_vertex_ratio+p2_*(1-2*Move_from_vertex_ratio);
                
                
                % We find optimal location by checking the maximum of
                % the minimal of CH face distances. Each face is
                % checked: normal direction created, and then distance
                % from the origin is checked.
                pos_range = linspace(p1_,p2_);
                torq_range = linspace(modified_EGW(3,1), modified_EGW(3,2));
                d_Full = zeros(size(W_CH,2),numel(pos_range));
                points_alonge_the_line = [repmat(n(:),1,100);torq_range];
                
                V1 = W_CH(:,end)-points_alonge_the_line;
                V2 = W_CH(:,1)-points_alonge_the_line;
                d_Full(size(W_CH,2),:) = dot(-points_alonge_the_line, cross(V1,V2)/norm(cross(V1,V2)));
                
                for f_i = 1:size(W_CH,2)-1
                    V1 = W_CH(:,f_i)-points_alonge_the_line;
                    V2 = W_CH(:,f_i+1)-points_alonge_the_line;
                    d_Full(f_i,:) = dot(-points_alonge_the_line, cross(V1,V2)/norm(cross(V1,V2)));
                end
                d = min(d_Full,[],1);
                optimal_pos = pos_range(find(d == max(d),1));
                contact_point = p.point_from_edgePosition(e_i,optimal_pos);
                f =  ContactVector(contact_point,...
                    p.find_normal_at_point(contact_point),finger_d,p_i);
                
                if DEBUG
                    figure(16)
                    f.plot_contact('b')
                    % DT = delaunayTriangulation([1e-3 1e-3 0; -1e-3 2.5e-3 1e-4 ; w]);
                    % tetramesh(DT,'FaceAlpha',0.1,'FaceColor','y');
                    % quiver3(0*w(:,1),0*w(:,1),0*w(:,1),w(:,1),w(:,2),w(:,3),'AutoScale','off')
%                     drawnow
%                     pause(0.2)
                end
                
                Fingers = [Fingers;table(p_i,...
                    max([1 1+max(Fingers.ContactGroup(Fingers.PolygonNum==p_i))]),...
                    e_i, [p1_ p2_], optimal_pos, 0, f,...
                    'VariableNames', TableVarNames)] %#ok<NOPTS>
                % This adds one finger that can immobilize the body
            end
        end
        clearvars n e1 e2 w modified_EGW
    else
        % 2 vectors (or 1?. But I believe that I always will have at
        % least 2, since I'm not doing vertex2edge contact.)
        norm_WCH = unit_vector(cross(W_min,W_max));
        
        while 0
            search for edge that can span both sides of the plane
            If no such edge found, search for plane couples.
        end
        % Try searching 1 edge
        for e_i = 1:p.N_e
            n = p.Inner_normals(e_i,:);
            e1 = (1-Move_from_vertex_ratio)*p.Edges(e_i,:) +...
                Move_from_vertex_ratio * p.Edges(e_i+1,:);
            e2 = Move_from_vertex_ratio*p.Edges(e_i,:) +...
                (1-Move_from_vertex_ratio) * p.Edges(e_i+1,:);
            w = [n, cross2d(e1-p.Center,n)/p.Area^0.5;
                n, cross2d(e2-p.Center,n)/p.Area^0.5]';
            if DEBUG
                figure(101); clf
                quiver3(0*W(1,:),0*W(1,:),0*W(1,:),...
                    W(1,:), W(2,:), W(3,:),'k-','AutoScale','off')
                hold on;axis equal; grid on;
                quiver3(0*w(1,:),0*w(1,:),0*w(1,:),w(1,:),w(2,:),w(3,:),'b')
%                 drawnow
%                 pause(0.2)

            end
            proj_1 = dot(w(:,1),norm_WCH);
            proj_2 = dot(w(:,2),norm_WCH);
            CH_temp = [W_CH,w];
            try
                K = convhull([CH_temp'; [0 0 0]]);
            catch
                % The EGW does not form a convex hull with given
                % contacts. Continue to next edge.
                continue
            end
            if ~ismember(size(CH_temp,2)+1,K)
                % the origin is not inside
                % the edge can form a fixation with given contacts
                m = -proj_2/(proj_1-proj_2);
                mid = Move_from_vertex_ratio+m*(1-2*Move_from_vertex_ratio);
                f1 =  ContactVector(e1, n, finger_d, p_i);
                f2 =  ContactVector(e2, n, finger_d, p_i);
                contact_group_new_index = max([1 1+max(Fingers.ContactGroup(Fingers.PolygonNum==p_i))]);
                Fingers = [Fingers;table(p_i,contact_group_new_index,...
                    e_i, [Move_from_vertex_ratio 1-mid-Move_from_vertex_ratio],...
                    Move_from_vertex_ratio, 0, f1,...
                    'VariableNames', TableVarNames)];
                Fingers = [Fingers;table(p_i,contact_group_new_index,...
                    e_i, [1-mid+Move_from_vertex_ratio 1-Move_from_vertex_ratio],...
                    1-Move_from_vertex_ratio, 0, f2,...
                    'VariableNames', TableVarNames)];
                
            end
        end
    end
    % Try searching edge combinations
    if sum(Fingers.PolygonNum==p_i) == 0
        fprintf("No single edge found to grasp polygon #%d.\nTrying 2-contact combinations\n",p_i);
        
        % Iterate over all variants or over all edges?
        % Going for variants:
        for var_i = p.Variants2'
            n1 = p.Inner_normals(var_i(1),:);
            n2 = p.Inner_normals(var_i(2),:);
            polygon = [W_CH(1:2,:), n1(:), n2(:)];
            try
                K = convhull(polygon');
            catch
                % The edge normals do not form a convex hull with given
                % contacts. Continue to next variant.
                continue
            end
            %                     K = convhull(polygon');
            polygon = polygon(:,K);
            [in,on] = inpolygon(Qx(:), Qy(:), polygon(1,:), polygon(2,:));
            
            if in
                % origin is inside
                if 1 % set to use Mode1
                    
                    d = -ones(10,10);
                    pos_range = linspace(Move_from_vertex_ratio,1-Move_from_vertex_ratio,10);
                    
                    for pos1_ind = 1:10
                        for pos2_ind = 1:10
                            p1 = p.point_from_edgePosition(var_i(1),pos_range(pos1_ind));
                            p2 = p.point_from_edgePosition(var_i(2),pos_range(pos2_ind));
                            W_combined = [W_CH,...
                                [n1(:);cross2d(p1-p.Center,n1)./sqrt(p.Area)] ,...
                                [n2(:);cross2d(p2-p.Center,n2)./sqrt(p.Area)]];
                            d(pos1_ind,pos2_ind) = GQM_from_W(W_combined');
                        end
                    end
                    max_d = max(d(:));
                    [i,j] = find(d==max_d,1);
                    pos1 = pos_range(i);
                    pos2 = pos_range(j);
                    contact_group_new_index = max([1 1+max(Fingers.ContactGroup(Fingers.PolygonNum==p_i))]);
                    f1 = ContactVector(p.point_from_edgePosition(var_i(1),pos1), n1,finger_d/2, p_i);
                    f2 = ContactVector(p.point_from_edgePosition(var_i(2),pos2), n2,finger_d/2, p_i);
                    Fingers = [Fingers;table(p_i,contact_group_new_index,...
                        var_i(1), [0 0], pos1, max_d, f1, 'VariableNames', TableVarNames)];
                    Fingers = [Fingers;table(p_i,contact_group_new_index,...
                        var_i(2), [0 0], pos2, max_d, f2, 'VariableNames', TableVarNames)];
                end
                
                
                if 0 % Set to use Mode2
                    n1 = p.Inner_normals( var_i(1),:); %#ok<UNRCH>
                    e11 = (1-Move_from_vertex_ratio)*p.Edges( var_i(1),:) +...
                        Move_from_vertex_ratio * p.Edges( var_i(1)+1,:);
                    e21 = Move_from_vertex_ratio*p.Edges( var_i(1),:) +...
                        (1-Move_from_vertex_ratio) * p.Edges( var_i(1)+1,:);
                    % Suggested EGW
                    EGW_1 = [n1, cross2d(e11-p.Center,n1)/p.Area^0.5;
                        n1, cross2d(e21-p.Center,n1)/p.Area^0.5]';
                    
                    n2 = p.Inner_normals( var_i(2),:);
                    e12 = (1-Move_from_vertex_ratio)*p.Edges( var_i(2),:) +...
                        Move_from_vertex_ratio * p.Edges( var_i(2)+1,:);
                    e22 = Move_from_vertex_ratio*p.Edges( var_i(2),:) +...
                        (1-Move_from_vertex_ratio) * p.Edges( var_i(2)+1,:);
                    % Suggested EGW
                    EGW_2 = [n2, cross2d(e12-p.Center,n2)/p.Area^0.5;
                        n2, cross2d(e22-p.Center,n2)/p.Area^0.5]';
                    
                    
                    if DEBUG
                        fprintf("Edges %d,%d can possibly complete the grasp:\n",...
                            var_i(1),var_i(2));
                        fprintf("\tNormal directions are: [%2.2f,%2.2f] and [%2.2f,%2.2f].\n",...
                            n1(1),n1(2),n2(1),n2(2));
                        figure(99);clf
                        quiver(0*W_CH(1,:), 0*W_CH(1,:), ...
                            W_CH(1,:), W_CH(2,:), 'k')
                        hold on; axis equal; grid on;
                        quiver([0 0], [0 0], [n1(1) n2(1)], [n1(2) n2(2)], 'b')
                        title(sprintf("Edges %d,%d",  var_i(1),var_i(2)));
                        
                        figure(101); clf
                        quiver3(0*W(1,:),0*W(1,:),0*W(1,:),...
                            W(1,:), W(2,:), W(3,:),'k-','AutoScale','off');
                        hold on; axis equal; grid on;
                        quiver3(0*W(1,:),0*W(1,:),0*W(1,:),...
                            -W(1,:), -W(2,:), -W(3,:),'Color',[.6 .6 .6],'AutoScale','off');
                        EGW_N = 10;
                        EGW_1_plot = [linspace(EGW_1(1,1),EGW_1(1,2),EGW_N);
                            linspace(EGW_1(2,1),EGW_1(2,2),EGW_N);
                            linspace(EGW_1(3,1),EGW_1(3,2),EGW_N)];
                        EGW_2_plot = [linspace(EGW_2(1,1),EGW_2(1,2),EGW_N);
                            linspace(EGW_2(2,1),EGW_2(2,2),EGW_N);
                            linspace(EGW_2(3,1),EGW_2(3,2),EGW_N)];
                        quiver3(0*EGW_1_plot(1,:),0*EGW_1_plot(1,:),0*EGW_1_plot(1,:),...
                            EGW_1_plot(1,:), EGW_1_plot(2,:), EGW_1_plot(3,:),...
                            'Color',[255 150 0]/255,'AutoScale','off')
                        quiver3(0*EGW_2_plot(1,:),0*EGW_2_plot(1,:),0*EGW_2_plot(1,:),...
                            EGW_2_plot(1,:), EGW_2_plot(2,:), EGW_2_plot(3,:),...
                            'Color',[255 69 0]/255,'AutoScale','off')
                        zlabel('\tau_z','FontSize',20)
                        xlabel('f_x','FontSize',20)
                        ylabel('f_y','FontSize',20)
                    end
                    
                end
            end
        end
        
        
    end
end

if sum(Fingers.PolygonNum==p_i) == 0
    % Meaning no fingers found to complete the grasp
    fprintf('No fingers for polygon %d\n',p_i)
    warning("YOSSI:NoEdgeOppositeToCH",...
        "No edge to complete the grasp.\nIterating over a polygon #%d in the PolyList to find an edge with normal direction to complete the ConvexHull formed by given contacts",p_i)
end
