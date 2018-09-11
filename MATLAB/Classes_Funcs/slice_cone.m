function cone2d_new = slice_cone(cone3d,cone2d)
%SLICE_CONE Gives the intersection of a 3D cone with 2D cone
%   Both variables dimensions: 3-by-mpts
cone3d = reorder_W(cone3d);
global DEBUG
cone2d_new = [];
N = size(cone3d,2); % Number of vertices of the 3D cone

n_cone2d = cross(cone2d(:,1), cone2d(:,2) );

% Inward normals of cone3d
n_cone3d = zeros(3,N);
n_cone3d(:,1:N-1) = cross(cone3d(:,1:N-1), cone3d(:,2:N) );
n_cone3d(:,N) = cross(cone3d(:,N), cone3d(:,1) );
n_cone3d = n_cone3d./vecnorm(n_cone3d);
Vert = [-cone3d cone2d];
try
    XY_CH_ind = convhull(Vert(1,:)',Vert(2,:)');
    XY_CH = polyshape(Vert(1:2,XY_CH_ind)');
    XY_CH_II = Polygon_mkII(Vert(1:2,XY_CH_ind)');
catch
    return
end
dx = 1e-2;
Qx = dx*cosd([0;120;240]);
Qy = dx*sind([0;120;240]);
%% DEBUG
if DEBUG
    figure(17); clf
    quiver3(0*cone3d(1,:),0*cone3d(1,:),0*cone3d(1,:), ...
        cone3d(1,:), cone3d(2,:), cone3d(3,:),'k-','AutoScale','off')
    hold on;
    %     quiver3(0*cone3d(1,:),0*cone3d(1,:),0*cone3d(1,:), ...
    %         -cone3d(1,:), -cone3d(2,:), -cone3d(3,:),'-','Color',[.5 .5 .5], 'AutoScale','off')
    % XY_CH.plot()
    quiver3(0*cone2d(1,:),0*cone2d(1,:),0*cone2d(1,:), ...
        cone2d(1,:), cone2d(2,:), cone2d(3,:),'r-','AutoScale','off')
    quiver3(0*n_cone3d(1,:),0*n_cone3d(1,:),0*n_cone3d(1,:), ...
        n_cone3d(1,:), n_cone3d(2,:), n_cone3d(3,:),'b-','AutoScale','off')
    % quiver3(0*n_cone3d(1,:),0*n_cone3d(1,:),0*n_cone3d(1,:), ...
    %     -n_cone3d(1,:), -n_cone3d(2,:), -n_cone3d(3,:),'b-','AutoScale','off')
    
    TR = delaunayTriangulation([cone3d'; 0 0 0]);
    if size(TR) > 0
        tetramesh(TR,'FaceAlpha',.2)
    end
    axis equal
    grid on
    zlabel('\tau_z','FontSize',20)
    xlabel('f_x','FontSize',20)
    ylabel('f_y','FontSize',20)
end
% DEBUG END
%%
[in,~] = inpolygon(Qx(:), Qy(:), Vert(1,XY_CH_ind), Vert(2,XY_CH_ind));
if in
    % The normal indeed spans the origin along with other wrenches
    
    w_top = cone2d(:,1);
    w_bot = cone2d(:,2);
    new_top = w_top;
    new_bot = w_bot;
    % Projection of 2D cone vectors on different face normals
    proj_top = dot(repmat(w_top,1,N), n_cone3d);
    proj_bot = dot(repmat(w_bot,1,N), n_cone3d);
    
    
    ind_top_min = find(proj_top == min(proj_top));
    ind_bot_min = find(proj_bot == min(proj_bot));
    if  ind_bot_min == ind_top_min 
        % Commented this to prevent  skipping when EGW is above a edge
        % of CH
        % The EGW is outside the inverse cone
        if (proj_top < 1e-4) & (proj_bot < 1e-4 )
            return
        end
    end
    
    if proj_top(ind_top_min(1))<0
        new_top = cross(n_cone3d(:,ind_top_min(1)),n_cone2d);
        new_top = new_top./norm(new_top(1:2));
    end
    
    if proj_bot(ind_bot_min(1))<0
        new_bot = cross(n_cone3d(:,ind_bot_min(1)),-n_cone2d);
        new_bot = new_bot./norm(new_bot(1:2));
    end
    cone2d_new=[new_top,new_bot];
    if DEBUG
        quiver3(0,0,0,new_top(1), new_top(2),new_top(3),'Color','m')
        quiver3(0,0,0,new_bot(1), new_bot(2),new_bot(3),'Color','m')
    end
else
    return
end


end

