function min_dist = GQM_from_W(W)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global DEBUG
if ~ismember(3,size(W))
    error("Convex hull dimensions are not Mx3")
end
if size(W,2) ~= 3
    W = W';
end

full = delaunayTriangulation(W);


id =  pointLocation(full,[0 0 0]);
min_dist = NaN;
if isnan(id)
    warning("YOSSI:OriginNotInCH",...
        "Given Convex Hull does not contain the origin.")
    return
end

[F,P] = freeBoundary(full);
for t_i = F'
    v1 = P(t_i(1),:);
    v2 = P(t_i(2),:);
    v3 = P(t_i(3),:);
    n = unit_vector(cross(v3-v1,v2-v1));
    min_dist = min([min_dist abs(dot(n,v1))]);
end

if DEBUG && 0 
    figure(98); clf
    quiver3(0*W(:,1),0*W(:,1),0*W(:,1),W(:,1),W(:,2),W(:,3),'Color','b','AutoScale','off')
    hold on; axis equal; grid on
    trisurf(F,P(:,1),P(:,2),P(:,3), ...
        'FaceColor','cyan','FaceAlpha',0.8);
end

end

