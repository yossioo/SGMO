function ordered_W = reorder_W(W)
%REORDER_W Rearranges given array of wrench space vectors
%   3D Vectors that form Wrench cone are rearranged to be in
%   (counter-)clockwise order around the axis of the cone.
[r,c] = size(W);
if r~=3
    W = W';
end
ax = mean(W,2);
hor = cross([0; 0; 1],ax(:));
ver = cross(ax(:),hor);
Wh = dot(repmat(hor,1,c),W);
Wv = dot(repmat(ver,1,c),W);
ang = atan2(Wv,Wh);
[~,i] = sort(ang);
ordered_W = W(:,i);
end

