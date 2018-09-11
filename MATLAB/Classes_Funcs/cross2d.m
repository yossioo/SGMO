function c = cross2d(a,b)
%CROSS2D Cross product of 2D vectors
% If these are matrices, one of the dimensions have to be 2
if ~ismember(2,size(a)) || ~ismember(2,size(b))
    error("Not 2D vectors")
end
if size(a,2) ~= 2
    a = a.';
end
if size(b,2) ~= 2
    b = b.';
end

c = a(:,1).*b(:,2)-a(:,2).*b(:,1);
end

