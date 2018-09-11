function varargout = plot_arrows3(array_3byM,color)
%PLOT_ARROWS3 Draws quiver3 for given array
%   Color specification is optional
[~,c] = size(array_3byM);
if nargin > 1
    q = quiver3(zeros(1,c),zeros(1,c),zeros(1,c),...
        array_3byM(1,:), array_3byM(2,:), array_3byM(3,:), 'Color',color);
else
    q = quiver3(zeros(1,c),zeros(1,c),zeros(1,c),...
        array_3byM(1,:), array_3byM(2,:), array_3byM(3,:));
end
if nargout > 0
    varargout{1} = q;
end
end

