function angle = angle_of_2_vec(v1,v2)
%ANGLE2VEC Summary of this function goes here
%   Given 2 2D vectors returns angle in degrees wrapped to 180 (-180:180).
% It can be seen as counterclockwise steering angle to change from V1 to V2
if nargin < 2
    error("Two vectors needed")
end
V1 = v1./vecnorm(v1.').';
V2 = v2./vecnorm(v2.');

c = acosd(round(dot(V1,V2),5));
s = asind(round(cross2d(V1,V2),5));

if s == 0
    angle = c;
else
    angle = c *sign(s);
end
end

