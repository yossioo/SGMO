clc; clear;

syms nx ny x1 x2 y1 y2 s wix wiy wiz wjx wjy wjz

wk = [nx; ny; (x1*ny - y1*nx) + s*(x2*ny + y2*nx -x1*ny -y1*nx)];
wi = [wix; wiy; wiz];
wj = [wjx; wjy; wjz];
d = -dot(wk, cross(wi-wk, wj-wk));

% latex(simplify(d))

ex = cross2d( (1-s)*[x1;y1]+s*[x2;y2], [nx;ny])
latex(simplify(cross2d( (1-s)*[x1;y1]+s*[x2;y2], [nx;ny])))
