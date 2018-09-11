clc; clear; %close all;

v1 = [1 0 -1];
v2 = [0 -1 1];
v3 = [0 1 1];

v1 = v1/norm(v1);
v2 = v2/norm(v2);
v3 = v3/norm(v3);

V = [v1; v2; v3];
DT = delaunayTriangulation([0 0 0; V]);
iDT = delaunayTriangulation(-[0 0 0; V]);

Ax = v1+v2+v3;

f = figure(1); clf;
qa = quiver3(0,0,0,V(1,1),V(1,2),V(1,3),...
    'AutoScale','off','LineWidth',2,'MaxHeadSize',2,...
    'DisplayName','$w_a$');% 'Color',[0 .3 0],
hold on; grid on;

quiver3(0,0,0,V(2,1),V(2,2),V(2,3),...
    'AutoScale','off','LineWidth',2,'MaxHeadSize',2);% 'Color',[0 .3 0],
quiver3(0,0,0,V(3,1),V(3,2),V(3,3),...
    'AutoScale','off','LineWidth',2,'MaxHeadSize',2);% 'Color',[0 .3 0],

% % In case you want all vectors to be draw together
% quiver3(0*V(:,1),0*V(:,1),0*V(:,1),...
%     V(:,1),V(:,2),V(:,3),...
%     'Color',[0 .3 0],'AutoScale','off','LineWidth',2)
zlabel('\tau_z','FontSize',20)
xlabel('f_x','FontSize',20)
ylabel('f_y','FontSize',20)
tetramesh(DT,'FaceAlpha',0.3,'FaceColor','g');
tetramesh(iDT,'FaceAlpha',0.1,'FaceColor','b');

plot3(Ax(1)*[-1 1], Ax(2)*[-1 1], Ax(3)*[-1 1], '--k','LineWidth',1.5)
view(-35,25)
legend({'$w_a$','$w_b$','$w_c$',...
    'Convex Cone','Inverted C.Cone','Cone Axis'},'Location','northeast','Interpreter','latex')
saveas(gcf,'../LyX/images/example_inv_cone_01.svg')

% N = 100;
% [A,E] = meshgrid(linspace(-180,180,N), linspace(-90,90,N));
d = 3;
[A,E] = meshgrid(-90:d:270, -90:d:90);
% N = size(A

p1 =  dot([reshape(cosd(A).*cosd(E),1,[]);
    reshape(sind(A).*cosd(E),1,[]);
    reshape(sind(E),1,[])] , repmat(cross(v1,v2)',1,numel(A)));
p1_fix = reshape(p1,size(A));

p2 =  dot([reshape(cosd(A).*cosd(E),1,[]);
    reshape(sind(A).*cosd(E),1,[]);
    reshape(sind(E),1,[])] , repmat(cross(v2,v3)',1,numel(A)));
p2_fix = reshape(p2,size(A));

p3 =  dot([reshape(cosd(A).*cosd(E),1,[]);
    reshape(sind(A).*cosd(E),1,[]);
    reshape(sind(E),1,[])] , repmat(cross(v3,v1)',1,numel(A)));
p3_fix = reshape(p3,size(A));

D_full = reshape(min([p3;p2;p1]),size(A));
% D_full = reshape(min([p3;p2;p1]),size(A));
% pos = find(d)
V = v1+v2+v3;
V = V/norm(V);
A_center = [atan2d(V(2),V(1)), atan2d(V(3),norm(V(1:2)))]

% figure(2)
% clf
% mesh(A,E,D_full,'EdgeAlpha',0, 'FaceColor','flat')
% % contour(A,E,D_full)
% hold on
% % surf(A,E,p1_fix,'EdgeAlpha',0.5)
% % surf(A,E,p2_fix,'EdgeAlpha',0.5)
% % surf(A,E,p3_fix,'EdgeAlpha',0.5)
% plot3(A_center(1)*[1 1], A_center(2)*[1 1],...
%     2.1*minmax(D_full(:)'),'r-','LineWidth',2)
% plot3(180+A_center(1)*[1 1], -A_center(2)*[1 1],...
%     2.1*minmax(D_full(:)'),'r-','LineWidth',2)
% % plot3(A_center(1)*[1 1]-180, -A_center(2)*[1 1]-180,...
% %     1.1*minmax(D_full(:)'),'r-','LineWidth',2)
% % plot3(A_center(1)*[1 1], A_center(2)*[1 1]-180,...
% %     1.1*minmax(D_full(:)'),'r-','LineWidth',2)
% xlabel('Azimuth from X [deg] ')
% ylabel('Elevation from XY [deg] ')
