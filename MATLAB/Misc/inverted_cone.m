clc; clear;


% p = Polygon_mkII([-100 -50 ; 100 -50; 0 100]);
% First case is the triangle, we can see what region is suitable
p = Polygon_mkII([-2 -2 ; -2 2; 2 2; 2 -2]);
figure(1)
clf
p.plot()
hold on; grid on;
axis equal

C1 = p.find_contacts_for_positions(3/16);
C2 = p.find_contacts_for_positions(5/16);
C3 = p.find_contacts_for_positions(13/16);

c1 = Vector(C1(1:2),C1(3:4),1);
c2 = Vector(C2(1:2),C2(3:4),1);
c3 = Vector(C3(1:2),C3(3:4),1);

c1.plot_contact();
c2.plot_contact();
c3.plot_contact();

W1 = [c1.direction_vector(:);  cross2d(c1.point_on_the_line,c1.direction_vector)];
W2 = [c2.direction_vector(:);  cross2d(c2.point_on_the_line,c2.direction_vector)]; %/sqrt(p.Area)
W3 = [c3.direction_vector(:);  cross2d(c3.point_on_the_line,c3.direction_vector)];

WC = [W1, W2, W3];
M = -mean(WC,2)*1;

V = [W1, W2, W3]';
DT = delaunayTriangulation([0 0 0; V]);
iDT = delaunayTriangulation(-[0 0 0; V]);

%%
f = figure(2); clf;
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


view(-35,25)
legend({'$w_a$','$w_b$','$w_c$',...
    'Convex Cone','Inverted C.Cone','Cone Axis'},'Location','northeast','Interpreter','latex')
saveas(gcf,'../LyX/images/example_inv_cone_01.svg')

%%
figure(3); clf;
cr = ceil(sqrt(p.N_e)); %cols-rows
% edge_names = 'Edge ' + string((1:p.N_e)');
edge_names = {'West edge','North edge','East edge','South edge'};
for i = 1:p.N_e
    subplot(cr,cr,i)
    tetramesh(iDT,'FaceAlpha',0.1,'FaceColor','b');
    hold on; grid on; axis equal;
%     a = gca;
%     a.Position = [0.1+(1-mod(i,2))*0.5, 0.6-0.475*(i>2), .35, .35];
    view(-35,15)
    zlabel('\tau_z','FontSize',20)
    xlabel('f_x','FontSize',20)
    ylabel('f_y','FontSize',20)
    title(edge_names{i})
    n = p.Inner_normals(i,:);
    e1 = p.Edges(i,:);
    e2 = p.Edges(i+1,:);
    w = [n, cross2d(e1,n);
        n, cross2d(e2,n)];
    DT_GEW = delaunayTriangulation([1e-3 1e-3 0; -1e-3 2.5e-3 1e-4 ; w]);
    tetramesh(DT_GEW,'FaceAlpha',0.1,'FaceColor','y');
    quiver3(0*w(:,1),0*w(:,1),0*w(:,1),w(:,1),w(:,2),w(:,3),'AutoScale','off')

end
saveas(gcf,'../LyX/images/example_inv_cone_02.svg')

%%
figure(4); clf;
i = 3;
nx = p.Inner_normals(i,1);
ny = p.Inner_normals(i,2);
x1 = p.Edges(i,1);
y1 = p.Edges(i,2);
x2 = p.Edges(i+1,1);
y2 = p.Edges(i+1,2);

s = linspace(0,1,100);
wk = [nx;ny;0]+0*s + cross([(1-s)*x1+s*x2;(1-s)*y1+s*y2;0*s] ,[nx;ny;0]+0*s);

% Make it work with any number of planes. But what if the cone is not
% convex? Make it convex first!
d1 = dot(wk, cross(W1-wk, W2-wk));
d2 = dot(wk, cross(W2-wk, W3-wk));
d3 = dot(wk, cross(W3-wk, W1-wk));
d = min([d1;d2;d3]);
mx = find(d==max(d));
pos = s(mx(1))*p.E_norm_lengths(i) + p.Norm_cum_len(i) - eps;

plot(s,d1,'-b',s,d2,'-g',s,d3,'-.r','LineWidth',1)
grid on;hold on
ylabel('Inscribed sphere radius');
xlabel('s')
a = gca;
x_pos = a.Position(1) + s(mx(1)) * a.Position(3);
y_pos = a.Position(2) + ...
    (d(mx(1))-a.YLim(1))/diff(a.YLim) * a.Position(4);

annotation('textarrow',x_pos*[0.8 1],y_pos*[0.9 1],...
    'String','Selected location',...
    'FontSize',20)
plot(s,d,'-k','LineWidth',1.5)
legend({'$d_1$','$d_2$','$d_3$','$\min{d}$'},'Location','northwest','Interpreter','latex')
saveas(gcf,'../LyX/images/example_inv_cone_03.svg')


figure(1)
C = p.find_contacts_for_positions(pos);
c = Vector(C(1:2),C(3:4),1);
c.plot_contact('b');
saveas(gcf,'../LyX/images/example_inv_cone_04.svg')
