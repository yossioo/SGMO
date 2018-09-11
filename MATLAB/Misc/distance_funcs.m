clc; clear;


p = Polygon_mkII([-100 -50 ; 100 -50; 0 100]);
figure(1)
clf
p.plot()
hold on
axis equal

C1 = p.find_contacts_for_positions(0.45);
C2 = p.find_contacts_for_positions(0.60);
C3 = p.find_contacts_for_positions(0.10);

c1 = Vector(C1(1:2)-p.Center,C1(3:4),100);
c2 = Vector(C2(1:2),C2(3:4),100);
c3 = Vector(C3(1:2),C3(3:4),100);

c1.plot_contact();
c2.plot_contact();
c3.plot_contact();

W1 = [c1.direction_vector(:);  cross2d(c1.point_on_the_line,c1.direction_vector)/sqrt(p.Area)];
W2 = [c2.direction_vector(:);  cross2d(c2.point_on_the_line,c2.direction_vector)/sqrt(p.Area)];
W3 = [c3.direction_vector(:);  cross2d(c3.point_on_the_line,c3.direction_vector)/sqrt(p.Area)];

WC = [W1, W2, W3];
M = -mean(WC,2)*1;
figure(2)
clf
quiver3(0*WC(1,:),0*WC(1,:),0*WC(1,:),WC(1,:),WC(2,:),WC(3,:),'AutoScale','off')
hold on
DT = delaunayTriangulation([WC, [0;0;0]]');
tetramesh(DT,'FaceAlpha',0.3);
DT = delaunayTriangulation([-WC, [0;0;0]]');
tetramesh(DT,'FaceAlpha',0.2,'FaceColor','b');

quiver3(0,0,0,M(1),M(2),M(3), 'Color',[0 .3 0],'AutoScale','off','LineWidth',2)
axis equal
grid on
zlabel('\tau_z','FontSize',20)
xlabel('f_x','FontSize',20)
ylabel('f_y','FontSize',20)

for i = 1:p.N_e
    n = p.Inner_normals(i,:);
    e1 = p.Edges(i,:);
    e2 = p.Edges(i+1,:);
    w = [n, cross2d(e1,n)/sqrt(p.Area);
        n, cross2d(e2,n)/sqrt(p.Area)];
    quiver3(0*w(:,1),0*w(:,1),0*w(:,1),w(:,1),w(:,2),w(:,3),'AutoScale','off')
end
