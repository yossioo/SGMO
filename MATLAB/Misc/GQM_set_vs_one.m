clc; clear;


p1 = Polygon_mkII([-2 0 ; 0 2; -3 2; -3.5 0; -3 -2; 0 -2 ]);
p2 = Polygon_mkII([-2 0 ; 0 2;  3 2;  3.5 0;  3 -2; 0 -2 ]);

C1 = p1.find_contacts_for_positions(7/16);
C2 = p1.find_contacts_for_positions(9/16);
C3 = p2.find_contacts_for_positions(7/16);
C4 = p2.find_contacts_for_positions(9/16);

c1 = Vector(C1(1:2),C1(3:4),1);
c2 = Vector(C2(1:2),C2(3:4),1);
c3 = Vector(C3(1:2),C3(3:4),1);
c4 = Vector(C4(1:2),C4(3:4),1);

IC1 = p1.find_contacts_for_positions(p1.Norm_cum_len(1));
IC2 = p1.find_contacts_for_positions(p1.Norm_cum_len(2));
IC3 = p1.find_contacts_for_positions(p1.Norm_cum_len(end));
IC4 = p1.find_contacts_for_positions(p1.Norm_cum_len(end-1));
ic1 = Vector(IC1(1:2),p1.Inner_normals(1,:),1);%[-1 1]/sqrt(2)
ic2 = Vector(IC2(1:2),p1.Inner_normals(1,:),1);
ic3 = Vector(IC3(1:2),p1.Inner_normals(end,:),1);
ic4 = Vector(IC4(1:2),p1.Inner_normals(end,:),1);

figure(5); clf
p1.plot('y')
hold on; axis equal; grid on 
p2.plot('b')
c1.plot_contact();
c2.plot_contact();
c3.plot_contact();
c4.plot_contact();
ic1.plot();
ic2.plot();
ic3.plot();
ic4.plot();

n_B = 1+0*1/sqrt(p1.Area+p2.Area);
n_L = 1+0*1/sqrt(p2.Area);
n_R = 1+0*1/sqrt(p1.Area);

W1 = [c1.direction_vector(:);  cross2d(c1.point_on_the_line,c1.direction_vector)*n_B];
W2 = [c2.direction_vector(:);  cross2d(c2.point_on_the_line,c2.direction_vector)*n_B];
W3 = [c3.direction_vector(:);  cross2d(c3.point_on_the_line,c3.direction_vector)*n_B];
W4 = [c4.direction_vector(:);  cross2d(c4.point_on_the_line,c4.direction_vector)*n_B];
CH_full= [W1, W2, W3, W4]';
DT_CH_full = delaunayTriangulation(CH_full);


W1 = [c1.direction_vector(:);  cross2d(c1.point_on_the_line,c1.direction_vector)*n_L];
W2 = [c2.direction_vector(:);  cross2d(c2.point_on_the_line,c2.direction_vector)*n_L];
Wl3 = [ic1.direction_vector(:);  cross2d(ic1.point_on_the_line,ic1.direction_vector)*n_L];
Wl4 = [ic2.direction_vector(:);  cross2d(ic2.point_on_the_line,ic2.direction_vector)*n_L];
Wl5 = [ic3.direction_vector(:);  cross2d(ic3.point_on_the_line,ic3.direction_vector)*n_L];
Wl6 = [ic4.direction_vector(:);  cross2d(ic4.point_on_the_line,ic4.direction_vector)*n_L];
CH_left = [W1, W2, Wl3, Wl4, Wl5, Wl6]';
DT_CH_left = delaunayTriangulation(CH_left);


W3 = [c3.direction_vector(:);  cross2d(c3.point_on_the_line,c3.direction_vector)*n_R];
W4 = [c4.direction_vector(:);  cross2d(c4.point_on_the_line,c4.direction_vector)*n_R];
Wl3 = -[ic1.direction_vector(:);  cross2d(ic1.point_on_the_line,ic1.direction_vector)*n_R];
Wl4 = -[ic2.direction_vector(:);  cross2d(ic2.point_on_the_line,ic2.direction_vector)*n_R];
Wl5 = -[ic3.direction_vector(:);  cross2d(ic3.point_on_the_line,ic3.direction_vector)*n_R];
Wl6 = -[ic4.direction_vector(:);  cross2d(ic4.point_on_the_line,ic4.direction_vector)*n_R];
CH_right = [W3, W4, Wl3, Wl4, Wl5, Wl6]';
DT_CH_right = delaunayTriangulation(CH_right);



figure(6); clf
tetramesh(DT_CH_full,'FaceAlpha',0.9,'FaceColor','r');
hold on; grid on; axis equal;
zlabel('\tau_z','FontSize',20)
xlabel('f_x','FontSize',20)
ylabel('f_y','FontSize',20)


tetramesh(DT_CH_left,'FaceAlpha',0.3,'FaceColor','y');
tetramesh(DT_CH_right,'FaceAlpha',0.3,'FaceColor','c');

view(-40,10)
f = gcf;
f.Position(3:4) = [490 500];
saveas(gcf,'../LyX/images/GQM_set_vs_one.svg')
% legend({'Both objects','Left','Right','5','6','7','8','9'},'Location','northeast','Interpreter','latex')
