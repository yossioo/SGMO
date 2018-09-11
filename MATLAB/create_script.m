warning('off', 'MATLAB:polyshape:repairedBySimplify')

Noise = 0;
%% N Random M-shapes
% 
N = 6 ; % number of objects
M = 3; % - number of vertices in polygons
t = 0: 2*pi/M :(2*pi-eps*1e5);
vert = [cos(t(:)), sin(t(:))];
colors = jet(N);
R = normrnd(20, Noise, numel(t),N);
for i=1:N
    name = "p" + char(65+floor(i/26)) +char(64+mod(i,26));
    P{i} = Polygon_mkII(R(:,i).*vert,name, colors(i,:)); %#ok<*SAGROW>
end

%% 4 shapes 3 triang and 1 convex - used in paper
% R = 30;
% colors = jet(4);
% P{1} = Polygon_mkII([0 0; R*[cosd(30) sind(30)]; R*1 0],"Tri_A", colors(1,:));
% P{2} = Polygon_mkII([0 0; R*[cosd(30) sind(30)]; R*1 0],"Tri_B", colors(2,:));
% P{3} = Polygon_mkII([0 0; R*[cosd(30) sind(30)]; R*1 0],"Tri_C", colors(3,:));
% P{4} = Polygon_mkII([0 0; 30 30; 25 -5; -25 -10; -30 30],"Concave",  colors(4,:));

%% Random objects A with concave angles
% P{1} = Polygon_mkII([0 0; 1 0 ;  .8 .1; .1 .8;0 1],"ConcaveA");
% P{2} = Polygon_mkII([0 0; 1 0 ; .3 .4; 0 1],"ConcaveB");
% P{2}.translate([1 1]);
% P{3} = Polygon_mkII([0 0; 1 0 ; 1 1; 0 1],"Convex");
% P{3}.translate([2 0])
% P{3}.rotate(90,[0 0 ]);

%% 3 Hex, - angles > 90
% t = 0: pi/3 :2*pi;
% vert = [cos(t(:)), sin(t(:))];
% P{1} = Polygon_mkII(vert,"Hex_A");
% P{2} = Polygon_mkII(vert,"Hex_B");
% P{3} = Polygon_mkII(vert,"Hex_3");

%% More shapes without concave vertices and less than 90 deg
% P{1} = Polygon_mkII([-1 -1; 1 -1; -.5 1],"Triangle A");
% P{2} = Polygon_mkII([-1 -1; 1 -1; 0.5 .5],"Poly B");

%% 6 equilateral triangles 
% t = 0: pi/1.5 :2*pi;
% vert = [cos(t(:)), sin(t(:))];
% P{1} = Polygon_mkII(vert,"Tri_A");
% P{2} = Polygon_mkII(vert,"Tri_B");
% P{3} = Polygon_mkII(vert,"Tri_C");
% P{4} = Polygon_mkII(vert,"Tri_D");
% P{5} = Polygon_mkII(vert,"Tri_E");
% P{6} = Polygon_mkII(vert,"Tri_F");

%% 1 Octagon Object with dents 
% M = 8; N = 1;
% t = 0: 2*pi/M :(2*pi-eps*1e5);
% vert = [cos(t(:)), sin(t(:))];
% colors = jet(N);
% R = normrnd(1, Noise, numel(t),N);
% for i=1:N
%     name = "Tr-" + char(65+floor(i/26)) +char(64+mod(i,26));
%     P{i} = Polygon_mkII(R(:,i).*vert,name, colors(i,:));
% end
% P{1}.Shape.Vertices(6,:) = [-0.4 0.4];
% P{1}.Shape.Vertices(7,:) = [0 0.4];
% P{1}.update_from_shape();
% % P{1} = Polygon_mkII([0 0; 1 0; 0.5 0.5*sqrt(3)],"Obj");
% Filtered_Contacts{1} = ContactVector(P{1}.point_from_edgePosition(3,0.7),...
%     P{1}.find_normal_at_point(P{1}.point_from_edgePosition(3,0.7)), 1, 1);
% Filtered_Contacts{2} = ContactVector(P{1}.point_from_edgePosition(3,0.3),...
%     P{1}.find_normal_at_point(P{1}.point_from_edgePosition(3,0.3)), 1, 1);
% Filtered_Contacts{3} = ContactVector(P{1}.point_from_edgePosition(8,0.5),...
%     P{1}.find_normal_at_point(P{1}.point_from_edgePosition(8,0.5)), 1, 1);
% Filtered_Contacts_poly_ind=[1 1 1];

% %% Triangle with 3 contacts
% 
% M = 3; N = 1;
% t = 0: 2*pi/M :(2*pi-eps*1e5);
% vert = [cos(t(:)), sin(t(:))];
% colors = jet(N);
% R = normrnd(1, Noise, numel(t),N);
% for i=1:N
%     name = "Tr-" + char(65+floor(i/26)) +char(64+mod(i,26));
%     P{i} = Polygon_mkII(R(:,i).*vert,name, colors(i,:));
% end
% 
% Filtered_Contacts{1} = ContactVector(P{1}.point_from_edgePosition(3,0.7),...
%     P{1}.find_normal_at_point(P{1}.point_from_edgePosition(3,0.7)), 1, 1);
% Filtered_Contacts{2} = ContactVector(P{1}.point_from_edgePosition(3,0.3),...
%     P{1}.find_normal_at_point(P{1}.point_from_edgePosition(3,0.3)), 1, 1);
% Filtered_Contacts{3} = ContactVector(P{1}.point_from_edgePosition(2,0.5),...
%     P{1}.find_normal_at_point(P{1}.point_from_edgePosition(2,0.5)), 1, 1);
% Filtered_Contacts_poly_ind=[1 1 1];

%% Triangle with 2 contacts
% 
% M = 3; N = 1;
% t = 0: 2*pi/M :(2*pi-eps*1e5);
% vert = 30*[cos(t(:)), sin(t(:))];
% colors = jet(N);
% R = normrnd(1, Noise, numel(t),N);
% for i=1:N
%     name = "Tr-" + char(65+floor(i/26)) +char(64+mod(i,26));
%     P{i} = Polygon_mkII(R(:,i).*vert,name, colors(i,:));
% end
% rot = [cosd(10), -sind(10)
%     sind(10), cosd(10)];
% Filtered_Contacts{1} = ContactVector(P{1}.point_from_edgePosition(3,0.2),...
%     rot*P{1}.find_normal_at_point(P{1}.point_from_edgePosition(3,0.1))', 10, 1);
% Filtered_Contacts{2} = ContactVector(P{1}.point_from_edgePosition(2,1),...
%     rot*P{1}.find_normal_at_point(P{1}.point_from_edgePosition(2,0.9))',10, 1);
% Filtered_Contacts_poly_ind=[1 1];
% 
% 


