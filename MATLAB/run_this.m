%% Demo 1
clc; clear all;
%#ok<*CLALL>
%#ok<*UNRCH>
tic
% warning on verbose
warning off YOSSI:NoEdgeOppositeToCH
addpath('Classes_Funcs')

global DEBUG
DEBUG = true;
finger_d = 10;
%% Create objects
% The script creates a cell array P, which contains the polygonal
% objects.

create_script

%% Rearrange the objects

rearrange_script

%% Find inter-object contacts

find_i_contacts_script
%% Find fingers

find_fingers_script

%% Evaluate grasp and Move fingers (if needed)

evaluate_script

%% Display results
%
figure(19);
clf;

PolyList{1}.plot(PolyList{1}.Polygon_color,0); hold on; axis equal; grid on;
for i = 2:numel(PolyList)
    PolyList{i}.plot(PolyList{i}.Polygon_color,0)
    %     text(PolyList{i}.Center(1)-5,PolyList{i}.Center(2),num2str(i))
end
axis manual
finger_d = 5;
t = linspace(0,2*pi);
x = finger_d/2*cos(t); y = finger_d/2*sin(t);
for p_i = 1:numel(PolyList) %#ok<*UNRCH>
        Contacts = Filtered_Contacts(Filtered_Contacts_poly_ind == p_i);
        for c_i = 1:numel(Contacts)

            Contacts{c_i}.plot('r')
            %             Contacts{f}.plot_contact()
        end
end

Fingers.ContactVector(4).point_on_the_line = PolyList{Fingers.PolygonNum(4)}.point_from_edgePosition(Fingers.EdgeNum(4),0.9);
if 1 %Show all results
    finger_centers = zeros(numel(Fingers(:,6)),2);
    for i = 1:numel(Fingers(:,7))
        f = Fingers{i,7};
%         f.plot_contact('b')
        allowed_reg = [PolyList{Fingers.PolygonNum(i)}.point_from_edgePosition(...
            Fingers.EdgeNum(i), Fingers.EdgeRange(i,1));
            PolyList{Fingers.PolygonNum(i)}.point_from_edgePosition(...
            Fingers.EdgeNum(i), Fingers.EdgeRange(i,2))];
        plot(allowed_reg(:,1),allowed_reg(:,2),'c','LineWidth',2);
        finger_centers(i,:) = f.get_finger_center(finger_d);
        f_c = fill(finger_centers(i,1)+x,finger_centers(i,2)+y,'b','FaceAlpha',.8);
        
        f.plot('b')
        % s = scatter(p(1),p(2),100,'b', 'filled','MarkerFaceAlpha',0.2);
%         f.draw_inf_line('k')
    end
else % show only selected results
    finger_centers = zeros(numel(BestFingersGroup(:,6)),2);
    for i = 1:numel(BestFingersGroup(:,7))
        f = BestFingersGroup{i,7};
        f.plot('b')
        allowed_reg = [PolyList{BestFingersGroup.PolygonNum(i)}.point_from_edgePosition(...
            BestFingersGroup.EdgeNum(i), BestFingersGroup.EdgeRange(i,1));
            PolyList{BestFingersGroup.PolygonNum(i)}.point_from_edgePosition(...
            BestFingersGroup.EdgeNum(i), BestFingersGroup.EdgeRange(i,2))];
        plot(allowed_reg(:,1),allowed_reg(:,2),'c','LineWidth',2);
        finger_centers(i,:) = f.get_finger_center(finger_d);
        f_c = fill(finger_centers(i,1)+x,finger_centers(i,2)+y,'b','FaceAlpha',.8);
        % s = scatter(p(1),p(2),100,'b', 'filled','MarkerFaceAlpha',0.2);
%         f.draw_inf_line('k')
    end
end

% Can draw filtered contacts
if DEBUG
    
end
if 0 % Draw selected contact group
    for f = Fingers{Fingers.ContactGroup == 2,7}'
        finger_centers(i,:) = f.get_finger_center(finger_d);
        f_c = fill(finger_centers(i,1)+x,finger_centers(i,2)+y,'b','FaceAlpha',.8);
        
        f.plot('b')
    end
end
%%
toc
