clc; clear all;
global DEBUG
DEBUG = 1
finger_d = 5;

PolyList{1} = Polygon_mkII([-1 -1; 1 -1; 1 1; -1 1]*10,"Square", 'y');
e = [1 2 3 4]';
pos = [.25;.75;.25;.25];
TableVarNames = {'PolygonNum','ContactGroup','EdgeNum','EdgeRange','OptimalPosition','Group_GQM','ContactVector'};
f1 = ContactVector(PolyList{1}.point_from_edgePosition(e(1),pos(1)),PolyList{1}.Inner_normals(e(1),:),finger_d,1);
f2 = ContactVector(PolyList{1}.point_from_edgePosition(e(2),pos(2)),PolyList{1}.Inner_normals(e(2),:),finger_d,1);
f3 = ContactVector(PolyList{1}.point_from_edgePosition(e(3),pos(3)),PolyList{1}.Inner_normals(e(3),:),finger_d,1);
f4 = ContactVector(PolyList{1}.point_from_edgePosition(e(4),pos(4)),PolyList{1}.Inner_normals(e(4),:),finger_d,1);
Fingers = table([1;1;1;1],[1;1;1;1],e,repmat([0 1],4,1),pos, [0;0;0;0], [f1;f2;f3;f4], 'VariableNames', TableVarNames)


%%

evaluate_script
%%

figure(19);
clf;

PolyList{1}.plot(); hold on; axis equal; grid on;
for i = 2:numel(PolyList)
    PolyList{i}.plot()
    %     text(PolyList{i}.Center(1)-5,PolyList{i}.Center(2),num2str(i))
end
t = linspace(0,2*pi);
x = finger_d/2*cos(t); y = finger_d/2*sin(t);
axis manual

finger_centers = zeros(numel(BestFingersGroup(:,6)),2);
for i = 1:numel(BestFingersGroup(:,7))
    f = BestFingersGroup{i,7};
    f.plot_contact('b')
    allowed_reg = [PolyList{BestFingersGroup.PolygonNum(i)}.point_from_edgePosition(...
        BestFingersGroup.EdgeNum(i), BestFingersGroup.EdgeRange(i,1));
        PolyList{BestFingersGroup.PolygonNum(i)}.point_from_edgePosition(...
        BestFingersGroup.EdgeNum(i), BestFingersGroup.EdgeRange(i,2))];
    plot(allowed_reg(:,1),allowed_reg(:,2),'c','LineWidth',2);
    finger_centers(i,:) = f.get_finger_center(finger_d);
    f_c = fill(finger_centers(i,1)+x,finger_centers(i,2)+y,'b','FaceAlpha',.2);
    % s = scatter(p(1),p(2),100,'b', 'filled','MarkerFaceAlpha',0.2);
    f.draw_inf_line('k')
end
