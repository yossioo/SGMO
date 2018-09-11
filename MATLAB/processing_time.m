clc;

x = [2 3 5 10 20];
num_fing = [4 6; 6 7; 6 9; 8 14; 13 24];
time = [1.64 1.93; 1.87 2.69; 2.246 4.47; 3.16 5.63; 5.71 8.69];

figure(90)
bar(num_fing)
grid
a = gca;
legend("Triangles",'Squares','Location','northwest')
xlabel("Number of objects")
ylabel("Number of fingers requiered")
a.XTickLabel = {'2', '3','5','10','20'};

figure(91)
bar(time)
grid
a = gca;
legend("Triangles",'Squares','Location','northwest')
xlabel("Number of objects")
ylabel("Computation time [s]")

a.XTickLabel = {'2', '3','5','10','20'};
