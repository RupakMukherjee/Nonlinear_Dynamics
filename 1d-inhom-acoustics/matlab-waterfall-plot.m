clc;clearvars;
L = 10*pi;
x = linspace(0,L,1023);
t = linspace(0,100,5);
[X,T]=meshgrid(x,t);
A = load('waterfalls.dat');
A(:,1) = NaN;
A(:,65:1023) = NaN;
h=waterfall(X,T,A);
set(h,'LineWidth',5);
colormap([0  0  0]);
xlabel('z','FontSize',20);
ylabel('Time','FontSize',20);
zlabel('b_y','FontSize',20);
xlim([0 2])
