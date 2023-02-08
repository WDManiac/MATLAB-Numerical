clc
close all
clear

x = logspace(0.0001,0.1,101);
y = x.*exp(-x);
xq = 1.0005:.006:1.2585;
yq = spline(x,y,xq);

figure, subplot(1,2,1)
% This is to create 1 row, 2 columns of subplots and draw within the 1st one
plot(x,y)
hold on, plot(xq,yq,'r+')

% Test of the composed function spline3:
[G,d,M,s]=spline3(x,y,xq);

subplot(1,2,2), plot(x,y)
hold on, plot(xq,s,'ro')