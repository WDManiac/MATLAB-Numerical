%% Test code for qr_sol
clc
clear
close all

% Specify the size of linear system:
N = 50;

% Let's generate random coefficient matrix and bias vector
A = 10*rand(N);
b = 5*rand(N,1);

% Then with MATLAB built-in method, we can solve for solution:
x_ref = A\b;

% Also we could solve with built-in QR function:
[Q_ref,R_ref] = qr(A);

% Then what our function produces are:
[Q,R,x] = qr_sol(A,b);

% Compare Q~Q_ref, R~R_ref and x~x_ref, accuracy of qr_sol could be testified

% Plot: solution comparison
figure
plot(x_ref,x,'r+','linewidth',1.5)
range = [min(x),max(x)];
hold on, plot(range,range,'k--','linewidth',.5)

% add a legend to indicate which curve is what
% Note: 'location' allows you to put the legend at specific corner of the figure
legend('Comparison','1:1','location','northwest')

grid on % enable grid on the figure
set(gca,'fontsize',14) % set font-size to be 14 for all labels

% here, "'interpreter','latex'" would show the axis label in LaTeX format
% this is recommended for mathematical expression
% Note: to enable LaTaX format, MUST have the 2 $ in the label string
xlabel('$x_{ref}$','interpreter','latex')
ylabel('$x_{calc}$','interpreter','latex') 
title 'Validation of Linear System Solving'