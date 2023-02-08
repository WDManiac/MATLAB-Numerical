%% CFD HW Assignment 4 Problem 2
%  Ran Li
%  N-S Equation Solver, using psi-vorticity approach


clear, clc, close all
%% Part A : Generating mesh for the problem
X = 101; %    number of nodes along horrizontal direction
Y = 101; %    number of nodes along vertical direction
gamma = 1; %     ratio of length to height: gamma = length/height
h = 1; %     height of physical domain
[dx,dy,mat_ini] = blockmesh(X,Y,gamma,h); 
%   This give an initial matrix of ALL ZEROs
W = mat_ini; % RHS of the equation
Re = 1000;


%% Part B : Applying Boundary condition
%   Configuration:
%
%           b1
%    ----------------
%    |              |
%  b4|              |b2
%    |              |
%    ----------------
%           b3
%   Insert BCs described in problem and calculate exact values
x = linspace(0,1,X);
y = linspace(0,1,Y);
b1 = 0;%-dy;     
b2 = 0;     
b3 = 0;     
b4 = 0;
%   Insert Dirichlet BCs (For Psi)
        mat_ini(:,X) = b2*ones(Y,1);
        mat_ini(Y,:) = b3*ones(1,X);
        mat_ini(:,1) = b4*ones(Y,1);
        mat_ini(1,:) = b1*ones(1,X);

%% Part C : Solving
lim = 1000000; % Max possible iteration number
sf = 0.95; % Scale Factor for time iteration
[Psi,Ome,K,E1,E2] = falseTrans2DNS(mat_ini,W,Re,dx,dy,lim,sf);

% Part D : Visualization

[Xmat, Ymat] = meshgrid(x,flip(y));

figure('units','normalized','outerposition',[0,0,1,1]),
% title('Re = 1')
subplot(1,2,1)
contour(x,flip(y),Psi)
title('\psi')
grid on, 
set(gca,'FontSize',18), 
colorbar, axis equal

subplot(1,2,2)
contour(x,flip(y),Ome)
title('\omega')
grid on, 
set(gca,'FontSize',18), 
colorbar, axis equal

saveas(gcf,['Re',num2str(Re),'sol_',num2str(X),'x',num2str(Y)],'fig')
saveas(gcf,['Re',num2str(Re),'sol_',num2str(X),'x',num2str(Y)],'jpg')

load('Psi_w.mat');
figure, subplot(1,2,1),plot(Ome((X+1)/2,:),'bo--')
title('\omega(0.5,y^{*}) Ran')
grid on,set(gca,'fontsize',14)
subplot(1,2,2), plot(flip(w_Re100(:,(X+1)/2)),'r*-')
title('\omega(0.5,y^{*}) Reference')
grid on,set(gca,'fontsize',14)
saveas(gcf,['OmeVertRe',num2str(Re),'sol_',num2str(X),'x',num2str(Y)],'fig')
saveas(gcf,['OmeVertRe',num2str(Re),'sol_',num2str(X),'x',num2str(Y)],'jpg')

figure, subplot(1,2,1),plot(Ome(:,(X+1)/2),'bo--')
title('\omega(x^{*},0.5) Ran')
grid on,set(gca,'fontsize',14)
subplot(1,2,2), plot(flip(w_Re100((X+1)/2,:)),'r*-')
title('\omega(x^{*},0.5) Reference')
grid on,set(gca,'fontsize',14)
saveas(gcf,['OmeHoriRe',num2str(Re),'sol_',num2str(X),'x',num2str(Y)],'fig')
saveas(gcf,['OmeHoriRe',num2str(Re),'sol_',num2str(X),'x',num2str(Y)],'jpg')