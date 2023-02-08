%% CFD HW Assignment 4 Problem 1
%  Ran Li
%  Poisson Equation Solver

clc, close all

%% Part A : Generating mesh for the problem
X = 201; %    number of nodes along horrizontal direction
Y = 201; %    number of nodes along vertical direction
gamma = 1; %     ratio of length to height: gamma = length/height
h = 1; %     height of physical domain
[dx,dy,mat_ini] = blockmesh(X,Y,gamma,h); 
%   This give an initial matrix of ALL ZEROs
W = mat_ini; % RHS of the equation
for i = 1:Y
    W(i,:) = (2+3*dy*(Y-i))*ones(1,X);
end

% Analytical Solution:
Phi = zeros(X,Y);
yf = flip(y);
for i = 1:X
    for j = 1:Y
        Phi(j,i) = x(i)^2 + 3*x(i)*yf(j) + 0.5*yf(j)^3;
    end
end


%% Part B : Applying Boundary condition----Dirichlet
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
b1 = x.^2 + 3*x +0.5;     
b2 = 1 + 3*y + 0.5*y.^3;     
b3 = x.^2;     
b4 = 0.5*y.^3;
%   Define intensity for each boundary condition, assume they are constant
        mat_ini(:,X) = flip(b2').*ones(Y,1);
        mat_ini(Y,:) = b3.*ones(1,X);
        mat_ini(:,1) = flip(b4').*ones(Y,1);
        mat_ini(1,:) = b1.*ones(1,X);

%   Apply Gauss Seidel iteration method for 2D problem

%% Part C : Solving
lim = 1000000; % Max possible iteration number
sf = 0.95; % Scale Factor for time iteration
[Mat_fin,K,e,Ek] = falseTrans2D(mat_ini,W,dx,dy,lim,sf,Phi);

%% Part D : Visualization

Em = abs(Phi - Mat_fin);
Error = max(max(Em));
[eidx, eidy] = find(abs(Phi - Mat_fin) == Error);

% Store error at each step
switch X
    case 201
        Estat(3,:) = [X,Error];
    case 101
        Estat(2,:) = [X,Error];
    case 51
        Estat(1,:) = [X,Error];
end

fy = flip(y);

% Part 1: Plot error field
figure,
contour(x,fy,Em)
hold on, plot(x(eidy),fy(eidx),'ro','LineWidth',2)
title(['Error with grid size: [',num2str(X),'\times',num2str(Y),']'])
grid on, 
set(gca,'FontSize',18), 
colorbar, axis equal
saveas(gcf,['Error_',num2str(X),'x',num2str(Y)],'fig')
saveas(gcf,['Error_',num2str(X),'x',num2str(Y)],'jpg')

% Part 2: Plot Simulation Result and Anlaytical Solution
figure('units','normalized','outerposition',[0,0,1,1]),
subplot(1,2,1)
contour(x,fy,Mat_fin)
hold on, plot(x(eidy),fy(eidx),'ro','LineWidth',2)
title(['Simulation [',num2str(X),'\times',num2str(Y),']'])
grid on, 
set(gca,'FontSize',18), 
colorbar, axis equal

subplot(1,2,2)
contour(x,flip(y),Phi)
title('Analytical')
grid on, 
set(gca,'FontSize',18), 
colorbar, axis equal

saveas(gcf,['Simulation_',num2str(X) ,'x',num2str(Y)],'fig')
saveas(gcf,['Simulation_',num2str(X),'x',num2str(Y)],'jpg')

% Part 3: Error against iteration number
figure,
plot(log(Ek),'r*-')
title(['Logrithm Error log(\epsilon) w.r.t. Iteration k --- grid: ',num2str(X),'x',num2str(Y)])
ylabel('log(\epsilon)'), xlabel('k')
grid on, 
set(gca,'FontSize',14) 

saveas(gcf,['ErrorK_',num2str(X),'x',num2str(Y)],'fig')
saveas(gcf,['ErrorK_',num2str(X),'x',num2str(Y)],'jpg')