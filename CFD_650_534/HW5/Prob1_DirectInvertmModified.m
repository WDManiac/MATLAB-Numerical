clc, clear
%% CFD HW 5: Problem 1 Solution Code
% Author: Ran Li
% In this problem, RHS : R = -8*pi^2*cos(2*pi*x)*sin(2*pi*y);

% Set up number of node along 2 directions
M = 201; % # of row (y)
N = 201; % # of col (x)

% Calculate Spatial Step Length
dx = 1/(N-1);
dy = 1/(M-1);

% Calculate R matrix with given discretization
R = zeros(M,N); Pana = R;
for i = 1:M
    for j = 1:N
        x = (j-1)*dx;
        y = (M-i)*dy;
        R(i,j) = -8*pi^2*cos(2*pi*x)*cos(2*pi*y);
        Pana(i,j) = cos(2*pi*x)*cos(2*pi*y);
    end
end

% Intialize 1st Order Diff. Matrix
D1x = zeros(M);
D1y = zeros(N); 

% Calculate 1st Order Diff. Matrix
D1x(2:M-1,2:M-1) = full(gallery('tridiag',M-2,-1,0,1));
D1x(2,1) = -1; D1x(M-1,M) = 1;
D1x(1,1:3) = [-3 4 -1]; D1x(M,M-2:M) = [1 -4 3];

D1y(2:N-1,2:N-1) = full(gallery('tridiag',N-2,-1,0,1));
D1y(2,1) = -1; D1y(N-1,N) = 1;
D1y(1,1:3) = [-3 4 -1]; D1y(N,N-2:N) = [1 -4 3];

% Derive 2nd Order Diff. Matrix
D1x = D1x/(2*dy);
D1y = D1y/(2*dx);
D2x = D1x*D1x;
D2y = D1y*D1y;

% Modify 2nd Order Diff. Matrix and keep only interior portion
A = D2x(2:M-1,2:M-1);
A(1,1) = D2x(2,2) + (4/3)*D2x(2,1); A(1,2) = D2x(2,3) - (1/3)*D2x(2,1);
A(2,1) = D2x(3,2) + (4/3)*D2x(3,1); A(2,2) = D2x(3,3) - (1/3)*D2x(3,1);
A(M-3,M-3) = D2x(M-2,M-2) - (1/3)*D2x(M-2,M); A(M-3,M-2) = D2x(M-2,M-1) + (4/3)*D2x(M-2,M); 
A(M-2,M-3) = D2x(M-1,M-2) - (1/3)*D2x(M-1,M); A(M-2,M-2) = D2x(M-1,M-1) + (4/3)*D2x(M-1,M); 

B = D2y(2:N-1,2:N-1);
B(1,1) = D2y(2,2) + (4/3)*D2y(2,1); B(1,2) = D2y(2,3) - (1/3)*D2y(2,1);
B(2,1) = D2y(3,2) + (4/3)*D2y(3,1); B(2,2) = D2y(3,3) - (1/3)*D2y(3,1);
B(N-3,N-3) = D2y(N-2,N-2) - (1/3)*D2y(N-2,N); B(N-3,N-2) = D2y(N-2,N-1) + (4/3)*D2y(N-2,N); 
B(N-2,N-3) = D2y(N-1,N-2) - (1/3)*D2y(N-1,N); B(N-2,N-2) = D2y(N-1,N-1) + (4/3)*D2y(N-1,N); 

% Transpose D2y Matrix
Bt = B';

% Characteristical Decomposition of Diff. Matrices
[P,Dx] = eig(A);
[Q,Dy] = eig(Bt);
% Examine whether there's any complex eigenvalues
% [x1,y1] = find(imag(Dx)~=0)
% [x2,y2] = find(imag(Dy)~=0)
Lx = diag(Dx);
Ly = diag(Dy);
PI = inv(P);
QI = inv(Q);

Rhat = PI*R(2:M-1,2:N-1)*Q;

% Calculate Pseudo Pressure Phat
for i = 1:length(Lx)
    for j = 1:length(Ly)
        L(i,j) = (Lx(i) + Ly(j));
        Phat(i,j) = Rhat(i,j)/L(i,j);
    end
end
Pressure(2:M-1,2:N-1) = P*Phat*QI;

% Compute Boundary Nodes' Values
Pressure(M,:) = (4/3)*Pressure(M-1,:) - (1/3)*Pressure(M-2,:);
Pressure(:,1) = (4/3)*Pressure(:,2) - (1/3)*Pressure(:,3);
Pressure(:,N) = (4/3)*Pressure(:,N-1) - (1/3)*Pressure(:,N-2);
Pressure(1,:) = (4/3)*Pressure(2,:) - (1/3)*Pressure(2,:);

% Compare Numerical and Analytical Solutions
figure, subplot(1,2,1),
contourf(0:dx:1,0:dy:1,Pressure),title('Numerical'),set(gca,'fontsize',18),colorbar,axis equal
subplot(1,2,2),
contourf(0:dx:1,0:dy:1,Pana),title('Analytical'),set(gca,'fontsize',18),colorbar,axis equal

% Compaire 2 Solutions along x = 0.5 
figure,subplot(1,2,1),
plot(0:dy:1,Pressure(:,(N-1)/2),'b*-')
hold on, plot(0:dy:1,Pana(:,(N-1)/2),'ro-')
title('Comparison in Y-direction (x = 0.5)'),legend('Num','Ana')
set(gca,'fontsize',18)
subplot(1,2,2),
plot(0:dx:1,Pressure((M-1)/2,:),'b*-')
hold on, plot(0:dx:1,Pana((M-1)/2,:),'ro-')
title('Comparison in X-direction (y = 0.5)'),legend('Num','Ana')
set(gca,'fontsize',18)