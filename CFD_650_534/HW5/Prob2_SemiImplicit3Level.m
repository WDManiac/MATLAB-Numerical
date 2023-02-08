%% CFD HW5 Problem 2: 
% Author: Ran Li
% Semi-implicit Three-level Method for N-S Equation Solving
clc, clear, close all

% Initialization of Simulation:
% Mesh Size
M = 21;
N = 21;
% Reynolds Number involved
Re = 100;
% Due to this nondimension, driving velocity V and size of cavity 
% domain L could be both set as 1

% Generate Mesh
[dx,dy,mat_ini] = blockmesh(N,M,1,1);
U = mat_ini; % Horizontal Velocity Field
V = mat_ini; % Vertical Velocity Field
Phi = mat_ini; % Pressure Field (Pseudo)
Us = U; Vs = V;

% Setup Time Iteration
tmax = 5; % Max Time for Iteration
t = 0; % Initial Time
dt = 0.1; % Time Interval

% BCs for Velocity Field:
Vtop = 0.5;
U(1,:) = Vtop*ones(1,N);

% 1st Order Diff. Operator Matrix
D1x = zeros(N-2,N); D1y = zeros(M-2,M);
D1x(1,1) = -1; D1x(N-2,N) = 1;
D1x(:,2:N-1) = full(gallery('tridiag',N-2,-1,0,1));
D1y(1,1) = -1; D1y(M-2,M) = 1;
D1y(:,2:M-1) = full(gallery('tridiag',M-2,-1,0,1));
D1x = D1x/(2*dx); D1y = D1y/(2*dy);

% Discretization Coefficients
alpha = dt/(2*Re*dx^2);
beta = dt/(2*Re*dy^2);

% 2nd Order Diff. Operator Matrix (Linear Term)
D2x = zeros(N-2,N); D2y = zeros(M-2,M);
D2x(:,2:N-1) = full(gallery('tridiag',N-2,alpha,-2*alpha,alpha));
D2x(1,1) = alpha; D2x(N-2,N) = alpha;
D2y(:,2:M-1) = full(gallery('tridiag',M-2,beta,-2*beta,beta));
D2y(1,1) = beta; D2y(M-2,M) = beta;
D2x = 2*D2x/3; D2y = 2*D2y/3;

% 2nd Order Diff. Operator Matrix (LHS Semi-Implicit Solving)
Ax = zeros(N-2,N); Ay = zeros(M-2,M);
Ax(:,2:N-1) = full(gallery('tridiag',N-2,-alpha/3,1+2*alpha/3,-alpha/3));
Ax(1,1) = -alpha/3; Ax(N-2,N) = -alpha/3;
Ay(:,2:M-1) = full(gallery('tridiag',M-2,-beta/3,1+2*beta/3,-beta/3));
Ay(1,1) = -beta/3; Ay(M-2,M) = -beta/3;

% For first step, setup Nonlinear Term
UU = U.*U;
UV = U.*V;
VV = V.*V;
Nx0 = -UU(2:M-1,:)*D1x' - D1y*UV(:,2:N-1);
Ny0 = -UV(2:M-1,:)*D1x' - D1y*VV(:,2:N-1);

sequence = [];
Utop = [];

% Time Iteration Starts Here
while t <= tmax
    
    % Compute Non-Linear Terms
    UU = U.*U;
    UV = U.*V;
    VV = V.*V;
    
    % Compute Non-linear Convection Term
    Nx = - (UU(2:M-1,:)*D1x') - (D1y*UV(:,2:N-1));
    Ny = - (UV(2:M-1,:)*D1x') - (D1y*VV(:,2:N-1));
    
    % Compute Linear Viscous Term
    Lx = D2y*U(:,2:N-1) + U(2:M-1,:)*D2x';
    Ly = D2y*V(:,2:N-1) + V(2:M-1,:)*D2x';
    
    % Compose RHS R for Advancing Process
    Rx = (dt/3)*(3*Nx - Nx0) + Lx;
    Ry = (dt/3)*(3*Ny - Ny0) + Ly;
    Nx0 = Nx; Ny0 = Ny;
    
    % Advance In Time
    dT = dt/3;
    % Advancing Step 1
    for i = 1:M-2
        % Insert Boundary Condition
%         Rx(i,1) = Rx(i,1) - Ax(1,1)*U(i+1,1);
%         Rx(i,N-2) = Rx(i,N-2)- Ax(N-2,N)*U(i+1,N); 
%         Ry(i,1) = Ry(i,1) - Ax(1,1)*V(i+1,1);
%         Ry(i,N-2) = Ry(i,N-2)- Ax(N-2,N)*V(i+1,N);

%         GP1(1,i+1) = (Phi(i,1) - Phi(i+2,1))/(2*dy);
% 
%         Psi1 = dt*(Phi(i,1) - Phi(i+2,1))/(2*dy) - ;
%         Psi2

        Rx(i,1) = Rx(i,1);
        Rx(i,N-2) = Rx(i,N-2);
        Ry(i,1) = Ry(i,1) - Ax(1,1)*dt*(-Phi(i,1) + Phi(i+2,1))/(2*dy);
        Ry(i,N-2) = Ry(i,N-2)- Ax(N-2,N)*dt*(-Phi(i,N) + Phi(i+2,N))/(2*dy);
        % Thomas Algorithm Solver
        [U1] = ThomasSolve(Ax(:,2:N-1),Rx(i,:)');
        [V1] = ThomasSolve(Ax(:,2:N-1),Ry(i,:)');
        % Generate U1/3 and V1/3
        Us(i+1,2:N-1) = U1'; 
        Vs(i+1,2:N-1) = V1'; 
    end
    % Advancing Step 2
    for i = 1:N-2
        % Insert Boundary Condition
%         Us(2,i+1) = Us(2,i+1) - Ay(1,1)*Us(1,i+1);
%         Us(M-1,i+1) = Us(M-1,i+1) - Ay(M-2,M)*Us(M,i+1);
%         Vs(2,i+1) = Vs(2,i+1) - Ay(1,1)*Vs(1,i+1);
%         Vs(M-1,i+1) = Vs(M-1,i+1) - Ay(M-2,M)*Vs(M,i+1);
        Us(2,i+1) = Us(2,i+1) - Ay(1,1)*dt*( - Phi(1,i) + Phi(1,i+2))/(2*dx);
        Us(M-1,i+1) = Us(M-1,i+1) - Ay(M-2,M)*dt*( - Phi(M,i) + Phi(M,i+2))/(2*dx);
        Vs(2,i+1) = Vs(2,i+1);
        Vs(M-1,i+1) = Vs(M-1,i+1);
        % Thomas Algorithm Solver
        [U2] = ThomasSolve(Ay(:,2:M-1),Us(2:M-1,i+1));
        [V2] = ThomasSolve(Ay(:,2:M-1),Vs(2:M-1,i+1));
        % Generate U2/3 and V2/3
        U(2:M-1,i+1) = U(2:M-1,i+1) + U2; 
        V(2:M-1,i+1) = V(2:M-1,i+1) + V2; 
    end
    
    U(M,2:N-1) = dt*( - Phi(M,1:N-2) + Phi(M,3:N))/(6*dx);
    U(1,2:N-1) = U(1,2:N-1) + dt*( - Phi(1,1:N-2) + Phi(1,3:N))/(6*dx);
    U(1,1) = U(1,1) + dt*( -3*Phi(1,1) + 4*Phi(1,2) - Phi(1,3))/(6*dx);
    U(1,N) = U(1,N) + dt*( 3*Phi(1,N) - 4*Phi(1,N-1) + Phi(1,N-2))/(6*dx);
    U(M,1) = dt*( -3*Phi(M,1) + 4*Phi(M,2) - Phi(M,3))/(6*dx);
    U(M,N) = dt*( 3*Phi(M,N) - 4*Phi(M,N-1) + Phi(M,N-2))/(6*dx);
    
    V(2:M-1,1) = dt*(Phi(3:M,1) - Phi(1:M-2,1))/(6*dy);
    V(2:M-1,N) = dt*(Phi(3:M,N) - Phi(1:M-2,N))/(6*dy);
    V(1,1) = dt*( -3*Phi(1,1) + 4*Phi(2,1) - Phi(3,1))/(6*dy);
    V(1,N) = dt*( -3*Phi(1,N) + 4*Phi(2,N) - Phi(3,N))/(6*dy);
    V(M,1) = dt*( 3*Phi(M,1) - 4*Phi(M-1,1) + Phi(M-2,1))/(6*dy);
    V(M,N) = dt*( 3*Phi(M,N) - 4*Phi(M-1,N) + Phi(M-2,N))/(6*dy);
    
    
    
    % Solve for Pseudo Pressure Phi 
    %       (Using Direct Inverse Method as in Problem 1)
    DIV = zeros(M,N);
    DIV(2:M-1,2:N-1) = U(2:M-1,:)*D1x' + D1y*V(:,2:N-1);
    
%     DIV(1,:) = -(V(2,:)-V(1,:))/dy; 
%     DIV(M,:) = (V(M,:)-V(M-1,:))/dy;
%     DIV(:,1) = -(U(:,2)-U(:,1))/dx; 
%     DIV(:,N) = (U(:,N)-U(:,N-1))/dx;

    DIV(1,2:N-1) = (-3*V(1,2:N-1)+4*V(2,2:N-1)-V(3,2:N-1))/(2*dy) + (U(1,3:N) - U(1,1:N-2))/(2*dx); 
    DIV(M,2:N-1) = (3*V(M,2:N-1)-4*V(M-1,2:N-1)+V(M-2,2:N-1))/(2*dy) + (U(M,3:N) - U(M,1:N-2))/(2*dx); 
    DIV(2:M-1,1) = (-3*U(2:M-1,1)+4*U(2:M-1,2)-U(2:M-1,3))/(2*dx) + (U(3:M,1) - U(1:M-2,1))/(2*dy); 
    DIV(2:M-1,N) = (3*U(2:M-1,N)-4*U(2:M-1,N-1)+U(2:M-1,N-2))/(2*dx) + (U(3:M,N) - U(1:M-2,N))/(2*dy); 
    
    DIV(1,1) = (-3*V(1,1)+4*V(2,1)-V(3,1))/(2*dy) + (-3*U(1,1)+4*U(1,2)-U(1,3))/(2*dx);
    DIV(1,N) = (-3*V(1,N)+4*V(2,N)-V(3,N))/(2*dy) + (3*U(1,N)-4*U(1,N-1)+U(1,N-2))/(2*dx);
    DIV(M,1) = (3*V(M,1)-4*V(M-1,1)+V(M-2,1))/(2*dy) + (-3*U(M,1)+4*U(M,2)-U(M,3))/(2*dx);
    DIV(M,N) = (3*V(M,N)-4*V(M-1,N)+V(M-2,N))/(2*dy) + (3*U(M,N)-4*U(M,N-1)+U(M,N-2))/(2*dx);
    
    RHS = 3*DIV/dt;
    Phi = DirectInv(M,N,dx,dy,RHS);
    
    % Correct U,V with Pressure Corrector
    U(2:M-1,2:N-1) = U(2:M-1,2:N-1) - (dt/3)*Phi(2:M-1,:)*D1x';
    V(2:M-1,2:N-1) = V(2:M-1,2:N-1) - (dt/3)*D1y*Phi(:,2:N-1);
    
%     % Correct U,V with Pressure Corrector
%     U(2:M-1,2:N-1) = U(2:M-1,2:N-1) - dt*Phi(2:M-1,:)*D1x';
%     V(2:M-1,2:N-1) = V(2:M-1,2:N-1) - dt*D1y*Phi(:,2:N-1);
    
    % Take note w.r.t. Time
    sequence = [sequence,U((M-1)/2,(N-1)/2)];
    Utop = [Utop,mean(U(1,:))];
    
%     if mod(t/dt,5)
% %        Uout = 
% %        Vout = 
% %        Pout =
%         
%     end
    
    % Time Counter
    t = t+dt;
end
figure, plot(sequence,'o-')

figure,quiver(flip(U),flip(V))
hold on, contour(flip(Phi))