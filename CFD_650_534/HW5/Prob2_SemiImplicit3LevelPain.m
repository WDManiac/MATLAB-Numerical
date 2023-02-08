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
tmax = 0.5; % Max Time for Iteration
t = 0; % Initial Time
dt = 0.1; % Time Interval

% BCs for Velocity Field:
Vtop = 0.5;
U(1,:) = Vtop*ones(1,N);

% 1st Order Diff. Operator Matrix
D1x = zeros(N,N); D1y = zeros(M,M);
D1x(1,1) = -3; D1x(1,2) = 4; D1x(1,3) = -1; 
D1x(N,N) = 3; D1x(N,N-1) = -4; D1x(N,N-2) = 1;
D1x(2,1) = -1; D1x(N-1,N) = 1;
D1x(2:N-1,2:N-1) = full(gallery('tridiag',N-2,-1,0,1));
D1y(1,1) = -3; D1y(1,2) = 4; D1y(1,3) = -1; 
D1y(M,M) = 3; D1y(M,M-1) = -4; D1y(M,M-2) = 1;
D1y(2,1) = -1; D1y(M-1,M) = 1;
D1y(2:M-1,2:M-1) = full(gallery('tridiag',M-2,-1,0,1));
D1x = D1x/(2*dx); D1y = D1y/(2*dy);

% Discretization Coefficients
alpha = dt/(2*Re*dx^2);
beta = dt/(2*Re*dy^2);

% 2nd Order Diff. Operator Matrix (Linear Term)
% D2x = dt*D1x*D1x/(3*Re);
% D2y = dt*D1y*D1y/(3*Re);
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
Nx0 = -UU(2:M-1,:)*D1x(2:N-1,:)' - D1y(2:M-1,:)*UV(:,2:N-1);
Ny0 = -UV(2:M-1,:)*D1x(2:N-1,:)' - D1y(2:M-1,:)*VV(:,2:N-1);

sequence = [];
Utop = [];

% Time Iteration Starts Here
while t <= tmax
    while error >= 10^-10
    
    % Compute Non-Linear Terms
    UU = U.*U;
    UV = U.*V;
    VV = V.*V;
    
    % Compute Non-linear Convection Term
%     Nx = - (UU(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*UV(:,2:N-1));
%     Ny = - (UV(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*VV(:,2:N-1));
    
    Nx = - U(2:M-1,2:N-1).*(U(2:M-1,:)*D1x(2:N-1,:)') - V(2:M-1,2:N-1).*(D1y(2:M-1,:)*U(:,2:N-1));
    Ny = - U(2:M-1,2:N-1).*(V(2:M-1,:)*D1x(2:N-1,:)') - V(2:M-1,2:N-1).*(D1y(2:M-1,:)*V(:,2:N-1));
    
    % Compute Linear Viscous Term
%     Lx = D2y(2:M-1,:)*U(:,2:N-1) + U(2:M-1,:)*D2x(2:N-1,:)';
%     Ly = D2y(2:M-1,:)*V(:,2:N-1) + V(2:M-1,:)*D2x(2:N-1,:)';
    Lx = D2y*U(:,2:N-1) + U(2:M-1,:)*D2x';
    Ly = D2y*V(:,2:N-1) + V(2:M-1,:)*D2x';
    
    % Compose RHS R for Advancing Process
    Rx = (dt/3)*(3*Nx - Nx0) + Lx;
    Ry = (dt/3)*(3*Ny - Ny0) + Ly;
    Nx0 = Nx; Ny0 = Ny;
    
    % Advance In Time
    dT = dt/3;
    % Advancing Step 1
    GP1(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,1)/3 - V(:,1));
    GP2(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,N)/3 - V(:,N));
%     GP1(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,1)/3);
%     GP2(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,N)/3);
    GP3(1,:) = -Ay*U(:,1);
    GP4(1,:) = -Ay*U(:,N);
    for i = 1:M-2
        % Insert Boundary Condition
%         Psi1 = GP1(i+1) - D2y(i+1,i:i+2)*GP1(i:i+2)';
%         Psi2 = GP2(i+1) - D2y(i+1,i:i+2)*GP2(i:i+2)';
        Psi1 = GP1(i);% - D2y(i,i:i+2)*GP1(i:i+2)';
        Psi2 = GP2(i);% - D2y(i,i:i+2)*GP2(i:i+2)';

%         Rx(i,1) = Rx(i,1) - GP3(i);
%         Rx(i,N-2) = Rx(i,N-2) - GP4(i);
        Ry(i,1) = Ry(i,1) - Ax(1,1)*Psi1;
        Ry(i,N-2) = Ry(i,N-2)- Ax(N-2,N)*Psi2;
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
%         Us(2,i+1) = Us(2,i+1) - Ay(1,1)*(dt*( Phi(1,i) - Phi(1,i+2))/(2*dx)/3 + Vtop);
%         Us(M-1,i+1) = Us(M-1,i+1) - Ay(M-2,M)*(dt*( Phi(M,i) - Phi(M,i+2))/(2*dx)/3);
        Us(2,i+1) = Us(2,i+1) - Ay(1,1)*(dt*( Phi(1,i) - Phi(1,i+2))/(2*dx)/3);
        Us(M-1,i+1) = Us(M-1,i+1) - Ay(M-2,M)*(dt*( Phi(M,i) - Phi(M,i+2))/(2*dx)/3);
%         Vs(2,i+1) = Vs(2,i+1) + Ay(1,1)*V(1,i+1);
%         Vs(M-1,i+1) = Vs(M-1,i+1) + Ay(M-2,M)*V(M,i+1);
        % Thomas Algorithm Solver
        [U2] = ThomasSolve(Ay(:,2:M-1),Us(2:M-1,i+1));
        [V2] = ThomasSolve(Ay(:,2:M-1),Vs(2:M-1,i+1));
        % Generate U2/3 and V2/3
        U(2:M-1,i+1) = U(2:M-1,i+1) + U2; 
        V(2:M-1,i+1) = V(2:M-1,i+1) + V2; 
    end
    U(M,:) = (dt/3)*(Phi(M,:)*D1x');
    V(:,N) = (dt/3)*(D1y*Phi(:,N));
    V(:,1) = (dt/3)*(D1y*Phi(:,1));
    U(1,:) = Vtop*ones(1,N) + (dt/3)*(Phi(1,:)*D1x');
    
    % Solve for Pseudo Pressure Phi 
    %       (Using Direct Inverse Method as in Problem 1)
    DIV = zeros(M,N);
    DIV = U*D1x' + D1y*V;
    RHS = 3*DIV/dt;
    Phi = DirectInv(M,N,dx,dy,RHS);
    
    % Correct U,V with Pressure Corrector
    U = U - (dt/3)*Phi*D1x';
    V = V - (dt/3)*D1y*Phi;
    
    error1 = mean(mean(abs(Une)))
        
        
        error = max([Max1,Max2,Max3,Max4]);
        disp(['Step No. ',num2str(counter),': Err = ', num2str(error)])
    
    counter = 1;
    end
    
    counter = 0;
    
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

figure, plot(Utop,'*-')

figure,quiver(flip(U),flip(V))
hold on, contour(flip(Phi))