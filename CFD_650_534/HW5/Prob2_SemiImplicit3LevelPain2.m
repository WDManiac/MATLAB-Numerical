%% CFD HW5 Problem 2: 
% Author: Ran Li
% Semi-implicit Three-level Method for N-S Equation Solving
clc, clear, close all

% Initialization of Simulation:
% Mesh Size
M = 41;
N = 41;
% Reynolds Number involved
Re = 100;
% Due to this nondimension, driving velocity V and size of cavity 
% domain L could be both set as 1

% Generate Mesh
[dx,dy,mat_ini] = blockmesh(N,M,1,1);
U = mat_ini; % Horizontal Velocity Field
V = mat_ini; % Vertical Velocity Field
Phi = mat_ini; % Pressure Field (Pseudo)


% Setup Time Iteration
tmax = 0.8; % Max Time for Iteration
t = 0; % Initial Time
dt = 0.005; % Time Interval

% BCs for Velocity Field:
Vtop = 1;
U(1,:) = Vtop*ones(1,N);

Us = U; Vs = V;
Unew = U; Vnew = V;
Ustar = U; Vstar = V;

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

% Record first step data:
tn = 1;
Urec(:,:,tn) = U;
Vrec(:,:,tn) = V;
Prec(:,:,tn) = Phi;

sequence = [];
Utop = [];

% Time Iteration Starts Here
while t <= tmax
    
    error = 1;
    counter = 0;
    % Iteration to get each time interval converge
    while error >= 10^-5
        
        % Compute Non-Linear Terms
        UU = U.*U;
        UV = U.*V;
        VV = V.*V;
        
        % Compute Non-linear Convection Term
%         Nx = - (UU(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*UV(:,2:N-1));
%         Ny = - (UV(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*VV(:,2:N-1));
        
                Nx = - U(2:M-1,2:N-1).*(U(2:M-1,:)*D1x(2:N-1,:)') - V(2:M-1,2:N-1).*(D1y(2:M-1,:)*U(:,2:N-1));
                Ny = - U(2:M-1,2:N-1).*(V(2:M-1,:)*D1x(2:N-1,:)') - V(2:M-1,2:N-1).*(D1y(2:M-1,:)*V(:,2:N-1));
        
        % Compute Linear Viscous Term
%             Lx = D2y*U(:,2:N-1) + U(2:M-1,:)*D2x';
%             Ly = D2y*V(:,2:N-1) + V(2:M-1,:)*D2x';
        
        % Compute un Term on RHS
        for i=1:M % (I + (dt/2Re)*Ay)
            for j=1:N-2
                L1x(i,j)=beta*U(i,j+2)/3 + (1 - 2*beta)*U(i,j+1)/3 + beta*U(i,j)/3;
                L1y(i,j)=beta*V(i,j+2)/3 + (1 - 2*beta)*V(i,j+1)/3 + beta*V(i,j)/3;
            end
        end
        for i=1:M-2 % (I + (dt/2Re)*Ax)
            for j=1:N-2
                Lx(i,j)=alpha*L1x(i+2,j) + (1 - 2*alpha)*L1x(i+1,j) + alpha*L1x(i,j);
                Ly(i,j)=alpha*L1y(i+2,j) + (1 - 2*alpha)*L1y(i+1,j) + alpha*L1y(i,j);
            end
        end
        
        % Compose RHS R for Advancing Process
        Rx = (dt/6)*(3*Nx - Nx0) + Lx;
        Ry = (dt/6)*(3*Ny - Ny0) + Ly;
        
        % Advancing Step 1
%         GP1(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,1)/3);
%         GP2(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,N)/3);
        GP1(1,:) = Ay*(dt*D1y*Phi(:,1)/3);
        GP2(1,:) = Ay*(dt*D1y*Phi(:,N)/3);
        for i = 1:M-2
            % Insert Boundary Condition
            Psi1 = GP1(i);% - D2y(i,i:i+2)*GP1(i:i+2)';
            Psi2 = GP2(i);% - D2y(i,i:i+2)*GP2(i:i+2)';
            Rx(i,1) = Rx(i,1);
            Rx(i,N-2) = Rx(i,N-2);
            Ry(i,1) = Ry(i,1) - Ax(1,1)*Psi1;
            Ry(i,N-2) = Ry(i,N-2) - Ax(N-2,N)*Psi2;
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
            Us(2,i+1) = Us(2,i+1) - Ay(1,1)*(Vtop + dt*(- Phi(1,i) + Phi(1,i+2))/(2*dx)/3);
            Us(M-1,i+1) = Us(M-1,i+1) - Ay(M-2,M)*(dt*(- Phi(M,i) + Phi(M,i+2))/(2*dx)/3);
            Vs(2,i+1) = Vs(2,i+1);
            Vs(M-1,i+1) = Vs(M-1,i+1);
            % Thomas Algorithm Solver
            [U2] = ThomasSolve(Ay(:,2:M-1),Us(2:M-1,i+1));
            [V2] = ThomasSolve(Ay(:,2:M-1),Vs(2:M-1,i+1));
            % Generate U2/3 and V2/3
            Ustar(2:M-1,i+1) = U2;
            Vstar(2:M-1,i+1) = V2;
        end
        Ustar(M,:) = (dt/3)*(Phi(M,:)*D1x');
        Vstar(:,N) = (dt/3)*(D1y*Phi(:,N));
        Vstar(:,1) = (dt/3)*(D1y*Phi(:,1));
        Ustar(1,:) = Vtop*ones(1,N) + (dt/3)*(Phi(1,:)*D1x');
        
        % Solve for Pseudo Pressure Phi
        %       (Using Direct Inverse Method as in Problem 1)
%         NxStar = Ustar.*(Ustar*D1x') + Vstar.*(D1y*Ustar);
%         NyStar = Ustar.*(Vstar*D1x') + Vstar.*(D1y*Vstar);
%         Nstar = NxStar*D1x' + D1y*NyStar;
        DIV = zeros(M,N);
        DIV = Ustar*D1x' + D1y*Vstar ;
        RHS = 3*(DIV)/dt;%+Nstar
        Phi = DirectInv(M,N,dx,dy,RHS);
        
        % Correct U,V with Pressure Corrector
        Unew = Ustar - (dt/3)*Phi*D1x';
        Vnew = Vstar - (dt/3)*D1y*Phi;
        
        U = Unew;
        V = Vnew;
        
        Max1 = max(abs(Unew(M,:)));
        Max2 = max(abs(Vnew(:,N)));
        Max3 = max(abs(Vnew(:,1)));
        Max4 = max(abs(Unew(1,:)-Vtop*ones(1,N)));
        
        
        error = max([Max1,Max2,Max3,Max4]);
        disp(['Step No. ',num2str(counter),': Err = ', num2str(error)])
        
        counter = counter + 1;
    end
    
    % Time Counter
    t = t+dt;
    tn = tn+1;
    
    Urec(:,:,tn) = Unew;
    Vrec(:,:,tn) = Vnew;
    Prec(:,:,tn) = Phi;
    
    U = Unew; V = Vnew;
    
    Nx0 = Nx; Ny0 = Ny;
    
    % Take note w.r.t. Time
    sequence = [sequence,U((M-1)/2,(N-1)/2)];
    Utop = [Utop,mean(U(1,:))];
    
    
end
figure, plot(sequence,'o-')

figure, plot(Utop,'*-')

figure,quiver(flip(U),flip(V))
hold on, contour(flip(Phi))