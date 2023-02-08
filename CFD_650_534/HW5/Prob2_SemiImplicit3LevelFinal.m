%% CFD HW5 Problem 2: 
% Author: Ran Li
% Semi-implicit Three-level Method for N-S Equation Solving
clc, clear, close all
 
% Initialization of Simulation:
% Mesh Size
M = 101;
N = 101;
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
tmax = 10; % Max Time for Iteration
t = 0; % Initial Time
dt = 0.001; % Time Interval
 
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
 
sequence1 = [];
sequence2 = [];
Utop = [];
 
% Time Iteration Starts Here
while t <= tmax
     
    error = 1;
    
    % Steady State Simulation Error Marin:  ERROR = 1;
    ERROR = 1; % Transient Simulation Error Marin
    
    counter = 0;
    % Iteration to get each time interval converge
    while ERROR <= 100 % Transient Simulation Error Marin
        % Steady State Simulation Error Marin: ERROR >= 10^-6
         
        % Compute Non-Linear Terms
        UU = U.*U;
        UV = U.*V;
        VV = V.*V;
        
        % Compute Non-linear Convection Term
        Nx = - (UU(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*UV(:,2:N-1));
        Ny = - (UV(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*VV(:,2:N-1));
        
        % Compute Linear Viscous Term
        Lx = D2y*U(:,2:N-1) + U(2:M-1,:)*D2x';
        Ly = D2y*V(:,2:N-1) + V(2:M-1,:)*D2x';
        
        
        % Compose RHS R for Advancing Process
        Rx = (dt/6)*(3*Nx - Nx0) + Lx;
        Ry = (dt/6)*(3*Ny - Ny0) + Ly;
        Nx0 = Nx; Ny0 = Ny;
         
        % Advancing Step 1
        for i = 1:M-2
            % Thomas Algorithm Solver
            [U1] = ThomasSolve(Ax(:,2:N-1),Rx(i,:)');
            [V1] = ThomasSolve(Ax(:,2:N-1),Ry(i,:)');
            % Generate U1/3 and V1/3
            Us(i+1,2:N-1) = U1';
            Vs(i+1,2:N-1) = V1';
        end
        % Advancing Step 2
        for i = 1:N-2
            % Thomas Algorithm Solver
            [U2] = ThomasSolve(Ay(:,2:M-1),Us(2:M-1,i+1));
            [V2] = ThomasSolve(Ay(:,2:M-1),Vs(2:M-1,i+1));
            % Generate U2/3 and V2/3
            Ustar(2:M-1,i+1) = U(2:M-1,i+1) + U2;
            Vstar(2:M-1,i+1) = V(2:M-1,i+1) + V2;
        end
         
         
        % Solve for Pseudo Pressure Phi
        %       (Using Direct Inverse Method as in Problem 1)
        DIV = zeros(M,N);
        DIV = Ustar*D1x' + D1y*Vstar ;
        RHS = 3*(DIV)/dt;%+Nstar
        Phi = DirectInv(M,N,dx,dy,RHS);
         
        % Correct U,V with Pressure Corrector
        Unew = Ustar - (dt/3)*Phi*D1x';
        Vnew = Vstar - (dt/3)*D1y*Phi;
         
        % 1> Steady State Error Margin
%         ERROR1 = mean(mean(abs(Unew - U)));
%         ERROR2 = mean(mean(abs(Vnew - V)));
%         ERROR = max([error1,error2]);
        
        % 2> Transient Error Margin
        ERROR = ERROR+1;
         
        U = Unew;
        V = Vnew;
        disp(['Time ',num2str(t),'-- Step No. ',num2str(counter),': Err = ', num2str(ERROR)])         
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
    sequence1 = [sequence1,U((M+1)/2,(N+1)/2)];
    sequence2 = [sequence2,V((M+1)/2,(N+1)/2)];
    Utop = [Utop,mean(U(1,:))];
     
     
end
 
%% Post Process
figure, plot(sequence1,'o-')
figure, plot(sequence2,'o-')
 
figure, plot(Utop,'*-')
 
figure,quiver(flip(U),flip(V))
hold on, contour(flip(Phi))
axis equal
 
Vorticity = D1y*U-V*D1x';
X=0:dx:1;
Y =0:dy:1;
[QuivX,QuivY] = meshgrid(0:dx:2,0:dy:1);
figure,contour(X,Y,flip(Vorticity))
hold on,quiver(X, Y, flip(U),-flip(V))
title('Vorticity, Re=100'),set(gca,'fontsize',14)
axis equal, grid on
% saveas(gcf,'Vorticity','jpg')
% saveas(gcf,'Vorticity','fig')
 

% steady state data processing
 
figure,plot(Y,flip(Vorticity(:,26)),'b')
hold on,plot(Y,flip(Vorticity(:,51)),'r','linewidth',2)
hold on,plot(Y,flip(Vorticity(:,76)),'c--')
legend('x=0.25','x=0.5','x=0.75')
set(gca,'fontsize',14)
grid on
title('w(y) from Simulation')
ylabel('w'),xlabel('y')
% saveas(gcf,'RanWsteady','jpg')
% saveas(gcf,'RanWsteady','fig')
 
figure,plot(Y,flip(U(:,26)),'b')
hold on,plot(Y,flip(U(:,51)),'r','linewidth',2)
hold on,plot(Y,flip(U(:,76)),'c--')
legend('x=0.25','x=0.5','x=0.75')
set(gca,'fontsize',14)
grid on
title('u(y) from Simulation')
ylabel('u'),xlabel('y')
% saveas(gcf,'RanUsteady','jpg')
% saveas(gcf,'RanUsteady','fig')
%  
figure,plot(Y,flip(-V(:,26)),'b')
hold on,plot(Y,flip(-V(:,51)),'r','linewidth',2)
hold on,plot(Y,flip(-V(:,76)),'c--')
legend('x=0.25','x=0.5','x=0.75')
set(gca,'fontsize',14)
grid on
title('v(y) from Simulation')
ylabel('v'),xlabel('y')
% saveas(gcf,'RanVsteady','jpg')
% saveas(gcf,'RanVsteady','fig')
 
%% Rescaling; L = 0.01 m; U = 1e-2 m/s; T = L/U = 1 s; nu = 1e-6 m^2/s
x = X*100;
y = Y*100;
t = (0:tmax/dt)*dt*10*3;
u = sequence1*1e2;
v = sequence2*1e2;
p0 = 1e3*1e4*1e4;
 
figure(1); clf;
plot(t,u,'--');grid on;
set(gca, 'FontSize', 14);
xlabel t; ylabel u;
title 'u(0.5, 0.5, t) for Re=1';
 
figure(2); clf;
plot(t,-v,'--');grid on;
set(gca, 'FontSize', 14);
xlabel t; ylabel v;
title 'v(0.5, 0.5, t) for Re=1';

%% Plotting with Reference Solution:
% x = X*100;
% y = Y*100;
% t = (0:tmax/dt)*dt*10*3;
% u = sequence1*1e2;
% v = sequence2*1e2;
% p0 = 1e3*1e4*1e4;
% figure(9)
% hold on,plot(t(1:334),0.01*u(1:334),'r--');grid on;
% legend('U_{ref}(0.5,0.5,t)','U_{sim}(0.5,0.5,t)')
% title 'u(0.5, 0.5, t) for Re=100';
% figure(10)
% hold on,plot(t(1:334),-0.01*v(1:334),'r--');grid on;
% legend('V_{ref}(0.5,0.5,t)','V_{sim}(0.5,0.5,t)')
% title 'v(0.5, 0.5, t) for Re=100';

% Pressure Analysis

% P2s = Prec(:,:,67);
% P4s = Prec(:,:,133);
% P6s = Prec(:,:,200);
% P8s = Prec(:,:,267);
% P10s = Prec(:,:,334);
% 
% % P2s = P2s-0.5*dt/Re
% 
% figure,subplot(1,2,1),contour(flip(P2s)),axis equal,subplot(1,2,2),contour(p_2s'),axis equal
% figure,subplot(1,2,1),contour(flip(P4s)),axis equal,subplot(1,2,2),contour(p_4s'),axis equal
% figure,subplot(1,2,1),contour(flip(P6s)),axis equal,subplot(1,2,2),contour(p_6s'),axis equal
% figure,subplot(1,2,1),contour(flip(P8s)),axis equal,subplot(1,2,2),contour(p_8s'),axis equal
% figure,subplot(1,2,1),contour(flip(P10s)),axis equal,subplot(1,2,2),contour(p_10s'),axis equal