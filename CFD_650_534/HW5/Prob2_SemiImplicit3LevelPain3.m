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

sequence1 = [];
sequence2 = [];
Utop = [];

% Time Iteration Starts Here
while t <= tmax
    
    error = 1;
    ERROR = 1;
%     ERROR0 = 0;
    counter = 0;
    % Iteration to get each time interval converge
    while ERROR <=100 %error >= 10^-6
        
        % Compute Non-Linear Terms
        UU = U.*U;
        UV = U.*V;
        VV = V.*V;
        
        % Compute Non-linear Convection Term
        
%         NxS = - U.*(U*D1x') - V.*(D1y*U);
%         NyS = - U.*(V*D1x') - V.*(D1y*V);
%         
        Nx = - (UU(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*UV(:,2:N-1));
        Ny = - (UV(2:M-1,:)*D1x(2:N-1,:)') - (D1y(2:M-1,:)*VV(:,2:N-1));
        
%         Nx = - U(2:M-1,2:N-1).*(U(2:M-1,:)*D1x(2:N-1,:)') - V(2:M-1,2:N-1).*(D1y(2:M-1,:)*U(:,2:N-1));
%         Ny = - U(2:M-1,2:N-1).*(V(2:M-1,:)*D1x(2:N-1,:)') - V(2:M-1,2:N-1).*(D1y(2:M-1,:)*V(:,2:N-1));
        
        % Compute Linear Viscous Term
            Lx = D2y*U(:,2:N-1) + U(2:M-1,:)*D2x';
            Ly = D2y*V(:,2:N-1) + V(2:M-1,:)*D2x';
        
        
        % Compose RHS R for Advancing Process
        Rx = (dt/6)*(3*Nx - Nx0) + Lx;
        Ry = (dt/6)*(3*Ny - Ny0) + Ly;
        Nx0 = Nx; Ny0 = Ny;
        
        % Advancing Step 1
        GP1(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,1)/3 - V(:,1));
        GP2(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,N)/3 - V(:,N));
%             GP1(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,1)/3 );
%             GP2(1,:) = Ay*(dt*D1y(:,2:M-1)*Phi(2:M-1,N)/3);
        GP3(1,:) = -Ay*U(:,1);
        GP4(1,:) = -Ay*U(:,N);
        for i = 1:M-2
            % Insert Boundary Condition
            %         Psi1 = GP1(i+1) - D2y(i+1,i:i+2)*GP1(i:i+2)';
            %         Psi2 = GP2(i+1) - D2y(i+1,i:i+2)*GP2(i:i+2)';
            Psi1 = GP1(i);% - D2y(i,i:i+2)*GP1(i:i+2)';
            Psi2 = GP2(i);% - D2y(i,i:i+2)*GP2(i:i+2)';
            
%             Rx(i,1) = Rx(i,1) - Ax(1,1)*GP3(i);
%             Rx(i,N-2) = Rx(i,N-2) - Ax(N-2,N)*GP4(i);
%             Ry(i,1) = Ry(i,1) - Ax(1,1)*Psi1;
%             Ry(i,N-2) = Ry(i,N-2)- Ax(N-2,N)*Psi2;
            % Thomas Algorithm Solver
            [U1] = ThomasSolve(Ax(:,2:N-1),Rx(i,:)');
            [V1] = ThomasSolve(Ax(:,2:N-1),Ry(i,:)');
            % Generate U1/3 and V1/3
            Us(i+1,2:N-1) = U1';
            Vs(i+1,2:N-1) = V1';
        end
        % Advancing Step 2
        for i = 1:N-2
%             % Insert Boundary Condition
%             Us(2,i+1) = Us(2,i+1) - Ay(1,1)*(Vtop - U(1,i+1) );
%             Us(M-1,i+1) = Us(M-1,i+1) - Ay(M-2,M)*( - U(M,i+1));
%             Vs(2,i+1) = Vs(2,i+1) - Ay(1,1)*V(1,i+1);
%             Vs(M-1,i+1) = Vs(M-1,i+1) - Ay(M-2,M)*V(M,i+1);
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
        
%         error1 = mean(mean(abs(Unew - U)));
%         error2 = mean(mean(abs(Vnew - V)));
%         error = max([error1,error2]);

%         ERROR1 = max(max(abs(Phi*D1x')));
%         ERROR2 = max(max(abs(D1y*Phi)));%abs(max(max(abs(Phi*D1x' + D1y*Phi))) - ERROR0) ;
%         ERROR = max([ERROR1, ERROR2]);

        ERROR = ERROR+1;
%         ERROR0 = ERROR;
        
        U = Unew;
        V = Vnew;
        
%         Max1 = max(abs(Unew(M,:)));
%         Max2 = max(abs(Vnew(:,N)));
%         Max3 = max(abs(Vnew(:,1)));
%         Max4 = max(abs(Unew(1,:)-Vtop*ones(1,N)));
        
        
%         error = max([Max1,Max2,Max3,Max4]);
%         disp(['Step No. ',num2str(counter),': Err = ', num2str(error)])
        disp(['Step No. ',num2str(counter),': Err = ', num2str(ERROR)])
        
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
    sequence1 = [sequence1,U((M-1)/2,(N-1)/2)];
    sequence2 = [sequence2,V((M-1)/2,(N-1)/2)];
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
saveas(gcf,'Vorticity','jpg')
saveas(gcf,'Vorticity','fig')

StreamX = dx/2:dx:(N-1)*dx-dx/2;
StreamY = dy/2:dy:(N-1)*dy-dy/2;
Uhat = zeros(M,N-1); Vhat = zeros(M-1,N);
for i = 1:M
    for j = 1:N
        if i~=M && j~=N
            Uhat(i,j) = (U(i,j) + U(i,j+1))/2;
            Vhat(i,j) = -(V(i,j) + V(i+1,j))/2;
        elseif i == M && j ~= N
            Uhat(i,j) = (U(i,j) + U(i,j+1))/2;
        elseif j == N && i ~= M
            Vhat(i,j) = -(V(i,j) + V(i+1,j))/2;
        end
    end
end
for i = 1:M-1
    for j = 1:N-1
        Stream(i,j) = (dy*((Uhat(i,j)) + (Uhat(i+1,j)))/2 - dx*((Vhat(i,j)) + (Vhat(i,j+1)))/2)*dx*dy;
    end
end
figure,contour(StreamX,StreamY,flip(Stream))
hold on,quiver(X, Y, flip(U),-flip(V))
axis equal, grid on

% steady state data

figure,plot(Y,flip(Vorticity(:,13)),'b')
hold on,plot(Y,flip(Vorticity(:,26)),'r','linewidth',2)
hold on,plot(Y,flip(Vorticity(:,48)),'c--')
legend('x=0.25','x=0.5','x=0.75')
set(gca,'fontsize',14)
grid on
title('w(y) from Simulation')
ylabel('w'),xlabel('y')
saveas(gcf,'RanWsteady','jpg')
saveas(gcf,'RanWsteady','fig')

figure,plot(Y,flip(U(:,13)),'b')
hold on,plot(Y,flip(U(:,26)),'r','linewidth',2)
hold on,plot(Y,flip(U(:,48)),'c--')
legend('x=0.25','x=0.5','x=0.75')
set(gca,'fontsize',14)
grid on
title('u(y) from Simulation')
ylabel('u'),xlabel('y')
saveas(gcf,'RanUsteady','jpg')
saveas(gcf,'RanUsteady','fig')

figure,plot(Y,flip(-V(:,13)),'b')
hold on,plot(Y,flip(-V(:,26)),'r','linewidth',2)
hold on,plot(Y,flip(-V(:,48)),'c--')
legend('x=0.25','x=0.5','x=0.75')
set(gca,'fontsize',14)
grid on
title('v(y) from Simulation')
ylabel('v'),xlabel('y')
saveas(gcf,'RanVsteady','jpg')
saveas(gcf,'RanVsteady','fig')

%% Rescaling; L = 0.01 m; U = 1e-2 m/s; T = L/U = 1 s; nu = 1e-6 m^2/s
x = X*100;
y = Y*100;
t = (0:tmax/dt)*dt;
u = sequence1*1e2;
v = sequence2*1e2;
p0 = 1e3*1e4*1e4;

figure(1); clf;
plot(t,u,'--');grid on;
set(gca, 'FontSize', 14);
xlabel t; ylabel u;
title 'u(0.5, 0.5, t) for Re=1';

figure(2); clf;
plot(t,v,'--');grid on;
set(gca, 'FontSize', 14);
xlabel t; ylabel v;
title 'v(0.5, 0.5, t) for Re=1';