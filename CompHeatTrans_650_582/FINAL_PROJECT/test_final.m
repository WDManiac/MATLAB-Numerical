function [T,Fo,aa,bb] = test_final(dt,dx,U,DetT,duration)
%   FTCS iteration method for FINAL PROJECT of CHT
%   dt ---- time step for iteration
%   dx ---- grid size of the mesh
%   DetT ---- time inteval
%   U ---- translation speed of the cylinder
%   duration ---- total time of our calculation

%   Fo ---- preset Fourier number
%   T ---- output matrix of result
%   Ta ---- ambient tempreture of atmosphere
%   T0 ---- boundary tempreture
%   aa ---- hD/k
%   bb ---- UD/a

a = 5*10^(-3);  %    Diffusion coeff
k = 25;  %    Conductivity coeff
h = 50; %    Convection film coeff
r = 0.001; % Radius of the rod
rhoC = k/a;

P = 2*pi*r;
A = pi*r^2;

aa=(h*r*2)/k;
bb=(U*r*2)/a;

T0 = 783.15; % initial temeprature : 500 degree Celsius assumed
Ta = 298.15; % ambient temeprature : 25 degree Celsius assumed

Fo= a*dt/(dx^2);

T=zeros(duration/dt,duration*U/dx);
T(:,1) = T0;% T(3,1)=T0;T(3,2)=T0;T(3,3)=T0;

L=zeros(duration/DetT,1);
N=zeros(duration/DetT,1);

lim = duration/DetT;

X0 = DetT*U/dx;
S = DetT/dt;

% T(1,1:X0)=T0;
N(1) = X0;

K=1;

while K<=lim
    
    L(K) = K*DetT*U; % Length of extruded rod as defined
    N(K) = K*X0; % Total number of nodes at this time step
    
    M1 =(K-1)*S+1;
    M2 =K*S; 
        
    T(M1,1:X0)=T0; 
        
    if K~=1
          T(M1,(X0+1):N(K))=T(M1-1,1:N(K-1));
    end
    

    
    for i=(M1+1):M2
    % calculate at each time step
        
        for j=2:N(K)
        % FTCS finite difference method for each time step
        
            if j~= N(K)
                T(i,j)= (Fo - U*dt/(2*dx))*T(i-1,j+1)+(1-2*Fo-h*P*dt/(rhoC*A))*T(i-1,j)+(Fo + U*dt/(2*dx))*T(i,j-1) + h*P*dt*Ta/(rhoC*A);
                % - U*dt/(2*dx)
                % + U*dt/(2*dx)
          
            elseif j==N(K)
                T(i,j) = (k*T(i,j-1)-h*dx*Ta)/(k-h*dx);%h*dx*(T(i-1,j-1)-Ta)/k + T(i,j-1);%(k*T(i,j-1)+dx*h*Ta)/(k+dx*h);
            end

        end
        
    end   
    
    K=K+1;
    
end

for i=1:duration/dt
    for j=1:duration*U/dx
        if T(i,j)~=0;
        T = (T-Ta)/(T0-Ta);
        end
    end
end

    
end