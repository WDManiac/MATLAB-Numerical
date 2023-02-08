%% Extrusion Heat Transfer Problem
%% Iteration parameters
dt = 0.001;% Time interval
U = 2;
dx=U*dt;

duration=10; % Total time length 10s

%% Physical parameters of the rod
D = 0.001; % Diameter of the rod

a = 0.001;  %    Diffusion coeff
k = 25;  %    Conductivity coeff
rhoC = k/a;

P = pi*D;
A = 0.25*pi*D^2;

Fo = a*dt/(dx^2);

%% Enviromental condition
h = 2500; %    Convection film coeff

T0 = 783.15; % initial temeprature : 500 degree Celsius assumed
Ta = 298.15; % ambient temeprature : 25 degree Celsius assumed

%% Indicators as shown in prof. Jaluria's case
C1=(h*D)/k;
C2=(U*D)/a;

%% Initialization
T=zeros(duration/dt,duration/dt);
T(:,1) = T0;% T(3,1)=T0;T(3,2)=T0;T(3,3)=T0;

[m,n]=size(T);

%% ITERATION START

for i=2:m
    for j=2:i
        if i==2 % first time interval
            T(i,j) = T(i,j-1) - (T(i-1,j-1)-Ta)*(h*dx/k);
        elseif j==i % last node of the rode
            T(i,j) = T(i,j-1) - (T(i-1,j-1)-Ta)*(h*dx/k);
        else % general case
            T(i,j)= (Fo )*T(i-1,j)+(1-2*Fo - h*P*dt/(rhoC*A))*T(i-1,j-1)+(Fo )*T(i,j-1) + h*P*dt*Ta/(rhoC*A);
            %- U*dt/(2*dx)
            %+ U*dt/(2*dx)
        end
    end
end
