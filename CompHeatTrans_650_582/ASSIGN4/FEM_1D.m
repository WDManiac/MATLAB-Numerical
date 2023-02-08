function [T,M,G] = FEM_1D(dx,L,Ta,h,T0,k,r)

N=L/dx;

M=zeros(N+1);

% Initializing output vector
T=zeros(N+1,1);
%  T(1)=T0; T(N+1)=T0;

% Cross section coefficients
A=pi*r^2; P=2*pi*r;

%% Initializing rigidity matrix
for i=1:N+1
    
    if i==1
        M(i,i)=dx/k;M(i,i+1)=1;
    elseif i~=N+1
        M(i,i-1)=k/(dx^2);
        M(i,i)=-2*k/(dx^2);
        M(i,i+1)=k/(dx^2);
    else
        M(i,i)=-dx/k;M(i,i-1)=1;
    end
end

% Eliminating nodes on Dirichlet BCs
%   M2=M(2:N,2:N);
 
%% Initializing convection terms

I=eye(N+1);

I=I*(h*P/A);
  I(1,1)=0;I(N+1,N+1)=0;

%% Combine into final rigidity matrix

G=M-I;

%% Bias vector

b=ones(N+1,1);
b=b*(h*P*Ta/A);
b(1)=T0;b(N+1)=T0;

%% Solve

T=G\b;%(2:N)
T(1,1)=T0;T(N+1,N+1)=T0;

