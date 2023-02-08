%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An exemplary 1D PPE solver w/   %%
%% Homogeneous Neumann BC          %%
%% 534/Psolver                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50; 
dx = 2*pi/N;
x = [0:dx:2*pi]';

%% Form 1st order difference matrix D
D = zeros(N+1,N+1);
D(1,1) = -3; D(1,2) = 4; D(1,3) = -1;
D(N+1,N+1) = 3; D(N+1,N) = -4; D(N+1,N-1) = 1;
for i= 2:N
    D(i,i-1) = -1;
    D(i,i) = 0;
    D(i,i+1) = 1;
end

D = D/2/dx;
D2 = D*D;

D2(2:3,2) = D2(2:3,2)+(4/3)*D2(2:3,1);
D2(2:3,3) = D2(2:3,3)-(1/3)*D2(2:3,1);
D2(N-1:N,N) = D2(N-1:N,N)+(4/3)*D2(N-1:N,N+1)
D2(N-1:N,N-1) = D2(N-1:N,N-1)-(1/3)*D2(N-1:N,N+1)

D2 = D2(2:N,2:N);

% Test, g = -cos(x); Exact solution p = cos(x)-1
g = -cos(x);
p = zeros(size(x));

[V, L] = eig(D2);
R = inv(V)*g(2:N);
L = diag(L);
for i = 2:N
    if abs(L(i-1))>1e-8
       p(i) = R(i-1)/L(i-1);
    else 
       i
       p(i) = 0;
    end
end

p(2:N) = V*p(2:N);

p(1) = (4/3)*p(2)-(1/3)*p(3); 
p(N+1) = (-1/3)*p(N-1)+(4/3)*p(N);
% Offset constant
p = p - mean(p)-1;


figure(1);
clf;
plot(x,cos(x)-1,'r');
hold on; 
plot(x,p,'ro');
set(gca,'FontSize',14);
legend('cos(x)-1','p','Location','SouthEast');
grid on;
xlabel 'x'; ylabel 'Solution'



