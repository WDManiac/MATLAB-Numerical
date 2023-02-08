function [G,d,M,s]=spline3(X,Y,xp)
% Get size of input first
n = length(X);
m = length(Y);
t = length(xp);

% if n == m % Natural Cubic Spline Interpolation
% Boundary condition: M0 = 0 = Mn

G = zeros(n); % Intialize Coefficient Matrix
h = zeros(n-1,1);
mu = zeros(n-2,1);
lambda = zeros(n-2,1);
d = zeros(n-2,1);

for i=1:n-1
    h(i)=X(i+1)-X(i);
end

% Then acquire parameters for Coefficient Matrix G
for i = 1:n-2
    mu(i) = h(i)/(h(i)+h(i+1));
    lambda(i) = 1 - mu(i);
end

[G] = CoeffOrg(mu,lambda);

% Then compute bias vector d:
for i = 2:n-1
    d(i-1)=(6/(h(i)+h(i-1)))*((Y(i+1) - Y(i))/h(i)-(Y(i) - Y(i-1))/h(i-1));
end

M = G\d;
M = [0; M; 0];

% Then interpolate with xp:
for k=1:t
    for i=1:n-1
        if (xp(k)<=X(i+1))&&(xp(k)>=X(i))
            p1=M(i,1)*(X(i+1)-xp(k))^3/(6*h(i));
            p2=M(i+1,1)*(xp(k)-X(i))^3/(6*h(i));
            p3=(Y(i)-M(i,1)/6*(h(i))^2)*(X(i+1)-xp(k))/h(i);
            p4=(Y(i+1)-M(i+1,1)/6*(h(i))^2)*(xp(k)-X(i))/h(i);
            s(k)=p1+p2+p3+p4;
            break;
        else
            s(k)=0;
        end
    end
end

% elseif n == m-2
%
% end
end

% Sub-function: Coefficient Matrix Organizer
function [G] = CoeffOrg(m,l)
L = length(m); L = L;

G = zeros(L);

for i = 1:L
    G(i,i) = 2;
    if i == 1
        G(i,2) = l(i);
        G(i,L) = m(i);
    elseif i == L
        G(i,1) = l(i);
        G(i,L-1) = m(i);
    else
        G(i,i+1) = l(i);
        G(i,i-1) = m(i);
    end
end

end
