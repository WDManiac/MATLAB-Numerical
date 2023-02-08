function [Mat_fin, MatW, K, Error1, Error2] = falseTrans2DNS(mat,W,Re,dx,dy,lim,sf)
%   False Transient Method for 2D domain
%   mat ---- initialized mesh
%   dx, dy ---- grid size of the mesh
%   lim ---- fail safe of iteration(max operation number)
%   sf ---- scale factor for stability

beta = dx/dy;
rho = sf * (2 * (1 + beta^2))^-1;

%   Grab size of input matrix for initialization
[m,n] = size(mat);

%   Initialization
Mat_fin = mat;%      Result matrix
Mat = mat;%     Intermediate matrix
MatW = W;
E1 = zeros(m-1,n-1);%    Error matrix
E2 = E1; 
Error1 = 1;%     Convergence indicator for psi
Error2 = 1;%     Convergence indicator for vorticity
K = 1;%     Counter of Iteration

U = zeros(m-2,n-2); V = U; % Initialize velocity matrices
N = U; % Initialize tail term in voritcity equation
%   Iteration start!
while (Error1 >=10^-8 || Error2 >=10^-8)
    
    % stop if iteration number exceeds limit
    if K >=lim
        break
    end
    
    % ---- Part 1: Solve for psi (stream function)
    for i = 2:(m-1)
        for j = 2:(n-1)
            
            Mat_fin(i,j) = Mat(i,j) + rho*(Mat(i+1,j) - 2*Mat(i,j) + Mat(i-1,j)...
                                         + beta^2*(Mat(i,j+1) - 2*Mat(i,j) + Mat(i,j-1)) ...
                                         + dx^2*W(i,j));
            
            E1(i,j) = abs(Mat_fin(i,j) - Mat(i,j));%  Store the error
            Mat(i,j) = Mat_fin(i,j);
        end
    end
    
    Error1 = max(max(E1));
    %error = sum(sum(E));%     Average error for this iteration
    
    % ---- Part 2: Insert BCs for Vorticity Field
    W(1,:) = (7*Mat_fin(1,:) - 8*Mat_fin(2,:) + Mat_fin(3,:))/(2*dy^2) - 3/dy;
    W(m,:) = (7*Mat_fin(m,:) - 8*Mat_fin(m-1,:) + Mat_fin(m-2,:))/(2*dy^2);
    W(:,1) = (7*Mat_fin(:,1) - 8*Mat_fin(:,2) + Mat_fin(:,3))/(2*dx^2);
    W(:,n) = (7*Mat_fin(:,n) - 8*Mat_fin(:,n-1) + Mat_fin(:,n-2))/(2*dx^2);
    
    % ---- Part 3: Acquire Velocity Field and Tail Term
    for i = 2:m-1
        for j = 2:n-1
            U(i-1,j-1) = (Mat_fin(i,j+1) - Mat_fin(i,j-1))/(2*dy);
            V(i-1,j-1) = -(Mat_fin(i+1,j) - Mat_fin(i-1,j))/(2*dx);
            Wx = (W(i+1,j) - W(i-1,j))/(2*dx);
            Wy = (W(i,j+1) - W(i,j-1))/(2*dy);
            N(i-1,j-1) = U(i-1,j-1)*Wx + V(i-1,j-1)*Wy;
        end
    end
    
    % ---- Part 4: Solve for Vorticity
    for i = 2:(m-1)
        for j = 2:(n-1)
            
            MatW(i,j) = W(i,j) + rho*((W(i+1,j) - 2*W(i,j) + W(i-1,j)...
                       + beta^2*(W(i,j+1) - 2*W(i,j) + W(i,j-1)))/Re ...
                       - dx^2*N(i-1,j-1));
            
            E2(i,j) = abs(MatW(i,j) - W(i,j));%  Store the error
            W(i,j) = MatW(i,j);
        end
    end
    
    Error2 = max(max(E2));
    %error = sum(sum(E));%     Average error for this iteration
    
    K = K+1;
    if mod(K,20) == 0
        disp(['Error of psi @ Step',num2str(K),' = ',num2str(Error1)])
        disp(['Error of vorticity @ Step',num2str(K),' = ',num2str(Error2)])
    end
end    