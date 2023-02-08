function [Mat_fin, K, Error, Enote] = falseTrans2D(mat,W,dx,dy,lim,sf,sol)
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
E = zeros(m-1,n-1);%    Error matrix
Error = 1;%     Convergence indicator
K = 1;

%   Iteration start!
while Error >=10^-8
    
    % stop if iteration number exceeds limit
    if K >=lim
        break
    end
    
    for i = 2:(m-1)
        for j = 2:(n-1)
            
            Mat_fin(i,j) = Mat(i,j) + rho*(Mat(i+1,j) - 2*Mat(i,j) + Mat(i-1,j)...
                                         + beta^2*(Mat(i,j+1) - 2*Mat(i,j) + Mat(i,j-1)) ...
                                         + dx^2*W(i,j));
            
            E(i,j) = abs(Mat_fin(i,j) - Mat(i,j));%  Store the error
            Mat(i,j) = Mat_fin(i,j);
        end
    end
    
    Error = max(max(E));
    %error = sum(sum(E));%     Average error for this iteration
    K = K+1;
    Enote(K) = max(max(abs(Mat_fin - sol)));
    if mod(K,200) == 0
        disp(['Error @ Step',num2str(K),' = ',num2str(Error)])
    end
end    