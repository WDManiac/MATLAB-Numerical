%% QR method for square real matrix eigenvalue and linear system solving

% Main Function
% ================================== %
function [Q,R,x] = qr_sol(A,b)
% nargin --- Number of input variables
[m,n] = size(A);

% This function only take care of square matrix, thus:
if m~=n
    disp 'Error: Oops, This function only deals with square matrix.'
    
    % Here is where actual computation begins:
else
    % ======== Initialize variables needed:
    M = A; % Using a dummy matrix M in case value of A is overwrote
    T = eye(m,n); % Dummy matrix T to record Q-matrix per iteration
    
    if nargin == 2
        % Initialize solution vector x:
        x = zeros(m,1);
        
        for i = 1:m-1
            u = M(:,i);
            
            % v-vector to be subtracted from u
            v = zeros(m,1);
            v(1:i) = u(1:i); % designate 1 to 1st ~ ith element
            % Note: sign of v(i) shall be opposite to that of A(i,i)
            if M(i,i) > 0
                sign = -1;
            else
                sign = 1;
            end
            
            % Include updated value of v(i)
            % Note: "sum" function take sum of all element in a vector
            %           "M.^2" means taking power of 2 over each element in M
            v(i) = sign*sqrt(sum(u(i:m).^2));
            
            % Householder matrix derived
            [H] = householder(u-v);
            
            % Recursively update dummy matrices M and T
            M = H*M;
            T = H*T;
            
        end
        R = M; % Upper triangle matrix R
        Q = T'; % Isotropic matrix Q, note that transpose is necessary!
        
        
        % Then solve for solution of linear system Ax=b:
        bT = Q'*b;
        % Here go row-by-row to compute each and every number in solution x:
        for i = m:-1:1
            if i == m
                x(i) = bT(i)/R(i,i);
            else
                x(i) = (bT(i) - sum(R(i,:)*x))/R(i,i);
            end
        end
        
    else
        disp 'Error: No or insuficient input variables!'
    end
end

end
% Here ends the main function
% ================================== %

% Sub-function - 1: Householder-matrix generator
function [H] = householder(u)
m = length(u);

% Normalize the vector
u = u/sqrt(sum(u.^2));

I = eye(m); % Generate identity matrix I of m dimension
H = I - 2*u*u';
end

