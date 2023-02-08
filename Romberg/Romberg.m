function [R,err,Y] = Romberg(func,a,b,tol,margin)
% func --- single variable function to be integrated, use function handle
% a --- lower bound of integrating domain
% b --- upper bound of integrating domain
% tol --- error tolerence to be designated by the user
% margin --- max number of  iterations step (avoid dead loop) to be designated by the user


if nargin ~=5 && nargin ~=4
    disp 'Error: Needs 4 or 5 inputs, check function code for details.'
else
    if a >= b
        disp 'Error: 2nd input must be smaller than the 3rd'
    else
        % open a figure window for error visualization
        figure
        
        if nargin == 5 % User designated maximum iteration
            halt = margin;
        elseif nargin == 4 % Use default maximum iteration of 1000
            disp 'Warning: Maximum iteration number un-specified, use 1000 by default.'
            halt = 1000;
        end
        
        % Size of each integrating sub-domain
        h = (b-a);        
        
        % For composite trapezoidal integration: T(h) = h(0.5*S0+S1)
        %      where: S0 = f(a)+f(b)
        %                  S1 = f(x_1) + f(x_2) + ... + f(x_n-1),   x_i = a + 0.5*(2*i-1)h,   i = 1,2, ... , n-1
        
        % Get the S0 with end points function values:
        Y = zeros(3,1);
        Y(1,1) = h*(func(a)+func(b))/2;
        
        % Initialize error to be 1 as starter
        err = 1.0;
        iter = 0; % Mark iteration step
        N = 1; % Intially take the whole margin
        
        while (err>=tol)%&&(iter<N)
            % Update iteration step
            iter = iter + 1;       
            h = h/2; % divide sub-division by half
            
            % Get the S1 per iteration:
            S1 = 0;
            for i = 1:N
                x = a + (2*i-1)*h;
                S1 =S1 + func(x);
            end
            Y(iter+1,1) = Y(iter,1)/2 + h*S1;
            
            N = 2*N;
            for k = 1:iter
                Y(iter+1,k+1) = Y(iter+1, k)+(Y(iter+1,k)-Y(iter,k))/(4^k-1);
            end
            
            err = abs(Y(iter,iter)-Y(iter+1,iter+1));
            eps(iter) = err;
            hold on, plot(iter,err,'ro','linewidth',1.5)
        end
        
        R = Y(iter+1,iter+1);
        
        hold on, plot(eps,'k--','linewidth',.5)
        grid on
        set(gca,'fontsize',14)
        xlabel 'Iteration Step'
        ylabel 'Error'
        title 'Romberg Iteration Error'
        
    end
    
end

end


