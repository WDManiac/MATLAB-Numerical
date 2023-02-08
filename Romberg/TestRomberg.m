clear 
clc
close all

%% Test Romberg integration code
a = 0;
b = 100;
func = @(x)x.^2.*exp(0.271*sqrt(x));
% error restriction
epsilon = 0.0001;

% Computed from the code
[R,e,Y]=Romberg(func,a,b,epsilon);

% Reference result:
R_ref = integral(func,a,b);

disp(['Referencial integral with built-in MATLAB code: ',num2str(R_ref)])
disp(['Integral value computed: ',num2str(R)])
disp(['Iterative error: ',num2str(e)])
disp(['Absolute error: ',num2str(abs(R-R_ref))])
