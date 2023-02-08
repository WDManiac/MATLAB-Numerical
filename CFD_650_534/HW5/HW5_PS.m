%% MAE 534 S15 HW5 Postscript
%% Data sets: Time_split_Re1 or Time_split_Re100

clear all;
load Time_split_Re100;

%% Rescaling; L = 0.01 m; U = 1e-2 m/s; T = L/U = 1 s; nu = 1e-6 m^2/s
x = x*100;
y = y*100;
t = t*1;
u = u05*1e2;
v = v05*1e2;
p0 = 1e3*1e4*1e4;

figure(1); clf;
plot(t,u);grid on;
set(gca, 'FontSize', 14);
xlabel t; ylabel u;
title 'u(0.5, 0.5, t) for Re=1';

figure(2); clf;
plot(t,v);grid on;
set(gca, 'FontSize', 14);
xlabel t; ylabel v;
title 'v(0.5, 0.5, t) for Re=1';