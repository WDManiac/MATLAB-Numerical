%% MAE 534 HW4 Post-Calcuation
%% Based on psi-w 
%% 4/17/2015
psi = Psi_Re100';
w = w_Re100';
contour(x0,y0,psi);
figure(1);set(gca,'FontSize',14);
xlabel x; ylabel y;title('\psi for Re=1');

%% Velocities
u = zeros(size(y0,2),size(x0,2));
v = zeros(size(y0,2),size(x0,2));
dx = x0(2)-x0(1);
dy = y0(2)-y0(1);
% u = dpsi/dx
u(2:end-1,:) = psi(3:end,:)-psi(1:end-2,:);
u(1,:)= -3*psi(1,:) + 4*psi(2,:) - psi(3,:);
u(end,:) = 3*psi(end,:) - 4*psi(end-1,:) + psi(end-2,:);
u = u/2/dy;
% v = -dpsi/dx
v(:,2:end-1) = psi(:,3:end)-psi(:,1:end-2);
v(:,1)= -3*psi(:,1) + 4*psi(:,2) - psi(:,3);
v(end,:) = 3*psi(:,end) - 4*psi(:,end-1) + psi(:,end-2);
v=-v/2/dx;

% figure(2);
% xp = [1:4:101];
% quiver(u(xp,xp),v(xp,xp));

figure(2);clf
plot(y0, u(:,(length(x0)-1)/4)*1e5);hold on;
plot(y0, u(:,(length(x0)-1)/2)*1e5,'LineWidth',2);
plot(y0, u(:,(length(x0)-1)/4*3)*1e5,'--');
legend('x=0.25','x=0.5','x=0.75','Location','NorthWest')
set(gca,'FontSize',14); grid on;
xlabel 'y'; ylabel 'u'; title 'u(y) at various locations for Re=100';

figure(3);clf
plot(y0, v(:,(length(x0)-1)/4)*1e5);hold on;
plot(y0, v(:,(length(x0)-1)/2)*1e5,'LineWidth',2);
plot(y0, v(:,(length(x0)-1)/4*3)*1e5,'--');
legend('x=0.25','x=0.5','x=0.75','Location','NorthWest')
set(gca,'FontSize',14); grid on;
xlabel 'y'; ylabel 'v'; title 'v(y) at various locations for Re=100';

figure(4);clf
plot(y0, w(:,(length(x0)-1)/4)*1e5);hold on;
plot(y0, w(:,(length(x0)-1)/2)*1e5,'LineWidth',2);
plot(y0, w(:,(length(x0)-1)/4*3)*1e5,'--');
legend('x=0.25','x=0.5','x=0.75','Location','NorthEast')
set(gca,'FontSize',14); grid on;
xlabel 'y'; ylabel '\omega'; title '\omega(y) at various locations for Re=100';