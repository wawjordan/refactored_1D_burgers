%% spatial_reconstruction class testing
clc; clear; close all;

x = linspace(-4,4,129);
grid = grid1D(x);

order = 4;
S.stencil_size = 5;
LS = spatial_reconstruction(grid,S,order);

S.stencil_size = 9;
LS2 = spatial_reconstruction(grid,S,order);
% U   =  sin(x');
% dU1 =  cos(x');
% dU2 = -sin(x');

U   =  sin(x').^2;
dU1 =  2*sin(x').*cos(x');
dU2 = -2*sin(x').^2+2*cos(x').^2; 

[u,du1,du2] = LS.eval(U);
[v,dv1,dv2] = LS2.eval(U);
% hold on;
% plot(x,U,'k')
% plot(x,u,'r')
% hold off

hold on;
plot(x,abs(U-u),'b')
plot(x,abs(dU1-du1),'k')
plot(x,abs(dU2-du2),'r')

plot(x,abs(U-v),'b--')
plot(x,abs(dU1-dv1),'k--')
plot(x,abs(dU2-dv2),'r--')
hold off
set(gca,'yscale','log')