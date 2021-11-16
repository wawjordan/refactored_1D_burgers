%% temporal_reconstruction class testing
clc; clear; close all;

U   = @(x,t) sin(x(:)-t);
dU1 = @(x,t) -cos(x(:)-t);


x = linspace(-4,4,129);
grid = grid1D(x);

order = 4;
S.stencil_size = 5;
S.dt = 0.1;
t = 0:S.dt:S.dt*(S.stencil_size-1);
sten = soln_stencil( S.stencil_size, length(x), 2, 'max_length', 10);
sten2 = soln_stencil( S.stencil_size, length(x), 2, 'max_length', 10);
for i = 1:S.stencil_size
    u_new = U(x,t(i));
    sten = sten.push(u_new,t(i));
    du_new = dU1(x,t(i));
    sten2 = sten.push(du_new,t(i));
end
% plot(x,sten.U(:,:,1))





T1 = temporal_reconstruction(grid,S,order,'method','default');
T2 = temporal_reconstruction(grid,S,order,'method','sgolay');

time = t(end);
[u,du1] = T1.eval(sten,time,1);
[v,dv1] = T2.eval(sten,time,1);

Uf = sten.U(:,end,end);
dUf = sten2.U(:,end,end);
% hold on;
% plot(x,U,'k')
% plot(x,u,'r')
% hold off
% hold on;
% plot(v(1:5))
% plot(U(1:5))
hold on;
% plot(x,abs(Uf-u),'b')
% plot(x,abs(dUf-du1),'r--')

% plot(x,abs(Uf(:,end)-v),'b--')
% plot(x,abs(dUf-dv1),'k--')
plot(x,du1-dv1,'k')
hold off
% set(gca,'yscale','log')