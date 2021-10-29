%% Testing New BDF2 type
clc; clear; close all;
x = linspace(-4,4,2049);
grid = grid1D(x);
soln = scalar_soln1D(grid);
S.dt = 0.05;
S.t0 = -2;
S.nu = 0.0001;
BD = BDF2_type(grid,soln,S);

LHS_BC1 = @(u,i) [0,1,0];
LHS_BC2 = @(u,i) [0,1,0];
RHS_BC1 = @(u,i) 0*(rand(1)-0.5);%u(grid.N);
RHS_BC2 = @(u,i) 0;%u(1);
% RHS_BC1 = @(u,i) u(i);
% RHS_BC2 = @(u,i) -u(i);

res = @(u) -ss_residual(u,grid.dx,S.nu,grid.N);
jac = @(u) jacobian(u,grid.dx,S.nu,grid.N);

u_old = soln.U;
% BD.um1 = sin(grid.x).^2;
% BD.um1 = rand(grid.N,1);
% BD.um1(1) = 0; BD.um1(grid.N) = 0;
% BD.um2 = BD.um1;
% BD.um2 = rand(grid.N,1);
% BD.um2(1) = 0; BD.um2(grid.N) = 0;
% u_old = rand(grid.N,1);
% u_old(1) = 0; u_old(grid.N) = 0;
% BD.um1 = sin(grid.x+S.dt);
% BD.um2 = sin(grid.x+2*S.dt);
% BD.max_newton_iter = 1;
% u_old = sin(grid.x);

figure(1);
% hold on;
for i = 1:10000
fprintf('%d\n',i)
n = i;
RHS_BC1 = @(u,i) u(i)-4*(sin(n/20)*cos(n/100));%u(grid.N);
RHS_BC2 = @(u,i) 0;%u(i)-3.999*(sin(n/100)*cos(n/20));%u(grid.N);
% RHS_BC1 = @(u,i) u(i)-4*(sin(n/20)*cos(n/100));%u(grid.N);
% RHS_BC2 = @(u,i) u(i)-3.999*(sin(n/100)*cos(n/20));%u(grid.N);
[u_new,R,S,BD] = BD.step(u_old,S,res,jac,LHS_BC1,LHS_BC2,RHS_BC1,RHS_BC2);
plot(grid.x,u_new,'linewidth',0.0001)
ylim([-4,4])
drawnow
u_old = u_new;
end
axis off;
% exportgraphics(gcf,'dots.png','Resolution',1200)

function val = jacobian(u,dx,nu,N)
val = zeros(N,3);
for ii = 2:N-1
    val(ii,1) = (nu/dx(ii-1)^2) + 0.5*u(ii-1)/dx(ii-1);
    val(ii,2) = -2*(nu/dx(ii)^2);
    val(ii,3) = (nu/dx(ii+1)^2) - 0.5*u(ii+1)/dx(ii+1);
end
end
