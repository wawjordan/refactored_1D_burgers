%% Primal solver testing
clc; clear; close all;
% N = 257;
% t0 = -2;
% tf = 2;
% dt = 0.4/8;
% ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);
N = 129;
t0 = 0.1;
tf = 0.6;
dt = 0.025;
ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
x = linspace(ex_soln.xmin,ex_soln.xmax,N);
grid = grid1D(x);
soln = scalar_soln1D(grid);

S = struct();
S.ex_soln = ex_soln;
S.nu = ex_soln.nu;
S.N = N;
S.t0 = t0;
S.tf = tf;
S.dt = dt;
S.dx = mean(grid.dx);
S.U_out_interval = 1;
S.Uex_out_interval = 4;
S.R_out_interval = 0;
S.E_out_interval = 2;
S.out_iters = [];
S.Niters = 1;
S.stencil_size = 5;
S.integrator = BDF2_type(grid,soln,S);

S.L_BC1 = @(~,~) [0,1,0];
S.L_BC2 = @(~,~) [0,1,0];
S.R_BC1 = @(u,i,u1) u(i)-u1;
S.R_BC2 = @(u,i,u2) u(i)-u2;
S.RHS = @(u) -ss_residual(u,grid.dx,S.nu,grid.N);
S.LHS = @(u) jacobian(u,grid.dx,S.nu,grid.N);

% [soln,OUT,S] = primal_solver(grid,soln,S);
[soln,OUT,S,stencil] = primal_solver_w_stencil(grid,soln,S);

function val = jacobian(u,dx,nu,N)
val = zeros(N,3);
for ii = 2:N-1
    val(ii,1) = (nu/dx(ii-1)^2) + 0.5*u(ii-1)/dx(ii-1);
    val(ii,2) = -2*(nu/dx(ii)^2);
    val(ii,3) = (nu/dx(ii+1)^2) - 0.5*u(ii+1)/dx(ii+1);
end
end