%% ETE solver testing
clc; clear; close all;

order = 4;

N = 8193;
t0 = -2;
tf = 2;
dt = 0.4/256;
ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);

% N = 129;
% t0 = 0.1;
% tf = 0.6;
% dt = 0.025;
% ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
x = linspace(ex_soln.xmin,ex_soln.xmax,N);
grid = grid1D(x);
soln = scalar_soln1D(grid);  % primal solution
Esoln = scalar_soln1D(grid); % ETE solution

S = struct();
S.ex_soln = ex_soln;
S.nu = ex_soln.nu;
S.N = N;
S.t0 = t0;
S.tf = tf;
S.dt = dt;
S.dx = mean(grid.dx);
S.U_out_interval = 1;
S.Uex_out_interval = 1;
S.R_out_interval = 1;
S.E_out_interval = 1;
S.out_iters = [];
S.Niters = 1;
% S.stencil_size = 11;
S.stencil_size = order+1+mod(order,2);
S.integrator    = BDF2_type(grid,soln,S);
S.ETEintegrator = BDF2_type(grid,Esoln,S);
S.LS_S = spatial_reconstruction(grid,S,order);
S.LS_T = temporal_reconstruction(grid,S,order,'method','svd2');

S.L_BC1 = @(~,~) [0,1,0];
S.L_BC2 = @(~,~) [0,1,0];
S.R_BC1 = @(u,i,u1) u(i)-u1;
S.R_BC2 = @(u,i,u2) u(i)-u2;

S.RHS = @(u) -ss_residual(u,grid.dx,S.nu,grid.N);
S.LHS = @(u) jacobian(u,grid.dx,S.nu,grid.N);

S.ETE_RHS = @(u,e,Ru,TE) ETE_residual(u,e,Ru,TE,grid.dx,S.nu,grid.N);
S.ETE_LHS = @(u,e) ETE_jacobian(u,e,grid.dx,S.nu,grid.N);


% [Esoln,soln,OUT,S,stencil] = ete_solver(grid,Esoln,soln,S);
[~,~,OUT,~,stencil] = ete_solver(grid,Esoln,soln,S);
% reset variables
[~,~,OUT2,~,stencil2] = ete_solver_w_lag(grid,Esoln,soln,S);
% [soln,OUT,S] = primal_solver(grid,soln,S);
% [soln,OUT,S,stencil] = primal_solver_w_stencil(grid,soln,S);

function val = jacobian(u,dx,nu,N)
val = zeros(N,3);
for ii = 2:N-1
    val(ii,1) = (nu/dx(ii-1)^2) + 0.5*u(ii-1)/dx(ii-1);
    val(ii,2) = -2*(nu/dx(ii)^2);
    val(ii,3) = (nu/dx(ii+1)^2) - 0.5*u(ii+1)/dx(ii+1);
end
end

function val = ETE_jacobian(u,e,dx,nu,N)
val = zeros(N,3);
ue = u-e;
for ii = 2:N-1
    val(ii,1) = (nu/dx(ii-1)^2) + 0.5*ue(ii-1)/dx(ii-1);
    val(ii,2) = -2*(nu/dx(ii)^2);
    val(ii,3) = (nu/dx(ii+1)^2) - 0.5*ue(ii+1)/dx(ii+1);
end
end

function val = ETE_residual(u,e,Ru,TE,dx,nu,N)
val = ss_residual(u-e,dx,nu,N)-TE-Ru;
end