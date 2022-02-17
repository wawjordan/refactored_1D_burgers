%% Completed Richardson Extrapolation solver testing
clc; clear; close all;
% N = 257;
% t0 = 0.1;
% tf = 5.1;
% dt = 0.025/2;
% ex_soln = burgers_exact_soln('#1',64,[-4,4]);

% N = 129;
% t0 = -2;
% tf = 4;
% dt = 0.4/8;
% ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);

% N = 65;
% t0 = 0.1;
% tf = 1.1;
% dt = 0.05;
% ex_soln = burgers_exact_soln('#1',64,[-4,4]);
N = 129;
t0 = -5;
tf = 5;
dt = 0.1;
ex_soln = burgers_exact_soln('move_shock',64,[-10,10]);


N2 = (N-1)*2 + 1;
dt2 = dt/2;
rng('default')
x2 = linspace(ex_soln.xmin,ex_soln.xmax,N2);
x2(2:end-1) = x2(2:end-1)+0*0.001*(0.5-rand(1,N2-2));
Fgrid = grid1D(x2);
Fsoln = scalar_soln1D(Fgrid);

x  = x2(1:2:end);
Cgrid = grid1D(x);
Csoln = scalar_soln1D(Cgrid);



% x  = linspace(ex_soln.xmin,ex_soln.xmax,N);
% Cgrid = grid1D(x);
% Csoln = scalar_soln1D(Cgrid);
% 
% x2 = linspace(ex_soln.xmin,ex_soln.xmax,N2);
% Fgrid = grid1D(x2);
% Fsoln = scalar_soln1D(Fgrid);

SC = struct();
SC.ex_soln = ex_soln;
SC.nu = ex_soln.nu;
SC.N = N;
SC.t0 = t0;
SC.tf = tf;
SC.dt = dt;
SC.dx = mean(Cgrid.dx);
SC.U_out_interval = 1;
SC.Uex_out_interval = 1;
SC.R_out_interval = 1;
SC.E_out_interval = 1;
SC.integrator = BDF2_type(Cgrid,Csoln,SC,'method','newton');
SC.L_BC1 = @(~,~) [0,1,0];
SC.L_BC2 = @(~,~) [0,1,0];
SC.R_BC1 = @(u,i,u1) u(i)-u1;
SC.R_BC2 = @(u,i,u2) u(i)-u2;
SC.RHS = @(u) -ss_residual(u,Cgrid.dx,SC.nu,Cgrid.N);
SC.LHS = @(u) jacobian(u,Cgrid.dx,SC.nu,Cgrid.N);

SF = struct();
SF.ex_soln = ex_soln;
SF.nu = ex_soln.nu;
SF.N = N2;
SF.t0 = t0;
SF.tf = tf;
SF.dt = dt2;
SF.dx = mean(Fgrid.dx);
SF.U_out_interval = 1;
SF.Uex_out_interval = 1;
SF.R_out_interval = 1;
SF.E_out_interval = 1;
SF.integrator = BDF2_type(Fgrid,Fsoln,SF,'method','newton');
SF.L_BC1 = @(~,~) [0,1,0];
SF.L_BC2 = @(~,~) [0,1,0];
SF.R_BC1 = @(u,i,u1) u(i)-u1;
SF.R_BC2 = @(u,i,u2) u(i)-u2;
SF.RHS = @(u) -ss_residual(u,Fgrid.dx,SF.nu,Fgrid.N);
SF.LHS = @(u) jacobian(u,Fgrid.dx,SF.nu,Fgrid.N);


[~,~,OUT,SC,SF] = CRE_solver(Cgrid,Fgrid,Csoln,Fsoln,SC,SF);


function val = jacobian(u,dx,nu,N)
val = zeros(N,3);
for ii = 2:N-1
    val(ii,1) = (nu/dx(ii-1)^2) + 0.5*u(ii-1)/dx(ii-1);
    val(ii,2) = -2*(nu/dx(ii)^2);
    val(ii,3) = (nu/dx(ii+1)^2) - 0.5*u(ii+1)/dx(ii+1);
end
end