function [grid,soln,Esoln,EsolnIC,S] = setup_problem_v1(IN)

x = linspace(IN.ex_soln.xmin,IN.ex_soln.xmax,IN.N);
grid = grid1D(x);
soln = scalar_soln1D(grid);    % primal solution
Esoln = scalar_soln1D(grid);   % ETE solution
EsolnIC = scalar_soln1D(grid); % ETE IC solution

S = struct();
S.ex_soln = IN.ex_soln;
S.nu = IN.ex_soln.nu;
S.N  = IN.N;
S.t0 = IN.t0;
S.tf = IN.tf;
S.dt = IN.dt;
S.dx = mean(grid.dx);
S.U_out_interval   = IN.U_out;
S.Uex_out_interval = IN.UE_out;
S.R_out_interval   = IN.R_out;
S.E_out_interval   = IN.E_out;
S.Niters = IN.N_IC;
S.out_iters = 1:S.Niters;
S.stencil_size = IN.order+1+mod(IN.order,2);
S.integrator      = BDF2_type(grid,soln,S);
S.ETEintegrator   = BDF2_type(grid,Esoln,S);
S.ETEintegratorIC = BDF2_type(grid,EsolnIC,S);
S.LS_S = spatial_reconstruction(grid,S,IN.order);
S.LS_T = temporal_reconstruction(grid,S,IN.order,'method',IN.Tmethod);

S.L_BC1 = @(~,~) [0,1,0];
S.L_BC2 = @(~,~) [0,1,0];
S.R_BC1 = @(u,i,u1) u(i)-u1;
S.R_BC2 = @(u,i,u2) u(i)-u2;

S.RHS = @(u) -ss_residual(u,grid.dx,S.nu,grid.N);
S.LHS = @(u) jacobian(u,grid.dx,S.nu,grid.N);

S.ETE_RHS = @(u,e,Ru,TE) ETE_residual(u,e,Ru,TE,grid.dx,S.nu,grid.N);
S.ETE_LHS = @(u,e) ETE_jacobian(u,e,grid.dx,S.nu,grid.N);

end
%% Supplementary Functions
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