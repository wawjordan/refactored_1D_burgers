%% Grid Refinement, Primal, Different Nonlinear Solvers
clc; clear; close all;
IN = struct();
% IN.t0      = 0.1;
% IN.tf      = 1.1;
% Ns  = 2.^(6:10)+1;
% dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);
IN.t0      = -2;
IN.tf      = 1;
Ns  = 2.^(6:10)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);
% S.integrator = BDF2_type(grid,soln,S,'method','broyden');
% S.integrator = BDF2_type(grid,soln,S,'method','broyden-lin');
% S.integrator = BDF2_type(grid,soln,S,'method','newton');
% S.integrator = BDF2_type(grid,soln,S,'method','newton_mod');
methods = {'newton','newton_mod','broyden','broyden-lin','broyden-I'};
K = length(methods);
M = length(Ns);
OUT = struct();
OUT.tstart = IN.t0;
OUT.tstop  = IN.tf;
OUT.Nx = Ns;
OUT.dx = zeros(M,1);
OUT.dt = zeros(M,1);
OUT.Local_Error_P = struct();
IN.U_out = 0;
IN.UE_out = 0;
IN.R_out = 0;
IN.E_out = 0;

primal_time = struct();

primal_time.SolverData = struct();
primal_time.Efinal = zeros(M,K,3);
primal_time.Espacetime = zeros(M,K,3);

primal_time.T = zeros(M,K);
for j = 1:M
    for i = 1:K
        IN.N = Ns(j);
        IN.dt = dts(j);
        [grid,soln,S] = setup_primal_problem_v1(IN);
        S.integrator = BDF2_type(grid,soln,S,'method',methods{i});
        S.integrator.max_newton_iter = 100;
        f = @() primal_solver(grid,soln,S);
        tmp = timeit(f,3);
        primal_time.T(j,i) = tmp;
        [soln,OUTPUT,S] = primal_solver(grid,soln,S);
        iters  = OUTPUT.PRI.Rnorm(2:end);
        counts = cellfun(@length,iters);
        primal_time.SolverData(j,i).iters  = iters;
        primal_time.SolverData(j,i).mean   = mean(counts);
        primal_time.SolverData(j,i).median = median(counts);
        primal_time.SolverData(j,i).max    = max(counts);
        primal_time.SolverData(j,i).min    = min(counts);
%         OUT.Local_Error_P(j,i).E = OUTPUT.PRI.E;
%         OUT.Local_Error_P(j,i).u = OUTPUT.PRI.U;
%         OUT.Local_Error_P(j,i).x = grid.x;
%         OUT.Local_Error_P(j,i).t = OUTPUT.t;
        
        L = length(OUTPUT.t)-1;
        primal_time.Efinal(j,i,:) = OUTPUT.PRI.EnormX(end,:);
        primal_time.Espacetime(j,i,1) = norm(OUTPUT.PRI.EnormX(2:L+1,1),1)/L;
        primal_time.Espacetime(j,i,2) = norm(OUTPUT.PRI.EnormX(2:L+1,2),2)/sqrt(L);
        primal_time.Espacetime(j,i,3) = norm(OUTPUT.PRI.EnormX(2:L+1,3),Inf);
    end
end
fname = [...
    'C:\Users\Will\Desktop\MATH5485_Project',...
    '\shock_coalesce_profile'];
save(fname,'primal_time');