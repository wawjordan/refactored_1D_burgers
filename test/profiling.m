%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 10;
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('move_shock',64,[-4,4]);

% IN.t0      = 0.05;
% IN.tf      = 1.05;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);
% IN.t0      = 0.1;
% IN.tf      = 0.6;
% Ns  = 2.^(7:14)+1;
% dts = 0.0125./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
IN.t0      = -2;
IN.tf      = 1;
Ns  = 2.^(7:14)+1;
dts = 0.05./(2.^(0:length(Ns)-1));
% Ns  = 2.^(13)+1;
% dts = 0.05./(2.^6);
IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);
IN.Tmethod = 'svd';
M = length(Ns);
IN.U_out = 1;
IN.UE_out = 0;
IN.R_out = 0;
IN.E_out = 0;

error_time = struct();
tStart = tic;
error_time.Efinal = zeros(M,IN.N_IC+1,3);
error_time.Espacetime = zeros(M,IN.N_IC+1,3);
error_time.T = zeros(M,1);
for i = 1:M
    IN.N = Ns(i);
    IN.dt = dts(i);
    [grid,soln,Esoln,EsolnIC,S] = setup_problem_v1(IN);
    S.ETEintegrator.max_newton_iter = 1;
    S.ETEintegratorIC.max_newton_iter = 1;
    tic
%     [Esoln,EsolnIC,soln,OUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
%     [Esoln,EsolnIC,soln,OUT,S,stencil] = ETEIC_solver_alg1_exact_startup(grid,Esoln,EsolnIC,soln,S);
    [Esoln,EsolnIC,soln,OUT,S,stencil] = ete_solver_w_IC_v4(grid,Esoln,EsolnIC,soln,S);
    error_time.T(i)= toc;
    L = length(OUT.t)-1;
    error_time.Efinal(i,:,:) = OUT.ERR.EnormX(end,:,:);
    for j = 1:IN.N_IC+1
        error_time.Espacetime(i,j,1) = norm(OUT.ERR.EnormX(2:L+1,j,1),1)/L;
        error_time.Espacetime(i,j,2) = norm(OUT.ERR.EnormX(2:L+1,j,2),2)/sqrt(L);
        error_time.Espacetime(i,j,3) = norm(OUT.ERR.EnormX(2:L+1,j,3),Inf);
    end
end
tEnd = toc(tStart);
error_time.Ttotal = tEnd-tStart;
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\shock_coalesce_profiles',...
    '\error_profile_ETEIC_single_newton_alg4_svd'];
save(fname,'error_time');