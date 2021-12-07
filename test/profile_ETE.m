%% Grid Refinement, ETE (No iterative correction)
clc; clear; close all;
IN = struct();
IN.order   = 4;
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% % dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('move_shock',64,[-4,4]);
% IN.t0      = 0.05;
% IN.tf      = 1.05;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);
IN.t0      = -2;
IN.tf      = 1;
Ns  = 2.^(7:12)+1;
dts = 0.05./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);
IN.Tmethod = 'svd';
M = length(Ns);
IN.U_out = 1;
IN.UE_out = 0;
IN.R_out = 0;
IN.E_out = 0;

error_time = struct();
tStart = tic;
error_time.Efinal = zeros(M,3);
error_time.Espacetime = zeros(M,3);
error_time.T = zeros(M,1);
for i = 1:M
    IN.N = Ns(i);
    IN.dt = dts(i);
    [grid,soln,Esoln,S] = setup_ETE_problem_v1(IN);
    tic
    [~,~,OUT,~,stencil] = ete_solver(grid,Esoln,soln,S);
    error_time.T(i)= toc;
    L = length(OUT.t)-1;
    error_time.Efinal(i,:) = OUT.ERR.EnormX(end,1,:);
    error_time.Espacetime(i,1) = norm(OUT.ERR.EnormX(2:L+1,1,1),1)/L;
    error_time.Espacetime(i,2) = norm(OUT.ERR.EnormX(2:L+1,1,2),2)/sqrt(L);
    error_time.Espacetime(i,3) = norm(OUT.ERR.EnormX(2:L+1,1,3),Inf);
end
tEnd = toc(tStart);
error_time.Ttotal = tEnd-tStart;
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\error_profile_ETE_svd'];
save(fname,'error_time');