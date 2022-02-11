%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 2;

% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('move_shock',64,[-4,4]);

IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(7:10)+1;
dts = 0.025./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

% IN.t0      = 0.1;
% IN.tf      = 0.6;
% Ns  = 2.^(7:14)+1;
% dts = 0.0125./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);

% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:10)+1;
% dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);


IN.Tmethod = 'svd';
IN.U_out  = 0;
IN.UE_out = 0;
IN.R_out  = 0;
IN.E_out  = 0;

M = length(Ns);
OUT = struct();
OUT.tstart = IN.t0;
OUT.tstop  = IN.tf;
OUT.Nx = Ns;
OUT.dx = zeros(M,1);
OUT.x  =  cell(M,1);
OUT.t  =  cell(M,1);
OUT.dt = dts;
OUT.Primal_History = struct();
OUT.ETE_History = struct();
OUT.Primal_Norms = zeros(M,3);
OUT.ETE_Norms = zeros(M,3);
OUT.ETE_IC_Norms = zeros(M,3);
for i = 1:M
    IN.N = Ns(i);
    IN.dt = dts(i);
    [grid,soln,Esoln,EsolnIC,S] = setup_problem_v1(IN);
%     S.integrator.max_newton_iter = 1;
%     [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1_exact_startup(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
P1 = zeros(M-1,3);
%%
tmp1 = OUT.Primal_Norms(1:end-1,:)./OUT.Primal_Norms(2:end,:);
tmp2 =    OUT.ETE_Norms(1:end-1,:)./   OUT.ETE_Norms(2:end,:);
tmp3 = OUT.ETE_IC_Norms(1:end-1,:)./OUT.ETE_IC_Norms(2:end,:);
PP = log(tmp1)/log(2);
EP = log(tmp2)/log(2);
IP = log(tmp3)/log(2);
% fname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\new\',...
%     '\post_processing\expand_exact_start_alg1_allsave'];
% save(fname,'OUT');
% save(fname, 'OUT', '-v7.3')