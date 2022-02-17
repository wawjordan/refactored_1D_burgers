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

% IN.t0      = 0.1;
% IN.tf      = 1.1;
% Ns  = 2.^(7:8)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);


IN.t0      = 0.2;
IN.tf      = 1.2;
Ns  = 2.^(5:12)+1;
dts = 0.05./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);



% IN.t0      = 0.1;
% IN.tf      = 0.6;
% Ns  = 2.^(7:14)+1;
% dts = 0.0125./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:12)+1;
% dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);


IN.Tmethod = 'svd';
% IN.U_out  = 1; IN.UE_out = 1; IN.R_out = 1; IN.E_out = 1;
M = length(Ns);
OUT = struct();
OUT.tstart = IN.t0;
OUT.tstop  = IN.tf;
OUT.Nx = Ns;
OUT.dx = zeros(M,1);
OUT.x  =  cell(M,1);
OUT.t  =  cell(M,1);
OUT.dt = dts;
OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,3);
OUT.Local_Error_E = struct();
OUT.Error_Norms_E = struct();
OUT.Final_Enorm_E = zeros(M,IN.N_IC+1,3);
intervals = uint16((OUT.tstop-OUT.tstart)./OUT.dt);
out_interval = 1;
for i = 1:40
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;
IN.U_out = 1;
IN.UE_out = 1;
IN.R_out = 1;
IN.E_out = 1;
for i = 1:M
    IN.N = Ns(i);
    IN.dt = dts(i);
    [grid,soln,Esoln,EsolnIC,S] = setup_problem_v1(IN);
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
%     [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ete_solver_w_IC_v4(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    OUT.Local_Error_P(i).E = OUTPUT.PRI.E;
    OUT.Local_Error_P(i).u = OUTPUT.PRI.U;
    OUT.Error_Norms_P(i).E = OUTPUT.PRI.EnormX;
    OUT.Final_Enorm_P(i,:) = OUTPUT.PRI.Enorm;
    
    OUT.Local_Error_E(i).Ee = OUTPUT.ERR.EE;
    OUT.Local_Error_E(i).e = OUTPUT.ERR.E;
    OUT.Error_Norms_E(i).E = OUTPUT.ERR.EnormX;
    OUT.Final_Enorm_E(i,:,:) = OUTPUT.ERR.Enorm;
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS10_full'];
save(fname,'OUT');
save(fname, 'OUT', '-v7.3')
% fname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\new\',...
%     '\post_processing\expand_exact_start_alg1_allsave'];
% save(fname,'OUT');
% save(fname, 'OUT', '-v7.3')