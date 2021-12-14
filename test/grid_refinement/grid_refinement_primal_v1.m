%% Grid Refinement, Primal
clc; clear; close all;
IN = struct();
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
IN.t0      = 0.1;
IN.tf      = 0.6;
Ns  = 2.^(7:14)+1;
dts = 0.0125./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);

IN.U_out  = 1; IN.UE_out = 1; IN.R_out = 1; IN.E_out = 1;
M = length(Ns);
OUT = struct();
OUT.tstart = IN.t0;
OUT.tstop  = IN.tf;
OUT.Nx = Ns;
OUT.dx = zeros(M,1);
OUT.dt = zeros(M,1);
OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,M,3);
OUT.Final_Enorm2_P = zeros(M,M,3);
IN.U_out = 0;
IN.UE_out = 0;
IN.R_out = 0;
IN.E_out = 0;
for i = 1:M
    for j = 1:M
        IN.N = Ns(i);
        IN.dt = dts(j);
        [grid,soln,S] = setup_primal_problem_v1(IN);
        [soln,OUTPUT,S] = primal_solver(grid,soln,S);
        OUT.dx(i) = S.dx;
        OUT.dt(j) = S.dt;
        L = length(OUTPUT.t)-1;
            Ef(1,1) = norm(OUTPUT.PRI.EnormX(2:L+1,1),1)/L;
            Ef(1,2) = norm(OUTPUT.PRI.EnormX(2:L+1,2),2)/sqrt(L);
            Ef(1,3) = norm(OUTPUT.PRI.EnormX(2:L+1,3),Inf);
        OUT.Local_Error_P(i,j).E = OUTPUT.PRI.E;
        OUT.Local_Error_P(i,j).u = OUTPUT.PRI.U;
        OUT.Local_Error_P(i,j).x = grid.x;
        OUT.Local_Error_P(i,j).t = OUTPUT.t;
        OUT.Error_Norms_P(i,j).E = OUTPUT.PRI.EnormX;
        OUT.Error_Norms_P(i,j).t = OUTPUT.t;
        OUT.Final_Enorm_P(i,j,:) = Ef;
        OUT.Final_Enorm2_P(i,j,:) = OUTPUT.PRI.EnormX(end,:);
    end
end
% fname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\new\',...
%     '\post_processing\expand_primal_1124'];
% save(fname,'OUT');