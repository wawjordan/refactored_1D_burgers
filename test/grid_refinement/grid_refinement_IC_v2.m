%% Grid Refinement, ETE w/ IC. testing (v2)
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 10;
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('move_shock',64,[-4,4]);

IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(4:9)+1;
dts = 0.2./(2.^(0:length(Ns)-1));
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
IN.U_out  = 1; IN.UE_out = 1; IN.R_out = 1; IN.E_out = 1;
M = length(Ns);
OUT = struct();
OUT.tstart = IN.t0;
OUT.tstop  = IN.tf;
OUT.Nx = Ns;
OUT.dx = zeros(M,1);
OUT.dt = dts;
OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,3);
OUT.Final_Enorm2_P = zeros(M,3);
OUT.Local_Error_E = struct();
OUT.Error_Norms_E = struct();
OUT.Final_Enorm_E = cell(M,1);
for i = 1:M
    IN.N = Ns(i);
    IN.dt = dts(i);
    [grid,soln,Esoln,EsolnIC,S] = setup_problem_v1(IN);
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_exact_startup_v2(grid,Esoln,EsolnIC,soln,S);
%     L = length(OUTPUT.t)-1;
%     P_L1  = sum(abs(([OUT.PRI.U{:}]-[OUT.PRI.Uex{:}])))/grid.N;
%     E1_L1 = sum(abs([OUT.ERR.E{:,1}]-([OUT.PRI.U{:}]-[OUT.PRI.Uex{:}])))/grid.N;
%     E2_L1 = sum(abs([OUT.ERR.E{:,2}]-([OUT.PRI.U{:}]-[OUT.PRI.Uex{:}])))/grid.N;
%     for j = 1:IN.N_IC+1
%         Ef(j,1) = norm(OUTPUT.ERR.EnormX(2:L+1,j,1),1)/L;
%         Ef(j,2) = norm(OUTPUT.ERR.EnormX(2:L+1,j,2),2)/sqrt(L);
%         Ef(j,3) = norm(OUTPUT.ERR.EnormX(2:L+1,j,3),Inf);
%     end
    OUT.dx(i) = S.dx;
    OUT.Local_Error_P(i).E = cellfun(@(a,b) a-b,OUTPUT.PRI.U,OUTPUT.PRI.Uex,'UniformOutput',false);
    OUT.Local_Error_P(i).u = OUTPUT.PRI.U;
    OUT.Local_Error_P(i).x = grid.x;
    OUT.Local_Error_P(i).t = OUTPUT.t;
    OUT.Error_Norms_P(i).E(:,1) = sum(abs(([OUT.Local_Error_P(i).E{:}])))/grid.N;
    OUT.Error_Norms_P(i).E(:,2) = sqrt(sum(([OUT.Local_Error_P(i).E{:}]).^2)/grid.N);
    OUT.Error_Norms_P(i).E(:,3) = max(abs(([OUT.Local_Error_P(i).E{:}])));
    OUT.Error_Norms_P(i).t = OUTPUT.t;
    
    OUT.Local_Error_E(i).EE(:,1) = cellfun(@(a,b) a-b,OUTPUT.ERR.E(:,1),OUT.Local_Error_P(i).E,'UniformOutput',false);
    OUT.Local_Error_E(i).EE(:,2) = cellfun(@(a,b) a-b,OUTPUT.ERR.E(:,2),OUT.Local_Error_P(i).E,'UniformOutput',false);
    OUT.Local_Error_E(i).E = OUTPUT.ERR.E;
    OUT.Local_Error_E(i).x = grid.x;
    OUT.Local_Error_E(i).t = OUTPUT.t;
    OUT.Error_Norms_E(i).E(:,1,1) = sum(abs(([OUT.Local_Error_E(i).EE{:,1}])))/grid.N;
    OUT.Error_Norms_E(i).E(:,1,2) = sqrt(sum(([OUT.Local_Error_E(i).EE{:,1}]).^2)/grid.N);
    OUT.Error_Norms_E(i).E(:,1,3) = max(abs(([OUT.Local_Error_E(i).EE{:,1}])));
    OUT.Error_Norms_E(i).E(:,2,1) = sum(abs(([OUT.Local_Error_E(i).EE{:,2}])))/grid.N;
    OUT.Error_Norms_E(i).E(:,2,2) = sqrt(sum(([OUT.Local_Error_E(i).EE{:,2}]).^2)/grid.N);
    OUT.Error_Norms_E(i).E(:,2,3) = max(abs(([OUT.Local_Error_E(i).EE{:,2}])));
    OUT.Error_Norms_E(i).t = OUTPUT.t;
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\results\12_16\',...
    'expand_exact_start'];
% save(fname,'OUT');
save(fname, 'OUT', '-v7.3')
%%
hold on;
for i = 1:6
    plot(OUT.Error_Norms_P(i).t,OUT.Error_Norms_P(i).E(:,1),'k');
    plot(OUT.Error_Norms_E(i).t,OUT.Error_Norms_E(i).E(:,:,1));
end
hold off;
set(gca,'yscale','log')