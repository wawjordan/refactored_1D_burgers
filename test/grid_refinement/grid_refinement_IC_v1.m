%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 10;
% IN.t0      = 0.1;
% IN.tf      = 0.6;
% Ns  = 2.^(7:9);
% dts = 0.025./(2.^(0:3));
% IN.ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
IN.t0      = -2;
IN.tf      = 2;
Ns  = 2.^(5:11)+1;
dts = 0.4./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);


IN.Tmethod = 'default';
IN.U_out  = 1; IN.UE_out = 1; IN.R_out = 0; IN.E_out = 0;
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
intervals = uint16((OUT.tstop-OUT.tstart)./OUT.dt);
out_interval = 1;
for i = 1:40
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;
for i = 1:M
    
    IN.N = Ns(i);
    IN.dt = dts(i);
    [grid,soln,Esoln,EsolnIC,S] = setup_problem_v1(IN);
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ete_solver_w_IC(grid,Esoln,EsolnIC,soln,S);
    L = length(OUTPUT.t)-1;
    for j = 1:IN.N_IC+1
        Ef(j,1) = norm(OUTPUT.ERR.EnormX(2:L,j,1),1)/L;
        Ef(j,2) = norm(OUTPUT.ERR.EnormX(2:L,j,1),2)/sqrt(L);
        Ef(j,3) = norm(OUTPUT.ERR.EnormX(2:L,j,1),Inf);
    end
    OUT.dx(i) = S.dx;
    OUT.Local_Error_P(i).E = OUTPUT.PRI.E;
    OUT.Local_Error_P(i).u = OUTPUT.PRI.U;
    OUT.Local_Error_P(i).x = grid.x;
    OUT.Local_Error_P(i).t = OUTPUT.t;
    OUT.Error_Norms_P(i).E = OUTPUT.PRI.EnormX;
    OUT.Error_Norms_P(i).t = OUTPUT.t;
    OUT.Final_Enorm_P(i,:) = OUTPUT.PRI.EnormX(end,:);
    
    OUT.Local_Error_E(i).Ee = OUTPUT.ERR.EE;
    OUT.Local_Error_E(i).e = OUTPUT.ERR.E;
    OUT.Local_Error_E(i).x = grid.x;
    OUT.Local_Error_E(i).t = OUTPUT.t;
    OUT.Error_Norms_E(i).E = OUTPUT.ERR.EnormX;
    OUT.Error_Norms_E(i).t = OUTPUT.t;
    OUT.Final_Enorm_E{i} = Ef;
%     OUT.Final_Enorm_E{i} = OUTPUT.ERR.EnormX(end,:,:);
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\post_processing\test3'];
save(fname,'OUT');