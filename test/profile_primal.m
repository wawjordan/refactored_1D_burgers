%% Grid Refinement, Primal
clc; clear; close all;
IN = struct();
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.025./(2.^(0:length(Ns)-1));
% % dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('move_shock',64,[-4,4]);
IN.t0      = 0.05;
IN.tf      = 1.05;
Ns  = 2.^(7:9)+1;
dts = 0.025./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);
% IN.t0      = -2;
% IN.tf      = 1;
% Ns  = 2.^(7:14)+1;
% dts = 0.05./(2.^(0:length(Ns)-1));
% IN.ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);

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

primal_time.Erich = zeros(M,3);
primal_time.Efinal = zeros(M,3);
primal_time.Espacetime = zeros(M,3);

% tStart = tic;
primal_time.T = zeros(M,1);
for j = 1:M
    IN.N = Ns(j);
    IN.dt = dts(j);
    [grid,soln,S] = setup_primal_problem_v1(IN);
    if j == 1 || j == 2
        f = @() primal_solver(grid,soln,S);
        tmp = timeit(f,3);
        primal_out.T(j) = tmp;
        [soln,OUTPUT,S] = primal_solver(grid,soln,S);
    else
        tic
        [soln,OUTPUT,S] = primal_solver(grid,soln,S);
        primal_time.T(j) = toc;
    end
    OUT.Local_Error_P(j).E = OUTPUT.PRI.E;
    OUT.Local_Error_P(j).u = OUTPUT.PRI.U;
    OUT.Local_Error_P(j).x = grid.x;
    OUT.Local_Error_P(j).t = OUTPUT.t;
    
    L = length(OUTPUT.t)-1;
    primal_time.Efinal(j,:) = OUTPUT.PRI.EnormX(end,:);
    primal_time.Espacetime(j,1) = norm(OUTPUT.PRI.EnormX(2:L+1,1),1)/L;
    primal_time.Espacetime(j,2) = norm(OUTPUT.PRI.EnormX(2:L+1,2),2)/sqrt(L);
    primal_time.Espacetime(j,3) = norm(OUTPUT.PRI.EnormX(2:L+1,3),Inf);
end

for j = 2:M
    E_rich = -(OUT.Local_Error_P(j).u{end}(1:2:end)-OUT.Local_Error_P(j-1).u{end})/(2^2-1);
    EE_rich = E_rich-OUT.Local_Error_P(j).E{end}(1:2:end);
    primal_time.Erich(j,1) = sum(abs(EE_rich))/Ns(j-1);
    primal_time.Erich(j,2) = sqrt(sum(EE_rich.^2)/Ns(j-1));
    primal_time.Erich(j,3) = max(abs(EE_rich));
end
% tEnd = toc(tStart);
% primal_time.Ttotal = tEnd-tStart;
% fname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\new\expansion_fan_profiles',...
%     '\primal_profile'];
% save(fname,'primal_time');