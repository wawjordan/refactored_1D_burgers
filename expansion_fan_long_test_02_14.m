%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_ES02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 2;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_ES02'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_ES04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 4;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_ES04'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_ES06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 6;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_ES06'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_ES08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 8;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_ES08'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_ES10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 10;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_ES10'];
save(fname,'OUT');

pause(10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_GS00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 0;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS00'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_GS02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 2;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS02'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_GS04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 4;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS04'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_GS06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 6;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS06'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_GS08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 8;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS08'];
save(fname,'OUT');

pause(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expansion_fan_GS10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grid Refinement, ETE w/ IC. testing
clc; clear; close all;
IN = struct();
IN.order   = 4;
IN.N_IC    = 10;
for fakeloop=1:1
IN.t0      = 0.1;
IN.tf      = 1.1;
Ns  = 2.^(5:12)+1;
dts = 0.1./(2.^(0:length(Ns)-1));
IN.ex_soln = burgers_exact_soln('#1',64,[-4,4]);

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
    [Esoln,EsolnIC,soln,OUTPUT,S,stencil] = ETEIC_solver_alg1(grid,Esoln,EsolnIC,soln,S);
    OUT.dx(i) = S.dx;
    OUT.x{i}  = grid.x;
    OUT.t{i}  = OUTPUT.t;
    
    OUT.Primal_History(i).E = OUTPUT.PRI.EnormX;
    OUT.Primal_Norms(i,:)   = OUTPUT.PRI.Enorm;

    OUT.ETE_History(i).E    = OUTPUT.ERR.EnormX;
    OUT.ETE_Norms(i,:)      = OUTPUT.ERR.Enorm(1,:);
    OUT.ETE_IC_Norms(i,:)   = OUTPUT.ERR.Enorm(end,:);
end
end
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\results\02_11\',...
    'expansion_fan_GS10'];
save(fname,'OUT');