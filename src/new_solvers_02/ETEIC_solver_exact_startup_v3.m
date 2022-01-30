function [Esoln,EsolnIC,soln,OUT,S,stencil] = ETEIC_solver_exact_startup_v3(grid,Esoln,EsolnIC,soln,S)
OUT = struct();
OUT.dx          = S.dx(1);  % equally spaced
OUT.dt          = S.dt;
OUT.t0          = S.t0;
OUT.tf          = S.tf;
OUT.t           = (S.t0:S.dt:S.tf)';
OUT.Nt          = length(OUT.t);
S.max_steps     = OUT.Nt;

OUT.PRI        = struct();
OUT.PRI.t      =  nan(S.max_steps,1);
OUT.PRI.U      = cell(S.max_steps,1);
OUT.PRI.Uex    = cell(S.max_steps,1);
OUT.PRI.R      = cell(S.max_steps,1);
OUT.PRI.Rt     =  nan(S.max_steps,1);

OUT.ERR        = struct();
OUT.ERR.t      =  nan(S.max_steps,1);
OUT.ERR.Rt     =  nan(S.max_steps,1);
if S.Niters > 0
    OUT.ERR.E      = cell(S.max_steps,2);
    OUT.ERR.EE     = cell(S.max_steps,2);
    OUT.ERR.R      = cell(S.max_steps,2);
else
    OUT.ERR.E      = cell(S.max_steps,1);
    OUT.ERR.EE     = cell(S.max_steps,1);
    OUT.ERR.R      = cell(S.max_steps,1);
end

%% Preprocessing steps
S = dt_options(S);
fprintf(S.string_fmt,OUT.t(1),0,S.max_steps-1);
S.count = 1;
if S.Niters > 0
    stencil = soln_stencil( S.stencil_size, S.N, 2 );
else
    stencil = soln_stencil( S.stencil_size, S.N, 1 );
end
tmp = stencil.U(:,:,1);
for i = 1:S.stencil_size
    time = OUT.t(i);
    Uex = S.ex_soln.eval(grid.x,time);
    stencil = stencil.push(Uex,time);
end
for i = S.stencil_size-1:-1:0
    time = OUT.t(1)-i*S.dt;
    [tmp(:,i+1),~] = S.LS_T.eval(stencil,time,1);
end
for i = S.stencil_size-1:-1:0
    time = OUT.t(1)-i*S.dt;
    stencil = stencil.push(tmp(:,i+1),time);
end
% for i = S.stencil_size-1:-1:0
%     time = OUT.t(1)-i*S.dt;
%     Uex = S.ex_soln.eval(grid.x,time);
%     stencil = stencil.push(Uex,time);
% end

S.isBDF2 = (isa(S.integrator,'BDF2_type'));
if (S.isBDF2)
    S.integrator.um2 = stencil.U(:,S.stencil_size-1);
    S.integrator.um1 = stencil.U(:,S.stencil_size);
end
soln.U = stencil.U(:,S.stencil_size);
OUT.PRI.U{1}   = soln.U;
OUT.PRI.Uex{1} = soln.U;
OUT.PRI.t(1)   = OUT.t(1);
OUT.PRI.R{1}   = 0;
OUT.PRI.Rt(1)  = OUT.t(1);

OUT.ERR.E{1,1}   = zeros(grid.N,1);
OUT.ERR.t(1)   = OUT.t(1);
OUT.ERR.R{1,1}   = 0;
OUT.ERR.Rt(1)  = OUT.t(1);

if S.Niters > 0
    OUT.ERR.E{1,2}   = zeros(grid.N,1);
    OUT.ERR.R{1,2}   = 0;
end

S.count = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advance Solution
if S.Niters > 0
    for i = S.count:S.max_steps
        time = OUT.t(i);
        fprintf(S.string_fmt,time,i-1,S.max_steps-1);
        [soln,stencil,resnorm,S] = step_primal(grid,soln,stencil,S,time);
        OUT = output_solution(OUT,S,soln.U,resnorm,i);
        [Esoln,stencil,resnorm,S] = step_ete(grid,Esoln,stencil,S);
        OUT = output_ete_solution(OUT,S,Esoln.U,resnorm,i);
        [EsolnIC,stencil,resnorm,S] = step_iterates(grid,EsolnIC,stencil,S);
        OUT = output_ete_iter_solution(OUT,S,EsolnIC.U,resnorm,i);
    end
else
    for i = S.count:S.max_steps
        time = OUT.t(i);
        fprintf(S.string_fmt,time,i-1,S.max_steps-1);
        [soln,stencil,resnorm,S] = step_primal(grid,soln,stencil,S,time);
        OUT = output_solution(OUT,S,soln.U,resnorm,i);
        [Esoln,stencil,resnorm,S] = step_ete(grid,Esoln,stencil,S);
        OUT = output_ete_solution(OUT,S,Esoln.U,resnorm,i);
    end
end
for i = S.count:S.max_steps
    time = OUT.t(i);
    if mod(i-1,S.U_out_interval) == 0 || i == S.max_steps
        OUT.PRI.Uex{i} = S.ex_soln.eval(grid.x,time);
    end
end
OUT = cleanup(OUT);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supplementary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = dt_options(S)
L1 = length(char(regexp(string(S.dt),'(?<=\.)\d*','match')));
L2 = length(char(string(S.max_steps)));
S.string_fmt = ...
    sprintf('t = %%#%d.%df (timestep: %%%dd/%%d)\\n',12,L1,L2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,stencil,resnorm,S] = step_primal(grid,soln,stencil,S,time)
    N = grid.N;
    Uex1 = S.ex_soln.eval(grid.x(1),time);
    UexN = S.ex_soln.eval(grid.x(N),time);
    L_BC1 = S.L_BC1;                         % LHS BC, left boundary
    L_BC2 = S.L_BC2;                         % LHS BC, right boundary
    R_BC1 = @(u,i) S.R_BC1(u,i,Uex1);        % RHS BC, left boundary
    R_BC2 = @(u,i) S.R_BC2(u,i,UexN);        % RHS BC, right boundary
    [soln.U,resnorm,S,S.integrator] = S.integrator.step(...
           soln.U, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
     stencil = stencil.push(soln.U,time); % update stencil
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,resnorm,S] = step_ete(grid,Esoln,stencil,S)
    time = stencil.t(stencil.queue_length);
    U = stencil.U(:,stencil.queue_length,1);
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,i) S.R_BC1(u,i,0);
    R_BC2 = @(u,i) S.R_BC2(u,i,0);
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    [u,dudt] = S.LS_T.eval(stencil,time,1);
    TE = ss_residual_cont(S,S.LS_S,u);
    TE = TE + dudt;
    RHS = @(e) S.ETE_RHS(U,e,RU,TE);
    LHS = @(e) S.ETE_LHS(U,e);
    [Esoln.U,resnorm,S,S.ETEintegrator] = S.ETEintegrator.step(...
              Esoln.U, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    if S.Niters > 0
        stencil.U(:,stencil.queue_length,2) = U - Esoln.U;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,resnorm,S] = step_iterates(grid,Esoln,stencil,S)
i = stencil.queue_length;
time = stencil.t(i);
for j = 1:S.Niters   % iterative correction steps
    S.ETEintegratorIC = S.ETEintegratorIC.integrator_reset();
    Esoln.U = zeros(grid.N,1);
    U = stencil.U(:,i,2); % corrected solution stencil
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,n) S.R_BC1(u,n,0); % dirichlet BC
    R_BC2 = @(u,n) S.R_BC2(u,n,0); % dirichlet BC
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    [u,dudt] = S.LS_T.eval(stencil,time,2);
    TE = ss_residual_cont(S,S.LS_S,u);
    TE = TE + dudt;
    RHS = @(e) S.ETE_RHS(U,e,RU,TE);
    LHS = @(e) S.ETE_LHS(U,e);
    [Esoln.U,resnorm,S,S.ETEintegratorIC] = S.ETEintegratorIC.step(...
        Esoln.U, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    UE = U - Esoln.U; % corrected solution
    stencil.U(:,i,2) = UE; % update stencil with corrected solution
end
Esoln.U = stencil.U(:,i,1) - stencil.U(:,i,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_solution(OUT,S,U,resnorm,count)
    if mod(count-1,S.R_out_interval) == 0 || count == S.max_steps
        OUT.PRI.R{count}  = resnorm;
        OUT.PRI.Rt(count) = OUT.t(count);
    end
    if mod(count-1,S.U_out_interval) == 0 || count == S.max_steps
        OUT.PRI.U{count} = U;
        OUT.PRI.t(count) = OUT.t(count);
    end
end
function OUT = output_ete_solution(OUT,S,E,resnorm,count)
    if mod(count-1,S.R_out_interval) == 0 || count == S.max_steps
        OUT.ERR.R{count,1}  = resnorm;
        OUT.ERR.Rt(count) = OUT.t(count);
    end
    if mod(count-1,S.E_out_interval) == 0 || count == S.max_steps
        OUT.ERR.E{count,1} = E;
        OUT.ERR.t(count) = OUT.t(count);
    end
end
function OUT = output_ete_iter_solution(OUT,S,E,resnorm,count)
    if mod(count-1,S.R_out_interval) == 0 || count == S.max_steps
        OUT.ERR.R{count,2}  = resnorm;
    end
    if mod(count-1,S.E_out_interval) == 0 || count == S.max_steps
        OUT.ERR.E{count,2} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = cleanup(OUT)
OUT.PRI.t   = OUT.PRI.t(~isnan(OUT.PRI.t));
OUT.PRI.U   = OUT.PRI.U(  ~cellfun('isempty',OUT.PRI.U  ));
OUT.PRI.Uex = OUT.PRI.Uex(~cellfun('isempty',OUT.PRI.Uex));
OUT.PRI.Rt  = OUT.PRI.Rt(~isnan(OUT.PRI.Rt));
OUT.PRI.R   = OUT.PRI.R(  ~cellfun('isempty',OUT.PRI.R  ));

OUT.ERR.t   = OUT.ERR.t(~isnan(OUT.ERR.t));
OUT.ERR.E   = OUT.ERR.E(~cellfun('isempty',OUT.ERR.E(:,1)),:);
OUT.ERR.Rt  = OUT.ERR.Rt(~isnan(OUT.ERR.Rt));
OUT.ERR.R   = OUT.ERR.R(~cellfun('isempty',OUT.ERR.R(:,1)),:);
end