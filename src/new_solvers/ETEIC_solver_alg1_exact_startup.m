function [Esoln,EsolnIC,soln,OUT,S,stencil] = ETEIC_solver_alg1_exact_startup(grid,Esoln,EsolnIC,soln,S)
%%
for fakeloop = 1:1
S.max_steps       = ceil((S.tf-S.t0)/S.dt)+1;
OUT = struct();

if    S.U_out_interval < 1; OUT.U_out_steps = 0;
else; OUT.U_out_steps = ceil(S.max_steps/S.U_out_interval);
end

if    S.Uex_out_interval < 1; OUT.Uex_out_steps = 0;
else; OUT.Uex_out_steps = ceil(S.max_steps/S.Uex_out_interval);
end

if    S.R_out_interval < 1; OUT.R_out_steps = 0;
else; OUT.R_out_steps = ceil(S.max_steps/S.R_out_interval);
end

if    S.E_out_interval < 1; OUT.E_out_steps = 0;
else; OUT.E_out_steps = ceil(S.max_steps/S.E_out_interval);
end

OUT.iters       = S.out_iters;
OUT.I_out       = length(S.out_iters);
OUT.dx          = S.dx;
OUT.dt          = S.dt;
OUT.t0          = S.t0;
OUT.tf          = S.tf;
OUT.t           = (S.t0:S.dt:S.tf)';

OUT.PRI        = struct();
OUT.PRI.Rnorm  = cell(S.max_steps,1);
OUT.PRI.EnormX =  nan(S.max_steps,3);
OUT.PRI.U      = cell(S.max_steps,1);
OUT.PRI.Uex    = cell(S.max_steps,1);
OUT.PRI.R      = cell(S.max_steps,1);
OUT.PRI.E      = cell(S.max_steps,1);

OUT.ERR        = struct();
OUT.ERR.Rnorm  = cell(S.max_steps,S.Niters+1);
OUT.ERR.EnormX = nan(S.max_steps,S.Niters+1,3);
OUT.ERR.R      = cell(S.max_steps,S.Niters+1);
OUT.ERR.E      = cell(S.max_steps,S.Niters+1);
OUT.ERR.EE     = cell(S.max_steps,S.Niters+1);
end
%% Preprocessing steps
S = dt_options(S);
S.count = 1;
stencil = soln_stencil( S.stencil_size, S.N, 2 );
for i = S.stencil_size-1:-1:0
    time = OUT.t(1)-i*S.dt;
    Uex = S.ex_soln.eval(grid.x,time);
    stencil = stencil.push(Uex,time);
end

S.isBDF2 = (isa(S.integrator,'BDF2_type'));
if (S.isBDF2)
    S.integrator.um2 = stencil.U(:,S.stencil_size-1);
    S.integrator.um1 = stencil.U(:,S.stencil_size);
end
soln.U = stencil.U(:,S.stencil_size);
OUT = output_primal_stuff(OUT,S,soln.U,soln.U,soln.E,0,1);
OUT = output_ete_stuff(OUT,S,Esoln.U,0,1);
for i = 1:S.Niters
    OUT = output_ete_iter_stuff(OUT,S,0,soln.U,1,i);
end
count = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advance Solution
solving = true;
while solving
    time = OUT.t(count);
    fprintf(S.string_fmt,time,count-1,S.max_steps-1);
    solving = (count<S.max_steps);
    [soln,stencil,OUT,S] = step_primal(grid,soln,stencil,OUT,S,count);
    [Esoln,stencil,OUT,S] = step_ete(grid,Esoln,stencil,OUT,S,count);
    [EsolnIC,stencil,OUT,S] = step_iterates(grid,EsolnIC,stencil,OUT,S,count);
    count = count + 1;
    OUT = primal_cleanup(OUT);
end
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
function [soln,stencil,OUT,S] = populate_stencil(grid,soln,stencil,OUT,S)
for i = S.count:S.stencil_size % fill stencil with primal solution
    count = i;
    time = OUT.t(count);
    fprintf(S.string_fmt,time,count-1,S.max_steps-1);
    Uex = S.ex_soln.eval(grid.x,time);
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,i) S.R_BC1(u,i,Uex(1));
    R_BC2 = @(u,i) S.R_BC2(u,i,Uex(grid.N));
    u_old = soln.U;
    [u_new,resnorm,S,S.integrator] = S.integrator.step(...
                      u_old, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    soln.U = u_new;
    stencil = stencil.push(soln.U,time); % update stencil
    soln.E = soln.U - Uex;
    OUT = output_primal_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,count);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,OUT,S] = populate_ete(grid,Esoln,stencil,OUT,S)
for i = S.count:S.stencil_size % march forward with ETE solution
    time = stencil.t(i);
    U = stencil.U(:,i,1);
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,i) S.R_BC1(u,i,0);
    R_BC2 = @(u,i) S.R_BC2(u,i,0);
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    [u,dudt] = S.LS_T.eval(stencil,time,1); % Reconstruction from primal
    TE = ss_residual_cont(S,S.LS_S,u);
    TE = TE + dudt;
    RHS = @(e) S.ETE_RHS(U,e,RU,TE);
    LHS = @(e) S.ETE_LHS(U,e);
    e_old = Esoln.U;
    [e_new,resnorm,S,S.ETEintegrator] = S.ETEintegrator.step(...
              e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    Esoln.U = e_new;
    stencil.U(:,i,2) = U - e_new; % update next stencil with corrected solution
    OUT = output_ete_stuff(OUT,S,Esoln.U,resnorm,i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,OUT,S] = populate_iterates(grid,Esoln,stencil,OUT,S)
for j = 1:S.Niters   % iterative correction steps
    estError = zeros(grid.N,S.stencil_size); % Make array for estimated error
    if (S.isBDF2)
        S.ETEintegratorIC.um2 = zeros(grid.N,1); % zero out error estimates
        S.ETEintegratorIC.um1 = zeros(grid.N,1);
    end
    Esoln.U = zeros(grid.N,1);
    for i = S.count:S.stencil_size % march forward in time with ETE solution
        U = stencil.U(:,i,2); % corrected solution stencil
        %% set BCS
        L_BC1 = S.L_BC1;
        L_BC2 = S.L_BC2;
        R_BC1 = @(u,n) S.R_BC1(u,n,0); % dirichlet BC
        R_BC2 = @(u,n) S.R_BC2(u,n,0); % dirichlet BC
        %% calculate residual at current primal solution
        RU = ss_residual(U,grid.dx,S.nu,grid.N);
        %% Estimate truncation error
        time = stencil.t(i);
        [u,dudt] = S.LS_T.eval(stencil,time,2);
        TE = ss_residual_cont(S,S.LS_S,u);
        TE = TE + dudt;
        %% Define system of equations
        RHS = @(e) S.ETE_RHS(U,e,RU,TE);
        LHS = @(e) S.ETE_LHS(U,e);
        %% set initial guess
        e_old = Esoln.U;
        %% solve ETE
        [e_new,resnorm,S,S.ETEintegratorIC] = S.ETEintegratorIC.step(...
            e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
        %% save error estimate
        estError(:,i) = e_new;
        Esoln.U = e_new;
    end
    %% apply correction to stencil
    stencil.U(:,:,2) = stencil.U(:,:,2) - estError;
    %% output info
    for i = 2:S.stencil_size
        UE = stencil.U(:,i,2); % corrected solution stencil
        OUT = output_ete_iter_stuff(OUT,S,resnorm,UE,i,j);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,OUT,S] = populate_iterates_start(grid,Esoln,stencil,OUT,S)
for j = 1:S.Niters   % iterative correction steps
    time = stencil.t(2);
    U = stencil.U(:,2,2);
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,i) S.R_BC1(u,i,0);
    R_BC2 = @(u,i) S.R_BC2(u,i,0);
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    % Reconstruction from primal
    [u,dudt] = S.LS_T.eval(stencil,time,2);
    TE = ss_residual_cont(S,S.LS_S,u);
    TE = TE + dudt;
    RHS = @(e) S.ETE_RHS(U,e,RU,TE);
    LHS = @(e) S.ETE_LHS(U,e);
    
    e_old = Esoln.U;
    tmp_integrator = SDIRK2_type(grid,Esoln,S);
    [e_new,~,S,~] = tmp_integrator.step(...
        e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    Esoln.U = e_new;
    estError = zeros(grid.N,S.stencil_size); % Make array for estimated error
    estError(:,2) = Esoln.U;
    S.ETEintegratorIC.um1 = e_new;
%     S.ETEintegratorIC.um2 = zeros(grid.N,1); % zero out error estimates
    
    for i = S.count:S.stencil_size % march forward in time with ETE solution
        U = stencil.U(:,i,2); % corrected solution stencil
        %% set BCS
        L_BC1 = S.L_BC1;
        L_BC2 = S.L_BC2;
        R_BC1 = @(u,n) S.R_BC1(u,n,0); % dirichlet BC
        R_BC2 = @(u,n) S.R_BC2(u,n,0); % dirichlet BC
        %% calculate residual at current primal solution
        RU = ss_residual(U,grid.dx,S.nu,grid.N);
        %% Estimate truncation error
        time = stencil.t(i);
        [u,dudt] = S.LS_T.eval(stencil,time,2);
        TE = ss_residual_cont(S,S.LS_S,u);
        TE = TE + dudt;
        %% Define system of equations
        RHS = @(e) S.ETE_RHS(U,e,RU,TE);
        LHS = @(e) S.ETE_LHS(U,e);
        %% set initial guess
        e_old = Esoln.U;
        %% solve ETE
        [e_new,resnorm,S,S.ETEintegratorIC] = S.ETEintegratorIC.step(...
            e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
        %% save error estimate
        estError(:,i) = e_new;
        Esoln.U = e_new;
    end
    %% apply correction to stencil
    stencil.U(:,:,2) = stencil.U(:,:,2) - estError;
    for i = 2:S.stencil_size
        UE = stencil.U(:,i,2); % corrected solution stencil
        OUT = output_ete_iter_stuff(OUT,S,resnorm,UE,i,j);% output info
    end
end
end
function [Esoln,stencil,OUT,S] = step_iterates(grid,Esoln,stencil,OUT,S,count)
% 1st step is initial condition (error is 0)
i = stencil.queue_length;
for j = 1:S.Niters   % iterative correction steps
    if (S.isBDF2)
        S.ETEintegratorIC.um2 = 0*(stencil.U(:,i-2,1)-stencil.U(:,i-2,2));
        S.ETEintegratorIC.um1 = 0*(stencil.U(:,i-1,1)-stencil.U(:,i-1,2));
    end
    e_old = 0*(stencil.U(:,i-1,1)-stencil.U(:,i-1,2));
    time = stencil.t(i);
    U = stencil.U(:,i,2); % corrected solution stencil
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,n) S.R_BC1(u,n,0); % dirichlet BC
    R_BC2 = @(u,n) S.R_BC2(u,n,0); % dirichlet BC
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    [u,dudt] = S.LS_T.eval(stencil,time,2);
%     [ua,dudta] = S.LS_T.eval(stencil,time,1);
    TE = ss_residual_cont(S,S.LS_S,u);
    TE = TE + dudt;
    
%     TEa = ss_residual_cont(S,S.LS_S,ua);
%     TEa = TEa + dudta;
    
    RHS = @(e) S.ETE_RHS(U,e,RU,TE);
    LHS = @(e) S.ETE_LHS(U,e);
    [e_new,resnorm,S,S.ETEintegratorIC] = S.ETEintegratorIC.step(...
        e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
%     Esoln.U = e_new;
    UE = U - e_new; % corrected solution
    stencil.U(:,i,2) = UE; % update stencil with corrected solution
    OUT = output_ete_iter_stuff(OUT,S,resnorm,UE,count,j);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,OUT,S] = step_ete(grid,Esoln,stencil,OUT,S,count)
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
    
    e_old = Esoln.U;
    [e_new,resnorm,S,S.ETEintegrator] = S.ETEintegrator.step(...
              e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    Esoln.U = e_new;
    OUT = output_ete_stuff(OUT,S,Esoln.U,resnorm,count);
    stencil.U(:,stencil.queue_length,2) = U - e_new; % update next stencil with corrected solution
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,stencil,OUT,S] = step_primal(grid,soln,stencil,OUT,S,count)
    time = OUT.t(count);
    Uex = S.ex_soln.eval(grid.x,time);
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,i) S.R_BC1(u,i,Uex(1));
    R_BC2 = @(u,i) S.R_BC2(u,i,Uex(grid.N));
    u_old = soln.U;
    [u_new,resnorm,S,S.integrator] = S.integrator.step(...
                      u_old, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    soln.U = u_new;
    stencil = stencil.push(soln.U,time); % update stencil
    soln.E = soln.U - Uex;
    OUT = output_primal_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,count);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_ete_stuff(OUT,S,E,resnorm,count)
    OUT.ERR.Rnorm{count,1}  = resnorm;
    tmp = abs(E-OUT.PRI.E{count});
    OUT.ERR.EnormX(count,1,1) = sum(tmp)/S.N;          % L1
    OUT.ERR.EnormX(count,1,2) = sqrt(sum(tmp.^2)/S.N); % L2
    OUT.ERR.EnormX(count,1,3) = max(tmp);              % L3
    if (S.E_out_interval~=0)&&(mod(count-1,S.E_out_interval) == 0)
        OUT.ERR.E{count,1}  = E;
        OUT.ERR.EE{count,1} = E-OUT.PRI.E{count};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_ete_iter_stuff(OUT,S,resnorm,UE,count,iter)
    OUT.ERR.Rnorm{count,iter+1}  = resnorm;
    E2 = OUT.PRI.U{count}- UE;
    tmp = abs(E2-OUT.PRI.E{count});
    OUT.ERR.EnormX(count,iter+1,1) = sum(tmp)/S.N;          % L1
    OUT.ERR.EnormX(count,iter+1,2) = sqrt(sum(tmp.^2)/S.N); % L2
    OUT.ERR.EnormX(count,iter+1,3) = max(tmp);              % L3
    if any(OUT.iters==iter)
        if (S.E_out_interval~=0)&&(mod(count-1,S.E_out_interval) == 0)
            OUT.ERR.E{count,iter+1}  = E2;
            OUT.ERR.EE{count,iter+1} = E2-OUT.PRI.E{count};
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_primal_stuff(OUT,S,U,Uex,E,resnorm,count)
    OUT.PRI.Rnorm{count}  = resnorm;
    OUT.PRI.EnormX(count,1) = sum(abs(E))/S.N;     % L1
    OUT.PRI.EnormX(count,2) = sqrt(sum(E.^2)/S.N); % L2
    OUT.PRI.EnormX(count,3) = max(abs(E));         % L3
    OUT.PRI.U{count} = U;
    OUT.PRI.E{count} = E;
    if mod(count-1,S.Uex_out_interval) == 0
        OUT.PRI.Uex{count} = Uex;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = primal_cleanup(OUT)
OUT.PRI.U   = OUT.PRI.U(~cellfun('isempty',OUT.PRI.U));
OUT.PRI.Uex = OUT.PRI.Uex(~cellfun('isempty',OUT.PRI.Uex));
OUT.PRI.R   = OUT.PRI.R(~cellfun('isempty',OUT.PRI.R));
OUT.PRI.E   = OUT.PRI.E(~cellfun('isempty',OUT.PRI.E));
end