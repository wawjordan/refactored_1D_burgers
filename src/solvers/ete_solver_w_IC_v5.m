function [Esoln,EsolnIC,soln,OUT,S,stencil] = ete_solver_w_IC_v5(grid,Esoln,EsolnIC,soln,S)
%{
S   - input structure

N  - number of spatial nodes
t0 - problem start time
tf - problem stop time
dt - problem time step
dx - spatial grid fineness measure (constant for now)

U_out_interval   - output interval for primal solution vector
Uex_out_interval - output interval for exact solution vector
R_out_interval   - output interval for unsteady residual vector
E_out_interval   - output interval for discretization error (DE) vector
out_iters        - (IC only) iterations to be saved
stencil_size     - size of stencil for solution reconstruction

OUT - output structure
U_out_steps - time step indices for primal output
R_out_steps - time step indices for resdidual output
E_out_steps - time step indices for DE output
t           - time vector
PRI         - primal solution data structure
ERR         - ETE solution data structure

Rnorm  - cell array containing residual norm history
EnormX - cell array containing spatial DE norm time history
EnormT - array containing time average of error for each spatial node
U      - cell array containing saved primal solution vectors
Uex    - cell array containing saved exact solution vectors
R      - cell array containing saved residual vectors
E      - cell array containing saved DE vectors

Local Variables:
max_steps         - maximum number of time steps
max_out_U_steps   - number of time steps for primal output
max_out_Uex_steps - number of time steps for exact solution output
max_out_R_steps   - number of time steps for residual output
max_out_E_steps   - number of time steps for DE output
max_out_iters     - number of iterations to be saved
%}
S.max_steps       = ceil((S.tf-S.t0)/S.dt)+1;
if S.U_out_interval < 1
    max_out_U_steps = 0;
else
    max_out_U_steps = ceil(S.max_steps/S.U_out_interval);
end

if S.Uex_out_interval < 1
    max_out_Uex_steps = 0;
else
    max_out_Uex_steps = ceil(S.max_steps/S.Uex_out_interval);
end

if S.R_out_interval < 1
    max_out_R_steps = 0;
else
    max_out_R_steps = ceil(S.max_steps/S.R_out_interval);
end

if S.E_out_interval < 1
    max_out_E_steps = 0;
else
    max_out_E_steps = ceil(S.max_steps/S.E_out_interval);
end

OUT = struct();
OUT.U_out_steps   = max_out_U_steps;
OUT.Uex_out_steps = max_out_Uex_steps;
OUT.R_out_steps = max_out_R_steps;
OUT.E_out_steps = max_out_E_steps;
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
%% Preprocessing steps
S = dt_options(S);
stencil = soln_stencil( S.stencil_size, S.N, S.Niters+1 );
int_test = (isa(S.integrator,'BDF2_type'));
if (int_test)
    time = OUT.t(1)-OUT.dt;
    Uex = S.ex_soln.eval(grid.x,time);
    stencil = stencil.push(Uex,time);
    S.integrator.um2 = Uex;
end

time = OUT.t(1);
Uex = S.ex_soln.eval(grid.x,time);
stencil = stencil.push(Uex,time);
if (int_test)
    S.integrator.um1 = Uex;
end
soln.U = Uex;
OUT = output_primal_stuff(OUT,S,Uex,Uex,soln.E,0,1);
OUT = output_ete_stuff(OUT,S,Esoln.U,0,1);
for i = 1:S.Niters
    OUT = output_ete_iter_stuff(OUT,S,0,Uex,1,i);
end
%% Solve Primal to fill stencil
[soln,stencil,OUT,S]  = populate_stencil(grid,soln,stencil,OUT,S);
[Esoln,stencil,OUT,S] = populate_ete(grid,Esoln,stencil,OUT,S);
[EsolnIC,stencil,OUT,S] = populate_iterates(grid,EsolnIC,stencil,OUT,S);

count = S.stencil_size+1;
%% Advance Primal
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
function S = dt_options(S)
L1 = length(char(regexp(string(S.dt),'(?<=\.)\d*','match')));
L2 = length(char(string(S.max_steps)));
S.string_fmt = ...
    sprintf('t = %%#%d.%df (timestep: %%%dd/%%d)\\n',12,L1,L2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,stencil,OUT,S] = populate_stencil(grid,soln,stencil,OUT,S)
% 1st step is initial condition
for i = 2:S.stencil_size % fill stencil with primal solution
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
% 1st step is initial condition (error is 0)
for i = 2:S.stencil_size % march forward with ETE solution
    time = stencil.t(i);
    U = stencil.U(:,i,1);
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,i) S.R_BC1(u,i,0);
    R_BC2 = @(u,i) S.R_BC2(u,i,0);
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    % Reconstruction from primal
    [u,dudt] = S.LS_T.eval(stencil,time,1);
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
%     stencil = stencil.push(Esoln.U,time); % update stencil
% %     OUT = output_ete_stuff(OUT,S,Esoln.U,Uex,soln.E,resnorm,count);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,OUT,S] = populate_iterates(grid,Esoln,stencil,OUT,S)
for j = 1:S.Niters   % iterative correction steps
    %% Make array for estimated error
    estError = zeros(grid.N,S.stencil_size);
    %% zero out error estimates
    S.ETEintegratorIC.um2 = zeros(grid.N,1);
    S.ETEintegratorIC.um1 = zeros(grid.N,1);
    Esoln.U = zeros(grid.N,1);
    for i = 2:S.stencil_size % march forward in time with ETE solution
        U = stencil.U(:,i,j+1); % corrected solution stencil
        %% set BCS
        L_BC1 = S.L_BC1;
        L_BC2 = S.L_BC2;
        R_BC1 = @(u,n) S.R_BC1(u,n,0); % dirichlet BC
        R_BC2 = @(u,n) S.R_BC2(u,n,0); % dirichlet BC
        %% calculate residual at current primal solution
        RU = ss_residual(U,grid.dx,S.nu,grid.N);
        %% Estimate truncation error
        time = stencil.t(i);
        [u,dudt] = S.LS_T.eval(stencil,time,j+1);
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
    stencil.U(:,:,j+1) = stencil.U(:,:,j+1) - estError;
    %% output info
    for i = 2:S.stencil_size
        UE = stencil.U(:,i,j+1); % corrected solution stencil
        if (j<S.Niters) % propagate correction to next iteration(?)
            stencil.U(:,:,j+2) = stencil.U(:,:,j+1);
        end
        OUT = output_ete_iter_stuff(OUT,S,resnorm,UE,i,j);
    end
end
% for j = 1:S.Niters
%     stencil.U(:,:,j+1) = stencil.U(:,:,end);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [Esoln,stencil,OUT,S] = step_iterates(grid,Esoln,stencil,OUT,S,count)
% 1st step is initial condition (error is 0)
i = stencil.queue_length;
for j = 1:S.Niters   % iterative correction steps
    S.ETEintegratorIC.um2 = (stencil.U(:,i-2,j)-stencil.U(:,i-2,j+1));
    S.ETEintegratorIC.um1 = (stencil.U(:,i-1,j)-stencil.U(:,i-1,j+1));
%     S.ETEintegratorIC.um2 = (stencil.U(:,i-2,j)-stencil.U(:,i-2,j+1));
%     S.ETEintegratorIC.um1 = (stencil.U(:,i-1,j)-stencil.U(:,i-1,j+1));
    e_old = (stencil.U(:,i,j)-stencil.U(:,i,j+1));
    time = stencil.t(i);
    U1 = stencil.U(:,i,1); % uncorrected solution stencil
    U = stencil.U(:,i,j); % corrected solution stencil
    L_BC1 = S.L_BC1;
    L_BC2 = S.L_BC2;
    R_BC1 = @(u,n) S.R_BC1(u,n,0); % dirichlet BC
    R_BC2 = @(u,n) S.R_BC2(u,n,0); % dirichlet BC
    RU = ss_residual(U,grid.dx,S.nu,grid.N);
    [u,dudt] = S.LS_T.eval(stencil,time,j);
    TE = ss_residual_cont(S,S.LS_S,u);
    TE = TE + dudt;
    RHS = @(e) S.ETE_RHS(U,e,RU,TE);
    LHS = @(e) S.ETE_LHS(U,e);
    [e_new,resnorm,S,S.ETEintegratorIC] = S.ETEintegratorIC.step(...
        e_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    ee = e_new-e_old;
%     if dot(ee,ee) > dot(e_old,e_old)
%         e_new = e_old;
%         rel = 0.1;%*dot(e_old,e_old)/dot(ee,ee);
%     else
%         rel = 1;
%     end
    Esoln.U = e_new;
    UE = U - e_new; % corrected solution
    stencil.U(:,i,j+1) = UE; % update stencil with corrected solution
    if (j<S.Niters)
        stencil.U(:,i,j+2) = UE; % update stencil with corrected solution
    end
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