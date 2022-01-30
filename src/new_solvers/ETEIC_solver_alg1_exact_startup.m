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
OUT.PRI.Enorm  = zeros(1,3);
OUT.PRI.EnormX =  nan(S.max_steps,3);
OUT.PRI.U      = cell(S.max_steps,1);
OUT.PRI.Uex    = cell(S.max_steps,1);
OUT.PRI.E      = cell(S.max_steps,1);

OUT.ERR        = struct();
OUT.ERR.Rnorm  = cell(S.max_steps,S.Niters+1);
OUT.ERR.Enorm  = zeros(S.Niters+1,3);
OUT.ERR.EnormX = nan(S.max_steps,S.Niters+1,3);
OUT.ERR.E      = cell(S.max_steps,S.Niters+1);
OUT.ERR.EE     = cell(S.max_steps,S.Niters+1);
end
%% Preprocessing steps
S = dt_options(S);
S.count = 1;
stencil = soln_stencil( S.stencil_size, S.N, 2 );

%% Backwards Extrapolation
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
%% Exact Solution
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
OUT = output_primal_stuff(OUT,S,soln.U,soln.U,soln.E,0,1);
OUT = output_ete_stuff(OUT,S,Esoln.U,soln.E,0,1);
for i = 1:S.Niters
    OUT = output_ete_iter_stuff(OUT,S,soln,0,soln.U,1,i);
end
count = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advance Solution
solving = true;
while solving
    time = OUT.t(count);
    fprintf(S.string_fmt,time,count-1,S.max_steps-1);
    solving = (count<S.max_steps);
    [soln,   stencil,OUT,S] = step_primal(  grid,soln,        stencil,OUT,S,count);
    [Esoln,  stencil,OUT,S] = step_ete(     grid,soln,Esoln,  stencil,OUT,S,count);
    [EsolnIC,stencil,OUT,S] = step_iterates(grid,soln,EsolnIC,stencil,OUT,S,count);
    count = count + 1;
    OUT = primal_cleanup(OUT,S);
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
function [Esoln,stencil,OUT,S] = step_iterates(grid,soln,Esoln,stencil,OUT,S,count)
i = stencil.queue_length; % location in stencil
for j = 1:S.Niters        % iterative correction steps
    % Error must be initialized to 0 for stability
    if (S.isBDF2)
        S.ETEintegratorIC.um2 = 0*S.ETEintegratorIC.um2;
        S.ETEintegratorIC.um1 = 0*S.ETEintegratorIC.um1;
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
%     OUT = output_ete_iter_stuff(OUT,S,resnorm,UE,count,j);
    OUT = output_ete_iter_stuff(OUT,S,soln,resnorm,UE,count,j);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Esoln,stencil,OUT,S] = step_ete(grid,soln,Esoln,stencil,OUT,S,count)
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
    OUT = output_ete_stuff(OUT,S,Esoln.U,soln.E,resnorm,count);
    % update next stencil with corrected solution
    stencil.U(:,stencil.queue_length,2) = U - e_new;
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
function OUT = output_primal_stuff(OUT,S,U,Uex,E,resnorm,count)
    OUT.PRI.Rnorm{count}  = resnorm;
    % Running totals
    OUT.PRI.Enorm(1) = OUT.PRI.Enorm(1) + sum(abs(E));    % L_1
    OUT.PRI.Enorm(2) = OUT.PRI.Enorm(2) + sum(E.^2);      % L_2
    OUT.PRI.Enorm(3) = max(OUT.PRI.Enorm(3),max(abs(E))); % L_inf
    % Spatial norms
    OUT.PRI.EnormX(count,1) = sum(abs(E))/S.N;            % L_1
    OUT.PRI.EnormX(count,2) = sqrt(sum(E.^2)/S.N);        % L_2
    OUT.PRI.EnormX(count,3) = max(abs(E));                % L_inf
    % Vector outputs
    SOLN_OUT  = mod(count-1,S.U_out_interval) == 0;
    EXACT_OUT = mod(count-1,S.Uex_out_interval) == 0;
    ERR_OUT   = mod(count-1,S.E_out_interval) == 0;
    ENDS      = count == S.max_steps || count == 1;
    if SOLN_OUT || ENDS
        OUT.PRI.U{count} = U;
    end
    if EXACT_OUT || ENDS
        OUT.PRI.Uex{count} = Uex;
    end
    if ERR_OUT || ENDS
        OUT.PRI.E{count} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_ete_stuff(OUT,S,E,EP,resnorm,count)
    tmp = abs(E-EP);
    % Running totals
    OUT.ERR.Enorm(1,1) =     OUT.ERR.Enorm(1,1) + sum(tmp);   % L_1
    OUT.ERR.Enorm(1,2) =     OUT.ERR.Enorm(1,2) + sum(tmp.^2);% L_2
    OUT.ERR.Enorm(1,3) = max(OUT.ERR.Enorm(1,3),max(tmp));    % L_inf
    % Spatial norms
    OUT.ERR.EnormX(count,1,1) = sum(tmp)/S.N;                 % L_1
    OUT.ERR.EnormX(count,1,2) = sqrt(sum(tmp.^2)/S.N);        % L_2
    OUT.ERR.EnormX(count,1,3) = max(tmp);                     % L_inf
    % Vector outputs
    ERR_OUT = (S.E_out_interval~=0)&&(mod(count-1,S.E_out_interval) == 0);
    RES_OUT = (S.R_out_interval~=0)&&(mod(count-1,S.R_out_interval) == 0);
    ENDS    = count == S.max_steps || count == 1;
    if ERR_OUT || ENDS
        OUT.ERR.E{count,1}  = E; % estimated error
        OUT.ERR.EE{count,1} = E-EP; % error in error estimate
    end
    if RES_OUT || ENDS
        OUT.ERR.Rnorm{count,1}  = resnorm;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_ete_iter_stuff(OUT,S,soln,resnorm,UE,count,iter)
    E2 = soln.U- UE;
    tmp = abs(E2-soln.E);
    % Running totals
    OUT.ERR.Enorm(iter+1,1) =     OUT.ERR.Enorm(1,1) + sum(tmp);   % L_1
    OUT.ERR.Enorm(iter+1,2) =     OUT.ERR.Enorm(1,2) + sum(tmp.^2);% L_2
    OUT.ERR.Enorm(iter+1,3) = max(OUT.ERR.Enorm(1,3),max(tmp));    % L_inf
    % Spatial norms
    OUT.ERR.EnormX(count,iter+1,1) = sum(tmp)/S.N;                 % L_1
    OUT.ERR.EnormX(count,iter+1,2) = sqrt(sum(tmp.^2)/S.N);        % L_2
    OUT.ERR.EnormX(count,iter+1,3) = max(tmp);                     % L_inf
    % Vector outputs
    ERR_OUT = (S.E_out_interval~=0)&&(mod(count-1,S.E_out_interval) == 0);
    RES_OUT = (S.R_out_interval~=0)&&(mod(count-1,S.R_out_interval) == 0);
    ENDS    = count == S.max_steps || count == 1;
    if any(OUT.iters==iter)
        if ERR_OUT || ENDS
            OUT.ERR.E{count,iter+1}  = E2; % estimated error
            OUT.ERR.EE{count,iter+1} = E2-soln.E; % error in error estimate
        end
        if RES_OUT || ENDS
            OUT.ERR.Rnorm{count,iter+1}  = resnorm;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = primal_cleanup(OUT,S)
OUT.PRI.U   = OUT.PRI.U(~cellfun('isempty',OUT.PRI.U));
OUT.PRI.Uex = OUT.PRI.Uex(~cellfun('isempty',OUT.PRI.Uex));
OUT.PRI.E   = OUT.PRI.E(~cellfun('isempty',OUT.PRI.E));

OUT.ERR.E  = OUT.ERR.E( ~cellfun('isempty',OUT.ERR.E(:,1)),:);
OUT.ERR.EE = OUT.ERR.EE(~cellfun('isempty',OUT.ERR.EE(:,1)),:);
OUT.ERR.Rnorm = OUT.ERR.Rnorm(~cellfun('isempty',OUT.ERR.Rnorm(:,1)),:);

Nt = S.max_steps - 1; % number of time steps (excluding 1st time step)
Nx = S.N;             % number of spatial nodes
OUT.PRI.Enorm(1) = OUT.PRI.Enorm(1)/( Nt*Nx );
OUT.PRI.Enorm(2) = sqrt( OUT.PRI.Enorm(2)/( Nt*Nx ) );
OUT.ERR.Enorm(:,1) = OUT.ERR.Enorm(:,1)/( Nt*Nx );
OUT.ERR.Enorm(:,2) = sqrt( OUT.ERR.Enorm(:,2)/( Nt*Nx ) );
end