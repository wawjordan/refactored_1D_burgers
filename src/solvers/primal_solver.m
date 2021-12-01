function [soln,OUT,S] = primal_solver(grid,soln,S)
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
OUT.dx          = S.dx;
OUT.dt          = S.dt;
OUT.t0          = S.t0;
OUT.tf          = S.tf;
OUT.t           = (S.t0:S.dt:S.tf)';

OUT.PRI        = struct();
OUT.PRI.Rnorm  = cell(S.max_steps,1);
OUT.PRI.EnormX =  nan(S.max_steps,1);
OUT.PRI.U      = cell(S.max_steps,1);
OUT.PRI.Uex    = cell(S.max_steps,1);
OUT.PRI.R      = cell(S.max_steps,1);
OUT.PRI.E      = cell(S.max_steps,1);
%% Preprocessing steps
S = dt_options(S);
int_test = (isa(S.integrator,'BDF2_type'));
if S.BDF2_startup == 0
if (int_test)
    time = OUT.t(1)-OUT.dt;
    Uex = S.ex_soln.eval(grid.x,time);
    S.integrator.um2 = Uex;
end

time = OUT.t(1);
Uex = S.ex_soln.eval(grid.x,time);
if (int_test)
    S.integrator.um1 = Uex;
end
Uex = S.ex_soln.eval(grid.x,time);
if (int_test)
    S.integrator.um1 = Uex;
end
soln.U = Uex;
OUT = output_primal_stuff(OUT,S,Uex,Uex,soln.E,[],1);
count = 2;
%%%%%%%%%%%%%%%%%
else
time = OUT.t(1);
Uex = S.ex_soln.eval(grid.x,time);
soln.U = Uex;
if (int_test)
    S.integrator.um2 = Uex;
    [soln,OUT,S] = startup_primal(grid,soln,OUT,S);
    count = 3;
else
    OUT = output_primal_stuff(OUT,S,Uex,Uex,soln.E,[],1);
    count = 2;
end
end
%% Advance Primal
solving = true;
while solving
    time = OUT.t(count);
    fprintf(S.string_fmt,time,count-1,S.max_steps-1);
    solving = (count<S.max_steps);
    [soln,OUT,S] = step_primal(grid,soln,OUT,S,count);
    count = count + 1;
end
OUT = primal_cleanup(OUT);
end
function S = dt_options(S)
L1 = length(char(regexp(string(S.dt),'(?<=\.)\d*','match')));
L2 = length(char(string(S.max_steps)));
S.string_fmt = ...
    sprintf('t = %%#%d.%dg (timestep: %%%dd/%%d)\\n',min(10,L1+4),L1,L2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,OUT,S] = step_primal(grid,soln,OUT,S,count)
    time = OUT.t(count);
    Uex = S.ex_soln.eval(grid.x,time);
    L_BC1 = S.L_BC1;                         % LHS BC, left boundary
    L_BC2 = S.L_BC2;                         % LHS BC, right boundary
    R_BC1 = @(u,i) S.R_BC1(u,i,Uex(1));      % RHS BC, left boundary
    R_BC2 = @(u,i) S.R_BC2(u,i,Uex(grid.N)); % RHS BC, right boundary
    u_old = soln.U;
    [u_new,resnorm,S,S.integrator] = S.integrator.step(...
           u_old, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
    soln.U = u_new;
    soln.E = soln.U - Uex;
    OUT = output_primal_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,count);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,OUT,S] = startup_primal(grid,soln,OUT,S)
time = OUT.t(2);
Uex = S.ex_soln.eval(grid.x,time);
L_BC1 = S.L_BC1;                         % LHS BC, left boundary
L_BC2 = S.L_BC2;                         % LHS BC, right boundary
R_BC1 = @(u,i) S.R_BC1(u,i,Uex(1));      % RHS BC, left boundary
R_BC2 = @(u,i) S.R_BC2(u,i,Uex(grid.N)); % RHS BC, right boundary
u_old = soln.U;
tmp_integrator = SDIRK2_type(grid,soln,S);
[u_new,resnorm,S,~] = tmp_integrator.step(...
    u_old, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );

soln.U = u_new;
soln.E = soln.U - Uex;
OUT = output_primal_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,2);
S.integrator.um1 = u_new;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [soln,OUT,S] = startup_primal(grid,soln,OUT,S)
%     utmp = zeros(grid.N,2);
%     for i = 1:2
%         time = OUT.t(i+1);
%         Uex = S.ex_soln.eval(grid.x,time);
%         L_BC1 = S.L_BC1;                         % LHS BC, left boundary
%         L_BC2 = S.L_BC2;                         % LHS BC, right boundary
%         R_BC1 = @(u,i) S.R_BC1(u,i,Uex(1));      % RHS BC, left boundary
%         R_BC2 = @(u,i) S.R_BC2(u,i,Uex(grid.N)); % RHS BC, right boundary
%         u_old = soln.U;
%         tmp_integrator = SDIRK2_type(grid,soln,S);
%         [u_new,resnorm,S,~] = tmp_integrator.step(...
%             u_old, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
%         
%         soln.U = u_new;
%         utmp(:,i) = u_new;
%         soln.E = soln.U - Uex;
%         OUT = output_primal_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,i+1);
%     end
%     S.integrator.um2 = utmp(:,1);
%     S.integrator.um1 = utmp(:,2);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_primal_stuff(OUT,S,U,Uex,E,resnorm,count)
    OUT.PRI.Rnorm{count}  = resnorm;
    OUT.PRI.EnormX(count,1) = sum(abs(E))/S.N;     % L1
    OUT.PRI.EnormX(count,2) = sqrt(sum(E.^2)/S.N); % L2
    OUT.PRI.EnormX(count,3) = max(abs(E));         % L3
    if mod(count-1,S.U_out_interval) == 0
        OUT.PRI.U{count} = U;
    end
    if mod(count-1,S.Uex_out_interval) == 0
        OUT.PRI.Uex{count} = Uex;
    end
    if mod(count-1,S.E_out_interval) == 0
        OUT.PRI.E{count} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = primal_cleanup(OUT)
OUT.PRI.U   = OUT.PRI.U(~cellfun('isempty',OUT.PRI.U));
OUT.PRI.Uex = OUT.PRI.Uex(~cellfun('isempty',OUT.PRI.Uex));
OUT.PRI.R   = OUT.PRI.R(~cellfun('isempty',OUT.PRI.R));
OUT.PRI.E   = OUT.PRI.E(~cellfun('isempty',OUT.PRI.E));
end