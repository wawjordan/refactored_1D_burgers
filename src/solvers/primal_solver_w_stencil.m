function [soln,OUT,S,stencil] = primal_solver_w_stencil(grid,soln,S)
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
OUT.t           = (S.t0:S.dt:max(S.tf,S.max_steps*S.dt))';

OUT.PRI        = struct();
OUT.PRI.Rnorm  = cell(S.max_steps,1);
OUT.PRI.EnormX =  nan(S.max_steps,1);
OUT.PRI.U      = cell(S.max_steps,1);
OUT.PRI.Uex    = cell(S.max_steps,1);
OUT.PRI.R      = cell(S.max_steps,1);
OUT.PRI.E      = cell(S.max_steps,1);
%% Preprocessing steps
S = dt_options(S);
stencil = soln_stencil( S.stencil_size, S.N );
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
OUT = output_primal_stuff(OUT,S,Uex,Uex,soln.E,[],1);
%% Solve Primal to fill stencil
[soln,stencil,OUT,S] = populate_stencil(grid,soln,stencil,OUT,S);
count = S.stencil_size+1;

%% Advance Primal
solving = true;
while solving
    time = OUT.t(count);
    fprintf(S.string_fmt,time,count-1,S.max_steps-1);
    solving = (count<S.max_steps);
    [soln,stencil,OUT,S] = step_primal(grid,soln,stencil,OUT,S,count);
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
    OUT.PRI.EnormX(count) = sum(abs(E))/S.N;
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