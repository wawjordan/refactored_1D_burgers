function [soln,OUT,S] = primal_solver_v2(grid,soln,S)
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

%% Preprocessing steps
S = dt_options(S);
int_test = (isa(S.integrator,'BDF2_type'));
fprintf(S.string_fmt,OUT.t(1),0,S.max_steps-1);
if (int_test)                                  % BDF2
    if S.BDF2_startup == 0                     % Exact Startup
        S.integrator.um2 = S.ex_soln.eval(grid.x, OUT.t(1)-OUT.dt);
        S.integrator.um1 = S.ex_soln.eval(grid.x, OUT.t(1));
        soln.U = S.integrator.um1;
        OUT.PRI.U{1}   = soln.U;
        OUT.PRI.Uex{1} = soln.U;
        OUT.PRI.t(1)   = OUT.t(1);
        OUT.PRI.R{1}   = 0;
        OUT.PRI.Rt(1)  = OUT.t(1);
        count = 2;
    else                                       % SDIRK2 startup
        S.integrator.um2 = S.ex_soln.eval(grid.x, OUT.t(1));
        soln.U = S.integrator.um2;
        [soln,OUT,S] = BDF2_SDIRK2_startup(grid,soln,OUT,S);
        count = 3;
    end
else                                           % Other
    soln.U = S.ex_soln.eval(grid.x,OUT.t(1));
    OUT.PRI.U{1}   = soln.U;
    OUT.PRI.Uex{1} = soln.U;
    OUT.PRI.t(1)   = OUT.t(1);
    OUT.PRI.R{1}   = 0;
    OUT.PRI.Rt(1)  = OUT.t(1);
    count = 2;
end

%% Advance Primal
for i = count:S.max_steps
    time = OUT.t(i);
    fprintf(S.string_fmt,time,i-1,S.max_steps-1);
    [soln,resnorm,S] = step_primal(grid,soln,S,time);
    OUT = output_solution(OUT,S,soln.U,resnorm,i);
end
for i = count:S.max_steps
    time = OUT.t(i);
    if mod(i-1,S.U_out_interval) == 0 || i == S.max_steps
        OUT.PRI.Uex{i} = S.ex_soln.eval(grid.x,time);
    end
end
OUT = primal_cleanup(OUT);
end
function S = dt_options(S)
L1 = length(char(regexp(string(S.dt),'(?<=\.)\d*','match')));
L2 = length(char(string(S.max_steps)));
S.string_fmt = ...
    sprintf('t = %%#%d.%df (timestep: %%%dd/%%d)\\n',10,L1,L2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,OUT,S] = BDF2_SDIRK2_startup(grid,soln,OUT,S)
N = grid.N;

OUT.PRI.Rnorm{1}  = 0;
OUT.PRI.U{1} = soln.U;
OUT.PRI.Uex{1} = S.ex_soln.eval(grid.x,OUT.t(1));

% 2nd time step
Uex = S.ex_soln.eval(grid.x,OUT.t(2));       % exact solution for error calc.
L_BC1 = S.L_BC1;                         % LHS BC, left boundary
L_BC2 = S.L_BC2;                         % LHS BC, right boundary
R_BC1 = @(u,i) S.R_BC1(u,i,Uex(1));      % RHS BC, left boundary
R_BC2 = @(u,i) S.R_BC2(u,i,Uex(N));      % RHS BC, right boundary
tmp_integrator = SDIRK2_type(grid,soln,S);
[soln.U,resnorm,S,~] = tmp_integrator.step(...
    soln.U, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
S.integrator.um1 = soln.U;

OUT.PRI.R{2}  = resnorm;
OUT.PRI.U{2} = soln.U;
OUT.PRI.Uex{2} = Uex;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,resnorm,S] = step_primal(grid,soln,S,time)
    N = grid.N;
    Uex1 = S.ex_soln.eval(grid.x(1),time);
    UexN = S.ex_soln.eval(grid.x(N),time);
    L_BC1 = S.L_BC1;                         % LHS BC, left boundary
    L_BC2 = S.L_BC2;                         % LHS BC, right boundary
    R_BC1 = @(u,i) S.R_BC1(u,i,Uex1);        % RHS BC, left boundary
    R_BC2 = @(u,i) S.R_BC2(u,i,UexN);        % RHS BC, right boundary
    [soln.U,resnorm,S,S.integrator] = S.integrator.step(...
           soln.U, S, S.RHS, S.LHS, L_BC1, L_BC2, R_BC1, R_BC2 );
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = primal_cleanup(OUT)
OUT.PRI.t   = OUT.PRI.t(~isnan(OUT.PRI.t));
OUT.PRI.U   = OUT.PRI.U(  ~cellfun('isempty',OUT.PRI.U  ));
OUT.PRI.Uex = OUT.PRI.Uex(~cellfun('isempty',OUT.PRI.Uex));
OUT.PRI.Rt  = OUT.PRI.Rt(~isnan(OUT.PRI.Rt));
OUT.PRI.R   = OUT.PRI.R(  ~cellfun('isempty',OUT.PRI.R  ));
end