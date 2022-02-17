function [Csoln,Fsoln,OUT,SC,SF] = CRE_solver(Cgrid,Fgrid,Csoln,Fsoln,SC,SF)
SC.max_steps = ceil((SC.tf-SC.t0)/SC.dt)+1;
SF.max_steps = ceil((SF.tf-SF.t0)/SF.dt)+1;

OUT = struct();

OUT.dxC          = SC.dx;
OUT.dtC          = SC.dt;
OUT.dxF          = SF.dx;
OUT.dtF          = SF.dt;

mx = 2; % should be 2 for this case
mt = 2; % should be 2 for this case
% mx = OUT.dxC/OUT.dxF; % should be 2 for this case
% mt = OUT.dtC/OUT.dtF; % should be 2 for this case

px = 2; % formal spatial OOA
% pt = 2; % formal temporal OOA

OUT.t0           = SC.t0;
OUT.tf           = SC.tf;

OUT.tC           = (SC.t0:SC.dt:SC.tf)';
OUT.tF           = (SF.t0:SF.dt:SF.tf)';

OUT.RICH        = struct();
OUT.RICH.COARSE = struct();
OUT.RICH.FINE   = struct();
OUT.RICH.COARSE.Enorm = zeros(1,3);
OUT.RICH.COARSE.EnormX =  nan(SC.max_steps,1);
OUT.RICH.COARSE.U      = cell(SC.max_steps,1);
OUT.RICH.COARSE.E      = cell(SC.max_steps,1);

OUT.RICH.FINE.Enorm = zeros(1,3);
OUT.RICH.FINE.EnormX =  nan(SC.max_steps,1);
OUT.RICH.FINE.U      = cell(SC.max_steps,1);
OUT.RICH.FINE.E      = cell(SC.max_steps,1);


OUT.COARSE        = struct();
OUT.COARSE.Rnorm  = cell(SC.max_steps,1);
OUT.COARSE.Enorm  = zeros(1,3);
OUT.COARSE.EnormX =  nan(SC.max_steps,1);
OUT.COARSE.U      = cell(SC.max_steps,1);
OUT.COARSE.R      = cell(SC.max_steps,1);
OUT.COARSE.E      = cell(SC.max_steps,1);

OUT.FINE        = struct();
OUT.FINE.Rnorm  = cell(SF.max_steps,1);
OUT.FINE.Enorm  = zeros(1,3);
OUT.FINE.EnormX =  nan(SF.max_steps,1);
OUT.FINE.U      = cell(SF.max_steps,1);
OUT.FINE.R      = cell(SF.max_steps,1);
OUT.FINE.E      = cell(SF.max_steps,1);

OUT.Uex         = cell(SF.max_steps,1);
%% Preprocessing steps
SC = dt_options(SC);
RnormC = 0;
RnormF = 0;
%% Assumed BDF2 type
% int_test = (isa(S.Cintegrator,'BDF2_type'));

% coarse grid (step 1)
SC.integrator.um2 = SC.ex_soln.eval(Cgrid.x,OUT.tC(1)-OUT.dtC);
SC.integrator.um1 = SC.ex_soln.eval(Cgrid.x,OUT.tC(1));
Csoln.U = SC.integrator.um1;
OUT = output_coarse_stuff(OUT,SC,Csoln.U,Csoln.U,Csoln.E,RnormC,1);

% fine grid (step 1)
SF.integrator.um2 = SF.ex_soln.eval(Fgrid.x,OUT.tF(1)-OUT.dtF);
SF.integrator.um1 = SF.ex_soln.eval(Fgrid.x,OUT.tF(1));
Fsoln.U = SF.integrator.um1;
OUT = output_fine_stuff(OUT,SF,Fsoln.U,Fsoln.U,Fsoln.E,RnormF,1);

% % inner loop
% for nf = 2:mt
%     [Fsoln,OUT,SF] = step_fine(Fgrid,Fsoln,OUT,SF,nf);
% end
% 
% % calculate extrapolation
% UR = ( (mx^px)*Fsoln.U(1:mx:end) - Csoln.U )/( (mx^px) - 1);
% 
% % transfer solution
% % Csoln.U = UR;
% % Fsoln.U(1:mx:end) = UR;
count  = 2;
count2 = 2;

%% Continue
solving = true;
while solving
    timeC = OUT.tC(count);
    fprintf(SC.string_fmt,timeC,count-1,SC.max_steps-1);
    solving = (count<SC.max_steps);
    [Csoln,OUT,SC] = step_coarse(Cgrid,Csoln,OUT,SC,count);
    % inner loop
    for nf = 1:mt
        [Fsoln,OUT,SF] = step_fine(Fgrid,Fsoln,OUT,SF,count2);
        count2 = count2 + 1;
    end
    
    % calculate extrapolation
    UR = ( (mx^px)*Fsoln.U(1:mx:end) - Csoln.U )/( (mx^px) - 1);
    % transfer solution (coarse)
%     Csoln.U = UR;
    Uex = SC.ex_soln.eval(Cgrid.x,OUT.tC(count));
    OUT = output_coarse_extrap_stuff(OUT,SC,UR,UR-Uex,count);
%     OUT = output_coarse_stuff(OUT,SC,Csoln.U,Uex,Csoln.E,0,count);
    % transfer solution (fine)
%     Fsoln.U = interp1(Fgrid.x(1:2:end),UR,Fgrid.x);
    URF = Fsoln.U + interp1(Fgrid.x(1:2:end),UR-Fsoln.U(1:2:end),Fgrid.x);
    Uex = SF.ex_soln.eval(Fgrid.x,OUT.tF(count2-1));
    OUT = output_fine_extrap_stuff(OUT,SC,URF,URF-Uex,count);
%     OUT = output_fine_stuff(OUT,SF,Fsoln.U,Uex,Fsoln.E,0,count2-1);
    count = count + 1;
end
OUT = primal_cleanup(OUT,SC,SF);
end
function S = dt_options(S)
L1 = length(char(regexp(string(S.dt),'(?<=\.)\d*','match')));
L2 = length(char(string(S.max_steps)));
S.string_fmt = ...
    sprintf('t = %%#%d.%dg (timestep: %%%dd/%%d)\\n',min(10,L1+4),L1,L2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,OUT,S] = step_coarse(grid,soln,OUT,S,count)
    time = OUT.tC(count);
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
    OUT = output_coarse_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,count);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [soln,OUT,S] = step_fine(grid,soln,OUT,S,count)
    time = OUT.tF(count);
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
    OUT = output_fine_stuff(OUT,S,soln.U,Uex,soln.E,resnorm,count);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_coarse_stuff(OUT,S,U,Uex,E,resnorm,count)
    % Residual norms from time sub-iterations
    OUT.COARSE.Rnorm{count}  = resnorm;
    % Running totals
    OUT.COARSE.Enorm(1) = OUT.COARSE.Enorm(1) + sum(abs(E));    % L_1
    OUT.COARSE.Enorm(2) = OUT.COARSE.Enorm(2) + sum(E.^2);      % L_2
    OUT.COARSE.Enorm(3) = max(OUT.COARSE.Enorm(3),max(abs(E))); % L_inf
    % Spatial norms
    OUT.COARSE.EnormX(count,1) = sum(abs(E))/S.N;            % L_1
    OUT.COARSE.EnormX(count,2) = sqrt(sum(E.^2)/S.N);        % L_2
    OUT.COARSE.EnormX(count,3) = max(abs(E));                % L_inf
    % Vector outputs
    SOLN_OUT  = mod(count-1,S.U_out_interval) == 0;
    EXACT_OUT = mod(count-1,S.Uex_out_interval) == 0;
    ERR_OUT   = mod(count-1,S.E_out_interval) == 0;
    ENDS      = count == S.max_steps || count == 1;
    if SOLN_OUT || ENDS
        OUT.COARSE.U{count} = U;
    end
    if EXACT_OUT || ENDS
        OUT.COARSE.Uex{count} = Uex;
    end
    if ERR_OUT || ENDS
        OUT.COARSE.E{count} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_coarse_extrap_stuff(OUT,S,UR,E,count)
    % Running totals
    OUT.RICH.COARSE.Enorm(1) = OUT.RICH.COARSE.Enorm(1) + sum(abs(E));    % L_1
    OUT.RICH.COARSE.Enorm(2) = OUT.RICH.COARSE.Enorm(2) + sum(E.^2);      % L_2
    OUT.RICH.COARSE.Enorm(3) = max(OUT.RICH.COARSE.Enorm(3),max(abs(E))); % L_inf
    % Spatial norms
    OUT.RICH.COARSE.EnormX(count,1) = sum(abs(E))/S.N;            % L_1
    OUT.RICH.COARSE.EnormX(count,2) = sqrt(sum(E.^2)/S.N);        % L_2
    OUT.RICH.COARSE.EnormX(count,3) = max(abs(E));                % L_inf
    % Vector outputs
    SOLN_OUT  = mod(count-1,S.U_out_interval) == 0;
    ERR_OUT   = mod(count-1,S.E_out_interval) == 0;
    ENDS      = count == S.max_steps || count == 1;
    if SOLN_OUT || ENDS
        OUT.RICH.COARSE.U{count} = UR;
    end
    if ERR_OUT || ENDS
        OUT.RICH.COARSE.E{count} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_fine_extrap_stuff(OUT,S,UR,E,count)
    % Running totals
    OUT.RICH.FINE.Enorm(1) = OUT.RICH.FINE.Enorm(1) + sum(abs(E));    % L_1
    OUT.RICH.FINE.Enorm(2) = OUT.RICH.FINE.Enorm(2) + sum(E.^2);      % L_2
    OUT.RICH.FINE.Enorm(3) = max(OUT.RICH.FINE.Enorm(3),max(abs(E))); % L_inf
    % Spatial norms
    OUT.RICH.FINE.EnormX(count,1) = sum(abs(E))/S.N;            % L_1
    OUT.RICH.FINE.EnormX(count,2) = sqrt(sum(E.^2)/S.N);        % L_2
    OUT.RICH.FINE.EnormX(count,3) = max(abs(E));                % L_inf
    % Vector outputs
    SOLN_OUT  = mod(count-1,S.U_out_interval) == 0;
    ERR_OUT   = mod(count-1,S.E_out_interval) == 0;
    ENDS      = count == S.max_steps || count == 1;
    if SOLN_OUT || ENDS
        OUT.RICH.FINE.U{count} = UR;
    end
    if ERR_OUT || ENDS
        OUT.RICH.FINE.E{count} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = output_fine_stuff(OUT,S,U,Uex,E,resnorm,count)
    % Residual norms from time sub-iterations
    OUT.FINE.Rnorm{count}  = resnorm;
    % Running totals
    OUT.FINE.Enorm(1) = OUT.FINE.Enorm(1) + sum(abs(E));    % L_1
    OUT.FINE.Enorm(2) = OUT.FINE.Enorm(2) + sum(E.^2);      % L_2
    OUT.FINE.Enorm(3) = max(OUT.FINE.Enorm(3),max(abs(E))); % L_inf
    % Spatial norms
    OUT.FINE.EnormX(count,1) = sum(abs(E))/S.N;            % L_1
    OUT.FINE.EnormX(count,2) = sqrt(sum(E.^2)/S.N);        % L_2
    OUT.FINE.EnormX(count,3) = max(abs(E));                % L_inf
    % Vector outputs
    SOLN_OUT  = mod(count-1,S.U_out_interval) == 0;
    EXACT_OUT = mod(count-1,S.Uex_out_interval) == 0;
    ERR_OUT   = mod(count-1,S.E_out_interval) == 0;
    ENDS      = count == S.max_steps || count == 1;
    if SOLN_OUT || ENDS
        OUT.FINE.U{count} = U;
    end
    if EXACT_OUT || ENDS
        OUT.FINE.Uex{count} = Uex;
    end
    if ERR_OUT || ENDS
        OUT.FINE.E{count} = E;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = primal_cleanup(OUT,SC,SF)
%% Coarse
OUT.COARSE.U   = OUT.COARSE.U(~cellfun('isempty',OUT.COARSE.U));
OUT.COARSE.Uex = OUT.COARSE.Uex(~cellfun('isempty',OUT.COARSE.Uex));
OUT.COARSE.R   = OUT.COARSE.R(~cellfun('isempty',OUT.COARSE.R));
OUT.COARSE.E   = OUT.COARSE.E(~cellfun('isempty',OUT.COARSE.E));

Nt = SC.max_steps - 1; % number of time steps (excluding 1st time step)
Nx = SC.N;             % number of spatial nodes
OUT.COARSE.Enorm(1) = OUT.COARSE.Enorm(1)/( Nt*Nx );
OUT.COARSE.Enorm(2) = sqrt( OUT.COARSE.Enorm(2)/( Nt*Nx ) );

%% Richardson Extrapolation
OUT.RICH.COARSE.U   = OUT.RICH.COARSE.U(~cellfun('isempty',OUT.RICH.COARSE.U));
OUT.RICH.COARSE.E   = OUT.RICH.COARSE.E(~cellfun('isempty',OUT.RICH.COARSE.E));
OUT.RICH.FINE.U     = OUT.RICH.FINE.U(~cellfun('isempty',OUT.RICH.FINE.U));
OUT.RICH.FINE.E     = OUT.RICH.FINE.E(~cellfun('isempty',OUT.RICH.FINE.E));

Nt = SC.max_steps - 1; % number of time steps (excluding 1st time step)
Nx = SC.N;             % number of spatial nodes
OUT.RICH.COARSE.Enorm(1) = OUT.RICH.COARSE.Enorm(1)/( Nt*Nx );
OUT.RICH.COARSE.Enorm(2) = sqrt( OUT.RICH.COARSE.Enorm(2)/( Nt*Nx ) );

Nt = SC.max_steps - 1; % number of time steps (excluding 1st time step)
Nx = SF.N;             % number of spatial nodes
OUT.RICH.FINE.Enorm(1) = OUT.RICH.FINE.Enorm(1)/( Nt*Nx );
OUT.RICH.FINE.Enorm(2) = sqrt( OUT.RICH.FINE.Enorm(2)/( Nt*Nx ) );

%% Fine
OUT.FINE.U   = OUT.FINE.U(~cellfun('isempty',OUT.FINE.U));
OUT.FINE.Uex = OUT.FINE.Uex(~cellfun('isempty',OUT.FINE.Uex));
OUT.FINE.R   = OUT.FINE.R(~cellfun('isempty',OUT.FINE.R));
OUT.FINE.E   = OUT.FINE.E(~cellfun('isempty',OUT.FINE.E));

Nt = SF.max_steps - 1; % number of time steps (excluding 1st time step)
Nx = SF.N;             % number of spatial nodes
OUT.FINE.Enorm(1) = OUT.FINE.Enorm(1)/( Nt*Nx );
OUT.FINE.Enorm(2) = sqrt( OUT.FINE.Enorm(2)/( Nt*Nx ) );
end