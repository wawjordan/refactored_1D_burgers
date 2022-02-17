function [SE1,PE1,DEBUG] = parse_IC_OOA_new(filename)

load(filename,'OUT');

% get dimensions for output arrays
L = size(OUT.ETE_History,2);   % (h) # of grid levels (space+time varying)
M = size(OUT.ETE_History(1).E,1); % (i) # of time steps  (coarsest mesh)
O = 3;                      % (j) norms

t = OUT.t{1}(:);

nNodes = OUT.Nx(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRI_Ex = zeros(L,M,O);      % Primal DE  (space norms)
PRI_px = zeros(L,M,O);      % Primal OOA (space norms)
PRI_Eg = zeros(L,O);        % Primal global DE norms
PRI_pg = zeros(L,O);        % Primal global OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETE_Ex = zeros(L,M,O);      % ETE corrected DE  (space norms)
ETE_px = zeros(L,M,O);      % ETE corrected OOA (space norms)
ETE_Eg = zeros(L,O);        % ETE global corrected DE norms
ETE_pg = zeros(L,O);        % ETE global corrected OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC_Ex = zeros(L,M,O);       % ETE w/IC corrected DE  (space norms)
IC_px = zeros(L,M,O);       % ETE w/IC corrected OOA (space norms)
IC_Eg = zeros(L,O);         % ETE w/IC global corrected DE norms
IC_pg = zeros(L,O);         % ETE w/IC global corrected OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handles Ex
for h = 1:L           % grid levels
    for i = 1:M       % time steps
        [~,index] = min( abs( OUT.t{1}(i) - OUT.t{h} ) ); % closest time step
        for k = 1:O   % norms
            PRI_Ex(h,i,k) = OUT.Primal_History(h).E(index,k);  % primal solve
            ETE_Ex(h,i,k) = OUT.ETE_History(h).E(index,1,k);   % 1 ETE solve
            IC_Ex(h,i,k)  = OUT.ETE_History(h).E(index,end,k); % Last IC
        end
    end
end

% handles px
for h = 2:L           % grid levels
    for i = 1:M       % time steps
        for k = 1:O
            ee = PRI_Ex(h-1,i,k)/PRI_Ex(h,i,k);  % primal solve
            PRI_px(h,i,k) = log(ee)/log(2);
            ee = ETE_Ex(h-1,i,k)/ETE_Ex(h,i,k);  % 1 ETE solve
            ETE_px(h,i,k) = log(ee)/log(2);
            ee = IC_Ex(h-1,i,k)/IC_Ex(h,i,k);    % Last IC
            IC_px(h,i,k) = log(ee)/log(2);
        end
    end
end

% handles Eg
for h = 1:L
    for k = 1:O
        PRI_Eg(h,k) = OUT.Primal_Norms(h,k);    % primal solve
        ETE_Eg(h,k) = OUT.ETE_Norms(h,k);       % 1 ETE solve
        IC_Eg(h,k) = OUT.ETE_IC_Norms(h,k);     % Last IC
    end
end

% handles pg
for h = 2:L
    for k = 1:O
        ee = PRI_Eg(h-1,k)/PRI_Eg(h,k);         % primal solve
        PRI_pg(h,k) = log(ee)/log(2);
        ee = ETE_Eg(h-1,k)/ETE_Eg(h,k);         % 1 ETE solve
        ETE_pg(h,k) = log(ee)/log(2);
        ee = IC_Eg(h-1,k)/IC_Eg(h,k);           % Last IC
        IC_pg(h,k) = log(ee)/log(2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEBUG = struct();
DEBUG.L = L;
DEBUG.M = M;
DEBUG.O = O;
DEBUG.t = t;
DEBUG.nNodes = nNodes;
DEBUG.PRI_Ex = PRI_Ex;
DEBUG.PRI_Eg = PRI_Eg;
DEBUG.PRI_px = PRI_px;
DEBUG.ETE_pg = PRI_pg;
DEBUG.ETE_Ex = ETE_Ex;
DEBUG.ETE_Eg = ETE_Eg;
DEBUG.ETE_px = ETE_px;
DEBUG.ETE_pg = ETE_pg;
DEBUG.IC_Ex = IC_Ex;
DEBUG.IC_Eg = IC_Eg;
DEBUG.IC_px = IC_px;
DEBUG.IC_pg = IC_pg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE1 = struct();
SE1.L = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:O
SE1.L(ii).title='blah';
SE1.L(ii).variables = cell(1,4);
SE1.L(ii).zoneFmt='%0.2f';
SE1.L(ii).zoneVar=1;
SE1.L(ii).Nzones=1;
SE1.L(ii).dataFmt='%g';
SE1.L(ii).variables{1} = 'Nodes';
SE1.L(ii).variables{2} = 'Base DE: time + space';
SE1.L(ii).variables{3} = 'Corrected DE: time + space';
SE1.L(ii).variables{4} = 'Corrected DE: time + space (final correction step)';
SE1.L(ii).DATA.dat = [...
    nNodes(:),...
    PRI_Eg(:,ii),...
    ETE_Eg(:,ii),...
    IC_Eg(:,ii)...
    ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PE1 = struct();
PE1.L = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:O
SE1.L(ii).title='blah';
SE1.L(ii).variables = cell(1,4);
SE1.L(ii).zoneFmt='%0.2f';
SE1.L(ii).zoneVar=1;
SE1.L(ii).Nzones=1;
SE1.L(ii).dataFmt='%g';
SE1.L(ii).variables{1} = 'Nodes';
SE1.L(ii).variables{2} = 'Base DE: time + space';
SE1.L(ii).variables{3} = 'Corrected DE: time + space';
SE1.L(ii).variables{4} = 'Corrected DE: time + space (final correction step)';
SE1.L(ii).DATA.dat = [...
    nNodes(:),...
    PRI_pg(:,ii),...
    ETE_pg(:,ii),...
    IC_pg(:,ii)...
    ];
end
end