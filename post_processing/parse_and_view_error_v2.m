%% Parse & View Error V2 (02/11/2022)
clc; clear; %close all;
src_dirname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\results',...
    '\02_11\'];
src_fname = 'unsteady_shock_GS10';
load([src_dirname,src_fname]);

%%
norms = 3;       % L-norm (3=infinity)
times = 1;

% separate primal and ETE info for easier handling
E_error  = OUT.ETE_History;
E_primal = OUT.Primal_History;

% get dimensions for output arrays
L = size(E_error,2);        % (h) # of grid levels (space+time varying)
M = size(E_error(1).E,1);   % (i) # of time steps  (coarsest mesh)
O = 3;                      % (j) norms

M2 = length(times);         % time steps for plotting
O2 = length(norms);         % norms for plotting
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
dx = OUT.dx;                % grid spacing
dt = OUT.dt;                % time steps
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



% Plotting

figure(1);
subplot(1,2,1)
hold on
for i = 1:M2
    for k = 1:O2
        plot(dx,PRI_Ex(:,times(i),norms(k)),'k');
        plot(dx,ETE_Ex(:,times(i),norms(k)),'b');
        plot(dx, IC_Ex(:,times(i),norms(k)),'r');
    end
end
for k = 1:O2
    plot(dx,PRI_Eg(:,norms(k)),'m');
    plot(dx,ETE_Eg(:,norms(k)),'g');
    plot(dx, IC_Eg(:,norms(k)),'c');
end
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(1,2,2)
hold on
for i = 1:M2
    for k = 1:O2
        plot(dx(2:end),PRI_px(2:end,times(i),norms(k)),'k');
        plot(dx(2:end),ETE_px(2:end,times(i),norms(k)),'b');
        plot(dx(2:end), IC_px(2:end,times(i),norms(k)),'r');
    end
end
for k = 1:O2
    plot(dx(2:end),PRI_pg(2:end,norms(k)),'m');
    plot(dx(2:end),ETE_pg(2:end,norms(k)),'g');
    plot(dx(2:end), IC_pg(2:end,norms(k)),'c');
end
hold off
set(gca,'xscale','log')