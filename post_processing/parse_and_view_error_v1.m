%% IC OOA analysis (simultaneous refinement in space + time)
clc; clear; close all;
% src_dirname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\unsteady_iterative_correction_results',...
%     '\iterative_weighting_experiments\'];
src_dirname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\post_processing\'];
src_fname = 'move_shock_ETEIC_newalg_mod';
load([src_dirname,src_fname]);

% fname='G:\My Drive\MATLAB\VT_Research\2021\shock_long_svd';
% load(fname);
%%
norms = 1;       % L-norm (3=infinity)
iters = 0:1;    % 0 corresponds to initial ETE solve
times = 1;

% separate primal and ETE info for easier handling
E_error = OUT.Error_Norms_E;
E_primal = OUT.Error_Norms_P;

% get dimensions for output arrays
L = size(E_error,2);        % (h) grid levels (space+time varying)
M = size(E_error(1,1).E,1); % (i) time levels (coarsest)
N = size(E_error(1,1).E,2); % (j) iteration levels
O = 3;                      % (k) norms

M2 = length(times);         % time levels for plotting
N2 = length(iters);         % iteration levels for plotting
O2 = length(norms);         % norms for plotting

Local_Errors = struct();    % ETE local errors
Local_Errors.X = struct();  % Primal local errors
Local_Errors.est_err = struct(); % estimated DE
Local_Errors.E_err = struct();   % error in error estimate
Local_Errors.P_err = struct();   % 'exact' DE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETE_Ex = zeros(L,M,N,O);         % ETE corrected DE norms
ETE_px = zeros(L,M,N,O);         % ETE corrected OOA

ETE_Ef = zeros(L,N,O);           % ETE final corrected DE norms
ETE_pf = zeros(L,N,O);           % ETE final corrected OOA

ETE_Eg = zeros(L,N,O);           % ETE global corrected DE norms
ETE_pg = zeros(L,N,O);           % ETE global corrected OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRI_Ex = zeros(L,M,O);           % Primal DE norms
PRI_px = zeros(L,M,O);           % Primal OOA

PRI_Ef = zeros(L,O);             % Primal final DE norms
PRI_pf = zeros(L,O);             % Primal final OOA

PRI_Eg = zeros(L,O);             % Primal global DE norms
PRI_pg = zeros(L,O);             % Primal global OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note: need to fix accumulation of error in time step calc
% these lines find the index of error at the correct time step
% to match the output interval
dx = OUT.dx;                % grid spacing
dt = OUT.dt;                % time steps



% handles Ex for ETE & PRI
for h = 1:L
    for i = 1:M
        for j = 1:N
            [~,index] = min( abs( OUT.Error_Norms_E(1).t(i)...
                - OUT.Error_Norms_E(h).t ) );
            for k = 1:O
                if j == N
                    % Conditional statement to include 50 iteration case
                    j2 = size(OUT.Error_Norms_E(h).E,2);
                    ETE_Ex(h,i,j,k) = OUT.Error_Norms_E(h).E(index,j2,k);
                else     
                    ETE_Ex(h,i,j,k) = OUT.Error_Norms_E(h).E(index,j,k);
                end
            end
        end
        [~,index2] = min( abs( OUT.Error_Norms_P(1).t(i)...
            - OUT.Error_Norms_P(h).t ) );
        for k = 1:O
            PRI_Ex(h,i,k) = OUT.Error_Norms_P(h).E(index,k);
        end
    end
end

% handles px for ETE & PRI
for h = 2:L
    for i = 1:M
        for j = 1:N
            for k = 1:O
                ee = ETE_Ex(h-1,i,j,k)/ETE_Ex(h,i,j,k);
                ETE_px(h,i,j,k) = log(ee)/log(2);
            end
        end
        for k = 1:O
            ee = PRI_Ex(h-1,i,k)/PRI_Ex(h,i,k);
            PRI_px(h,i,k) = log(ee)/log(2);
        end
    end
end

% handles Ef & Eg for ETE & PRI
for h = 1:L
    for j = 1:N
        for k = 1:O
            ETE_Ef(h,j,k) = OUT.Error_Norms_E(h).E(end,j,k);
            ETE_Eg(h,j,k) = OUT.Final_Enorm_E{h,1}(j,k);
        end
    end
    for k = 1:O
        PRI_Ef(h,k) = OUT.Error_Norms_P(h).E(end,k);
        PRI_Eg(h,k) = OUT.Final_Enorm_P(h,k);
    end
end

% handles pf & pg for ETE & PRI
for h = 2:L
    for j = 1:N
        for k = 1:O
            ee1 = ETE_Ef(h-1,j,k)/ETE_Ef(h,j,k);
            ETE_pf(h,j,k) = log(ee1)/log(2);
            ee2 = ETE_Eg(h-1,j,k)/ETE_Eg(h,j,k);
            ETE_pg(h,j,k) = log(ee2)/log(2);
        end
    end
    for k = 1:O
        ee1 = PRI_Ef(h-1,k)/PRI_Ef(h,k);
        PRI_pf(h,k) = log(ee1)/log(2);
        ee2 = PRI_Eg(h-1,k)/PRI_Eg(h,k);
        PRI_pg(h,k) = log(ee2)/log(2);
    end
end



% Plotting

figure(1);
subplot(1,2,1)
hold on
for i = 1:M2
    for k = 1:O2
        plot(dx,PRI_Ex(:,times(i),norms(k)),'k');
    end
    for j = 1:N2
        for k = 1:O2
            plot(dx,ETE_Ex(:,times(i),iters(j)+1,norms(k)),'r');
        end
    end
end
for k = 1:O2
%     plot(dx,PRI_Ef(:,norms(k)),'g');
    plot(dx,PRI_Eg(:,norms(k)),'m');
end
for j = 1:N2
    for k = 1:O2
%         plot(dx,ETE_Ef(:,iters(j)+1,norms(k)),'c');
        plot(dx,ETE_Eg(:,iters(j)+1,norms(k)),'b');
    end
end
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(1,2,2)
hold on
for i = 1:M2
    for k = 1:O2
        plot(dx(2:end),PRI_px(2:end,times(i),norms(k)),'k');
    end
    for j = 1:N2
        for k = 1:O2
            plot(dx(2:end),ETE_px(2:end,times(i),iters(j)+1,norms(k)),'r');
        end
    end
end
for k = 1:O2
%     plot(dx(2:end),PRI_pf(2:end,norms(k)),'g');
    plot(dx(2:end),PRI_pg(2:end,norms(k)),'m');
end
for j = 1:N2
    for k = 1:O2
%         plot(dx(2:end),ETE_pf(2:end,iters(j)+1,norms(k)),'c');
        plot(dx(2:end),ETE_pg(2:end,iters(j)+1,norms(k)));
    end
end
hold off
set(gca,'xscale','log')