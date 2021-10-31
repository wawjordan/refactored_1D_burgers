function [SE1,SE2,PE1,PE2,DEBUG] = parse_IC_OOA_concise(filename)

load(filename,'OUT');

% get dimensions for output arrays
L = size(OUT.Error_Norms_E,2);        % (h) grid levels (space+time comb.)
M = size(OUT.Error_Norms_E(1,1).E,1); % (i) time levels (coarsest)
N = size(OUT.Error_Norms_E(1,1).E,2); % (j) iteration levels
O = 3;                                % (k) norms

t = OUT.Error_Norms_E(1,1).t(:);

nNodes = OUT.Nx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRI_Ex = zeros(L,M,O);           % Primal DE norms
PRI_px = zeros(L,M,O);           % Primal OOA

PRI_Ef = zeros(L,O);             % Primal final DE norms
PRI_pf = zeros(L,O);             % Primal final OOA

PRI_Eg = zeros(L,O);             % Primal global DE norms
PRI_pg = zeros(L,O);             % Primal global OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETE_Ex = zeros(L,M,N,O);         % ETE corrected DE norms
ETE_px = zeros(L,M,N,O);         % ETE corrected OOA

ETE_Ef = zeros(L,N,O);           % ETE final corrected DE norms
ETE_pf = zeros(L,N,O);           % ETE final corrected OOA

ETE_Eg = zeros(L,N,O);           % ETE global corrected DE norms
ETE_pg = zeros(L,N,O);           % ETE global corrected OOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% handles Ex for ETE & PRI
for h = 1:L
    for i = 1:M
        for j = 1:N
            [~,index] = min( abs( OUT.Error_Norms_E(1).t(i)...
                - OUT.Error_Norms_E(h).t ) );
            for k = 1:O
                if j == N
                    %% Conditional statement to include 50 iteration case
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
            PRI_Ex(h,i,k) = OUT.Error_Norms_P(h).E(index2,k);
        end
    end
end

%% handles px for ETE & PRI
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

%% handles Ef & Eg for ETE & PRI
for h = 1:L
    for j = 1:N
        for k = 1:O
            if j == N
                %% Conditional statement to include 50 iteration case
                j2 = size(OUT.Error_Norms_E(h).E,2);
                ETE_Ef(h,j,k) = OUT.Error_Norms_E(h).E(end,j2,k);
                ETE_Eg(h,j,k) = OUT.Final_Enorm_E{h,1}(j2,k);
                else     
            ETE_Ef(h,j,k) = OUT.Error_Norms_E(h).E(end,j,k);
            ETE_Eg(h,j,k) = OUT.Final_Enorm_E{h,1}(j,k);
            end
        end
    end
    for k = 1:O
        PRI_Ef(h,k) = OUT.Error_Norms_P(h).E(end,k);
        PRI_Eg(h,k) = OUT.Final_Enorm_P(h,k);
    end
end

%% handles pf & pg for ETE & PRI
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEBUG = struct();
DEBUG.L = L;
DEBUG.M = M;
DEBUG.N = N;
DEBUG.O = O;
DEBUG.t = t;
DEBUG.nNodes = nNodes;
DEBUG.PRI_Ex = PRI_Ex;
DEBUG.PRI_Ef = PRI_Ef;
DEBUG.PRI_Eg = PRI_Eg;
DEBUG.PRI_px = PRI_px;
DEBUG.PRI_pf = PRI_pf;
DEBUG.ETE_pg = PRI_pg;
DEBUG.ETE_Ex = ETE_Ex;
DEBUG.ETE_Ef = ETE_Ef;
DEBUG.ETE_Eg = ETE_Eg;
DEBUG.ETE_px = ETE_px;
DEBUG.ETE_pf = ETE_pf;
DEBUG.ETE_pg = ETE_pg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE1 = struct();
SE1.L = struct();

SE2 = struct();
SE2.L = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:O
SE1.L(ii).title='blah';
SE1.L(ii).variables = cell(1,3 + 2*N);
SE1.L(ii).zoneFmt='%0.2f';
SE1.L(ii).zoneVar=1;
SE1.L(ii).Nzones=1;
SE1.L(ii).dataFmt='%g';
SE1.L(ii).variables{1} = 'Nodes';
SE1.L(ii).variables{2} = 'Base DE: time + space';
SE1.L(ii).variables{3} = 'Base DE: final time step';
SE1.L(ii).variables{4} = 'Corrected DE: time + space';
SE1.L(ii).variables{5} = 'Corrected DE: final time step';
for j = 2:N
    k = j + 4;
    if j == 2
        SE1.L(ii).variables{k} = 'Corrected DE: time + space (step: #1)';
    elseif j == N
        SE1.L(ii).variables{k} = sprintf('Corrected DE: time + space (step: #%d)',N-1);
    else
        SE1.L(ii).variables{k} = sprintf('Corrected DE: time + space (steps: #%d - %d)',2,N-2);
    end
end
for j = 2:N
    k = j + 3 + N;
    if j == 2
        SE1.L(ii).variables{k} = 'Corrected DE: final time step (step: #1)';
    elseif j == N
        SE1.L(ii).variables{k} = sprintf('Corrected DE: final time step (step: #%d)',N-1);
    else
        SE1.L(ii).variables{k} = sprintf('Corrected DE: final time step (steps: #%d - %d)',2,N-2);
    end
end
SE1.L(ii).DATA.dat = [...
    nNodes(:),...
    PRI_Eg(:,ii),...
    PRI_Ef(:,ii),...
    ETE_Eg(:,1,ii),...
    ETE_Ef(:,1,ii),...
    ETE_Eg(:,2:end,ii),...
    ETE_Ef(:,2:end,ii)...
    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE2.L(ii).title='blah';
SE2.L(ii).variables = cell(1,2 + N);
SE2.L(ii).zoneFmt='%0.2f';
SE2.L(ii).zoneVar=round(100*t)/100;
SE2.L(ii).Nzones=M-1;
SE2.L(ii).dataFmt='%g';
SE2.L(ii).variables{1} = 'Nodes';
SE2.L(ii).variables{2} = 'Base DE';
SE2.L(ii).variables{3} = 'Corrected DE';
for j = 2:N
    k = j + 2;
    if j == 2
        SE2.L(ii).variables{k} = 'Corrected DE: (step: #1)';
    elseif j == N
        SE2.L(ii).variables{k} = sprintf('Corrected DE: (step: #%d)',N-1);
    else
        SE2.L(ii).variables{k} = sprintf('Corrected DE: (steps: #%d - %d)',2,N-2);
    end
end
for k = 2:M
    tmp1 = zeros(L,N-1);
    for jj = 1:N
        tmp1(:,jj) = ETE_Ex(:,k,jj,ii);
    end
    SE2.L(ii).DATA(k-1).dat = [...
        nNodes(:),...
        PRI_Ex(:,k,ii),...
        tmp1...
        ];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PE1 = struct();
PE1.L = struct();

PE2 = struct();
PE2.L = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:O
PE1.L(ii).title='blah';
PE1.L(ii).variables = cell(1,3 + 2*N);
PE1.L(ii).zoneFmt='%0.2f';
PE1.L(ii).zoneVar=1;
PE1.L(ii).Nzones=1;
PE1.L(ii).dataFmt='%g';
PE1.L(ii).variables{1} = 'Nodes';
PE1.L(ii).variables{2} = 'Base DE: time + space';
PE1.L(ii).variables{3} = 'Base DE: final time step';
PE1.L(ii).variables{4} = 'Corrected DE: time + space';
PE1.L(ii).variables{5} = 'Corrected DE: final time step';
for j = 2:N
    k = j + 4;
    if j == 2
        PE1.L(ii).variables{k} = 'Corrected DE: time + space (step: #1)';
    elseif j == N
        PE1.L(ii).variables{k} = sprintf('Corrected DE: time + space (step: #%d)',N-1);
    else
        PE1.L(ii).variables{k} = sprintf('Corrected DE: time + space (steps: #%d - %d)',2,N-2);
    end
end
for j = 2:N
    k = j + 3 + N;
    if j == 2
        PE1.L(ii).variables{k} = 'Corrected DE: final time step (step: #1)';
    elseif j == N
        PE1.L(ii).variables{k} = sprintf('Corrected DE: final time step (step: #%d)',N-1);
    else
        PE1.L(ii).variables{k} = sprintf('Corrected DE: final time step (steps: #%d - %d)',2,N-2);
    end
end
PE1.L(ii).DATA.dat = [...
    nNodes(:),...
    PRI_pg(:,ii),...
    PRI_pf(:,ii),...
    ETE_pg(:,1,ii),...
    ETE_pf(:,1,ii),...
    ETE_pg(:,2:end,ii),...
    ETE_pf(:,2:end,ii)...
    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PE2.L(ii).title='blah';
PE2.L(ii).variables = cell(1,2 + N);
PE2.L(ii).zoneFmt='%0.2f';
PE2.L(ii).zoneVar=round(100*t(2:end))/100;
PE2.L(ii).Nzones=M-1;
PE2.L(ii).dataFmt='%g';
PE2.L(ii).variables{1} = 'Nodes';
PE2.L(ii).variables{2} = 'Base DE';
PE2.L(ii).variables{3} = 'Corrected DE';
for j = 2:N
    k = j + 2;
    if j == 2
        PE2.L(ii).variables{k} = 'Corrected DE: (step: #1)';
    elseif j == N
        PE2.L(ii).variables{k} = sprintf('Corrected DE: (step: #%d)',N-1);
    else
        PE2.L(ii).variables{k} = sprintf('Corrected DE: (steps: #%d - %d)',2,N-2);
    end
end
for k = 2:M
    tmp1 = zeros(L,N-1);
    for jj = 1:N
        tmp1(:,jj) = ETE_px(:,k,jj,ii);
    end
    PE2.L(ii).DATA(k-1).dat = [...
        nNodes(:),...
        PRI_px(:,k,ii),...
        tmp1...
        ];
end

end



end