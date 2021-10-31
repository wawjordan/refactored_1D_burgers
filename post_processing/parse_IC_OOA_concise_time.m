function [SE,PE,DEBUG] = parse_IC_OOA_concise_time(filename)

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
                ETE_Ex(h,i,j,k) = OUT.Error_Norms_E(h).E(index,j,k);
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
            ETE_Ef(h,j,k) = OUT.Error_Norms_E(h).E(end,j,k);
            ETE_Eg(h,j,k) = OUT.Final_Enorm_E{h,1}(j,k);
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
SE = struct();
SE.L = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:O
SE.L(ii).title='blah';
SE.L(ii).variables = cell(1,2 + N);
SE.L(ii).zoneFmt='%0.2d';
SE.L(ii).zoneVar=nNodes;  % Zone matches mesh level
SE.L(ii).Nzones=L;
SE.L(ii).dataFmt='%g';
SE.L(ii).variables{1} = 'Time';
SE.L(ii).variables{2} = 'Base DE';
SE.L(ii).variables{3} = 'Corrected DE';
for j = 2:N
    k = j + 2;
    if j == 2
        SE.L(ii).variables{k} = 'Corrected DE: (step: #1)';
    elseif j == N
        SE.L(ii).variables{k} = sprintf('Corrected DE: (step: #%d)',N-1);
    else
        SE.L(ii).variables{k} = sprintf('Corrected DE: (steps: #%d - %d)',2,N-2);
    end
end
for h = 1:L
    tmp1 = zeros(M,N-1);
    for jj = 1:N
        tmp1(:,jj) = ETE_Ex(h,:,jj,ii);
    end
    SE.L(ii).DATA(h).dat = [...
        t(:),...
        PRI_Ex(h,:,ii)',...
        tmp1...
        ];
end
% for k = 2:M
%     tmp1 = zeros(L,N-1);
%     for jj = 1:N
%         tmp1(:,jj) = ETE_Ex(:,k,jj,ii);
%     end
%     SE.L(ii).DATA(k-1).dat = [...
%         nNodes(:),...
%         PRI_Ex(:,k,ii),...
%         tmp1...
%         ];
% end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PE = struct();
PE.L = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:O
PE.L(ii).title='blah';
PE.L(ii).variables = cell(1,2 + N);
PE.L(ii).zoneFmt='%0.2f';
PE.L(ii).zoneVar=round(100*t(2:end))/100;
PE.L(ii).Nzones=M-1;
PE.L(ii).dataFmt='%g';
PE.L(ii).variables{1} = 'Time';
PE.L(ii).variables{2} = 'Base DE';
PE.L(ii).variables{3} = 'Corrected DE';
for j = 2:N
    k = j + 2;
    if j == 2
        PE.L(ii).variables{k} = 'Corrected DE: (step: #1)';
    elseif j == N
        PE.L(ii).variables{k} = sprintf('Corrected DE: (step: #%d)',N-1);
    else
        PE.L(ii).variables{k} = sprintf('Corrected DE: (steps: #%d - %d)',2,N-2);
    end
end
for h = 1:L
    tmp1 = zeros(M,N-1);
    for jj = 1:N
        tmp1(:,jj) = ETE_px(h,:,jj,ii);
    end
    PE.L(ii).DATA(h).dat = [...
        t(:),...
        PRI_px(h,:,ii)',...
        tmp1...
        ];
end
% for k = 2:M
%     tmp1 = zeros(L,N-1);
%     for jj = 1:N
%         tmp1(:,jj) = ETE_px(:,k,jj,ii);
%     end
%     PE.L(ii).DATA(k-1).dat = [...
%         nNodes(:),...
%         PRI_px(:,k,ii),...
%         tmp1...
%         ];
% end

end



end