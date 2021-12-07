clc; clear; close all;
filename = 'shock_exact_start_alg1_allsave';
load(filename);

NX    = OUT.Nx;
N_grids = length(NX);
Nt = length(OUT.Local_Error_P(1).t);
E_RICH = struct();
for j = 2:N_grids
    E_RICH(j).E = cell(Nt,1);
    for i = 1:Nt
        s1 = 2^(j-2); % 1, 2, 4, ...
        s2 = 2*s1;    % 2, 4, 8, ...
        k1 = s1*(i-1)+1; % 1, 2, 3, 4, ...
        k2 = s2*(i-1)+1; % 1, 3, 5, 7, ...
%         fprintf('%d, %d, %d %d\n',i,j,k1,k2);
        E_RICH(j).E{i} = (OUT.Local_Error_P(j).u{k2}(1:2:end)-OUT.Local_Error_P(j-1).u{k1})/(2^2-1);
    end
end