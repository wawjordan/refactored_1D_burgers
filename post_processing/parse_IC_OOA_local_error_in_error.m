function [SE] = parse_IC_OOA_local_error_in_error(filename)

load(filename,'OUT');

% get dimensions for output arrays
L = length(OUT.Nx);                   % (h) grid levels (space+time comb.)
N = size(OUT.Error_Norms_E(1,1).E,2); % (j) iteration levels
nNodes = OUT.Nx';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE.title='blah';
SE.variables = cell(1,2 + N);
SE.zoneFmt='%0.2d';
SE.zoneVar=nNodes(2:L);  % Zone matches mesh level
SE.Nzones=L-1;
SE.dataFmt='%g';
SE.variables{1} = 'x';
SE.variables{2} = 'Richardson Extrapolation';
SE.variables{3} = 'Estimated Error (ETE)';
for j = 2:N
    k = j + 2;
    if j == 2
        SE.variables{k} = 'Estimated Error (step: #1)';
    elseif j == N
        SE.variables{k} = sprintf('Estimated Error (step: #%d)',N-1);
    else
        SE.variables{k} = sprintf('Estimated Error (steps: #%d - %d)',2,N-2);
    end
end
for h = 2:L
E_RICH = -( OUT.Local_Error_P(h).u{end}(1:2:end)-OUT.Local_Error_P(h-1).u{end})/(2^2-1);
E_RICH = E_RICH-OUT.Local_Error_P(h).E{end}(1:2:end);
E_RICH = interp1(OUT.Local_Error_P(h-1).x(:),E_RICH,OUT.Local_Error_P(end).x(:));
% E_EX = interp1(OUT.Local_Error_P(h).x(:),OUT.Local_Error_P(h).E{end}(:),OUT.Local_Error_P(end).x(:));
% E_RICH = E_RICH-E_EX;
tmp1 = zeros(nNodes(L),N);
for jj = 1:N
    tmp1(:,jj) = interp1(OUT.Local_Error_P(h).x(:),OUT.Local_Error_E(h).Ee{end,jj},OUT.Local_Error_P(end).x(:));
end

SE.DATA(h-1).dat = [...
    OUT.Local_Error_P(end).x(:),...
    E_RICH,...
    tmp1,...
    ];
end
end