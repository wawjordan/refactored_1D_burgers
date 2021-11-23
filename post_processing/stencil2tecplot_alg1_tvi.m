function DAT = stencil2tecplot_alg1_tvi(DAT,S,grid,stencil,p,iter,start)
Ncorr = S.Niters+1;
M = S.stencil_size;
if (start)
%     DAT = struct();
    DAT.title = 'blah';
    DAT.variables = cell(Ncorr+1,1);
    DAT.zoneFmt='T = %g';
    DAT.zoneVar=stencil.t(end);
    DAT.Nzones=M;
    DAT.dataFmt='%g';
    DAT.variables{1} = 't';
    DAT.variables{2} = 'Primal';
    for j = 2:Ncorr
        DAT.variables{j+1} = sprintf('Corrected: (step: #%d)',j-1);
    end
    DAT.dim = M;
    DAT.customlabels = cell(1,M);
    for i = 1:M
        if (i < M-1)
            DAT.customlabels{i} = sprintf('n-%d',(M-(i+1)));
        elseif (i == M-1)
            DAT.customlabels{i} = 'n';
        else
            DAT.customlabels{i} = 'n+1';
        end
    end
    DAT.aux = true;
    DAT.auxdata(1).name  = 'xloc';
    DAT.auxdata(1).value = sprintf('%0.4f',grid.x(p));
    DAT.auxdata(2).name  = 'nodeNum';
    DAT.auxdata(2).value = sprintf('%d',p);
    DAT.cust = true;
    for j = 1:M
        DAT.DATA.dat(j,1) = j;
    end
    Uex = zeros(M,1);
    for j = 1:M
        Uex(j) = S.ex_soln.eval(grid.x(p),stencil.t(j));
    end
    % primal solution stencil
    for K = 1:M % time in stencil
        DAT.DATA.dat(K,2) = abs(stencil.U(p,K,1)-Uex(K));
    end
else
    Uex = zeros(M,1);
    for j = 1:M
        Uex(j) = S.ex_soln.eval(grid.x(p),stencil.t(j));
    end
    % corrected solution stencil
    for K = 1:M % time in stencil
        DAT.DATA.dat(K,iter+1) = abs(stencil.U(p,K,2)-Uex(K));
    end
end

end