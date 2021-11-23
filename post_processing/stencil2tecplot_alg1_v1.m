function DAT = stencil2tecplot_alg1_v1(DAT,S,grid,stencil,p,iter,start)
N = grid.N;
Ncorr = S.Niters+1;
M = S.stencil_size;
if (start)
%     DAT = struct();
    DAT.title = 'blah';
    DAT.variables = cell(M+1,1);
    DAT.zoneFmt='T = %g';
    DAT.zoneVar=stencil.t(end);
    DAT.Nzones=M;
    DAT.dataFmt='%g';
    DAT.variables{1} = 'iter';
    for i = 1:M
        if (i < M-1)
            DAT.variables{1+i} = sprintf('n-%d',(M-(i+1)));
        elseif (i == M-1)
            DAT.variables{1+i} = 'n';
        else
            DAT.variables{1+i} = 'n+1';
        end
    end
    DAT.dim = [Ncorr];
    DAT.aux = true;
    DAT.auxdata(1).name  = 'xloc';
    DAT.auxdata(1).value = sprintf('%0.4f',grid.x(p));
    DAT.auxdata(2).name  = 'nodeNum';
    DAT.auxdata(2).value = sprintf('%d',p);
end
Uex = zeros(M,1);
for j = 1:M
    Uex(j) = S.ex_soln.eval(grid.x(p),stencil.t(j));
end
DAT.DATA.dat(iter,1) = iter-1;
for K = 1:M
    DAT.DATA.dat(iter,1+K) = abs(stencil.U(p,K,2)-Uex(K));
end

% for J = 1:Ncorr
%     DAT.DATA.dat(J,1) = J-1;
%     for K = 1:M
%         DAT.DATA.dat(J,1+K) = abs(stencil.U(p,K,J)-Uex(K));
%     end
% end

end