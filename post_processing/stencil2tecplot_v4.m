function DAT = stencil2tecplot_v4(S,grid,stencil)
N = grid.N;
Ncorr = S.Niters+1;
M = S.stencil_size;

DAT = struct();
DAT.title = 'blah';
DAT.variables = cell(M+2,1);
DAT.zoneFmt='T = %g';
DAT.zoneVar=stencil.t(end);
DAT.Nzones=M;
DAT.dataFmt='%g';
DAT.variables{1} = 'x';
DAT.variables{2} = 'iter';
for i = 1:M
    if (i < M-1)
        DAT.variables{2+i} = sprintf('n-%d',(M-(i+1)));
    elseif (i == M-1)
        DAT.variables{2+i} = 'n';
    else
        DAT.variables{2+i} = 'n+1';
    end
end
DAT.dim = [N,Ncorr];
Uex = zeros(N,M);
for j = 1:M
    Uex(:,j) = S.ex_soln.eval(grid.x,stencil.t(j));
end
for J = 1:Ncorr
    for I = 1:N
        DAT.DATA.dat(I,J,1) = grid.x(I);
        DAT.DATA.dat(I,J,2) = J-1;
        for K = 1:M
            DAT.DATA.dat(I,J,2+K) = abs(stencil.U(I,K,J)-Uex(I,K));
        end
    end
end

end