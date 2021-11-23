function DAT = stencil2tecplot_v2(S,grid,stencil)
N = grid.N;
Ncorr = S.Niters+1;
M = S.stencil_size;

DAT = struct();
DAT.title = 'blah';
DAT.variables = cell(1);
DAT.zoneFmt='%g';
DAT.zoneVar=stencil.t;
DAT.Nzones=M;
DAT.dataFmt='%g';
DAT.variables{1} = 'x';
DAT.variables{2} = 'iter';
DAT.variables{3} = 'u-u<sub>exact</sub>';
DAT.dim = [N,Ncorr];
Uex = zeros(N,M);
for j = 1:M
    Uex(:,j) = S.ex_soln.eval(grid.x,stencil.t(j));
end
for K = 1:M
    for J = 1:Ncorr
        for I = 1:N
            DAT.DATA(K).dat(I,J,1) = grid.x(I);
            DAT.DATA(K).dat(I,J,2) = J-1;
            DAT.DATA(K).dat(I,J,3) = abs(stencil.U(I,K,J)-Uex(I,K));
        end
    end
end

end