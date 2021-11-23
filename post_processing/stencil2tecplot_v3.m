function DAT = stencil2tecplot_v3(S,grid,stencil)
N = grid.N;
Ncorr = S.Niters+1;
M = S.stencil_size;

DAT = struct();
DAT.title = 'blah';
DAT.variables = cell(1);
DAT.zoneFmt='%g';
DAT.zoneVar=1;
DAT.Nzones=1;
DAT.dataFmt='%g';
DAT.variables{1} = 'x';
DAT.variables{2} = 't';
DAT.variables{3} = 'iter';
DAT.variables{4} = 'u-u<sub>exact</sub>';
DAT.dim = [N,M,Ncorr];
Uex = zeros(N,M);
for j = 1:M
    Uex(:,j) = S.ex_soln.eval(grid.x,stencil.t(j));
end
% for K = 1:Ncorr
%     for J = 1:M
%         for I = 1:N
%             DAT.DATA.dat(I,J,K,1) = grid.x(I);
%             DAT.DATA.dat(I,J,K,2) = stencil.t(J);
%             DAT.DATA.dat(I,J,K,3) = K-1;
%             DAT.DATA.dat(I,J,K,4) = abs(stencil.U(I,J,K)-Uex(I,J));
%         end
%     end
% end
for K = 1:Ncorr
    for J = 1:M
        for I = 1:N
            DAT.DATA.dat(I,J,K,1) = grid.x(I);
            DAT.DATA.dat(I,J,K,2) = stencil.t(J);
            DAT.DATA.dat(I,J,K,3) = K-1;
            DAT.DATA.dat(I,J,K,4) = abs(stencil.U(I,J,K)-Uex(I,J));
        end
    end
end

end