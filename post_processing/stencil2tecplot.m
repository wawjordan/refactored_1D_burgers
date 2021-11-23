function [UDAT,EDAT,EXDAT] = stencil2tecplot(S,grid,stencil)
N = grid.N;
Ncorr = S.Niters+1;
M = S.stencil_size;
t = stencil.t(:);
ind = 1:M;

UDAT = struct();
UDAT.title = 'blah';
UDAT.variables = cell(1);
UDAT.zoneFmt='%0.2d';
UDAT.zoneVar=ind;
UDAT.Nzones=M;
UDAT.dataFmt='%g';
UDAT.variables{1} = 'x';
UDAT.variables{2} = 'Exact Solution';
UDAT.variables{3} = 'Primal Solution';
for j = 2:Ncorr
    UDAT.variables{j+2} = sprintf('Corrected Solution: (step: #%d)',j-1);
end
for j = 1:M
    Uex = S.ex_soln.eval(grid.x,stencil.t(j));
    UDAT.DATA(j).dat = [grid.x(:),Uex(:),reshape(stencil.U(:,j,:),N,Ncorr)];
end

EDAT = struct();
EDAT.title = 'blah';
EDAT.variables = cell(1);
EDAT.zoneFmt='%0.2d';
EDAT.zoneVar=ind;
EDAT.Nzones=M;
EDAT.dataFmt='%g';
EDAT.variables{1} = 'x';
EDAT.variables{2} = 'Primal Solution';
for j = 2:Ncorr
    EDAT.variables{j+1} = sprintf('Corrected Solution: (step: #%d)',j-1);
end
for j = 1:M
    Uex = S.ex_soln.eval(grid.x,stencil.t(j));
    EDAT.DATA(j).dat = [grid.x(:),reshape(stencil.U(:,j,:),N,Ncorr)-Uex(:)];
end

EXDAT = struct();
EXDAT.title = 'blah';
EXDAT.variables = cell(1);
EXDAT.zoneFmt='%g';
EXDAT.zoneVar=grid.x;
EXDAT.Nzones=N;
EXDAT.dataFmt='%g';
EXDAT.variables{1} = 'iter';
for j = 1:M
    EXDAT.variables{j+1} = sprintf('Position #%d',j);
end
Uex = zeros(N,M);
for j = 1:M
    Uex(:,j) = S.ex_soln.eval(grid.x,stencil.t(j));
end
for j = 1:N
    EXDAT.DATA(j).dat = [(0:Ncorr-1)',abs(reshape(stencil.U(j,:,:),M,Ncorr)'-Uex(j,:))];
end



end