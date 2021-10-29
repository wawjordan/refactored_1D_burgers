function pd = ddpoly(c,x,nd,mu)
nx = length(x);

dxdx = 1/mu(2);

nc = length(c);
% x = x*mu(2) + mu(1);
x = (x-mu(1))/mu(2);
pd = zeros(nx,nd);
pd(:,1) = c(nc);
for i = nc-1:-1:1
    nnd = min(nd,nc+1-i);
    for j = nnd:-1:2
        pd(:,j) = pd(:,j).*x + pd(:,j-1);
    end
    pd(:,1) = pd(:,1).*x + c(i);
end
const = dxdx;
for i = 2:nd
    pd(:,i) = const*pd(:,i);
    const = const*i*dxdx;
end

% pd(:,2) = pd(:,2)*dxdx;
% const = 2*dxdx*dxdx;
% for i = 3:nd
%     pd(:,i) = const*pd(:,i);
%     const = const*i*dxdx;
% end

end