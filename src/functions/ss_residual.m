function R = ss_residual(u,dx,nu,N)
R = zeros(N,1);
for i = 2:N-1
    R(i) = (nu/(dx(i)^2))*(u(i+1) - 2*u(i) + u(i-1)) ...
           - (0.25/dx(i))*(u(i+1)^2 - u(i-1)^2);
end
end