function x = tridiag(a,b,c,d)
N = length(b);
x = zeros(N,1);
for i = 2:N
    w = a(i)/b(i-1);
    b(i) = b(i)-w*c(i-1);
    d(i) = d(i)-w*d(i-1);
end
x(N) = d(N)/b(N);
for i = N-1:-1:1
    x(i) = (d(i)-c(i)*x(i+1))/b(i);
end
% M = zeros(N);
% M(N+1:N+1:N*N-1) = a(2:N);
% M(1:N+1:N*N) = b;
% M(2:N+1:N*(N-1)) = c(1:N-1);
end