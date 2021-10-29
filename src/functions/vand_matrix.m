function [A,mu] = vand_matrix(x,n)
x = x(:);
m = length(x);
mu = [mean(x),std(x)];
xhat = (x-mu(1))/mu(2); % centering and scaling

A = ones(m,n+1);
for i = 1:n
    A(:,i+1) = A(:,i).*xhat;
end

end