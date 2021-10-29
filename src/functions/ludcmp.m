function [A,indx,d] = ludcmp(A)
% integer :: n, indx(n), NMAX
% real(prec) :: d, a(np,np), TINY
% integer :: i, imax, j, k
% real(prec) :: aamax, dum, sum, vv(NMAX)
n = size(A,1);
indx = zeros(n,1);
vv = zeros(n,1);
d = 1;
for i = 1:n
    aamax = 0;
    for j = 1:n
        if (abs(A(i,j)) > aamax)
            aamax = abs(A(i,j));
        end
    end
    if (aamax == 0)
        error('singular matrix in ludcmp')
    end
    vv(i)=1/aamax;
end
for j = 1:n
    for i = 1:j-1
        sum = A(i,j);
        for k = 1:i-1
            sum = sum - A(i,k)*A(k,j);
        end
        A(i,j) = sum;
    end
    aamax = 0;
    for i = j:n
        sum =  A(i,j);
        for k = 1:j-1
            sum = sum-A(i,k)*A(k,j);
        end
        A(i,j) = sum;
        dum = vv(i)*abs(sum);
        if (dum >= aamax)
            imax = i;
            aamax = dum;
        end
    end
    if (j ~= imax)
        for k = 1:n
            dum = A(imax,k);
            A(imax,k) = A(j,k);
            A(j,k) = dum;
        end
        d = -d;
        vv(imax) = vv(j);
    end
    indx(j) = imax;
    if (A(j,j) == 0)
        A(j,j) = 1.0e-12;
    end
    if (j ~= n)
        dum = 1/A(j,j);
        for i = j+1:n
            A(i,j) = A(i,j)*dum;
        end
    end
end
end