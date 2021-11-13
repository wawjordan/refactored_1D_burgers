function c = my_savgol(np,nl,nr,ld,m)
%{
% INTEGER ld,m,nl,np,nr,MMAX
% REAL c(np)
% PARAMETER (MMAX=6)
% INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
% REAL d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
% USES lubksb, ludcmp
  Returns in c(1:np), in wrap-around order (N.B.!) consistent with the 
argument respns in routine convlv, a set of Savitzky-Golay filter 
coefficients. 
  nl is the number of leftward (past) data points used, while 
nr is the number of rightward (future) data points, making the total 
number of data points used nl + nr + 1. 
  ld is the order of the derivative desired (e.g., ld = 0 for 
smoothed function). 
  m is the order of the smoothing polynomial, also equal to the highest
conserved moment; usual values are m = 2 or m = 4.
%}
c = zeros(np,1);
b = zeros(m+1,1);
A = zeros(m+1,m+1);
if (np < nl+nr+1)||(nl < 0)||(nr < 0)||(ld > m)||(nl+nr < m)
    error('bad args in savgol');
end
for ipj = 0:2*m % Set up normal equations of the desired least squares fit.
    sum = 0.0;
    if (ipj==0)
        sum = 1.0;
    end
    for k = 1:nr
        sum = sum + k^ipj;
    end
    for k = 1:nl
        sum = sum + (-k)^ipj;
    end
    mm = min(ipj,2*m-ipj);
    for imj = -mm:2:mm
        A(1+(ipj+imj)/2,1+(ipj-imj)/2) = sum;
    end
end
[A,indx,~] = ludcmp(A); % Solve them: LU decomposition.
for j = 1:m+1
    b(j) = 0.0;
end
% Right-hand side vector is unit vector, depending on which derivative.
b(ld+1)=1.0;
[b] = lubksb(A,indx,b); % Backsubstitute => one row of the inverse matrix.
%Zero the output array (it may be bigger than number of coefficients).
for kk = 1:np
    c(kk) = 0.0;
end
% Each Savitzky-Golay coefficient is the dot product of 
% powers of an integer with the inverse matrix row.
for k = -nl:nr
    sum = b(1);
    fac = 1.0;
    for mm = 1:m
        fac = fac*k;
        sum = sum + b(mm+1)*fac;
    end
    kk = mod(np-k,np)+1; % Store in wrap-around order.
    c(kk) = sum;
end

end