function g = my_savgol_end(m,np,dir)
if (mod(uint8(np-1),uint8(2))~=0)
    error('bad args in my_savgol_end');
end
nlr = (np-1)/2;
g = zeros(np,m+1,nlr);
if (dir == 1) % beginning
for i = 0:m
    for j = nlr+2:np
        nl = j-1;
        nr = np-j;
        c = my_savgol(np,nl,nr,i,m);
        c = circshift(c,nr);
        g(:,i+1,(np+1)-j) = c;
    end
end
elseif (dir == -1) % end
for i = 0:m
    for j = np:-1:nlr+2
        nl = np-j;
        nr = j-1;
        c = my_savgol(np,nl,nr,i,m);
        c = circshift(c,nr);
        g(:,i+1,j-(nlr+1)) = c;
    end
end
elseif (dir == 0) %center
    g = zeros(np,m+1);
    nl = (np-1)/2;
    nr = nl;
    for i = 0:m
        c = my_savgol(np,nl,nr,i,m);
        c = circshift(c,nl);
        g(:,i+1) = c;
    end
else
    error('dir must be -1,0, or 1')
end
% make same sign as MATLAB's sgolay
% g = g.*(-1).^(0:m);
end