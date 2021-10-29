function [b] = lubksb(A,indx,b)

n = size(A,1);
ii = 0;
for i = 1:n
  LL = indx(i);
  sum = b(LL);
  b(LL) = b(i);
  if (ii ~= 0)
      for j = ii:i-1
          sum = sum - A(i,j)*b(j);
      end
  elseif (sum ~= 0)
      ii = i;
  end
  b(i) = sum;
end
for i = n:-1:1
  sum = b(i);
  for j = i+1:n
      sum = sum - A(i,j)*b(j);
  end
  b(i) = sum/A(i,i);
end

end