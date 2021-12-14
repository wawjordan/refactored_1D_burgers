function [xk,Rks] = modified_newton_with_backtracking(x0,F,dF,gamma,res_tol,max_iter)
Rks = nan(max_iter,1);     % allocate residual norm history
converged = false;         % initialize convergence booleans
diverged  = false;

k = 1;                     % initialize iteration counter
xk = x0;                   % initial guess
Fk = F(xk);                % evaluate function at x_k
Rinit = norm (Fk,2);       % residual norm at first iteration
Rks(k) = 1;
dFk = dF(xk);              % evaluate jacobian at x_k

while ( (~converged)&&(~diverged) )
  k = k + 1;               % advance iteration counter
  eta = 1.0;               % initialize eta
  dxk = tridiag(dFk(:,1),dFk(:,2),dFk(:,3),Fk); % solve for update
  xkp1 = xk - eta*dxk;     % calculate new x_k
  Fkp1 = F(xkp1);          % evaluate function at new x_k
  counter = 0;
  stalled = false;
  while (dot(Fkp1,Fkp1) > dot(Fk,Fk))&&(~stalled) % backtracking loop
      eta = gamma*eta;     % calculate new eta
      xkp1 = xk - eta*dxk; % calculate new x_k
      Fkp1 = F(xkp1);      % evaluate function at new x_k
      counter = counter+1;
      stalled = (counter > max_iter);
  end
  Rks(k)= norm(Fkp1,2)/Rinit;   % calculate relative residual reduction
  converged = (Rks(k)<res_tol); % update convergence booleans
  diverged  = ( k > max_iter );
  xk = xkp1;                    % x_k <-- x_k+1
  Fk = Fkp1;                    % f_k <-- f_k+1
end
Rks = Rks(~isnan(Rks));         % remove NaNs
end
