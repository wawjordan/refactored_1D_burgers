function [xk,Rks] = broyden_with_backtracking(x0,F,J0,gamma,res_tol,max_iter)
Rks = nan(max_iter,1);     % allocate residual norm history
converged = false;         % initialize convergence booleans
diverged  = false;
N = length(x0);
k = 1;                     % initialize iteration counter
xk = x0;                   % initial guess
Fk = F(xk);                % evaluate function at x_k
Rinit = norm (Fk,2);       % residual norm at first iteration
Rks(k) = 1;
Bk = J0;                   % initialize approximate Jacobian

while ( (~converged)&&(~diverged) )
    k = k + 1;               % advance iteration counter
    eta = 1.0;               % initialize eta
    dxk = tridiag(Bk(:,1),Bk(:,2),Bk(:,3),Fk); % solve for update
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
    sk = xkp1 - xk;
    yk = Fkp1 - Fk;
    
    a        = yk       - Bk(:    ,2).*sk;
    a(2:N)   = a(2:N)   - Bk(2:N  ,1).*sk(1:N-1);
    a(1:N-1) = a(1:N-1) - Bk(1:N-1,3).*sk(2:N);
    
    den = [0;sk(1:end-1)].^2 + sk.^2 + [sk(2:end);0].^2;
    den(den<=eps(1)) = 0;
    den(den >eps(1)) = 1./den(den>eps(1));
    
    Bk(2:N-1,1) = Bk(2:N-1,1) + a(2:N-1).*sk(1:N-2).*den(2:N-1);
    Bk(2:N-1,2) = Bk(2:N-1,2) + a(2:N-1).*sk(2:N-1).*den(2:N-1);
    Bk(2:N-1,3) = Bk(2:N-1,3) + a(2:N-1).*sk(3:N  ).*den(2:N-1);
    
    Rks(k)= norm(Fkp1,2)/Rinit;   % calculate relative residual reduction
    converged = (Rks(k)<res_tol); % update convergence booleans
    diverged  = ( k > max_iter );
    xk = xkp1;                    % x_k <-- x_k+1
    Fk = Fkp1;                    % f_k <-- f_k+1
end
Rks = Rks(~isnan(Rks));         % remove NaNs
end
