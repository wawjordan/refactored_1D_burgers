classdef SDIRK2_type < time_integrator_type
    properties
        N
        alpha
    end
    methods
        function this = SDIRK2_type( grid, ~, S )
            this = this@time_integrator_type(S.dt);
            this.N   = grid.N;
        end
        function [u_new,R,S,this] = step( this, u_old, S, RHS, LHS,...
                                             L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.SDIRK2_res1( u, u_old,...
                                     this.dt, this.N, RHS, R_BC1, R_BC2 );
            dF = @(u) this.SDIRK2_jac1( u,...
                                     this.dt, this.N, LHS, L_BC1, L_BC2 );
            gamma = 0.5;
            [u_alpha,R1] = newton_with_backtracking( u_old, F, dF, gamma,...
                                  this.newton_tol, this.max_newton_iter );
                              
            F  = @(u) this.SDIRK2_res2( u, u_old, u_alpha,...
                                     this.dt, this.N, RHS, R_BC1, R_BC2 );
            dF = @(u) this.SDIRK2_jac2( u,...
                                     this.dt, this.N, LHS, L_BC1, L_BC2 );
            [u_new,R] = newton_with_backtracking( u_alpha, F, dF, gamma,...
                                  this.newton_tol, this.max_newton_iter );                  
        end
        
    end
    methods (Static)
        function val = SDIRK2_res1( u, u_old, dt, N, RHS, R_BC1, R_BC2 )
            alpha = (2-sqrt(2))/2;
            val    = u - u_old + alpha*dt*RHS(u);
%             val(1) = R_BC1(u,1); % apply RHS boundary conditions
%             val(N) = R_BC2(u,N);
        end
        function val = SDIRK2_res2( u, u_old, u_alpha, dt, N, RHS, R_BC1, R_BC2 )
            alpha = (2-sqrt(2))/2;
            val    = u - u_old + (1-alpha)*dt*RHS(u_alpha) + alpha*dt*RHS(u);
            val(1) = R_BC1(u,1); % apply RHS boundary conditions (?)
            val(N) = R_BC2(u,N);
        end
        function val = SDIRK2_jac1( u, dt, N, LHS, L_BC1, L_BC2 )
            alpha = (2-sqrt(2))/2;
            val      = -alpha*dt*LHS(u);  % change to SDIRK2 LHS
            val(:,2) = val(:,2) + 1;
            val(1,:) = L_BC1(u,1); % apply LHS boundary conditions
            val(N,:) = L_BC2(u,N);
        end
        function val = SDIRK2_jac2( u, dt, N, LHS, L_BC1, L_BC2 )
            alpha = (2-sqrt(2))/2;
            val      = -alpha*dt*LHS(u);  % change to SDIRK2 LHS
            val(:,2) = val(:,2) + 1;
            val(1,:) = L_BC1(u,1); % apply LHS boundary conditions
            val(N,:) = L_BC2(u,N);
        end
    end
end