classdef BDF2_type < time_integrator_type
    properties
        N
        um1, um2
    end
    methods
        function this = BDF2_type( grid, soln, S )
            this = this@time_integrator_type(S.dt);
            this.N   = grid.N;
            this.um1 = soln.U;
            this.um2 = soln.U;
        end
        function [u_new,R,S,this] = step( this, u_old, S, LHS, RHS,...
                                             L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.BDF2_res( u, this.um1, this.um2,...
                                     this.dt, this.N, LHS, R_BC1, R_BC2 );
            dF = @(u) this.BDF2_jac( u,...
                                     this.dt, this.N, RHS, L_BC1, L_BC2 );
            gamma = 0.5;
            [u_new,R] = newton_with_backtracking( u_old, F, dF, gamma,...
                                  this.newton_tol, this.max_newton_iter );
            this.um2 = this.um1;
            this.um1 = u_new;
        end
        
    end
    methods (Static)
        function val = BDF2_res( u, um1, um2, dt, N, RHS, R_BC1, R_BC2 )
            val    = u - (4/3)*um1 + (1/3)*um2 +(2/3)*dt*RHS(u);
            val(1) = R_BC1(u,1); % apply RHS boundary conditions
            val(N) = R_BC2(u,N);
        end
        function val = BDF2_jac( u, dt, N, LHS, L_BC1, L_BC2 )
            val      = -(2/3)*dt*LHS(u);  % change to BDF2 LHS
            val(:,2) = val(:,2) + 1;
            val(1,:) = L_BC1(u,1); % apply LHS boundary conditions
            val(N,:) = L_BC2(u,N);
        end
    end
end