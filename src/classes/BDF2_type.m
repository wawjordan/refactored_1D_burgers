classdef BDF2_type < time_integrator_type
    properties
        N
        um1, um2
        method
        gamma
    end
    methods
        function this = BDF2_type( grid, soln, S, varargin )
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x)&&isscalar(x)&&(x > 0);
            expectedMethods = {'newton','newton_mod','broyden',...
                               'broyden-lin','broyden-I'};
            defaultMethod = 'newton';
            defaultGamma  = 0.5;
            addParameter(p,'gamma',defaultGamma,validScalarPosNum);
            addOptional(p,'method',defaultMethod,...
                      @(x) any(validatestring(x,expectedMethods)));
            addRequired(p,'grid');
            addRequired(p,'soln');
            addRequired(p,'S');
            this = this@time_integrator_type(S.dt);
            parse(p,grid,soln,S,varargin{:});
            this.N   = grid.N;
            this.um1 = soln.U;
            this.um2 = soln.U;
            this.gamma = p.Results.gamma;
            switch(p.Results.method)
               case 'newton'
                   this.method = @newton_step;
               case 'newton_mod'
                   this.method = @newton_step_mod;
               case 'broyden'
                   this.method = @broyden_step;
               case 'broyden-lin'
                   this.method = @broyden_step_lin;
               case 'broyden-I'
                   this.method = @broyden_step_I;
            end
        end
        function [u_new,R,S,this] =...
               step( this, u_old, S, RHS, LHS, L_BC1, L_BC2, R_BC1, R_BC2 )
            [u_new,R,S,this] = this.method( this, u_old, S, RHS, LHS,...
                                            L_BC1, L_BC2, R_BC1, R_BC2 );
        end
        function [u_new,R,S,this] =...
            newton_step( this, u_old, S, RHS, LHS,...
                               L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.BDF2_res( u, this.um1, this.um2, this.dt,...
                                        this.N, RHS, R_BC1, R_BC2 );
            dF = @(u) this.BDF2_jac( u, this.dt, this.N,...
                                        LHS, L_BC1, L_BC2 );
            [u_new,R] = newton_with_backtracking(...
                              u_old, F, dF, this.gamma,...
                              this.newton_tol, this.max_newton_iter );
            this.um2 = this.um1;
            this.um1 = u_new;
        end
        function [u_new,R,S,this] =...
            newton_step_mod( this, u_old, S, RHS, LHS,...
                                   L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.BDF2_res( u, this.um1, this.um2, this.dt,...
                                        this.N, RHS, R_BC1, R_BC2 );
            dF = @(u) this.BDF2_jac( u, this.dt, this.N,...
                                        LHS, L_BC1, L_BC2 );
            [u_new,R] = modified_newton_with_backtracking(...
                                  u_old, F, dF, this.gamma,...
                                  this.newton_tol, this.max_newton_iter );
            this.um2 = this.um1;
            this.um1 = u_new;
        end
        function [u_new,R,S,this] =...
            broyden_step( this, u_old, S, RHS, LHS,...
                                L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.BDF2_res( u, this.um1, this.um2, this.dt,...
                                        this.N, RHS, R_BC1, R_BC2 );
            J0 = this.BDF2_jac( u_old,this.dt, this.N, LHS, L_BC1, L_BC2 );
            [u_new,R] = broyden_with_backtracking(...
                                  u_old, F, J0, this.gamma,...
                                  this.newton_tol, this.max_newton_iter );
            this.um2 = this.um1;
            this.um1 = u_new;
        end
        function [u_new,R,S,this] =...
            broyden_step_lin( this, u_old, S, RHS, ~ ,...
                                    L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.BDF2_res( u, this.um1, this.um2, this.dt,...
                                        this.N, RHS, R_BC1, R_BC2 );
            J0 = this.BDF2_lin_jac( u_old, this.dt, this.N,...
                                           S, L_BC1, L_BC2 );
            [u_new,R] = broyden_with_backtracking(...
                                  u_old, F, J0, this.gamma,...
                                  this.newton_tol, this.max_newton_iter );
            this.um2 = this.um1;
            this.um1 = u_new;
        end
        function [u_new,R,S,this] =...
            broyden_step_I( this, u_old, S, RHS, ~ ,...
                                    L_BC1, L_BC2, R_BC1, R_BC2 )
            F  = @(u) this.BDF2_res( u, this.um1, this.um2, this.dt,...
                                        this.N, RHS, R_BC1, R_BC2 );
            J0 = zeros(this.N,3); J0(:,2) = 1;
            J0(1,:) = L_BC1(u_old,1); % apply LHS boundary conditions
            J0(this.N,:) = L_BC2(u_old,this.N);
            [u_new,R] = broyden_with_backtracking(...
                                  u_old, F, J0, this.gamma,...
                                  this.newton_tol, this.max_newton_iter );
            this.um2 = this.um1;
            this.um1 = u_new;
        end
        function this = integrator_reset(this)
            this.um2 = zeros(this.N,1);
            this.um1 = this.um2;
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
        function val = BDF2_lin_jac( u, dt, N, S, L_BC1, L_BC2 )
            nu = S.ex_soln.nu;
            Uref = S.ex_soln.Uref;
            dx = S.dx;
            val = zeros(N,3);
            val(:,1) = -(2/3)*dt*( nu/dx^2 + 0.5*Uref/dx );
            val(:,2) = 1 -(2/3)*dt*( -2*nu/dx^2 );
            val(:,3) = -(2/3)*dt*( nu/dx^2 - 0.5*Uref/dx );
            val(1,:) = L_BC1(u,1); % apply LHS boundary conditions
            val(N,:) = L_BC2(u,N);
        end
    end
end