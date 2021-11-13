classdef spatial_reconstruction_v2
    properties
        N, M, M2, R
        X, x, mu
    end
    methods
        function this = spatial_reconstruction_v2(grid,S,order)
            this.N = grid.N;
            this.M  = (S.stencil_size-1)/2;
            this.M2 = S.stencil_size;
            this.R = order; % polynomial order
            this.x = grid.x;
            
            this.X = cell(3,1);
            this.X{1} = my_savgol_end(this.R,this.M2, 1);
            this.X{2} = my_savgol_end(this.R,this.M2, 0);
            [~,this.X{2}] = sgolay(this.R,this.M2);
            this.X{3} = my_savgol_end(this.R,this.M2,-1);
            this.mu = factorial(0:this.R-1)./(-grid.dx).^(0:this.R-1);
        end
        function [u,du1,du2] = eval(this,U)
            u = this.mu(:,1).*conv(U, this.X{2}(:,1), 'same');
            du1 = this.mu(:,2).*conv(U, this.X{2}(:,2), 'same');
            du2 = this.mu(:,3).*conv(U, this.X{2}(:,3), 'same');
            
            for i = 1:this.M
                i1 = 1:this.M2;
                u(i) = U(i1)'*this.mu(i,1)*this.X{1}(:,1,i);
                du1(i) = U(i1)'*this.mu(i,2)*this.X{1}(:,2,i);
                du2(i) = U(i1)'*this.mu(i,3)*this.X{1}(:,3,i);
            end
%             for i = this.N:-1:this.N-(this.M-1)
            for i = this.N-(this.M-1):this.N
%                 ioff = this.N+this.M-1-i;
                ioff = i-(this.N-this.M);
                i1 = this.N-(this.M2-1):this.N;
                u(i) = U(i1)'*this.mu(i,1)*this.X{3}(:,1,ioff);
                du1(i) = U(i1)'*this.mu(i,2)*this.X{3}(:,2,ioff);
                du2(i) = U(i1)'*this.mu(i,3)*this.X{3}(:,3,ioff);
            end
        end
    end
end