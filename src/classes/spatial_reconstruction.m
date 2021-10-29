classdef spatial_reconstruction
    properties
        N, M, M2, R
        X, x, mu
    end
    methods
        function this = spatial_reconstruction(grid,S,order)
            this.N = grid.N;
            this.M  = (S.stencil_size-1)/2;
            this.M2 = S.stencil_size;
            this.R = order; % polynomial order
            this.x = grid.x;
            
            this.X = cell(this.N,1);
            this.mu = zeros(this.N,2);
            for i = 1:this.M+1
                i1 = 1:this.M2;
                [A,mu] = vand_matrix(grid.x(i1),this.R);
                this.X{i} = decomposition(A);
                this.mu(i,:) = mu;
            end
            for i = this.M+2:this.N-(this.M+1)
                i1 = i-this.M:i+this.M;
                [A,mu] = vand_matrix(grid.x(i1),this.R);
                this.X{i} = decomposition(A);
                this.mu(i,:) = mu;
            end
            for i = this.N-this.M:this.N
                i1 = this.N-this.M2+1:this.N;
                [A,mu] = vand_matrix(grid.x(i1),this.R);
                this.X{i} = decomposition(A);
                this.mu(i,:) = mu;
            end
        end
        function [u,du1,du2] = eval(this,U)
            u   = zeros(this.N,1);
            du1 = zeros(this.N,1);
            du2 = zeros(this.N,1);
            
            for i = 1:this.M+1
                i1 = 1:this.M2;
                P = this.X{i}\U(i1,:);
                pd = ddpoly(P,this.x(i),3,this.mu(i,:));
                  u(i) = pd(:,1);
                du1(i) = pd(:,2);
                du2(i) = pd(:,3);
            end
            for i = this.M+2:this.N-(this.M+1)
                i1 = i-this.M:i+this.M;
                P = this.X{i}\U(i1,:);
                pd = ddpoly(P,this.x(i),3,this.mu(i,:));
                  u(i) = pd(:,1);
                du1(i) = pd(:,2);
                du2(i) = pd(:,3);
            end
            for i = this.N-this.M:this.N
                i1 = this.N-this.M2+1:this.N;
                P = this.X{i}\U(i1,:);
                pd = ddpoly(P,this.x(i),3,this.mu(i,:));
                  u(i) = pd(:,1);
                du1(i) = pd(:,2);
                du2(i) = pd(:,3);
            end
        end
    end
end