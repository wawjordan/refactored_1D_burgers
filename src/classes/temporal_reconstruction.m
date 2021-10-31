classdef temporal_reconstruction
    properties
        N, M, M2, R, method
    end
    methods
        function this = temporal_reconstruction(grid,S,order,varargin)
            this.N = grid.N;
            this.M  = (S.stencil_size-1)/2;
            this.M2 = S.stencil_size;
            this.R = order; % polynomial order
            defaultMethod = 'default';
            expectedMethods = {'default','svd'};
            
            p = inputParser;
            validPosInt = @(x) isnumeric(x) && isscalar(x) ...
                && (x > 0) && mod(x,1) == 0;
            addRequired(p,'grid');
            addRequired(p,'S');
            addRequired(p,'order',validPosInt);
            addParameter(p,'method',defaultMethod,...
                @(x) any(validatestring(x,expectedMethods)));
            parse(p,grid,S,order,varargin{:});
            this.method = p.Results.method;
            switch(p.Results.method)
                case 'default'
                    this.method = @eval_default;
                case 'svd'
                    this.method = @eval_svd;
            end
        end
        function [u0,du1] = eval(this,stencil,time,iter)
            [u0,du1] = this.method(this,stencil,time,iter);
        end
        function [u0,du1] = eval_default(this,stencil,time,iter)
            U = stencil.U(:,:,iter);
            [A,mu] = vand_matrix(stencil.t,this.R);
            A = decomposition(A);
            u0 = zeros(this.N,1);
            du1 = zeros(this.N,1);
            for i = 1:this.N
                P = A\(U(i,:)');
                pd = ddpoly(P,time,2,mu);
                u0(i) = pd(:,1);
                du1(i) = pd(:,2);
            end
        end
        function [u0,du1] = eval_svd(this,stencil,time,iter)
            [U,S,V]=svd(stencil.U(:,:,iter),0);
            [A,mu] = vand_matrix(stencil.t,this.R);
            V0 = zeros(this.M2,1);
            V1 = zeros(this.M2,1);
            for i = 1:this.M2
                P = A\V(:,i);
                pd = ddpoly(P,time,2,mu);
                V0(i) = pd(:,1);
                V1(i) = pd(:,2);
            end
            u0  = U(:,1:this.M2)*S*V0;
            du1 = U(:,1:this.M2)*S*V1;
        end
    end
end