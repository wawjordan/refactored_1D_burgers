classdef temporal_reconstruction
    properties
        N, M, M2, R
    end
    methods
        function this = temporal_reconstruction(grid,S,order)
            this.N = grid.N;
            this.M  = (S.stencil_size-1)/2;
            this.M2 = S.stencil_size;
            this.R = order; % polynomial order
        end
        function [u0,du1] = eval(this,stencil,time,iter)
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