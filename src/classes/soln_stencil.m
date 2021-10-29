classdef soln_stencil
    properties
        M, N, L
        max_length
        queue_length
        U, t
    end
    methods
        function this = soln_stencil( M, N, varargin )
            this.M = M;
            this.N = N;
            default_L = 1;
            default_max_length = 10;
            p = inputParser;
            validPosInt = @(x) isnumeric(x) && isscalar(x) ...
                                       && (x > 0) && mod(x,1) == 0;
            addRequired(p,'M');
            addRequired(p,'N');
            addOptional(p,'L',default_L,validPosInt);
            addParameter(p,'max_length',default_max_length,validPosInt);
            parse(p,M,N,varargin{:});
            this.max_length = max(M,p.Results.max_length);
            this.L = p.Results.L;
            this.U = zeros(N,M,this.L);
            this.t = zeros(1,M);
            this.queue_length = 1;
        end
        function this = push(this,u_new,t_new)
            this.U(:,1:end-1,:) = this.U(:,2:end,:); % shift left
            this.t(1:end-1) = this.t(:,2:end);
            this.U(:,end,:) = repmat(u_new,1,1,this.L); % add new vector
            this.t(end) = t_new;
            this.queue_length = min(this.M, this.queue_length+1);
        end
        function this = grow(this,k)
            if (this.M < this.max_length)
                if (this.M+k < this.max_length)
                    this.U = cat(2,zeros(this.N,k,this.L),this.U);
                    this.t = [zeros(1,k),this.t];
                    this.M = this.M+k;
                else
                    warning('specified k exceeds max_length')
                    k2 = this.max_length-this.M;
                    this.U = cat(2,zeros(this.N,k2,this.L),this.U);
                    this.t = [zeros(1,k2),this.t];
                    this.M = this.max_length;
                end
            else
                warning('stencil has max_length')
            end
        end
    end
end
            