classdef time_integrator_type
    properties
        t, dt
        max_newton_iter
        newton_tol
    end
    methods (Abstract)
        step
    end
    methods
        function this = time_integrator_type(dt,varargin)
            default_t   = 0;
            default_nmi = 10;
            default_nt  = 1e-10;
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x)...
                                  && isscalar(x) ...
                                  && (x > 0);
            validScalarPosInt = @(x) isnumeric(x)...
                                  && isscalar(x) ...
                                  && (x > 0)     ...
                                  && mod(x,1) == 0;
            addRequired(p,'dt');
            addParameter(p,'time',...
                default_t,validScalarPosNum);
            addParameter(p,'max_newton_iter',...
                default_nmi,validScalarPosInt);
            addParameter(p,'newton_tol',...
                default_nt,validScalarPosNum);
            parse(p,dt,varargin{:});
            this.dt = dt;
            this.t = p.Results.time;
            this.max_newton_iter = p.Results.max_newton_iter;
            this.newton_tol = p.Results.newton_tol;
        end
    end
end