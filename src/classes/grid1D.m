classdef grid1D
    properties
        x, dx, xmin, xmax, N, L {mustBeNumeric}
    end
    methods
        function this = grid1D( x )
            this.xmin = min(x);
            this.xmax = max(x);
            this.N = length(x);
            this.L = max(x)-min(x);
            this.x = x(:);
            xplus = [(2*x(1)-x(2));x(:);(2*x(end)-x(end-1))];
            this.dx = 0.5*(diff(xplus(1:end-1))+diff(xplus(2:end)));
        end
    end
end