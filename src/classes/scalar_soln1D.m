classdef scalar_soln1D
   properties
      U, R, E, Rnorm, Rinit, N {mustBeNumeric}
   end
   methods
       function this = scalar_soln1D(grid)
           imax = grid.N;
           this.U = zeros(imax,1);
           this.R = zeros(imax,1);
           this.E = zeros(imax,1);
           this.Rnorm = 0;
           this.Rinit = 0;
           this.N = imax;
       end
   end
end