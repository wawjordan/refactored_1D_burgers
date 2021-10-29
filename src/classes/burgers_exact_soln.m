classdef burgers_exact_soln
   properties
      soln_type
      xmin, xmax, L
      t0
      Re, nu, alpha, Uref {mustBeNumeric}
   end
   methods
       function this = burgers_exact_soln(soln_type,Re,xlim,varargin)
           expectedSolutions = {'steady_shock','unsteady_shock',...
               'pulse_plus','pulse_minus'};
           defaultUref  = 2;
           defaultNu    = 0.25;
           defaultAlpha = 1;
           defaultT0    = 0;
           p = inputParser;
           validRange = @(x) isnumeric(x) && numel(x)==2 && x(1)<x(2);
           validScalarNum = @(x) isnumeric(x) && isscalar(x);
           validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
           addRequired(p,'soln_type',...
               @(x) any(validatestring(x,expectedSolutions)));
           addRequired(p,'Re',validScalarPosNum);
           addRequired(p,'xlim',validRange);
           addParameter(p,'Uref',defaultUref,validScalarPosNum);
           addParameter(p,'nu',defaultNu,validScalarPosNum);
           addParameter(p,'alpha',defaultAlpha,validScalarPosNum);
           addParameter(p,'t0',defaultT0,validScalarNum);
           parse(p,soln_type,Re,xlim,varargin{:});
           this.Re   = p.Results.Re;
           this.Uref = p.Results.Uref;
           this.nu   = p.Results.nu;
           this.xmin = p.Results.xlim(1);
           this.xmax = p.Results.xlim(2);
           this.L = this.xmax - this.xmin;
           this.t0 = p.Results.t0;
           this.alpha = p.Results.alpha;
           switch(p.Results.soln_type)
               case 'steady_shock'
                   this.soln_type = @steady_shock;
                   this.nu = this.Uref*this.L/this.Re;
                   this.alpha = this.Re/2;
               case 'unsteady_shock'
                   this.soln_type = @shock_coalesce;
                   this.nu = this.Uref*this.L/this.Re;
                   this.alpha = this.Re/2;
               case 'pulse_plus'
                   this.soln_type = @pulse_decay_plus;
                   this.nu = this.Uref*this.L/this.Re;
                   this.alpha = this.Re/2;
               case 'pulse_minus'
                   this.soln_type = @pulse_decay_minus;
           end
       end
       function Uex = eval(this,x,t)
          Uex = this.soln_type(this,x,t);
       end
       function uex = shock_coalesce(this,x,t)
           x2 = x*0.5*this.Re/this.L;
           t2 = t*0.25*this.Re*this.Uref/this.L;
           uex = -2*sinh(x2)./(cosh(x2)+exp(-t2));
       end
       function uex = pulse_decay_plus(this,x,t)
           a = 0.5*this.Re/this.L;
           x2 = x*a;
           t2 = t*a^2*this.nu;
           uex = x2./t2./(1 + sqrt(t2).*exp(x2.^2./(4*t2)));
       end
       function uex = pulse_decay_minus(this,x,t)
           a = 0.5*this.Re/this.L;
           x2 = x*a;
           t2 = t*a^2*this.nu;
           uex = x2./t2./(1 - sqrt(t2).*exp(x2.^2./(4*t2)));
       end
       function uex = steady_shock(this,x,~)
           a1 = -this.Re*this.nu/this.L;
           x2 = x*0.5*this.Re/this.L;
           uex = a1*tanh(x2);
       end
   end
end