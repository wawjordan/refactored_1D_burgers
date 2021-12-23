classdef burgers_exact_soln
   properties
      soln_type
      xmin, xmax, L
      t0
      Re, nu, alpha, Uref {mustBeNumeric}
   end
   methods
       function this = burgers_exact_soln(soln_type,Re,xlim,varargin)
           expectedSolutions = {'#1','#1mod','steady_shock','unsteady_shock',...
               'move_shock','pulse_plus','pulse_minus'};
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
               case '#1'
                   this.soln_type = @soln1;
                   this.nu = this.Uref*this.L/this.Re;
                   this.alpha = this.Re/2;
               case '#1mod'
                   this.soln_type = @soln1mod;
                   this.nu = this.Uref*this.L/this.Re;
                   this.alpha = this.Re/2;
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
               case 'move_shock'
                   this.soln_type = @moving_shock;
                   this.nu = this.Uref*this.L/this.Re;
                   this.alpha = this.Re/2;
           end
       end
       function Uex = eval(this,x,t)
          Uex = this.soln_type(this,x,t);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function uex = soln1(this,x,t)
           if t <= 0
               uex(x<=0) = 0;
               uex(x>0) = 1;
%                uex(x<0) = 0;
%                uex(x>=0) = 1;
               uex = uex(:);
           else
               E  = exp((t-2*x)./(4*this.nu));
               tx = (t-x)./sqrt(4*this.nu*t);
               Etx = erf(tx);
               Ex = erf(x./sqrt(4*this.nu*t));
               num = E.*(1-Etx);
               den = E.*(1-Etx) + 1 - Ex;
               uex = num./den;
           end
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function uex = soln1mod(this,x,t)
%            this.Uref = 1;
             G = @(x,t) exp(t/this.Uref^2-x/this.Uref).*erfc((t/this.Uref-0.5*x)./sqrt(t));
%            G = @(x,t) 0.5*exp(t-x).*erfc((2*t-x)./(2*sqrt(t)));
           if t <= 0
               uex(x<0) = 0;
               uex(x>=0) = this.Uref;
               uex = uex(:);
           else
               a = this.alpha/this.L;
%                a = 0.25*this.Re/this.L;
               x2 = x*a;
               t2 = t*a^2*this.nu;
%                x2 = x;
%                t2 = t;
               tmp1 = G(x2,t2);
               tmp2 = G(-x2,t2);
%                tmp3 = 0.25*erfc(x2./(2*sqrt(t2)));
               tmp3 = erfc(x2./(2*sqrt(t2)));
               uex = this.Uref*(tmp1)./(tmp1+tmp3);
           end
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
       function uex = moving_shock(this,x,t)
           uR = this.Uref;
           uL = 0;
           Vs = 0.5*(uR+uL);
           uex = uR + 0.5*(uL - uR)*...
               (1-tanh(((uL-uR)*(x-Vs*t))/(4*this.nu)));
       end
           
           
   end
end