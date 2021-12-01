%% burgers_exac_soln class testing
clc; clear; close all;
% ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);
% x = linspace(-4,4,129)';
% hold on;
% for t = -2:0.4:2
% Uex = ex_soln.eval(x,t);
% plot(x,Uex)
% end
ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
f = @(x) 10*x./(1+(2/sqrt(5))*exp(20*x.^2));
% ex_soln1 = burgers_exact_soln('#1',64,[-4,4]);
% ex_soln2 = burgers_exact_soln('#1mod',64,[-4,4],'Uref',1);
% x = linspace(-4,4,1025)';
x = linspace(-2,2,1025)';
hold on
% plot(x,ex_soln.eval(x,0.1),'k')
plot(x,f(x)-ex_soln.eval(x,0.1),'r')
hold off;
hold on;
for t = 0:0.1:1
Uex1 = ex_soln1.eval(x,t);
Uex2 = ex_soln2.eval(x,t);
plot(x,Uex1,'k')
plot(x,Uex2)
end