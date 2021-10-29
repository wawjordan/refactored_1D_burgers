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
x = linspace(-2,2,129)';
hold on;
for t = 0.1:0.1:0.6
Uex = ex_soln.eval(x,t);
plot(x,Uex)
end