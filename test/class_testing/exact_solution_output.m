%% Exact Solution Output
clc; clear; close all;

dirname = 'C:\Users\Will\Desktop\';
filename1 = 'UEX_shock_coalesce.dat';
% filename1 = 'UEX_pulse_decay.dat';
% filename1 = 'UEX_moving_shock.dat';
% filename1 = 'UEX_expansion_fan.dat';

% x = linspace(-2,2,1025)';
% ex_soln = burgers_exact_soln('pulse_plus',64,[-2,2]);
% t = 0.1:0.1:0.6;

x = linspace(-4,4,1025)';
ex_soln = burgers_exact_soln('unsteady_shock',64,[-4,4]);
t = -2:0.01:1;
% ex_soln = burgers_exact_soln('move_shock',64,[-4,4]);
% t = -2:0.5:1;
% ex_soln = burgers_exact_soln('#1',64,[-4,4]);
% t = [0,0.1:0.01:1.1];
L = length(t);

name=[dirname,filename1];
fid=fopen(name,'wt');
fprintf(fid,'TITLE = solution\n');
str1 = 'variables="x","u"';
fprintf(fid,str1);
for k = 1:L
Uex = ex_soln.eval(x,t(k));
UX = [x,Uex];
fprintf(fid,'ZONE\n');
fprintf(fid,'T = "t = %0.2f"\n',round(100*t(k))/100);
fprintf(fid,'I=%d\n',length(Uex));
fprintf(fid,'ZONETYPE = Ordered\n');
fprintf(fid,'DT=(DOUBLE DOUBLE)\n');
fprintf(fid,'DATAPACKING=BLOCK\n');
[mm,nn]=size(UX);
for j = 1:nn
    for i = 1:mm
        fprintf(fid, '%g\n',UX(i,j));
    end
end
fprintf(fid, '\n');
end
fclose(fid);