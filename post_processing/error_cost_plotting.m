%% Plotting cost of ETE & IC vs error
clc; clear; close all;
% fnames = cell(4,1);
dirname1 = ['C:\Users\Will\Documents\MATLAB\',...
               'VT_Research\new\Himalia\shock_coalesce_profiles\'];
           
ET = struct();
fname1 = [dirname1,'primal_shock_coalesce_fan.mat'];
load(fname1)
ET(1).T = primal_time.T(:);
ET(1).E = primal_time.Efinal(:,1);

ET(2).T = primal_time.T(:);
ET(2).E = primal_time.Erich(:,1);

fname = [dirname1,'error_ETE_no_IC_newton10_alg1_exact_startup_shock_coalesce.mat'];
load(fname)
ET(3).T = [0;0;error_time.T(:,1)];
ET(3).E = [0;0;error_time.Efinal(:,1,1)];

fnames = [dirname1,'error_ETE_no_IC_newton10_alg1_shock_coalesce.mat'];
load(fname)
ET(4).T = [0;0;error_time.T(:,1)];
ET(4).E = [0;0;error_time.Efinal(:,1,1)];

fname = [dirname1,'error_ETEIC_newton10_alg4_shock_coalesce.mat'];
load(fname)

ET(5).T = [0;0;error_time.T(:,1)];
ET(5).E = [0;0;error_time.Efinal(:,3,1)];

fname = [dirname1,'error_ETEIC_newton10_alg1_exact_startup_shock_coalesce.mat'];
load(fname)

ET(6).T = [0;0;error_time.T(:,1)];
ET(6).E = [0;0;error_time.Efinal(:,end,1)];

fname = [dirname1,'error_ETEIC_newton10_alg1_shock_coalesce_fan.mat'];
load(fname)

ET(7).T = [0;0;error_time.T(:,1)];
ET(7).E = [0;0;error_time.Efinal(:,end,1)];



% fname1 = [dirname1,'primal_profile_expansion_fan.mat'];
% load(fname1)
% ET(1).T = primal_time.T(:);
% ET(1).E = primal_time.Efinal(:,1);
% 
% ET(2).T = primal_time.T(:);
% ET(2).E = primal_time.Erich(:,1);
% 
% % fname = [dirname1,'error_ETE_no_IC_newton1_alg1_exact_startup_expansion_fan.mat'];
% % load(fname)
% % ET(2).T = [error_time.T(:,1,1);0;0];
% % ET(2).E = [error_time.Efinal(:,1,1);0;0];
% 
% fname = [dirname1,'error_ETE_no_IC_newton10_alg1_exact_startup_expansion_fan.mat'];
% load(fname)
% ET(3).T = [0;0;error_time.T(:,1)];
% ET(3).E = [0;0;error_time.Efinal(:,1,1)];
% 
% fnames = [dirname1,'error_ETE_no_IC_newton10_alg1_expansion_fan.mat'];
% load(fname)
% ET(4).T = [0;0;error_time.T(:,1)];
% ET(4).E = [0;0;error_time.Efinal(:,1,1)];
% 
% fname = [dirname1,'error_ETEIC_newton10_alg4_svd_expansion_fan.mat'];
% load(fname)
% 
% ET(5).T = [0;0;error_time.T(:,1)];
% ET(5).E = [0;0;error_time.Efinal(:,3,1)];
% 
% fname = [dirname1,'error_ETEIC_newton10_alg1_exact_startup_expansion_fan.mat'];
% load(fname)
% 
% ET(6).T = [0;0;error_time.T(:,1)];
% ET(6).E = [0;0;error_time.Efinal(:,end,1)];
% 
% fname = [dirname1,'error_ETEIC_newton10_alg1_svd_expansion_fan.mat'];
% load(fname)
% 
% ET(7).T = [0;0;error_time.T(:,1)];
% ET(7).E = [0;0;error_time.Efinal(:,end,1)];



% dirname1 = ['C:\Users\Will\Documents\MATLAB\',...
%                'VT_Research\new\Himalia\shock_coalesce_profiles\'];
% fname1 = [dirname1,'primal_profile_shock_coalesce.mat'];
% fnames{1} = [dirname1,'error_ETE_no_IC_newton1_alg1_exact_startup_shock_coalesce.mat'];
% fnames{2} = [dirname1,'error_ETE_no_IC_newton10_alg1_exact_startup_shock_coalesce.mat'];
% fnames{3} = [dirname1,'error_ETE_no_IC_newton1_alg1_shock_coalesce.mat'];
% fnames{4} = [dirname1,'error_ETE_no_IC_newton10_alg1_shock_coalesce.mat'];
% fnames{5} = [dirname1,'error_ETEIC_newton1_alg4_shock_coalesce.mat'];
% fnames{6} = [dirname1,'error_ETEIC_newton10_alg4_shock_coalesce.mat'];
% fnames{7} = [dirname1,'error_ETEIC_newton10_alg1_exact_startup_shock_coalesce.mat'];
% fnames{8} = [dirname1,'error_ETEIC_newton10_alg1_shock_coalesce.mat'];



% load(fname1)
% PT = primal_time;
% ET = struct();
% for i = 1:length(fnames)
%     load(fnames{i})
%     ET(i).E = error_time.Efinal(:,1,1);
%     ET(i).T = error_time.T(:,1,1);
% end

% 
% hold on;
% plot(PT.T,PT.Efinal(:,1),'k');
% plot(PT.T,PT.Erich(:,1),'r');
% plot(ET(1).E.T,ET(1).E.Efinal(:,:,1),'b');
% plot(ET(2).E.T,ET(2).E.Efinal(:,:,1),'b--');
% plot(ET(3).E.T,ET(3).E.Efinal(:,:,1),'m');
% plot(ET(4).E.T,ET(4).E.Efinal(:,:,1),'m--');
% plot(ET(5).E.T,ET(5).E.Efinal(:,:,1),'c');
% plot(ET(6).E.T,ET(6).E.Efinal(:,:,1),'c--');
% plot(ET(7).E.T,ET(7).E.Efinal(:,:,1));
% hold off;
% set(gca,'yscale','log')
% set(gca,'xscale','log')




% get dimensions for output arrays
L = 9;  % (h) grid levels (space+time comb.)
N = 10; % (j) iteration levels
O = 3;  % (k) norms

nNodes  = 2.^(7:15)+1;

SE = struct();
SE.L = struct();
SE.L.title='blah';
SE.L.variables = cell(14);
SE.L.zoneFmt='%0.2f';
SE.L.zoneVar=1;
SE.L.Nzones=1;
SE.L.dataFmt='%g';
SE.L.variables{1} = 'wall-clock time 1';
SE.L.variables{2} = 'wall-clock time 2';
SE.L.variables{3} = 'wall-clock time 3';
SE.L.variables{4} = 'wall-clock time 4';
SE.L.variables{5} = 'wall-clock time 5';
SE.L.variables{6} = 'wall-clock time 6';
SE.L.variables{7} = 'wall-clock time 7';
SE.L.variables{8} = 'Primal';
SE.L.variables{9} = 'Richardson Extrapolation';
SE.L.variables{10} = 'ETE - exact startup';
SE.L.variables{11} = 'ETE - general startup';
SE.L.variables{12} = 'ETE - iterative correction at startup';
SE.L.variables{13} = 'ETE + 10 steps IC - exact startup';
SE.L.variables{14} = 'ETE + 10 steps IC - general startup';
SE.L.DATA.dat = [[ET(:).T],[ET(:).E]];
format_for_tecplot(dirname1,'L1_Error_cost_estimates.dat',SE.L)

