%% script for concise IC tecplot output
clc; clear;
src_dirname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\',...
    '\post_processing\'];
src_fname = 'shock_exact_start_alg1_allsave';
% src_fname = 'ETE-IC-BD-4_pulseplus_test';
% src_dirname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\unsteady_iterative_correction_results',...
%     '\iterative_weighting_experiments\'];
% src_fname = 'ETE-IC-BD-4_alg1_lin-weights2';
% src_fname = 'ETE-IC-BD-4_alg1_inv-weights1';
% src_fname = 'ETE-IC-BD-4_alg1_s2-weights1';
% src_fname = 'ETE-IC-BD-4_alg1_no-weights';

target_dirname = 'C:\Users\Will\Documents\MATLAB\VT_Research\new\post_processing\';
% target_dirname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\new',...
%     '\p',...
%     '\tecplot_output\'];


% [SE] = parse_IC_OOA_local_error(src_fname);
[SE] = parse_IC_OOA_local_error_in_error(src_fname);
format_for_tecplot(target_dirname,'Local_Error_in_Error_Shock_4097.dat',SE);

% [SE1,SE2,PE1,PE2,DEBUG] = parse_IC_OOA_concise([src_dirname,src_fname]);
% [SE,PE,DEBUG] = parse_IC_OOA_concise_time([src_dirname,src_fname]);
% [SE,PE,DEBUG] = parse_IC_OOA_concise_time_fine([src_dirname,src_fname]);
% src_fname = 'pulse_ETEIC_newalg';
% format_for_tecplot(target_dirname,[src_fname,'L1_Error_fine.dat'],SE.L(1));

% format_for_tecplot(target_dirname,[src_fname,'L1_Error.dat'],SE1.L(1));
% format_for_tecplot(target_dirname,[src_fname,'L1_OOA.dat'],PE1.L(1));

% format_for_tecplot(target_dirname,['time_',src_fname,'L1_OOA.dat'],PE2.L(1));