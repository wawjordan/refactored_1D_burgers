%% script for IC tecplot output (new formatting, 02/14)
clc; clear;
src_dirname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\results',...
    '\02_11\'];
src_fname = 'expansion_fan_GS10';
% src_fname = 'expansion_fan_GS102';

target_dirname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\results',...
    '\02_15\'];
% target_dirname = src_dirname;
[SE] = parse_IC_OOA_new(src_fname);
format_for_tecplot(target_dirname,[src_fname,'L1.dat'],SE.L(1));
% format_for_tecplot(target_dirname,[src_fname,'L2.dat'],SE.L(2));
% format_for_tecplot(target_dirname,[src_fname,'L3.dat'],SE.L(3));