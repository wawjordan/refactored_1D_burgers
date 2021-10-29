%% soln_stencil testing
clc; clear; close all;
M = 5;
N = 101;
L = 2;
sten = soln_stencil( M, N, L, 'max_length', 10);

x = linspace(0,pi,N)';
sten = sten.grow(2);
for t = 0:7
    u_new = sin(x+t/pi);
    sten = sten.push(u_new,t);
end

plot(x,sten.U(:,:,1))