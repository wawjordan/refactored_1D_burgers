function R = ss_residual_cont(S,LS,U)
[u,du1,du2] = LS.eval(U);
R = -S.nu*du2 + u.*du1;
end