clear; clc;

% Load data for Assignment 4, Q4
A4Q4_wlsse_data;   % this defines nfrom, nto, r, x, b, Pinj, Qinj, Pflow, Qflow, Vnode, etc.

toler   = 1e-4;
maxiter = 20;

[delta, V, N, time] = fdwlsse(nfrom, nto, r, x, b, Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter);

fprintf('=== WLS State Estimation Results ===\n');
fprintf('Iterations: %d, CPU time: %.6f s\n\n', N, time);
fprintf('Bus   Delta(rad)   Delta(deg)   V(p.u.)\n');
for i = 1:length(V)
    fprintf('%3d   %9.5f   %9.5f   %8.4f\n', ...
        i, delta(i), delta(i)*180/pi, V(i));
end
