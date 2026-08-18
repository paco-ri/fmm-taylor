%% GMRES system (A11 D = B_alpha) isolation on my_square_torus
%
% Isolates solveforD (TaylorState.solveforD, which calls
% TaylorState.mtxBalpha for the RHS and TaylorState.gmresA for the
% operator) from the flux/alpha bookkeeping in
% RefTaylorState.compute_sigma_alpha. Reports GMRES iteration counts and
% residuals for my_square_torus vs. a circular torus (prepare_torus) at
% comparable patch counts and tolerance.
%
% A large jump in iteration count, or a residual that stagnates above
% eps_gmres, on the square torus but not the circular torus would point
% at the sharp edges degrading the A11 operator's conditioning,
% independently of fluxalpha/fluxsigma.

clear
addpath('../examples')

tols = 1e-7; % [eps_gmres, eps_taylor, eps_laphelm]
zk = 0;

%% Square torus
n = 5;
m = 6;
nu = 3*m;
nv = 4*m;
dom = my_square_torus(n, nu, nv);
domparams = [n, nu, nv];

ts_sq = TaylorState(dom, domparams, zk, 0, tols);
ts_sq = ts_sq.get_quad_corr_laphelm();
ts_sq = ts_sq.get_quad_corr_taylor();

t1 = tic;
dfunc_sq = ts_sq.solveforD(true);
fprintf('square torus  (n=%d, nu=%d, nv=%d, npatches=%d): solveforD wall time = %.2f s\n', ...
    n, nu, nv, numel(dom.x), toc(t1));

%% Circular torus baseline at a comparable patch count
n2 = 5;
nv2 = round(sqrt(nu*nv/3));
nu2 = nv2*3;
domc = prepare_torus(n2,nu2,nv2,16,40);
domc = domc{1}; % prepare_torus always returns a cell of surfacemeshes

ts_c = TaylorState(domc, [n2, nu2, nv2], zk, 0, tols);
ts_c = ts_c.get_quad_corr_laphelm();
ts_c = ts_c.get_quad_corr_taylor();

t1 = tic;
dfunc_c = ts_c.solveforD(true);
fprintf('circular torus (n=%d, nu=%d, nv=%d, npatches=%d): solveforD wall time = %.2f s\n', ...
    n2, nu2, nv2, numel(domc.x), toc(t1));

fprintf(['\ntestgmres_square_torus: compare the printed GMRES iteration ' ...
    'counts/times (printed by solveforD itself). Re-run with the two ' ...
    'resolutions used by examples/conv_square_torus.m (m = 6 and m = 8) ' ...
    'to check whether the iteration count grows abnormally with ' ...
    'refinement for the square torus.\n']);

fprintf('||dfunc_sq{1,1}||_inf = %.6e\n', norm(dfunc_sq{1,1}, inf));
fprintf('||dfunc_c{1,1}||_inf  = %.6e\n', norm(dfunc_c{1,1}, inf));
