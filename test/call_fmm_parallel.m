clear 
close

% set up domain
n = 8;
nu = 8;
nv = 3 * nu;

dom = surfacemesh.stellarator(n, nu, nv);
C = dom.connectivity.elem2elem; % patch connectivity prior to refinement

amr_tol = 1e-6;
rmax = 5;
mode = 1;

domparams = [n, nu, nv];

% specify flux, Beltrami parameter, tolerance
flux = 1.0;
zk = 0.5;
tol = 1e-6;

ts = TaylorState(dom, domparams, zk, flux, tol);
ts = ts.get_quad_corr_laphelm();
% ts = ts.get_quad_corr_taylor();
