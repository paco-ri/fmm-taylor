clear

% set up domain
n = 8;
nu = 8;
nv = 3 * nu;

dom = surfacemesh.stellarator(n, nu, nv);
C = dom.connectivity.elem2elem; % patch connectivity prior to refinement

amr_tol = 1e-6;
rmax = 5;
mode = 1;

% 75 patches initially
marked = 1:10; % just mark first ten for refinement for demonstration
% TURN OFF AMR FOR NOW
% dom = surfacemesh.adap_ref(dom, amr_tol, rmax, mode, marked);
% plot(dom)

domparams = [n, nu, nv];

% specify flux, Beltrami parameter, tolerance
flux = 1.0;
zk = 0.5;
tol = 1e-6;

ts = TaylorState(dom, domparams, zk, flux, tol);
ts = ts.solve(true);
% ts.get_quad_corr_taylor()
t1 = tic;
B = ts.surface_B();
t2 = toc(t1);
fprintf("compute surface B: \t%.4f seconds\n", t2);

% export to vtk for visualization
fname = "adaptive_stellarator_B.vtk";
surfacemesh_to_vtk(dom, fname, B{1})
fname = "adaptive_stellarator_sigma.vtk";
surfacemesh_to_vtk(dom, fname, ts.sigma{1})
fname = "adaptive_stellarator_mH.vtk";
surfacemesh_to_vtk(dom, fname, norm(ts.domain.mH{1}))
