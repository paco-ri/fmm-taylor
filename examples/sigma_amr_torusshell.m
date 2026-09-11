% Sigma-driven iterative mesh refinement on a two-surface toroidal shell.
%
% Writes examples/sigma_amr_torusshell.csv and per-level VTK.

fprintf('=== sigma_amr_torusshell start %s\n', datestr(now));

% --- geometry (outer surface is dom{1}) ---
n  = 5;
nv = 3;
nu = 3*nv;
ao = 1.0;        % outer minor radius
ai = 0.6;        % inner minor radius
nr = 16; nt = 40; np = 40;   % cross-section quadrature
[dom0, qnodes, qweights] = prepare_torus(n, nu, nv, n, nu, nv, ao, ai, ...
    nr, nt, np);
domparams = [n, nu, nv];
fprintf('base meshes: %d + %d patches, n = %d\n', ...
    length(dom0{1}.x), length(dom0{2}.x), n);

% --- reference Taylor state ---
zk     = 0.5;
ntheta = 200;
rmin   = 2.0;
rmaj   = 2.0;
jmag   = 1.0;
tol    = 1e-6;

% Fluxes of B0 through the two cross-sections. These depend only on B0 and
% the cross-section quadrature, not on the surface mesh, so they are fixed
% across refinement levels. Signs follow run_n9nv8_twosurface.m: toroidal
% flux on surface 1 takes +B0(2), poloidal on surface 2 takes -B0(3).
t0 = tic;
flux = zeros(1,2);
for i = 1:nr*nt
    B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{1}(:,i));
    flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
end
for i = 1:nr*np
    B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{2}(:,i));
    flux(2) = flux(2) - B0eval(3)*qweights{2}(i);
end
fprintf('flux = [%.16e, %.16e]  (%.1f s)\n', flux(1), flux(2), toc(t0));

refstate = [];
refstate.xs_nodes   = qnodes;
refstate.xs_weights = qweights;
% B0 has to be re-evaluated on whatever mesh each level produces.
refstate.B0fun = @(d, s) reftaylorsurffun(d, n, ntheta, rmin, ...
    rmaj, jmag, zk);

% --- refinement controls ---
rmax = 2;
opts = [];
opts.marking   = 'dorfler';
opts.theta     = 0.6;
opts.eps_sigma = 1e-4;
opts.refstate  = refstate;
opts.vtkbase   = fullfile(fileparts(mfilename('fullpath')), ...
                          'sigma_amr_torusshell');
opts.verbose   = true;

out = sigma_amr(dom0, domparams, zk, flux, tol, rmax, opts);

% --- report ---
csv = fullfile(fileparts(mfilename('fullpath')), 'sigma_amr_torusshell.csv');
fid = fopen(csv, 'w');
fprintf(fid, ['level,npat_outer,npat_inner,npts,nmarked,max_eta,l2_eta,' ...
    'err_vs_B0,time_s\n']);
for k = 1:numel(out.history)
    h = out.history(k);
    fprintf(fid, '%d,%d,%d,%d,%d,%.6e,%.6e,%.6e,%.3f\n', ...
        h.level, h.npatches(1), h.npatches(2), sum(h.npts), h.nmarked, ...
        h.max_eta, h.l2_eta, h.err_vs_B0, h.time_s);
end
fclose(fid);
fprintf('SAVED %s\n', csv);

fprintf('\n level    patches      npts   marked      max eta    err vs B0\n');
for k = 1:numel(out.history)
    h = out.history(k);
    fprintf('%6d %5d+%-5d %9d %8d %12.3e %12.3e\n', h.level, ...
        h.npatches(1), h.npatches(2), sum(h.npts), h.nmarked, ...
        h.max_eta, h.err_vs_B0);
end

fprintf('\nstopped on: %s\n', out.reason);

fprintf('DONE_SENTINEL %s\n', datestr(now));
