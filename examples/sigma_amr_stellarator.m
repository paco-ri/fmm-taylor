% Sigma-driven iterative mesh refinement on a one-surface stellarator.
%
% Solves the Taylor state, refines, solves again, and uses the change in the
% density sigma between successive discretizations as the per-patch error
% indicator. Patches where the coarse discretization already resolved sigma
% are left alone. See sigma_amr.m.
%
% Correctness check: there is no closed-form reference Taylor state for this
% surface (the reftaylor current ring links the torus, and a ring radius
% valid for the circular torus is not guaranteed to clear the stellarator),
% so this script reports
%   - the sigma indicator, which must fall as the mesh refines, and
%   - fd_test, a finite-difference check of curl B = zk*B at an interior
%     point, which is an independent physical check of the solution.
% For a true error against a known field see sigma_amr_torusshell.m, which
% uses RefTaylorState on nested circular tori.
%
% Writes examples/sigma_amr_stellarator.csv and per-level VTK.

fprintf('=== sigma_amr_stellarator start %s\n', datestr(now));

% --- geometry ---
n  = 5;          % polynomial order + 1
nv = 3;          % patches in poloidal direction
nu = 3*nv;       % patches in toroidal direction
dom0 = prepare_stellarator(n, nu, nv, 16, 40);
dom0 = dom0{1};
domparams = [n, nu, nv];
fprintf('base mesh: %d patches, n = %d\n', length(dom0.x), n);

% --- physics ---
zk   = 0.5;      % Beltrami parameter
flux = 1.0;
tol  = 1e-6;

% --- refinement controls ---
rmax = 3;   % first completed run; raise once it lands
opts = [];
opts.marking   = 'dorfler';
opts.theta     = 0.5;
opts.eps_sigma = 1e-4;
opts.vtkbase   = fullfile(fileparts(mfilename('fullpath')), ...
                          'sigma_amr_stellarator');
opts.verbose   = true;

out = sigma_amr(dom0, domparams, zk, flux, tol, rmax, opts);

% --- independent physical check: curl B = zk*B at an interior point ---
% Interior point recipe follows examples/ux2.m for this geometry.
rmaj = 5.0;
phi = 4*pi/7;
h = 1e-6;
center = [rmaj*cos(phi)+.2 rmaj*sin(phi)-.1 .5];
[errB, curlB, kB] = out.ts.fd_test(center, h);
fprintf('\nfd_test at interior point: norm(curl B - zk*B) = %.6e\n', norm(errB));
fprintf('   (relative to norm(zk*B) = %.6e -> %.3e)\n', norm(kB), ...
    norm(errB)/norm(kB));

% --- report ---
csv = fullfile(fileparts(mfilename('fullpath')), 'sigma_amr_stellarator.csv');
fid = fopen(csv, 'w');
fprintf(fid, 'level,npatches,npts,nmarked,max_eta,l2_eta,re_alpha,im_alpha,time_s\n');
for k = 1:numel(out.history)
    h = out.history(k);
    fprintf(fid, '%d,%d,%d,%d,%.6e,%.6e,%.16e,%.16e,%.3f\n', ...
        h.level, sum(h.npatches), sum(h.npts), h.nmarked, ...
        h.max_eta, h.l2_eta, real(h.alpha(1)), imag(h.alpha(1)), h.time_s);
end
fclose(fid);
fprintf('SAVED %s\n', csv);

fprintf('\n level  patches      npts   marked      max eta      alpha (real)\n');
for k = 1:numel(out.history)
    h = out.history(k);
    fprintf('%6d %8d %9d %8d %12.3e %17.10e\n', h.level, sum(h.npatches), ...
        sum(h.npts), h.nmarked, h.max_eta, real(h.alpha(1)));
end

% Compare against the cost of refining everything to the same depth.
lev = out.history(end).level;
fprintf('\nfinal mesh %d patches; uniform refinement to level %d would be %d\n', ...
    sum(out.history(end).npatches), lev, length(dom0.x)*4^lev);

fprintf('stopped on: %s\n', out.reason);

fprintf('DONE_SENTINEL %s\n', datestr(now));
