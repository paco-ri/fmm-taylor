% Fundamental-form iterative mesh refinement on a one-surface stellarator.
%
% Every geometry and physics parameter below is held identical to
% sigma_amr_stellarator.m.
%
% Correctness check: there is no obvious closed-form Taylor state for 
% this surface, so this script reports
%   - the fundamental-form indicator, which must fall as the mesh refines,
%   - fd_test, a finite-difference check of curl B = zk*B at an interior
%     point.
% For a true error against a known field see ff_amr_torusshell.m.
%
% Writes examples/ff_amr_stellarator.csv, per-level sigma VTK, and per-level
% ff_amr_stellarator_lev<k>*.mat for render_ff_amr.sh.

fprintf('=== ff_amr_stellarator start %s\n', datestr(now));

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
% Indicator scale on this geometry, measured on the uniform sequence:
%   level 0  max eta = 2.53e-1   (min 9.40e-2, median 1.84e-1)
%   level 1  max eta = 1.95e-2
%   level 2  max eta = 8.14e-4
% Dorfler marking ignores amr_tol; it is set here for the record and takes
% effect only under marking = 'threshold'.
rmax = 3;
opts = [];
opts.marking = 'dorfler';
opts.theta   = 0.5;
opts.amr_tol = 1e-3;
opts.mode    = 1;
opts.vtkbase = fullfile(fileparts(mfilename('fullpath')), ...
                        'ff_amr_stellarator');
% B, nodes and patch depths per level, for render_ff_amr.sh.
opts.savebase = opts.vtkbase;
opts.verbose = true;

out = ff_amr(dom0, domparams, zk, flux, tol, rmax, opts);

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
csv = fullfile(fileparts(mfilename('fullpath')), 'ff_amr_stellarator.csv');
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
