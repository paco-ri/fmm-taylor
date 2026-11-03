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
% Runs once per fundamental form (mode 1: first; mode 2: first and second).
% Writes, per mode m, examples/ff_amr_stellarator_ff<m>.csv, per-level sigma
% VTK, and per-level ff_amr_stellarator_ff<m>_lev<k>.mat for render_ff_amr.sh.

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
tol  = [1e-6 1e-8 1e-8];

% --- refinement controls ---
% Indicator scale on this geometry, measured on the uniform sequence
% (max eta, with min and median):
%             mode 1                       mode 2
%   level 0   2.53e-1 (9.40e-2, 1.84e-1)   4.23e-1 (2.13e-1, 2.90e-1)
%   level 1   1.95e-2 (2.14e-3, 6.73e-3)   1.64e-1 (6.82e-3, 3.85e-2)
%   level 2   8.14e-4 (5.65e-5, 1.97e-4)   5.52e-2 (8.32e-5, 1.45e-3)
rmax = 3;
for ffmode = 1:2
    opts = [];
    opts.amr_tol = 1e-2;
    opts.mode    = ffmode;
    opts.vtkbase = fullfile(fileparts(mfilename('fullpath')), ...
                            sprintf('ff_amr_stellarator_ff%d', ffmode));
    % B, nodes and patch depths per level, for render_ff_amr.sh.
    opts.savebase = opts.vtkbase;
    opts.verbose = true;

    fprintf('\n=== fundamental-form mode %d\n', ffmode);
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
    csv = [opts.vtkbase '.csv'];
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
end

fprintf('DONE_SENTINEL %s\n', datestr(now));
