% Sigma-driven iterative mesh refinement on the one-surface torus of
% convoop_sweep.m, against a reference Taylor state whose current ring passes
% close to the surface.
%
% Same refinement setup as sigma_amr_stellarator.m. The differences:
%   - geometry is prepare_torus (one surface), as in convoop_sweep.m
%   - the reftaylor ring (in the plane x = 0, centred at (0, rmaj, 0)) links
%     the torus and has its radius chosen so that its minimum distance to the
%     surface is exactly `gap`. B0 then peaks sharply where the ring passes
%     closest, which gives the indicator something localized to find.
%   - RefTaylorState is used, so the true error against B0 is reported per
%     level (there is no fd_test).
%
% Variables read if set before running:
%   setup_only [false]      stop after the geometry check plot
%   nv         [3]          base-mesh patches poloidally (nu = 3*nv)
%   rmax       [3]          maximum refinement depth, and the number of
%                           refinements
%   frac       [0.3]        opts.frac of sigma_amr
%
% With BASE = sigma_amr_torusring_fraction, suffixed _nv<nv> for a
% non-default nv, _rmax<rmax> for a non-default rmax and _frac<100*frac> for
% a non-default frac, writes
% BASE_setup.{png,fig} (the geometry check plot), BASE.csv, per-level sigma
% VTK, per-level BASE_lev<k>.mat for make_amr_vtk.py / view_amr_levels.py,
% and BASE_levels_matlab.png showing |B| and the patch mesh at every level.

if ~exist('setup_only', 'var'); setup_only = false; end
if ~exist('nv', 'var'); nv = 3; end
if ~exist('rmax', 'var'); rmax = 3; end
if ~exist('frac', 'var'); frac = 0.3; end
here = fileparts(mfilename('fullpath'));
base = fullfile(here, 'sigma_amr_torusring_fraction');
if nv ~= 3; base = sprintf('%s_nv%d', base, nv); end
if rmax ~= 3; base = sprintf('%s_rmax%d', base, rmax); end
if frac ~= 0.3; base = sprintf('%s_frac%g', base, 100*frac); end
fprintf('=== sigma_amr_torusring start %s\n', datestr(now));

% --- geometry (convoop_sweep.m) ---
n  = 5;
nu = 3*nv;
nr = 12; nt = 50;   % cross-section quadrature for the flux
[dom0, qnodes, qweights] = prepare_torus(n, nu, nv, nr, nt);
dom0 = dom0{1}; qnodes = qnodes{1}; qweights = qweights{1};
domparams = [n, nu, nv];
fprintf('base mesh: %d patches, n = %d\n', length(dom0.x), n);

% --- current ring ---
% Distance from p to the ring is sqrt(px^2 + (rho - rmin)^2), with rho the
% in-plane distance of p from the ring centre. The ring links the torus
% through the u = pi/2 cross-section, which lies exactly in x = 0, so it
% clears the surface only if rmin exceeds rho all along that cross-section.
% Surface sampled on a fine mesh; its node spacing (~0.03) bounds how far the
% sampled minimum can overshoot the true one, to O(spacing^2).
gap  = 0.2;         % minimum ring-surface distance (minor radius ~1)
rmaj = 2.0;         % ring centre, on the torus's magnetic axis
jmag = 1.0;
zk   = 0.5;
tol  = 1e-6;

nuf = 36;           % multiple of 4, so u = pi/2 is a patch edge
domf = prepare_torus(16, nuf, 12, 2, 2);
domf = domf{1};
P = zeros(3, 0); slice = zeros(2, 0);
for i = 1:length(domf.x)
    P = [P, [domf.x{i}(:).'; domf.y{i}(:).'; domf.z{i}(:).']]; %#ok<AGROW>
    if max(abs(domf.x{i}(:,1))) < 1e-12 && domf.y{i}(1,1) > 0
        slice = [slice, [domf.y{i}(:,1).'; domf.z{i}(:,1).']]; %#ok<AGROW>
    end
end
rho = sqrt((P(2,:) - rmaj).^2 + P(3,:).^2);
rhoslice = sqrt((slice(1,:) - rmaj).^2 + slice(2,:).^2);
dmin = @(r) min(sqrt(P(1,:).^2 + (rho - r).^2));
rmin = fzero(@(r) dmin(r) - gap, [max(rhoslice), max(rhoslice) + 1]);
[dclose, iclose] = min(sqrt(P(1,:).^2 + (rho - rmin).^2));
pclose = P(:, iclose);
if rmaj + rmin >= 2*rmaj
    error('ring reaches the u = 3*pi/2 cross-section');
end

% The ring is summed with the trapezoidal rule, whose error at distance gap
% decays like exp(-2*pi*gap/h), h the ring point spacing. h <= gap/4 puts it
% below 1e-10.
ntheta = max(200, ceil(8*pi*rmin/gap));
fprintf(['ring: rmaj = %.4f, rmin = %.6f, ntheta = %d\n' ...
    '      min distance to surface = %.6f at (%.4f, %.4f, %.4f)\n' ...
    '      (u = pi/2 cross-section extends to rho = %.6f)\n'], ...
    rmaj, rmin, ntheta, dclose, pclose, max(rhoslice));

% --- geometry check plot ---
th = linspace(0, 2*pi, 400);
ring = [zeros(size(th)); rmaj + rmin*cos(th); rmin*sin(th)];
cring = [0.835 0.369 0.000];   % Okabe-Ito vermillion
cpt   = [0.000 0.447 0.698];   % Okabe-Ito blue
fig = figure('Visible', 'off', 'Position', [0 0 1500 650]);
theme(fig, 'light');

subplot(1, 2, 1); hold on
for i = 1:length(dom0.x)
    surf(dom0.x{i}, dom0.y{i}, dom0.z{i}, 'FaceColor', [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceLighting', 'gouraud');
    for e = {{1,':'}, {n,':'}, {':',1}, {':',n}}
        xe = dom0.x{i}(e{1}{:}); ye = dom0.y{i}(e{1}{:}); ze = dom0.z{i}(e{1}{:});
        plot3(xe, ye, ze, 'k', 'LineWidth', 0.6);
    end
end
plot3(ring(1,:), ring(2,:), ring(3,:), 'Color', cring, 'LineWidth', 2.5);
plot3(pclose(1), pclose(2), pclose(3), 'o', 'MarkerFaceColor', cpt, ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 8);
axis equal; view(60, 25); camlight headlight; grid on
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('base mesh (%d patches) and current ring', length(dom0.x)));

subplot(1, 2, 2); hold on
[~, ord] = sort(atan2(slice(2,:), slice(1,:) - rmaj));
plot(slice(1,[ord ord(1)]), slice(2,[ord ord(1)]), 'k', 'LineWidth', 1.5);
plot(ring(2,:), ring(3,:), 'Color', cring, 'LineWidth', 2);
plot(rmaj, 0, '+', 'Color', cring, 'MarkerSize', 10);
plot(pclose(2), pclose(3), 'o', 'MarkerFaceColor', cpt, ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 8);
axis equal; grid on
xlabel('y'); ylabel('z');
legend({'surface, u = \pi/2 cross-section (x = 0)', 'current ring', ...
    'ring centre', 'closest surface point'}, 'Location', 'southoutside');
title(sprintf('plane x = 0: r_{min} = %.4f, gap = %.4f', rmin, dclose));

exportgraphics(fig, [base '_setup.png'], 'Resolution', 150);
savefig(fig, [base '_setup.fig']);
fprintf('SAVED %s_setup.png / .fig\n', base);

if setup_only
    fprintf('DONE_SENTINEL %s\n', datestr(now));
    return
end

% --- reference Taylor state ---
t0 = tic;
flux = 0;
for i = 1:nr*nt
    B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes(:,i));
    flux = flux + B0eval(2)*qweights(i);
end
fprintf('flux = %.16e  (%.1f s)\n', flux, toc(t0));

refstate = [];
refstate.xs_nodes   = {qnodes};
refstate.xs_weights = {qweights};
refstate.B0fun = @(d, s) reftaylorsurffun(d, n, ntheta, rmin, rmaj, jmag, zk);

% --- refinement controls ---
opts = [];
opts.frac      = frac;
opts.refstate  = refstate;
opts.vtkbase   = base;
opts.savebase  = base;
opts.verbose   = true;

out = sigma_amr(dom0, domparams, zk, flux, tol, rmax, opts);
hist = out.history;
save([base '_hist.mat'], 'hist');

% --- report ---
csv = [base '.csv'];
fid = fopen(csv, 'w');
fprintf(fid, ['level,npatches,npts,nmarked,max_eta,l2_eta,err_vs_B0,' ...
    'time_s,l2err_vs_B0\n']);
for k = 1:numel(out.history)
    h = out.history(k);
    fprintf(fid, '%d,%d,%d,%d,%.6e,%.6e,%.6e,%.3f,%.6e\n', ...
        h.level, sum(h.npatches), sum(h.npts), h.nmarked, ...
        h.max_eta, h.l2_eta, h.err_vs_B0, h.time_s, h.l2err_vs_B0);
end
fclose(fid);
fprintf('SAVED %s\n', csv);

fprintf(['\n level  patches      npts   marked      max eta  inf err vs B0' ...
    '   L2 err vs B0\n']);
for k = 1:numel(out.history)
    h = out.history(k);
    fprintf('%6d %8d %9d %8d %12.3e %14.3e %14.3e\n', h.level, ...
        sum(h.npatches), sum(h.npts), h.nmarked, h.max_eta, h.err_vs_B0, ...
        h.l2err_vs_B0);
end

lev = out.history(end).level;
fprintf('\nfinal mesh %d patches; uniform refinement to level %d would be %d\n', ...
    sum(out.history(end).npatches), lev, length(dom0.x)*4^lev);

fprintf('stopped on: %s\n', out.reason);

% --- |B| and patch mesh at every level ---
% Patches are resampled from their n x n Chebyshev grid to M x M points so
% coarse patches do not render as facets. The torus is turned 135 degrees
% clockwise seen from above, which brings the ring's closest approach
% (azimuth 90 deg) round to azimuth -45 deg, facing the camera.
M = 13;
Im = barymat(linspace(-1, 1, M).', chebpts(n), chebtech2.barywts(n));
edges = {{1,':'}, {M,':'}, {':',1}, {':',M}};
Rz = [cosd(-135) -sind(-135); sind(-135) cosd(-135)];
lift = 0.004;   % patch lines sit just off the surface so it cannot hide them

Bmax = 0;
for k = 0:lev
    L = load(sprintf('%s_lev%d.mat', base, k));
    Bmax = max(Bmax, max(sqrt(sum(L.Bre.^2 + L.Bim.^2, 4)), [], 'all'));
end

% Axes are placed by hand with fixed inner positions, so every torus gets
% the same plot box; TILEDLAYOUT resized tiles unevenly.
ncol = ceil(sqrt(lev+1));
nrow = ceil((lev+1)/ncol);
w = 0.9/ncol; h = 0.94/nrow;
fig = figure('Visible', 'off', 'Position', [0 0 900*ncol 650*nrow]);
theme(fig, 'light');
ax = gobjects(lev+1, 1);
for k = 0:lev
    L = load(sprintf('%s_lev%d.mat', base, k));
    r = floor(k/ncol); c = mod(k, ncol);
    ax(k+1) = axes(fig, 'PositionConstraint', 'innerposition', ...
        'InnerPosition', [c*w + 0.005, 0.94 - (r+1)*h + 0.005, w - 0.01, h - 0.08]);
    hold on
    for i = 1:L.npatch
        X = Im*L.xx(:,:,i)*Im.'; Y = Im*L.yy(:,:,i)*Im.'; Z = Im*L.zz(:,:,i)*Im.';
        Bmag = zeros(M);
        for c = 1:3
            Bmag = Bmag + abs(Im*complex(L.Bre(:,:,i,c), L.Bim(:,:,i,c))*Im.').^2;
        end
        XY = Rz*[X(:).'; Y(:).'];
        X = reshape(XY(1,:), M, M); Y = reshape(XY(2,:), M, M);
        surf(X, Y, Z, sqrt(Bmag), 'EdgeColor', 'none', 'FaceColor', 'interp', ...
            'FaceLighting', 'gouraud', 'AmbientStrength', 0.5, ...
            'SpecularStrength', 0.1);
        % outward from the magnetic axis (the circle R = 2, z = 0)
        R = hypot(X, Y);
        nx = X - 2*X./R; ny = Y - 2*Y./R;
        nn = sqrt(nx.^2 + ny.^2 + Z.^2);
        for e = edges
            ix = e{1};
            plot3(X(ix{:}) + lift*nx(ix{:})./nn(ix{:}), ...
                  Y(ix{:}) + lift*ny(ix{:})./nn(ix{:}), ...
                  Z(ix{:}) + lift*Z(ix{:})./nn(ix{:}), ...
                  'Color', [0.1 0.1 0.12], 'LineWidth', 0.7);
        end
    end
    axis equal off
    view(45, 21.3);
    camlight headlight
    colormap(parula);
    clim([0 Bmax]);
    title({sprintf('level %d: %d patches', k, L.npatch), ...
        sprintf('error vs B_0: %.3e (L^\\infty), %.3e (L^2)', ...
        out.history(k+1).err_vs_B0, out.history(k+1).l2err_vs_B0)}, ...
        'FontSize', 14);
end
set(ax, 'XLim', ax(1).XLim, 'YLim', ax(1).YLim, 'ZLim', ax(1).ZLim);
cb = colorbar(ax(end), 'Position', [0.93 0.1 0.012 0.75]);
cb.Label.String = '|B|';
cb.Label.FontSize = 16;
ttl = sprintf('frac = %g', frac);
annotation(fig, 'textbox', [0 0.95 0.9 0.05], 'String', ttl, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 18);
png = sprintf('%s_levels_matlab.png', base);
exportgraphics(fig, png, 'Resolution', 120);
close(fig);
fprintf('SAVED %s\n', png);

fprintf('DONE_SENTINEL %s\n', datestr(now));
