% Sigma-driven iterative mesh refinement on the two-surface toroidal shell of
% sigma_amr_torusshell.m, against a reference Taylor state made by two
% current rings, one on each side of the shell.
%
% Same refinement setup as sigma_amr_torusring.m. The differences:
%   - geometry is the two-surface shell (outer minor radius 1.0, inner 0.6)
%   - B0 is the sum of two reftaylor-type rings, neither linking the tori:
%       exterior  horizontal, radius 0.5, beside the outboard side of the
%                 outer torus at azimuth 90 deg
%       interior  vertical (normal radial), inside the inner tube at
%                 azimuth 180 deg
%     Each ring passes about 0.2 from the nearer surface.
%   - fluxes are toroidal (surface 1 cross-section) and poloidal (surface 2)
%
% Variables read if set before running:
%   setup_only [false]      stop after the geometry check plot
%   nv         [3]          base-mesh patches poloidally (nu = 3*nv)
%   rmax       [3]          maximum refinement depth
%   frac       [0.3]        opts.frac of sigma_amr
%   figext     ['pdf']      format of the figures, 'pdf' or 'png'
%
% With BASE = sigma_amr_shellrings_fraction, suffixed _nv<nv> for a
% non-default nv and _frac<100*frac> for a non-default frac, writes
% BASE_setup.{<figext>,fig} (the geometry check plot), BASE_{outer,inner}_geometry.<figext>
% (surfaces and rings only, oriented as the |B| figures), BASE.csv, per-level sigma VTK, per-level BASE_lev<k>_surf<s>.mat, and
% BASE_{outer,inner}_levels_matlab.<figext> and
% BASE_{outer,inner}_errmap_levels_matlab.<figext> showing |B| and the error
% |B - B0| with the patch mesh on that surface at every level, and
% BASE_conv.<figext> with the L^inf and L^2 relative errors per level.

if ~exist('setup_only', 'var'); setup_only = false; end
if ~exist('nv', 'var'); nv = 3; end
if ~exist('rmax', 'var'); rmax = 3; end
if ~exist('frac', 'var'); frac = 0.3; end
if ~exist('figext', 'var'); figext = 'pdf'; end
here = fileparts(mfilename('fullpath'));
base = fullfile(here, 'sigma_amr_shellrings_fraction');
if nv ~= 3; base = sprintf('%s_nv%d', base, nv); end
if frac ~= 0.3; base = sprintf('%s_frac%g', base, 100*frac); end
fprintf('=== sigma_amr_shellrings start %s\n', datestr(now));

% --- geometry (sigma_amr_torusshell.m; outer surface is dom0{1}) ---
n  = 8;
nu = 3*nv;
ao = 1.0;
ai = 0.6;
nr = 16; nt = 40; np = 40;   % cross-section quadrature for the fluxes
[dom0, qnodes, qweights] = prepare_torus(n, nu, nv, n, nu, nv, ao, ai, ...
    nr, nt, np);
domparams = [n, nu, nv];
fprintf('base meshes: %d + %d patches, n = %d\n', ...
    length(dom0{1}.x), length(dom0{2}.x), n);

% --- current rings ---
% In a ring's own frame p' = Q.'*(p - c), the ring is the circle of radius a
% about the p'_1 axis, so the distance from p to it is
% sqrt(p'_1^2 + (rho - a)^2) with rho = hypot(p'_2, p'_3). The closest
% surface point is found on a fine sampling of the surfaces.
rmaj = 2.0;
jmag = 1.0;
zk   = 0.5;
tol  = 1e-6;

nuf = 36;           % multiple of 4, so u = pi/2 and u = pi are patch edges
domf = prepare_torus(16, nuf, 12, 16, nuf, 12, ao, ai, 2, 2, 2);
P = cell(1, 2); sl90 = cell(1, 2); sl180 = cell(1, 2);
for s = 1:2
    P{s} = zeros(3, 0); sl90{s} = zeros(2, 0); sl180{s} = zeros(2, 0);
    for i = 1:length(domf{s}.x)
        xi = domf{s}.x{i}; yi = domf{s}.y{i}; zi = domf{s}.z{i};
        P{s} = [P{s}, [xi(:).'; yi(:).'; zi(:).']];
        if max(abs(xi(:,1))) < 1e-12 && yi(1,1) > 0
            sl90{s} = [sl90{s}, [yi(:,1).'; zi(:,1).']];
        end
        if max(abs(yi(:,1))) < 1e-12 && xi(1,1) < 0
            sl180{s} = [sl180{s}, [-xi(:,1).'; zi(:,1).']];
        end
    end
    % order each cross-section by poloidal angle so it plots as a closed curve
    [~, o] = sort(atan2(sl90{s}(2,:), sl90{s}(1,:) - rmaj));
    sl90{s} = sl90{s}(:, [o o(1)]);
    [~, o] = sort(atan2(sl180{s}(2,:), sl180{s}(1,:) - rmaj));
    sl180{s} = sl180{s}(:, [o o(1)]);
end

% Q maps reftaylor's frame (ring normal e_x) to the ring's. rings(k) is
% closest to surface k, and neither ring links the tori.
%   rings(1), exterior: horizontal (normal e_z), radius ae, centred at
%     (0, Re, 0) beside the outboard side of the outer torus.
%   rings(2), interior: vertical (normal e_x, radial at azimuth 180 deg),
%     radius ai_ring, centred mid-tube in the inner u = pi cross-section at
%     (-xi, 0, 0).
ringdist = @(Q, c, a, X) sqrt((Q(:,1).'*(X - c)).^2 + ...
    (vecnorm(Q(:,2:3).'*(X - c)) - a).^2);
Qe = [0 1 0; 0 0 1; 1 0 0];
ae = 0.50;
Re = 3.59;
Qi = eye(3);
xi = 1.96;
ci = [-xi; 0; 0];
ai_ring = 0.54;
rings = struct('c', {[0; Re; 0], ci}, 'Q', {Qe, Qi}, 'a', {ae, ai_ring}, ...
               'ntheta', {0, 0}, 'jmag', {jmag, jmag}, 'pclose', {[], []});
for k = 1:2
    [dk, ik] = min(ringdist(rings(k).Q, rings(k).c, rings(k).a, P{k}));
    rings(k).pclose = P{k}(:, ik);
    % The trapezoidal rule along the ring has error ~exp(-2*pi*d/h) at
    % distance d, h the ring point spacing; h <= dk/4 puts it below 1e-10 on
    % the surface.
    rings(k).ntheta = max(200, ceil(8*pi*rings(k).a/dk));
    fprintf(['ring %d: a = %.6f, ntheta = %d, min distance to surface %d = ' ...
        '%.6f at (%.4f, %.4f, %.4f)\n'], k, rings(k).a, rings(k).ntheta, ...
        k, dk, rings(k).pclose);
end
if ~inpolygon(-ci(1), ai_ring, sl180{2}(1,:), sl180{2}(2,:))
    error('interior ring is not inside the inner tube');
end

% --- geometry check plot ---
th = linspace(0, 2*pi, 400);
rpts = cell(1, 2);
for k = 1:2
    rpts{k} = rings(k).c + rings(k).Q*[zeros(size(th)); ...
        rings(k).a*cos(th); rings(k).a*sin(th)];
end
cring = [0.835 0.369 0.000; 0.000 0.447 0.698];   % Okabe-Ito vermillion, blue
csurf = [0.55 0.55 0.55; 0.902 0.624 0.000];       % grey, Okabe-Ito orange
salpha = [0.12 0.35];
sname = {'outer', 'inner'};
fig = figure('Visible', 'off', 'Position', [0 0 1800 650]);
theme(fig, 'light');

subplot(1, 3, 1); hold on
for s = 1:2
    for i = 1:length(dom0{s}.x)
        surf(dom0{s}.x{i}, dom0{s}.y{i}, dom0{s}.z{i}, 'FaceColor', csurf(s,:), ...
            'FaceAlpha', salpha(s), 'EdgeColor', 'none');
    end
end
for k = 1:2
    plot3(rpts{k}(1,:), rpts{k}(2,:), rpts{k}(3,:), 'Color', cring(k,:), ...
        'LineWidth', 2.5);
    plot3(rings(k).pclose(1), rings(k).pclose(2), rings(k).pclose(3), 'o', ...
        'MarkerFaceColor', cring(k,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
end
axis equal; view(-30, 35); grid on
xlabel('x'); ylabel('y'); zlabel('z');
title('base meshes (outer grey, inner orange) and current rings');

subplot(1, 3, 2); hold on
plot(sl90{1}(1,:), sl90{1}(2,:), 'k-', 'LineWidth', 1.5);
plot(sl90{2}(1,:), sl90{2}(2,:), 'k--', 'LineWidth', 1.5);
plot(Re + ae*[-1 1], [0 0], 'x', 'Color', cring(1,:), 'MarkerSize', 12, ...
    'LineWidth', 2);
plot(rings(1).pclose(2), rings(1).pclose(3), 'o', 'MarkerFaceColor', ...
    cring(1,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
axis equal; grid on
xlabel('y'); ylabel('z');
legend({'outer surface', 'inner surface', 'exterior ring crossings', ...
    'closest outer point'}, 'Location', 'southoutside');
title(sprintf('plane x = 0 (u = \\pi/2): exterior ring a = %.2f, centre y = %.4f', ...
    ae, Re));

subplot(1, 3, 3); hold on
plot(sl180{1}(1,:), sl180{1}(2,:), 'k-', 'LineWidth', 1.5);
plot(sl180{2}(1,:), sl180{2}(2,:), 'k--', 'LineWidth', 1.5);
plot(-ci(1)*[1 1], ai_ring*[-1 1], 'x', 'Color', cring(2,:), ...
    'MarkerSize', 12, 'LineWidth', 2);
plot(hypot(rings(2).pclose(1), rings(2).pclose(2)), rings(2).pclose(3), 'o', ...
    'MarkerFaceColor', cring(2,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
axis equal; grid on
xlabel('-x'); ylabel('z');
legend({'outer surface', 'inner surface', 'interior ring crossings', ...
    'closest inner point (rotated into plane)'}, 'Location', 'southoutside');
title(sprintf('plane y = 0, x < 0 (u = \\pi): interior ring a = %.4f', ai_ring));

exportgraphics(fig, [base '_setup.' figext], 'Resolution', 150);
savefig(fig, [base '_setup.fig']);
fprintf('SAVED %s_setup.%s / .fig\n', base, figext);

% --- surfaces and rings only, in each orientation of the |B| figures ---
% The outer surface is see-through so the inner one and its ring show.
% Figure t turns everything about z so that ring t's closest approach faces
% the camera (azimuth -45 deg), as the |B| figures for surface t do.
for t = 1:2
    ang = -45 - atan2d(rings(t).pclose(2), rings(t).pclose(1));
    Rz = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [0 0 3.25 2.4]);
    theme(fig, 'light');
    axes(fig, 'Position', [0 0 1 1]); hold on
    for s = 1:2
        for i = 1:length(domf{s}.x)
            XY = Rz*[domf{s}.x{i}(:).'; domf{s}.y{i}(:).'];
            surf(reshape(XY(1,:), size(domf{s}.x{i})), ...
                reshape(XY(2,:), size(domf{s}.x{i})), domf{s}.z{i}, ...
                'FaceColor', csurf(s,:), 'FaceAlpha', salpha(s), ...
                'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
                'AmbientStrength', 0.4, 'SpecularStrength', 0.2);
        end
    end
    for k = 1:2
        XY = Rz*rpts{k}(1:2,:);
        plot3(XY(1,:), XY(2,:), rpts{k}(3,:), 'Color', cring(k,:), 'LineWidth', 1.5);
    end
    axis equal off
    view(15, 30);   % the |B| figures' view(45, 21.3) sees the rings edge-on
    camzoom(1.75);   % the default view angle leaves wide margins
    camlight headlight
    figfile = sprintf('%s_%s_geometry.%s', base, sname{t}, figext);
    if strcmp(figext, 'pdf')
        % vector output would drop the gouraud shading
        exportgraphics(fig, figfile, 'ContentType', 'image', 'Resolution', 600);
    else
        exportgraphics(fig, figfile, 'Resolution', 200);
    end
    close(fig);
    fprintf('SAVED %s\n', figfile);
end

if setup_only
    fprintf('DONE_SENTINEL %s\n', datestr(now));
    return
end

% --- reference Taylor state ---
% Signs follow sigma_amr_torusshell.m: toroidal flux on surface 1 takes
% +B0(2), poloidal on surface 2 takes -B0(3).
t0 = tic;
B0q1 = ringfield(rings, zk, qnodes{1});
B0q2 = ringfield(rings, zk, qnodes{2});
flux = [sum(B0q1(2,:).*qweights{1}), -sum(B0q2(3,:).*qweights{2})];
fprintf('flux = [%.16e, %.16e]  (%.1f s)\n', flux(1), flux(2), toc(t0));

refstate = [];
refstate.xs_nodes   = qnodes;
refstate.xs_weights = qweights;
refstate.B0fun = @(d, s) ringfield(rings, zk, d);

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

lev = out.history(end).level;
fprintf('\nfinal mesh %d patches; uniform refinement to level %d would be %d\n', ...
    sum(out.history(end).npatches), lev, ...
    (length(dom0{1}.x) + length(dom0{2}.x))*4^lev);

fprintf('stopped on: %s\n', out.reason);

% --- |B|, error vs B0 and patch mesh at every level, per surface ---
% Patches are resampled from their n x n Chebyshev grid to M x M points so
% coarse patches do not render as facets. Each surface is turned about z so
% that its ring's closest approach faces the camera (azimuth -45 deg).
M = 13;
Im = barymat(linspace(-1, 1, M).', chebpts(n), chebtech2.barywts(n));
edges = {{1,':'}, {M,':'}, {':',1}, {':',M}};
lift = 0.004;   % patch lines sit just off the surface so it cannot hide them

% dB = B - B0 at the saved nodes
Ls = cell(lev+1, 2);
for k = 0:lev
    for s = 1:2
        L = load(sprintf('%s_lev%d_surf%d.mat', base, k, s));
        L.B = complex(L.Bre, L.Bim);
        L.dB = L.B - reshape(ringfield(rings, zk, ...
            [L.xx(:).'; L.yy(:).'; L.zz(:).']).', size(L.B));
        Ls{k+1, s} = L;
    end
end
emax = max(cellfun(@(L) max(vecnorm(L.dB, 2, 4), [], 'all'), Ls), [], 'all');
eclim = ceil(log10(emax)) + [-6 0];

% Axes are placed by hand with fixed inner positions, so every torus gets
% the same plot box; TILEDLAYOUT resized tiles unevenly. Sizes are in inches
% for a 6.5 in text width, so fonts print at their nominal size.
ncol = ceil(sqrt(lev+1));
nrow = ceil((lev+1)/ncol);
tilew = 2.85; tileh = 2.0; cbw = 0.8;
figname = {'levels', 'errmap_levels'};
for s = 1:2
    ang = -45 - atan2d(rings(s).pclose(2), rings(s).pclose(1));
    Rz = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
    Bmax = max(cellfun(@(L) max(vecnorm(L.B, 2, 4), [], 'all'), Ls(:, s)));
    for q = 1:2   % 1: |B|, 2: log10 |B - B0|
        fig = figure('Visible', 'off', 'Units', 'inches', ...
            'Position', [0 0 ncol*tilew + cbw, nrow*tileh]);
        theme(fig, 'light');
        ax = gobjects(lev+1, 1);
        for k = 0:lev
            L = Ls{k+1, s};
            r = floor(k/ncol); c = mod(k, ncol);
            ax(k+1) = axes(fig, 'Units', 'inches', 'Clipping', 'off', ...
                'PositionConstraint', 'innerposition', ...
                'InnerPosition', [c*tilew, (nrow-r-1)*tileh, tilew, tileh - 0.2]);
            hold on
            for i = 1:L.npatch
                X = Im*L.xx(:,:,i)*Im.'; Y = Im*L.yy(:,:,i)*Im.'; Z = Im*L.zz(:,:,i)*Im.';
                if q == 1, F = L.B; else, F = L.dB; end
                C = zeros(M);
                for j = 1:3
                    C = C + abs(Im*F(:,:,i,j)*Im.').^2;
                end
                C = sqrt(C);
                if q == 2, C = log10(max(C, 10^eclim(1))); end
                XY = Rz*[X(:).'; Y(:).'];
                X = reshape(XY(1,:), M, M); Y = reshape(XY(2,:), M, M);
                surf(X, Y, Z, C, 'EdgeColor', 'none', 'FaceColor', 'interp', ...
                    'FaceLighting', 'gouraud', 'AmbientStrength', 0.5, ...
                    'SpecularStrength', 0.1);
                % outward from the magnetic axis (the circle R = rmaj, z = 0)
                R = hypot(X, Y);
                nx = X - rmaj*X./R; ny = Y - rmaj*Y./R;
                nn = sqrt(nx.^2 + ny.^2 + Z.^2);
                for e = edges
                    ix = e{1};
                    plot3(X(ix{:}) + lift*nx(ix{:})./nn(ix{:}), ...
                          Y(ix{:}) + lift*ny(ix{:})./nn(ix{:}), ...
                          Z(ix{:}) + lift*Z(ix{:})./nn(ix{:}), ...
                          'Color', [0.1 0.1 0.12], 'LineWidth', 0.5);
                end
            end
            axis equal off
            view(45, 21.3);
            camzoom(ax(k+1), 2.2);   % the default view angle leaves wide margins
            camlight headlight
            colormap(parula);
            if q == 1, clim([0 Bmax]); else, clim(eclim); end
            % TITLE would sit above the zoomed plot box, off the tile
            annotation(fig, 'textbox', 'Units', 'inches', ...
                'Position', [c*tilew, (nrow-r)*tileh - 0.2, tilew, 0.2], ...
                'String', sprintf('level %d', k), 'FontSize', 10, 'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        set(ax, 'XLim', ax(1).XLim, 'YLim', ax(1).YLim, 'ZLim', ax(1).ZLim);
        cb = colorbar(ax(end), 'Units', 'inches', 'FontSize', 9, ...
            'Position', [ncol*tilew + 0.1, 0.3, 0.12, nrow*tileh - 0.8]);
        if q == 1
            cb.Title.String = '|B|';
            cb.Title.FontSize = 10;
            cb.Title.FontWeight = 'normal';
        else
            cb.Ticks = eclim(1):eclim(2);
            cb.TickLabels = compose('10^{%d}', cb.Ticks);
            cb.Label.String = '|B - B_0|';
            cb.Label.FontSize = 10;
        end
        figfile = sprintf('%s_%s_%s_matlab.%s', base, sname{s}, figname{q}, figext);
        if strcmp(figext, 'pdf')
            % vector output would drop the gouraud shading
            exportgraphics(fig, figfile, 'ContentType', 'image', 'Resolution', 600);
        else
            exportgraphics(fig, figfile, 'Resolution', 200);
        end
        close(fig);
        fprintf('SAVED %s\n', figfile);
    end
end

% --- relative error vs B0 against the number of points ---
npts = arrayfun(@(h) sum(h.npts), out.history);
fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [0 0 3.25 2.4]);
theme(fig, 'light');
cblue = [0.000 0.447 0.698]; cverm = [0.835 0.369 0.000];   % Okabe-Ito
loglog(npts, [out.history.err_vs_B0], 'o-', 'Color', cblue, ...
    'MarkerFaceColor', cblue, 'MarkerSize', 4, 'LineWidth', 1);
hold on
loglog(npts, [out.history.l2err_vs_B0], 's-', 'Color', cverm, ...
    'MarkerFaceColor', cverm, 'MarkerSize', 4, 'LineWidth', 1);
grid on
set(gca, 'FontSize', 9, 'XTick', npts, 'XTickLabel', string(npts), ...
    'XMinorTick', 'off', 'XMinorGrid', 'off');
xlim(npts([1 end]).*[0.85 1.15]);
xlabel('number of points');
ylabel('relative error vs B_0');
legend({'L^\infty', 'L^2'}, 'Location', 'southwest');
figfile = sprintf('%s_conv.%s', base, figext);
if strcmp(figext, 'pdf')
    exportgraphics(fig, figfile, 'ContentType', 'vector');
else
    exportgraphics(fig, figfile, 'Resolution', 200);
end
close(fig);
fprintf('SAVED %s\n', figfile);

fprintf('DONE_SENTINEL %s\n', datestr(now));
