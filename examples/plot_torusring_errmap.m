% Per-patch error against B0 and per-patch indicator eta on one level of a
% sigma_amr_torusring.m run, plus eta against error per patch.
%
% Set the variables the run was made with (nv, rmax, frac),
% and level [last saved level]. Reruns sigma_amr_torusring.m with
% setup_only = true for the ring geometry and BASE, then writes
% BASE_errmap_lev<level>.png.
%
% The per-patch error is max |B - B0| over the patch's nodes and components,
% divided by max |B0| over the surface, so its maximum over patches is the
% err_vs_B0 that sigma_amr reports.

setup_only = true;
sigma_amr_torusring;

H = load([base '_hist.mat']);
if ~exist('level', 'var'); level = H.hist(end).level; end
L = load(sprintf('%s_lev%d.mat', base, level));
eta = H.hist(level+1).eta{1};
depth = L.depth;

B = complex(L.Bre, L.Bim);
B0 = zeros(size(B));
for i = 1:L.npatch
    for a = 1:L.n
        for b = 1:L.n
            B0(a,b,i,:) = reftaylor(ntheta, rmin, rmaj, jmag, zk, ...
                [L.xx(a,b,i); L.yy(a,b,i); L.zz(a,b,i)]);
        end
    end
end
err = reshape(max(abs(B - B0), [], [1 2 4]), [], 1) / max(abs(B0), [], 'all');
fprintf('max per-patch error = %.4e (CSV err_vs_B0 = %.4e)\n', max(err), ...
    H.hist(level+1).err_vs_B0);

[~, worst] = sort(err, 'descend');
[~, ordeta] = sort(eta, 'descend');
rk = zeros(size(eta)); rk(ordeta) = 1:numel(eta);
fprintf('\n  patch  depth   err vs B0        eta   eta rank\n');
for i = worst(1:8).'
    fprintf('  %5d  %5d  %10.3e  %9.3e  %4d/%d\n', i, depth(i), err(i), ...
        eta(i), rk(i), numel(eta));
end

% --- figure ---
M = 13;
Im = barymat(linspace(-1, 1, M).', chebpts(L.n), chebtech2.barywts(L.n));
edges = {{1,':'}, {M,':'}, {':',1}, {':',M}};
Rz = [cosd(-135) -sind(-135); sind(-135) cosd(-135)];
lift = 0.004;
okabeito = [0.000 0.447 0.698; 0.902 0.624 0.000; 0.000 0.620 0.451; ...
            0.835 0.369 0.000; 0.800 0.475 0.655];

fig = figure('Visible', 'off', 'Position', [0 0 2100 700]);
theme(fig, 'light');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
vals = {err, eta};
names = {'log_{10} error vs B_0 (per patch)', 'log_{10} \eta (tail, per patch)'};
for p = 1:2
    ax = nexttile(tl);
    hold(ax, 'on');
    for i = 1:L.npatch
        X = Im*L.xx(:,:,i)*Im.'; Y = Im*L.yy(:,:,i)*Im.'; Z = Im*L.zz(:,:,i)*Im.';
        XY = Rz*[X(:).'; Y(:).'];
        X = reshape(XY(1,:), M, M); Y = reshape(XY(2,:), M, M);
        surf(ax, X, Y, Z, log10(vals{p}(i))*ones(M), 'EdgeColor', 'none', ...
            'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
            'AmbientStrength', 0.6, 'SpecularStrength', 0.1);
        R = hypot(X, Y);
        nx = X - 2*X./R; ny = Y - 2*Y./R;
        nn = sqrt(nx.^2 + ny.^2 + Z.^2);
        for e = edges
            ix = e{1};
            plot3(ax, X(ix{:}) + lift*nx(ix{:})./nn(ix{:}), ...
                  Y(ix{:}) + lift*ny(ix{:})./nn(ix{:}), ...
                  Z(ix{:}) + lift*Z(ix{:})./nn(ix{:}), ...
                  'Color', [0.1 0.1 0.12], 'LineWidth', 0.5);
        end
    end
    axis(ax, 'equal', 'off');
    view(ax, 45, 21.3);
    camlight(ax, 'headlight');
    colormap(ax, parula);
    cb = colorbar(ax, 'southoutside');
    cb.FontSize = 12;
    title(ax, names{p}, 'FontSize', 16);
end

ax = nexttile(tl);
hold(ax, 'on');
for dd = unique(depth).'
    i = depth == dd;
    plot(ax, eta(i), err(i), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', okabeito(mod(dd, 5) + 1, :), ...
        'MarkerEdgeColor', 'none', 'DisplayName', sprintf('depth %d', dd));
end
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 13);
grid(ax, 'on');
xlabel(ax, '\eta (tail)');
ylabel(ax, 'error vs B_0');
legend(ax, 'Location', 'southeast');
title(ax, 'per patch', 'FontSize', 16);

title(tl, sprintf('frac = %g, level %d: %d patches, error vs B_0 = %.3e', ...
    frac, level, L.npatch, max(err)), 'FontSize', 18, ...
    'FontWeight', 'bold');
png = sprintf('%s_errmap_lev%d.png', base, level);
exportgraphics(fig, png, 'Resolution', 120);
close(fig);
fprintf('SAVED %s\nDONE_SENTINEL\n', png);
