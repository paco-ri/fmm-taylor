% Pointwise error of B against B0 on both surfaces at every level of a
% sigma_amr_shellrings.m run.
%
% Set the variables the run was made with (nv, frac).
% Reruns sigma_amr_shellrings.m with setup_only = true for the rings and
% BASE, then writes BASE_errmap.png: one row per level, outer surface left,
% inner right, on a shared colour scale.
%
% The error is max over components of |B - B0|, divided by max |B0| over
% both surfaces, so its maximum over the nodes is the err_vs_B0 that
% sigma_amr reports. Nodal errors are interpolated to the plotting grid.

setup_only = true;
sigma_amr_shellrings;

H = load([base '_hist.mat']);
nlev = numel(H.hist);
L = cell(nlev, 2); E = cell(nlev, 2);
for k = 1:nlev
    for s = 1:2
        L{k,s} = load(sprintf('%s_lev%d_surf%d.mat', base, k-1, s));
        X = [L{k,s}.xx(:).'; L{k,s}.yy(:).'; L{k,s}.zz(:).'];
        B0 = ringfield(rings, zk, X);
        B = reshape(complex(L{k,s}.Bre, L{k,s}.Bim), [], 3).';
        E{k,s} = reshape((B - B0).', size(L{k,s}.Bre));
        L{k,s}.B0max = max(abs(B0), [], 'all');
    end
end

M = 13;
Im = barymat(linspace(-1, 1, M).', chebpts(n), chebtech2.barywts(n));
edges = {{1,':'}, {M,':'}, {':',1}, {':',M}};
lift = 0.004;

errs = cell(nlev, 2);
for k = 1:nlev
    den = max(L{k,1}.B0max, L{k,2}.B0max);
    for s = 1:2
        errs{k,s} = zeros(M, M, L{k,s}.npatch);
        for i = 1:L{k,s}.npatch
            for c = 1:3
                errs{k,s}(:,:,i) = max(errs{k,s}(:,:,i), ...
                    abs(Im*E{k,s}(:,:,i,c)*Im.'));
            end
        end
        errs{k,s} = errs{k,s}/den;
    end
    nodemax = max(cellfun(@(e) max(abs(e), [], 'all'), E(k,:)))/den;
    fprintf('level %d: max nodal error %.4e (CSV err_vs_B0 = %.4e)\n', ...
        k-1, nodemax, H.hist(k).err_vs_B0);
end
hi = ceil(log10(max(cellfun(@(e) max(e, [], 'all'), errs), [], 'all')));
lo = hi - 5;

% Axes are placed by hand with one fixed camera per surface; tiledlayout
% resizes 3D axes unevenly.
fig = figure('Visible', 'off', 'Position', [0 0 1900 620*nlev]);
theme(fig, 'light');
hr = 0.94/nlev;
for k = 1:nlev
    for s = 1:2
        Ls = L{k,s};
        ang = -45 - atan2d(rings(s).pclose(2), rings(s).pclose(1));
        Rz = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
        % look from below when the closest approach is on the underside
        el = 21.3;
        if rings(s).pclose(3) < 0; el = -40; end
        pos = [0.01 + 0.45*(s-1), 0.98 - k*hr, 0.44, hr - 0.03];
        ax = axes(fig, 'Position', pos, 'Clipping', 'off');
        hold(ax, 'on');
        for i = 1:Ls.npatch
            X = Im*Ls.xx(:,:,i)*Im.'; Y = Im*Ls.yy(:,:,i)*Im.'; Z = Im*Ls.zz(:,:,i)*Im.';
            XY = Rz*[X(:).'; Y(:).'];
            X = reshape(XY(1,:), M, M); Y = reshape(XY(2,:), M, M);
            surf(ax, X, Y, Z, log10(max(errs{k,s}(:,:,i), 10^lo)), ...
                'EdgeColor', 'none', 'FaceColor', 'interp', ...
                'FaceLighting', 'gouraud', 'AmbientStrength', 0.6, ...
                'SpecularStrength', 0.1);
            % outward from the magnetic axis (the circle R = rmaj, z = 0)
            R = hypot(X, Y);
            nx = X - rmaj*X./R; ny = Y - rmaj*Y./R;
            nn = sqrt(nx.^2 + ny.^2 + Z.^2);
            for e = edges
                ix = e{1};
                plot3(ax, X(ix{:}) + lift*nx(ix{:})./nn(ix{:}), ...
                      Y(ix{:}) + lift*ny(ix{:})./nn(ix{:}), ...
                      Z(ix{:}) + lift*Z(ix{:})./nn(ix{:}), ...
                      'Color', [0.1 0.1 0.12], 'LineWidth', 0.5);
            end
        end
        rl = max(hypot(L{1,s}.xx, L{1,s}.yy), [], 'all');
        zl = max(abs(L{1,s}.zz), [], 'all');
        axis(ax, [-rl rl -rl rl -zl zl]);
        daspect(ax, [1 1 1]);
        axis(ax, 'off');
        view(ax, 45, el);
        camtarget(ax, [0 0 0]);
        camva(ax, 7);
        camlight(ax, 'headlight');
        colormap(ax, parula);
        clim(ax, [lo hi]);
        annotation(fig, 'textbox', [pos(1), pos(2) + pos(4), pos(3), 0.02], ...
            'String', sprintf('level %d, %s surface: %d patches, max error %.3e', ...
            k-1, sname{s}, Ls.npatch, max(errs{k,s}, [], 'all')), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight', 'bold');
    end
end
cb = colorbar(ax, 'Position', [0.925 0.1 0.012 0.8]);
ax.Position = pos;
cb.Label.String = 'log_{10} |B - B_0| / max |B_0|';
cb.Label.FontSize = 16;
cb.FontSize = 13;
annotation(fig, 'textbox', [0 0.975 1 0.025], 'String', sprintf(['frac = %g: ' ...
    'error vs B_0 per level (each surface turned to face its ring''s ' ...
    'closest approach)'], frac), 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');
png = [base '_errmap.png'];
exportgraphics(fig, png, 'Resolution', 90);
close(fig);
fprintf('SAVED %s\nDONE_SENTINEL\n', png);
