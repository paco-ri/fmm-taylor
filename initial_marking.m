function [marked, leg, tail] = initial_marking(dom, nmark, plotpatch)
%INITIAL_MARKING Choose the leaves to refine on the first SIGMA_AMR pass.
%   MARKED = INITIAL_MARKING(DOM, NMARK) takes the level-0 mesh (a cell, one
%   per surface) and marks the NMARK patches, over all surfaces, with the
%   largest tail norm norm(TAIL{s}(i,:)) (see below). MARKED is a cell, one
%   per surface, of sorted column vectors of patch indices into DOM{s}.
%   NMARK larger than the number of patches marks every patch.
%
%   [MARKED, LEG] = INITIAL_MARKING(...) also returns the Legendre
%   coefficients of the geometry: LEG{s}{i}(:,:,c) holds those of coordinate
%   c (x, y, z) on patch i, entry (k,l) multiplying P_{k-1}(v) P_{l-1}(u),
%   with P_k normalized so that P_k(1) = 1.
%
%   [MARKED, LEG, TAIL] = INITIAL_MARKING(...) also returns TAIL{s}(i,c),
%   the 2-norm of the coefficients in LEG{s}{i}(:,:,c) of total degree
%   (k-1) + (l-1) >= n-1, i.e. on and above the anti-diagonal.
%
%   INITIAL_MARKING(..., PLOTPATCH) with PLOTPATCH = [s, i] plots
%   log10|LEG{s}{i}| for each coordinate. The default, false, plots nothing.
%
%   See also SIGMA_AMR.

if nargin < 3, plotpatch = false; end
ns = numel(dom);

leg = cell(1, ns);
tail = cell(1, ns);
score = []; sid = []; pid = [];
for s = 1:ns
    n = size(dom{s}.x{1}, 1);
    % values at 2nd-kind Chebyshev points -> Legendre coefficients
    T = cheb2leg(chebtech2.vals2coeffs(eye(n)));
    % total degree (k-1) + (l-1) >= n-1: on and above the anti-diagonal
    [K, L] = ndgrid(1:n);
    intail = K + L >= n + 1;
    npat = length(dom{s}.x);
    leg{s} = cell(npat, 1);
    tail{s} = zeros(npat, 3);
    for i = 1:npat
        leg{s}{i} = cat(3, T*dom{s}.x{i}*T.', T*dom{s}.y{i}*T.', ...
            T*dom{s}.z{i}*T.');
        for c = 1:3
            C = leg{s}{i}(:,:,c);
            tail{s}(i, c) = norm(C(intail));
        end
    end
    score = [score; vecnorm(tail{s}, 2, 2)]; %#ok<AGROW>
    sid = [sid; s*ones(npat, 1)]; %#ok<AGROW>
    pid = [pid; (1:npat).']; %#ok<AGROW>
end

[~, ord] = sort(score, 'descend');
top = ord(1:min(nmark, end));
marked = cell(1, ns);
for s = 1:ns
    marked{s} = sort(pid(top(sid(top) == s)));
end

if ~isequal(plotpatch, false)
    s = plotpatch(1); i = plotpatch(2);
    C = leg{s}{i};
    n = size(C, 1);
    lc = log10(abs(C) + eps);
    figure;
    tl = tiledlayout(1, 3);
    names = {'x', 'y', 'z'};
    for c = 1:3
        nexttile;
        imagesc(0:n-1, 0:n-1, lc(:,:,c));
        axis image; set(gca, 'YDir', 'normal');
        xlabel('degree in u'); ylabel('degree in v');
        title(names{c});
        clim([min(lc, [], 'all'), max(lc, [], 'all')]);
    end
    colormap(parula);
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'log_{10} |coefficient|';
    title(tl, sprintf('Legendre coefficients, surface %d, patch %d', s, i));
end

end
