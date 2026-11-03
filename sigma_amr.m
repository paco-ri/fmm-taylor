function out = sigma_amr(dom0, domparams, zk, flux, tol, rmax, opts)
%SIGMA_AMR Iteratively refine a surface mesh until sigma is resolved.
%   OUT = SIGMA_AMR(DOM0, DOMPARAMS, ZK, FLUX, TOL, RMAX) solves the Taylor
%   state problem on a sequence of meshes, refining at each step the patches
%   on which the density sigma is least resolved, in the sense that the L2 norm 
%   of the part of sigma with Legendre total degree >= n-1 is large relative
%   to the L2 norm of sigma over the surface.
%
%   This is in contrast to SURFACEMESH/ADAP_REF, which refines the fundamental 
%   forms.
%
%   Arguments:
%     DOM0      a surfacemesh or a cell of one or two, the initial mesh(es)
%     DOMPARAMS [n, nu, nv]: Chebyshev nodes per patch in each direction,
%               and patches in the toroidal and poloidal directions of DOM0
%     ZK        Beltrami parameter lambda; 0 gives the vacuum field
%     FLUX      cross-sectional flux, one per surface
%     TOL       scalar, or [GMRES, grad/curl S_lambda quadrature,
%               S_lambda quadrature] tolerances
%     RMAX      maximum quadtree depth (each base patch can split into up to 
%               4^RMAX leaves)
%
%   OUT = SIGMA_AMR(..., OPTS) accepts a struct with fields (default values in 
%   brackets):
%
%     maxlevels [RMAX]        maximum number of refinements. Distinct from
%                             RMAX: a patch left unmarked for several levels
%                             can still be eligible to refine when the
%                             budget runs out.
%     frac      [0.3]         before every refinement, the CEIL(FRAC*N)
%                             leaves with the largest ETA are refined, N
%                             being the number of leaves shallower than
%                             RMAX
%     refstate  []            struct with fields xs_nodes and xs_weights
%                             (each a cell, one per surface) and B0fun, a
%                             handle @(dom, s) returning B0 on a mesh. When
%                             given, RefTaylorState is used and the true
%                             error against B0 is reported per level.
%     vtkbase   ''            if nonempty, write <base>_lev<k>_sigma.vtk
%     savebase  [false]       false, or a path prefix. When a prefix is
%                             given, B is evaluated on the surface at every
%                             level and saved, with the mesh nodes and the
%                             patch depths, to <base>_lev<k>.mat (one
%                             surface) or <base>_lev<k>_surf<s>.mat (two),
%                             in the layout examples/make_amr_vtk.py reads.
%     verbose   [true]        print progress
%
%   OUT is a struct with the final DOM, QF, P2Q, SIGMA and ALPHA, plus a
%   HISTORY struct array carrying, for each level, the patch counts, the
%   indicator statistics, the per-patch indicator ETA (a cell, one per
%   surface, indexed like that level's mesh), the number of patches marked,
%   and the solve time. With opts.refstate it also carries ERR_VS_B0,
%   max|B - B0| / max|B0| over the surface nodes, and L2ERR_VS_B0, the
%   relative L2 error ||B - B0|| / ||B0|| over the surface.
%   OUT.REASON says why the loop stopped:
%
%     'depth'     every leaf is already at depth RMAX
%     'maxlevels' the refinement budget ran out with leaves still eligible
%
%   See also TAYLORSTATE, SURFACEMESH/REFINE_LEAVES.

% -- options ---------------------------------------------------------------
if nargin < 7, opts = []; end
maxlevels = getopt(opts, 'maxlevels', rmax);
frac      = getopt(opts, 'frac', 0.3);
refstate  = getopt(opts, 'refstate', []);
vtkbase   = getopt(opts, 'vtkbase', '');
savebase  = getopt(opts, 'savebase', false);
verbose   = getopt(opts, 'verbose', true);

if ~isposint(rmax)
    error('SIGMA_AMR:rmax', 'RMAX must be a positive integer scalar.');
end
if ~isposint(maxlevels)
    error('SIGMA_AMR:maxlevels', ...
        'opts.maxlevels must be a positive integer scalar.');
end
if ~(isnumeric(frac) && isscalar(frac) && frac > 0 && frac <= 1)
    error('SIGMA_AMR:frac', 'opts.frac must be a scalar in (0, 1].');
end
if ~(isequal(savebase, false) || ...
        ((ischar(savebase) || isstring(savebase)) && strlength(savebase) > 0))
    error('SIGMA_AMR:savebase', ...
        'opts.savebase must be false or a nonempty path prefix.');
end
if ~isempty(refstate) && ~(isfield(refstate, 'B0fun') && ...
        isa(refstate.B0fun, 'function_handle'))
    error('SIGMA_AMR:refstate', ...
        'opts.refstate.B0fun is required and must be a function handle.');
end

if isa(dom0, 'surfacemesh')
    base = {dom0};
else
    base = dom0;
end
ns = numel(base);
% base = base mesh

% values at 2nd-kind Chebyshev points -> Legendre coefficients
n = size(base{1}.x{1}, 1);
T = cheb2leg(chebtech2.vals2coeffs(eye(n)));
% total degree (k-1) + (l-1) >= n-1: on and above the anti-diagonal
[K, L] = ndgrid(1:n);
intail = K + L >= n + 1;

p2q = cell(1, ns);
for s = 1:ns
    np = length(base{s}.x);
    p2q{s} = [(1:np).', zeros(np, 2)];
end
dom = cell(1, ns); qf = cell(1, ns); marked = cell(1, ns);

% reason for stopping refinement
reason = 'maxlevels';
for lev = 0:maxlevels
    % -- refine (at level 0 this only builds the mesh and quadforest) -------
    for s = 1:ns
        [dom{s}, qf{s}, p2q{s}] = surfacemesh.refine_leaves( ...
            base{s}, p2q{s}, marked{s}, rmax);
    end
    if verbose
        fprintf('[sigma_amr] level %d: %s patches\n', lev, strjoin( ...
            cellfun(@(d) sprintf('%d', length(d.x)), dom, ...
            'UniformOutput', false), '+'));
    end

    % -- solve -------------------------------------------------------------
    t0 = tic;
    D = Domain(dom, domparams, qf, p2q);
    if isempty(refstate)
        ts = TaylorState(D, domparams, zk, flux, tol);
    else
        B0 = cell(1, ns);
        for s = 1:ns
            B0{s} = refstate.B0fun(dom{s}, s);
        end
        ts = RefTaylorState(D, domparams, zk, flux, B0, ...
            refstate.xs_nodes, refstate.xs_weights, tol);
    end
    ts = ts.solve(false);
    sigma = ts.sigma;
    tsolve = toc(t0);

    % -- compute tail error measure -----------------------------------------
    eta = cell(1, ns);
    for s = 1:ns
        npat = length(dom{s}.x);
        e = zeros(npat, 1);
        g = 0;
        for j = 1:npat
            J = dom{s}.J{j}; % Jacobian on patch j
            v = sigma{s}.vals{j}; 
            C = T*v*T.'; % Legendre coefficients of sigma on patch j
            C(~intail) = 0; % zero out the low-degree coefficients
            e(j) = patchL2norm(T\C/T.', J); % L2 norm of the tail of sigma on patch j
            g = g + patchL2norm(v, J)^2;
        end
        eta{s} = e / sqrt(g);
    end

    % -- mark: of the N leaves at level <RMAX, the CEIL(FRAC*N) with the
    % largest ETA ----------------------------------------------------------
    alleta = vertcat(eta{:}); % error measure on all leaves
    dep = cellfun(@(p) p(:, 2), p2q(:), 'UniformOutput', false); % depths of all leaves
    eligible = vertcat(dep{:}) < rmax; % leaves at level <RMAX
    refinable = find(eligible); % leaves that can be refined
    [~, ord] = sort(alleta(refinable), 'descend'); % indices of refinable leaves sorted by decreasing ETA
    ismarked = false(size(alleta)); 
    ismarked(refinable(ord(1:ceil(frac*numel(refinable))))) = true; % boolean array of which leaves are marked for refinement 
    marked = cellfun(@find, mat2cell(ismarked, cellfun(@numel, eta)), ...
        'UniformOutput', false);

    % -- record ------------------------------------------------------------
    h = struct();
    h.level = lev;
    h.npatches = cellfun(@(d) length(d.x), dom);
    h.npts = ts.domain.nptspersurf;
    h.nmarked = nnz(ismarked);
    h.max_eta = max(alleta);
    h.l2_eta = norm(alleta);
    h.eta = eta;
    h.alpha = ts.alpha;
    h.time_s = tsolve;
    h.err_vs_B0 = nan;
    h.l2err_vs_B0 = nan;
    if ~isempty(refstate) || ~isequal(savebase, false)
        B = ts.surface_B();
    end
    if ~isempty(refstate)
        % RefTaylorState solves n.B = n.B0, so B should equal B0 exactly; the
        % residual is the discretization error.
        num = 0; den = 0; num2 = 0; den2 = 0;
        for s = 1:ns
            dB = ts.B0{s} - B{s};
            num = max(num, vecinfnorm(dB));
            den = max(den, vecinfnorm(ts.B0{s}));
            % norm of a surfacefunv is the pointwise magnitude; norm of that
            % is its L2 norm over the surface
            num2 = num2 + norm(norm(dB))^2;
            den2 = den2 + norm(norm(ts.B0{s}))^2;
        end
        h.err_vs_B0 = num/den;
        h.l2err_vs_B0 = sqrt(num2/den2);
    end
    if lev == 0
        history = h;
    else
        history(end+1) = h; %#ok<AGROW>
    end

    % -- write -------------------------------------------------------------
    for s = 1:ns
        tag = sprintf('_lev%d', lev);
        if ns > 1, tag = sprintf('%s_surf%d', tag, s); end
        if ~isempty(vtkbase)
            surfacemesh_to_vtk(dom{s}, sprintf('%s%s_sigma.vtk', vtkbase, tag), ...
                real(sigma{s}), 'Title', ...
                sprintf('sigma (real part), level %d', lev));
        end
        if ~isequal(savebase, false)
            % B, the mesh nodes and the patch depths for make_amr_vtk.py
            npatch = length(dom{s}.x);
            S = struct('n', n, 'npatch', npatch, 'level', lev, ...
                'depth', double(p2q{s}(:, 2)));
            S.xx = zeros(n, n, npatch); S.yy = S.xx; S.zz = S.xx;
            for i = 1:npatch
                S.xx(:,:,i) = dom{s}.x{i};
                S.yy(:,:,i) = dom{s}.y{i};
                S.zz(:,:,i) = dom{s}.z{i};
            end
            % surfacefun overloads subsref, so cat(3, f.vals{:}) would
            % silently return only the first patch. Index one patch at a time.
            Bv = zeros(n, n, npatch, 3);
            for c = 1:3
                comp = B{s}.components{c};
                for i = 1:npatch
                    Bv(:,:,i,c) = comp.vals{i};
                end
            end
            S.Bre = real(Bv);
            S.Bim = imag(Bv);
            save(sprintf('%s%s.mat', savebase, tag), '-struct', 'S', '-v7.3');
        end
    end

    if verbose
        es = '';
        if ~isempty(refstate)
            es = sprintf(', err vs B0 = %.3e (inf), %.3e (L2)', ...
                h.err_vs_B0, h.l2err_vs_B0);
        end
        fprintf('    max eta = %.3e, marked %d for next level%s\n', ...
            h.max_eta, h.nmarked, es);
    end

    if ~any(eligible)
        reason = 'depth';
        break
    end
end

if verbose
    switch reason
        case 'depth'
            fprintf('[sigma_amr] nothing left to refine; stopping.\n');
        case 'maxlevels'
            fprintf(['[sigma_amr] reached maxlevels = %d with %d patches ' ...
                'still eligible to refine; stopping.\n'], maxlevels, ...
                nnz(eligible));
    end
end

out = [];
out.dom = dom;
out.qf = qf;
out.p2q = p2q;
out.sigma = sigma;
out.alpha = ts.alpha;
out.ts = ts;
out.history = history;
out.reason = reason;
out.opts = struct('rmax', rmax, 'maxlevels', maxlevels, 'frac', frac);

end

% ==========================================================================

function tf = isposint(x)
tf = isnumeric(x) && isscalar(x) && isreal(x) && x >= 1 && x == round(x);
end

function v = getopt(opts, name, default)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
