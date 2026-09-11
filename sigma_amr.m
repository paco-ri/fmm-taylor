function out = sigma_amr(dom0, domparams, zk, flux, tol, rmax, varargin)
%SIGMA_AMR Iteratively refine a surface mesh until sigma is resolved.
%   OUT = SIGMA_AMR(DOM0, DOMPARAMS, ZK, FLUX, TOL, RMAX) solves the Taylor
%   state problem on a sequence of meshes, using the change in the density
%   sigma between successive discretizations as a per-patch error indicator.
%
%   The scheme is the usual solve-estimate-mark-refine loop. Starting from
%   DOM0, the Taylor state is solved, every patch is refined once, and the
%   problem is solved again. On each patch of the finer mesh the coarse
%   sigma (interpolated onto that patch) is compared against the fine
%   sigma. Where the two agree the coarse discretization already resolved
%   sigma and the patch is left alone; where they disagree the patch is
%   refined again and the comparison repeats one level down.
%
%   This is a different criterion from SURFACEMESH/ADAP_REF, which refines
%   until the *geometry* (the fundamental forms) is resolved and says
%   nothing about the solution.
%
%   RMAX is the maximum quadtree depth. It is required, and there is no
%   sensible default: it is the compute budget for the whole run, since a
%   mesh may grow by up to 4^RMAX patches, and it also sizes the quadforest
%   handed to TAYLORSTATE.INTACYC / INTBCYC.
%
%   OUT = SIGMA_AMR(..., OPTS) accepts a struct with fields:
%
%     maxlevels [RMAX]        maximum number of solve-estimate-mark-refine
%                             cycles. Distinct from RMAX: a patch left
%                             unmarked for several levels can still be
%                             eligible to refine when the budget runs out.
%     marking   ['threshold'] 'threshold' or 'dorfler'
%     eps_sigma [1e-6]        threshold on the relative per-patch indicator
%     theta     [0.5]         Dorfler bulk-marking fraction
%     p2q0      []            starting leaf set (cell, one per surface).
%                             Pass the P2Q returned by SURFACEMESH.ADAP_REF
%                             to start from a geometrically adapted mesh;
%                             leave empty to start from DOM0 uniformly.
%     refstate  []            struct with required fields xs_nodes and
%                             xs_weights (each a cell, one per surface) and
%                             B0fun, a handle @(dom, s) returning B0 on a
%                             mesh. B0fun is required because the mesh
%                             changes at every level, so B0 has to be
%                             re-evaluated on each one. An optional field
%                             B0 (a cell, one per surface) supplies B0 on
%                             the starting mesh, saving one evaluation;
%                             later levels always go through B0fun.
%                             When given, RefTaylorState is used and the
%                             true error against B0 is reported per level.
%     vtkbase   ''            if nonempty, write <base>_lev<k>_sigma.vtk
%     verbose   [true]        print progress
%
%   Arguments DOM0 (a surfacemesh or a cell of one or two), DOMPARAMS, ZK,
%   FLUX and TOL are as for TAYLORSTATE.
%
%   OUT is a struct with the final DOM, QF, P2Q, SIGMA and ALPHA, plus a
%   HISTORY struct array carrying, for each level, the patch counts, the
%   indicator statistics, the number of patches marked, and the timings.
%   OUT.REASON says why the loop stopped:
%
%     'resolved'  the indicator marked nothing, i.e. sigma is resolved
%     'depth'     everything still marked is already at depth RMAX
%     'maxlevels' the cycle budget ran out with patches still eligible
%
%   See also TAYLORSTATE, SURFACEMESH/REFINE_LEAVES,
%   SURFACEFUN/PROLONG_LEAVES.

% -- options ---------------------------------------------------------------
opts = [];
if ~isempty(varargin), opts = varargin{1}; end
if isstruct(opts) && isfield(opts, 'rmax')
    error('SIGMA_AMR:rmax', ...
        ['RMAX is now the sixth positional argument, not an option. ' ...
         'Call SIGMA_AMR(DOM0, DOMPARAMS, ZK, FLUX, TOL, RMAX, OPTS).']);
end
maxlevels = getopt(opts, 'maxlevels', []);
marking   = getopt(opts, 'marking', 'threshold');
eps_sigma = getopt(opts, 'eps_sigma', 1e-6);
theta     = getopt(opts, 'theta', 0.5);
p2q0      = getopt(opts, 'p2q0', []);
refstate  = getopt(opts, 'refstate', []);
vtkbase   = getopt(opts, 'vtkbase', '');
verbose   = getopt(opts, 'verbose', true);

if ~isposint(rmax)
    error('SIGMA_AMR:rmax', 'RMAX must be a positive integer scalar.');
end
if isempty(maxlevels)
    maxlevels = rmax;
elseif ~isposint(maxlevels)
    error('SIGMA_AMR:maxlevels', ...
        'opts.maxlevels must be a positive integer scalar.');
end

if ~ismember(marking, {'threshold','dorfler'})
    error('SIGMA_AMR:marking', ...
        'opts.marking must be ''threshold'' or ''dorfler''.');
end

if ~isempty(refstate) && ~(isfield(refstate, 'B0fun') && ...
        isa(refstate.B0fun, 'function_handle'))
    error('SIGMA_AMR:refstate', ...
        ['opts.refstate.B0fun is required and must be a function handle: ' ...
         'the mesh changes at every level, so B0 has to be re-evaluated ' ...
         'on each one.']);
end

% -- normalize the base mesh to a cell -------------------------------------
if isa(dom0, 'surfacemesh')
    base = {dom0};
else
    base = dom0;
end
ns = numel(base);

% -- starting leaf sets ----------------------------------------------------
p2q = cell(1, ns);
for s = 1:ns
    if isempty(p2q0)
        np = length(base{s}.x);
        p2q{s} = [(1:np).', zeros(np, 2)];
    else
        p2q{s} = p2q0{s};
        % A starting mesh deeper than RMAX would silently desynchronize the
        % quadforest's L_max from the leaves REFINE_LEAVES materializes.
        if max(p2q{s}(:, 2)) > rmax
            error('SIGMA_AMR:rmax', ...
                ['opts.p2q0 for surface %d reaches depth %d, deeper than ' ...
                 'RMAX = %d.'], s, max(p2q{s}(:, 2)), rmax);
        end
    end
end

% -- level 0 ----------------------------------------------------------------
domk = cell(1, ns); qfk = cell(1, ns); p2qk = cell(1, ns);
% domk, qfk, p2qk will eventually hold the output surfacemesh, quadforest, p2q, resp.
for s = 1:ns
    [domk{s}, qfk{s}, p2qk{s}] = surfacemesh.refine_leaves( ...
        base{s}, p2q{s}, [], rmax); 
    % no refining done; just returns surfacemesh and quadfores
end
if verbose
    fprintf('[sigma_amr] level 0: %s patches\n', patchstr(domk));
end
% A supplied refstate.B0 describes the starting mesh only, so it is used for
% this solve and nowhere else. Note that with opts.p2q0 the starting mesh is
% the adapted one, not DOM0; the patch counts below catch a B0 built on the
% wrong mesh.
B00 = getopt(refstate, 'B0', []);
if ~isempty(B00)
    if numel(B00) ~= ns
        error('SIGMA_AMR:refstate', ...
            'opts.refstate.B0 has %d entries but there are %d surfaces.', ...
            numel(B00), ns);
    end
    for s = 1:ns
        npB = numel(B00{s}.components{1}.vals);
        if npB ~= length(domk{s}.x)
            error('SIGMA_AMR:refstate', ...
                ['opts.refstate.B0{%d} covers %d patches but the starting ' ...
                 'mesh has %d. Omit B0 to have it built by B0fun.'], ...
                s, npB, length(domk{s}.x));
        end
    end
end

% carry out a Taylor state solve
t0 = tic;
[tsk, sigmak] = solve_level(domk, qfk, p2qk, domparams, zk, flux, tol, ...
    refstate, B00);
tk = toc(t0);

% record information about level 0
history = struct('level', {}, 'npatches', {}, 'npts', {}, 'nmarked', {}, ...
    'max_eta', {}, 'l2_eta', {}, 'alpha', {}, 'time_s', {}, 'err_vs_B0', {});
history(1) = record(0, domk, tsk, [], [], tk, refstate, sigmak);

if ~isempty(vtkbase)
    writevtk(vtkbase, 0, domk, sigmak);
end

% First pass refines everything, so that every patch gets compared once.
marked = cell(1, ns);
for s = 1:ns
    marked{s} = (1:size(p2qk{s}, 1)).';
end

% -- refinement loop --------------------------------------------------------
reason = 'maxlevels'; % reason why loop stopped
for lev = 1:maxlevels
    % `eligible` prevents marking a leaf that is already as deep as we allow
    [marked, nmark] = eligible(marked, p2qk, rmax);
    if nmark == 0
        reason = 'depth';
        if verbose
            fprintf('[sigma_amr] nothing left to refine; stopping.\n');
        end
        break
    end

    domf = cell(1, ns); qff = cell(1, ns); p2qf = cell(1, ns);
    for s = 1:ns
        % refine marked leaves and balance the quadtree
        [domf{s}, qff{s}, p2qf{s}] = surfacemesh.refine_leaves( ...
            base{s}, p2qk{s}, marked{s}, rmax);
    end
    if verbose
        fprintf('[sigma_amr] level %d: %s patches (marked %d)\n', ...
            lev, patchstr(domf), nmark);
    end

    % carry out Taylor state solve
    % this mesh is new, so we evaluate a new B0 if `refstate` provided
    t0 = tic;
    [tsf, sigmaf] = solve_level(domf, qff, p2qf, domparams, zk, flux, ...
        tol, refstate, []);
    tf = toc(t0);

    % -- estimate ----------------------------------------------------------
    eta = cell(1, ns);
    for s = 1:ns
        % interpolate coarse sigma on the refined mesh
        Psig = prolong_leaves(sigmak{s}, p2qk{s}, domf{s}, p2qf{s});
        npf = length(domf{s}.x); % number of patches on surface mesh s
        e = zeros(npf, 1);
        g = 0;
        for j = 1:npf
            J = domf{s}.J{j};
            % compute the L2 norm of the difference between the interpolated and refined sigma
            e(j) = patchL2norm(Psig.vals{j} - sigmaf{s}.vals{j}, J);
            g = g + patchL2norm(sigmaf{s}.vals{j}, J)^2;
        end
        eta{s} = e / sqrt(g);
    end

    % -- mark --------------------------------------------------------------
    newmarked = cell(1, ns);
    for s = 1:ns
        switch marking
            case 'threshold'
                idx = find(eta{s} > eps_sigma);
            case 'dorfler'
                [se, ord] = sort(eta{s}, 'descend');
                c = cumsum(se.^2);
                if c(end) == 0
                    idx = [];
                else
                    k = find(c >= theta*c(end), 1, 'first');
                    % k is the minumum integer such that the most 
                    %     k error-producing patches account for 
                    %     `theta` fraction of the total error
                    idx = ord(1:k);
                end
        end
        newmarked{s} = idx(:);
    end

    history(end+1) = record(lev, domf, tsf, eta, newmarked, tf, ...
        refstate, sigmaf); %#ok<AGROW>
    if ~isempty(vtkbase)
        writevtk(vtkbase, lev, domf, sigmaf);
    end
    if verbose
        if isnan(history(end).err_vs_B0)
            es = '';
        else
            es = sprintf(', err vs B0 = %.3e', history(end).err_vs_B0);
        end
        fprintf('    max eta = %.3e, marked %d for next level%s\n', ...
            history(end).max_eta, sum(cellfun(@numel, newmarked)), es);
    end

    domk = domf; qfk = qff; p2qk = p2qf; sigmak = sigmaf; tsk = tsf;
    marked = newmarked;

    if sum(cellfun(@numel, marked)) == 0
        reason = 'resolved';
        if verbose
            fprintf('[sigma_amr] sigma resolved everywhere; stopping.\n');
        end
        break
    end
end

% print correct output message if we stopped because of maxlevels but 
%     there are still eligible patches
if strcmp(reason, 'maxlevels')
    [~, nleft] = eligible(marked, p2qk, rmax);
    if nleft == 0
        reason = 'depth';
        if verbose
            fprintf('[sigma_amr] nothing left to refine; stopping.\n');
        end
    elseif verbose
        fprintf(['[sigma_amr] reached maxlevels = %d with %d patches ' ...
            'still eligible to refine; stopping.\n'], maxlevels, nleft);
    end
end

% fill in output struct
out = [];
out.dom = domk;
out.qf = qfk;
out.p2q = p2qk;
out.sigma = sigmak;
out.alpha = tsk.alpha;
out.ts = tsk;
out.history = history;
out.reason = reason;
out.opts = struct('rmax', rmax, 'maxlevels', maxlevels, 'marking', marking, ...
    'eps_sigma', eps_sigma, 'theta', theta);

end

% ==========================================================================

function tf = isposint(x)
tf = isnumeric(x) && isscalar(x) && isreal(x) && x >= 1 && x == round(x);
end

function [marked, n] = eligible(marked, p2q, rmax)
%ELIGIBLE Drop marked leaves that are already at the maximum depth.
n = 0;
for s = 1:numel(marked)
    marked{s} = marked{s}(p2q{s}(marked{s}, 2) < rmax);
    n = n + numel(marked{s});
end
end

function v = getopt(opts, name, default)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end

function [ts, sigma] = solve_level(dom, qf, p2q, domparams, zk, flux, ...
    tol, refstate, B0)
%SOLVE_LEVEL Build a Domain carrying the quadforest and solve on it.
%   B0 may be empty, in which case it is evaluated on DOM through B0FUN.
D = Domain(dom, domparams, qf, p2q);
if isempty(refstate)
    ts = TaylorState(D, domparams, zk, flux, tol);
else
    if isempty(B0)
        ns = numel(dom);
        B0 = cell(1, ns);
        for s = 1:ns
            B0{s} = refstate.B0fun(dom{s}, s);
        end
    end
    ts = RefTaylorState(D, domparams, zk, flux, B0, ...
        refstate.xs_nodes, refstate.xs_weights, tol);
end
ts = ts.solve(false);
sigma = ts.sigma;
end

function h = record(lev, dom, ts, eta, marked, t, refstate, sigma)
h.level = lev;
h.npatches = cellfun(@(d) length(d.x), dom);
h.npts = ts.domain.nptspersurf;
if isempty(eta)
    h.nmarked = nan;
    h.max_eta = nan;
    h.l2_eta = nan;
else
    h.nmarked = sum(cellfun(@numel, marked));
    h.max_eta = max(cellfun(@max, eta));
    h.l2_eta = sqrt(sum(cellfun(@(e) sum(e.^2), eta)));
end
h.alpha = ts.alpha;
h.time_s = t;
h.err_vs_B0 = nan;
if ~isempty(refstate)
    % RefTaylorState solves n.B = n.B0, so B should equal B0 exactly; the
    % residual is the discretization error.
    B = ts.surface_B();
    num = 0; den = 0;
    for s = 1:numel(B)
        num = max(num, vecinfnorm(ts.B0{s} - B{s}));
        den = max(den, vecinfnorm(ts.B0{s}));
    end
    h.err_vs_B0 = num/den;
end
end

function s = patchstr(dom)
s = strjoin(arrayfun(@(k) sprintf('%d', length(dom{k}.x)), ...
    1:numel(dom), 'UniformOutput', false), '+');
end

function writevtk(base, lev, dom, sigma)
for s = 1:numel(dom)
    if isscalar(dom)
        fname = sprintf('%s_lev%d_sigma.vtk', base, lev);
    else
        fname = sprintf('%s_lev%d_surf%d_sigma.vtk', base, lev, s);
    end
    surfacemesh_to_vtk(dom{s}, fname, real(sigma{s}), ...
        'Title', sprintf('sigma (real part), level %d', lev));
end
end
