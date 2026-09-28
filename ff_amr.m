function out = ff_amr(dom0, domparams, zk, flux, tol, rmax, varargin)
%FF_AMR Iteratively refine a surface mesh until the geometry is resolved.
%   OUT = FF_AMR(DOM0, DOMPARAMS, ZK, FLUX, TOL, RMAX) solves the Taylor
%   state problem on a sequence of meshes, using the fundamental-form error
%   of SURFACEMESH/FF_INDICATOR as the per-patch refinement indicator: a
%   patch is refined when interpolating its first fundamental form onto its
%   four children differs from recomputing the form there.
%
%   Unlike SIGMA_AMR the indicator needs no previous-level solve, so it is
%   available on the starting mesh and every level is marked from its own
%   geometry. The solve at each level is for reporting only.
%
%   RMAX is the maximum quadtree depth. It is required, and there is no
%   sensible default: it is the compute budget for the whole run, since a
%   mesh may grow by up to 4^RMAX patches, and it also sizes the quadforest
%   handed to TAYLORSTATE.INTACYC / INTBCYC.
%
%   OUT = FF_AMR(..., OPTS) accepts a struct with fields:
%
%     maxlevels [RMAX]        maximum number of solve-estimate-mark-refine
%                             cycles. Distinct from RMAX: a patch left
%                             unmarked for several levels can still be
%                             eligible to refine when the budget runs out.
%     marking   ['threshold'] 'threshold' or 'dorfler'
%     amr_tol   [1e-3]        threshold on the per-patch indicator. Unlike
%                             SIGMA_AMR's EPS_SIGMA this is an absolute
%                             quantity, in the level-0 patch's parameter
%                             units, so probe FF_INDICATOR on DOM0 before
%                             choosing it.
%     theta     [0.5]         Dorfler bulk-marking fraction
%     mode      [1]           1 or 2, selecting the fundamental form
%     p2q0      []            starting leaf set (cell, one per surface)
%     refstate  []            as for SIGMA_AMR: a struct with xs_nodes,
%                             xs_weights and B0fun. When given,
%                             RefTaylorState is used and the true error
%                             against B0 is reported per level.
%     vtkbase   ''            if nonempty, write <base>_lev<k>_sigma.vtk
%     savebase  [false]       false, or a path prefix. When a prefix is
%                             given, B is evaluated on the surface at every
%                             level and saved, with the mesh nodes and the
%                             patch depths, to <base>_lev<k>.mat (one
%                             surface) or <base>_lev<k>_surf<s>.mat (two),
%                             in the layout examples/make_amr_vtk.py reads.
%                             Costs one extra surface_B per level.
%     verbose   [true]        print progress
%
%   OUT is a struct with the final DOM, QF, P2Q, SIGMA and ALPHA, plus a
%   HISTORY struct array carrying, for each level, the patch counts, the
%   indicator statistics, the number of patches marked, and the timings.
%   MAX_ETA and L2_ETA hold the fundamental-form error, not a sigma
%   difference. OUT.REASON says why the loop stopped:
%
%     'resolved'  the indicator marked nothing
%     'depth'     everything still marked is already at depth RMAX
%     'maxlevels' the cycle budget ran out with patches still eligible
%
%   See also SIGMA_AMR, SURFACEMESH/FF_INDICATOR, SURFACEMESH/REFINE_LEAVES.

% -- options ---------------------------------------------------------------
opts = [];
if ~isempty(varargin), opts = varargin{1}; end
if isstruct(opts) && isfield(opts, 'eps_sigma')
    error('FF_AMR:amr_tol', ...
        ['FF_AMR thresholds an absolute fundamental-form error, not a ' ...
         'relative sigma indicator. Use opts.amr_tol, not opts.eps_sigma.']);
end
maxlevels = getopt(opts, 'maxlevels', []);
marking   = getopt(opts, 'marking', 'threshold');
amr_tol   = getopt(opts, 'amr_tol', 1e-3);
theta     = getopt(opts, 'theta', 0.5);
mode      = getopt(opts, 'mode', 1);
p2q0      = getopt(opts, 'p2q0', []);
refstate  = getopt(opts, 'refstate', []);
vtkbase   = getopt(opts, 'vtkbase', '');
savebase  = getopt(opts, 'savebase', false);
verbose   = getopt(opts, 'verbose', true);

if ~isposint(rmax)
    error('FF_AMR:rmax', 'RMAX must be a positive integer scalar.');
end
if isempty(maxlevels)
    maxlevels = rmax;
elseif ~isposint(maxlevels)
    error('FF_AMR:maxlevels', ...
        'opts.maxlevels must be a positive integer scalar.');
end

if ~ismember(marking, {'threshold','dorfler'})
    error('FF_AMR:marking', ...
        'opts.marking must be ''threshold'' or ''dorfler''.');
end

if ~ismember(mode, [1 2])
    error('FF_AMR:mode', 'opts.mode must be 1 or 2.');
end

if ~(isequal(savebase, false) || ...
        ((ischar(savebase) || isstring(savebase)) && strlength(savebase) > 0))
    error('FF_AMR:savebase', ...
        'opts.savebase must be false or a nonempty path prefix.');
end

if ~isempty(refstate) && ~(isfield(refstate, 'B0fun') && ...
        isa(refstate.B0fun, 'function_handle'))
    error('FF_AMR:refstate', ...
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
        if max(p2q{s}(:, 2)) > rmax
            error('FF_AMR:rmax', ...
                ['opts.p2q0 for surface %d reaches depth %d, deeper than ' ...
                 'RMAX = %d.'], s, max(p2q{s}(:, 2)), rmax);
        end
    end
end

% -- level 0 ----------------------------------------------------------------
domk = cell(1, ns); qfk = cell(1, ns); p2qk = cell(1, ns);
for s = 1:ns
    [domk{s}, qfk{s}, p2qk{s}] = surfacemesh.refine_leaves( ...
        base{s}, p2q{s}, [], rmax);
end
if verbose
    fprintf('[ff_amr] level 0: %s patches\n', patchstr(domk));
end
B00 = getopt(refstate, 'B0', []);
if ~isempty(B00)
    if numel(B00) ~= ns
        error('FF_AMR:refstate', ...
            'opts.refstate.B0 has %d entries but there are %d surfaces.', ...
            numel(B00), ns);
    end
    for s = 1:ns
        npB = numel(B00{s}.components{1}.vals);
        if npB ~= length(domk{s}.x)
            error('FF_AMR:refstate', ...
                ['opts.refstate.B0{%d} covers %d patches but the starting ' ...
                 'mesh has %d. Omit B0 to have it built by B0fun.'], ...
                s, npB, length(domk{s}.x));
        end
    end
end

t0 = tic;
[tsk, sigmak] = solve_level(domk, qfk, p2qk, domparams, zk, flux, tol, ...
    refstate, B00);
tk = toc(t0);

% The indicator is geometric, so it is available on this mesh already.
etak = get_eta(domk, p2qk, mode);
marked = mark(etak, marking, amr_tol, theta);

history = struct('level', {}, 'npatches', {}, 'npts', {}, 'nmarked', {}, ...
    'max_eta', {}, 'l2_eta', {}, 'alpha', {}, 'time_s', {}, 'err_vs_B0', {});
history(1) = record(0, domk, tsk, etak, marked, tk, refstate);

if ~isempty(vtkbase)
    writevtk(vtkbase, 0, domk, sigmak);
end
if ~isequal(savebase, false)
    savelevel(savebase, 0, domk, tsk, p2qk);
end
if verbose
    fprintf('    max eta = %.3e, marked %d for next level%s\n', ...
        history(end).max_eta, history(end).nmarked, errstr(history(end)));
end

% -- refinement loop --------------------------------------------------------
reason = 'maxlevels';
if sum(cellfun(@numel, marked)) == 0
    reason = 'resolved';
    if verbose
        fprintf('[ff_amr] geometry resolved everywhere; stopping.\n');
    end
end

if ~strcmp(reason, 'resolved')
for lev = 1:maxlevels
    [marked, nmark] = eligible(marked, p2qk, rmax);
    if nmark == 0
        reason = 'depth';
        if verbose
            fprintf('[ff_amr] nothing left to refine; stopping.\n');
        end
        break
    end

    domf = cell(1, ns); qff = cell(1, ns); p2qf = cell(1, ns);
    for s = 1:ns
        [domf{s}, qff{s}, p2qf{s}] = surfacemesh.refine_leaves( ...
            base{s}, p2qk{s}, marked{s}, rmax);
    end
    if verbose
        fprintf('[ff_amr] level %d: %s patches (marked %d)\n', ...
            lev, patchstr(domf), nmark);
    end

    t0 = tic;
    [tsf, sigmaf] = solve_level(domf, qff, p2qf, domparams, zk, flux, ...
        tol, refstate, []);
    tf = toc(t0);

    etaf = get_eta(domf, p2qf, mode);
    newmarked = mark(etaf, marking, amr_tol, theta);

    history(end+1) = record(lev, domf, tsf, etaf, newmarked, tf, ...
        refstate); %#ok<AGROW>
    if ~isempty(vtkbase)
        writevtk(vtkbase, lev, domf, sigmaf);
    end
    if ~isequal(savebase, false)
        savelevel(savebase, lev, domf, tsf, p2qf);
    end
    if verbose
        fprintf('    max eta = %.3e, marked %d for next level%s\n', ...
            history(end).max_eta, history(end).nmarked, errstr(history(end)));
    end

    domk = domf; qfk = qff; p2qk = p2qf; sigmak = sigmaf; tsk = tsf;
    marked = newmarked;

    if sum(cellfun(@numel, marked)) == 0
        reason = 'resolved';
        if verbose
            fprintf('[ff_amr] geometry resolved everywhere; stopping.\n');
        end
        break
    end
end
end

if strcmp(reason, 'maxlevels')
    [~, nleft] = eligible(marked, p2qk, rmax);
    if nleft == 0
        reason = 'depth';
        if verbose
            fprintf('[ff_amr] nothing left to refine; stopping.\n');
        end
    elseif verbose
        fprintf(['[ff_amr] reached maxlevels = %d with %d patches ' ...
            'still eligible to refine; stopping.\n'], maxlevels, nleft);
    end
end

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
    'amr_tol', amr_tol, 'theta', theta, 'mode', mode);

end

% ==========================================================================

function tf = isposint(x)
tf = isnumeric(x) && isscalar(x) && isreal(x) && x >= 1 && x == round(x);
end

function eta = get_eta(dom, p2q, mode)
eta = cell(1, numel(dom));
for s = 1:numel(dom)
    eta{s} = surfacemesh.ff_indicator(dom{s}, p2q{s}, mode);
end
end

function marked = mark(eta, marking, amr_tol, theta)
marked = cell(1, numel(eta));
for s = 1:numel(eta)
    switch marking
        case 'threshold'
            idx = find(eta{s} > amr_tol);
        case 'dorfler'
            [se, ord] = sort(eta{s}, 'descend');
            c = cumsum(se.^2);
            if c(end) == 0
                idx = [];
            else
                k = find(c >= theta*c(end), 1, 'first');
                idx = ord(1:k);
            end
    end
    marked{s} = idx(:);
end
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

function h = record(lev, dom, ts, eta, marked, t, refstate)
h.level = lev;
h.npatches = cellfun(@(d) length(d.x), dom);
h.npts = ts.domain.nptspersurf;
h.nmarked = sum(cellfun(@numel, marked));
h.max_eta = max(cellfun(@max, eta));
h.l2_eta = sqrt(sum(cellfun(@(e) sum(e.^2), eta)));
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

function s = errstr(h)
if isnan(h.err_vs_B0)
    s = '';
else
    s = sprintf(', err vs B0 = %.3e', h.err_vs_B0);
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

function savelevel(base, lev, dom, ts, p2q)
%SAVELEVEL Save B, the mesh nodes and the patch depths for make_amr_vtk.py.
B = ts.surface_B();
for s = 1:numel(dom)
    if isscalar(dom)
        fname = sprintf('%s_lev%d.mat', base, lev);
    else
        fname = sprintf('%s_lev%d_surf%d.mat', base, lev, s);
    end
    d = dom{s};
    n = size(d.x{1}, 1);
    npatch = length(d.x);
    xx = zeros(n, n, npatch); yy = xx; zz = xx;
    for i = 1:npatch
        xx(:,:,i) = d.x{i}; yy(:,:,i) = d.y{i}; zz(:,:,i) = d.z{i};
    end
    % surfacefun overloads subsref, so cat(3, f.vals{:}) would silently
    % return only the first patch. Index one patch at a time.
    Bv = zeros(n, n, npatch, 3);
    for c = 1:3
        comp = B{s}.components{c};
        for i = 1:npatch
            Bv(:,:,i,c) = comp.vals{i};
        end
    end
    Bre = real(Bv);
    Bim = imag(Bv);
    depth = double(p2q{s}(:, 2));
    level = lev;
    save(fname, 'xx', 'yy', 'zz', 'Bre', 'Bim', 'depth', 'n', 'npatch', ...
        'level', '-v7.3');
end
end
