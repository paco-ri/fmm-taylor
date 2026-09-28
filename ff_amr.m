function out = ff_amr(dom0, domparams, zk, flux, tol, rmax, varargin)
%FF_AMR Solve the Taylor state on each mesh of a fundamental-form refinement.
%   OUT = FF_AMR(DOM0, DOMPARAMS, ZK, FLUX, TOL, RMAX) refines each surface
%   of DOM0 with SURFACEMESH/ADAP_REF, which marks patches by the
%   fundamental-form error of SURFACEMESH/FF_INDICATOR, and solves the
%   Taylor state problem on every intermediate mesh. The refinement uses
%   only the geometry; the solves are for reporting.
%
%   RMAX is the maximum quadtree depth and the number of refinement passes.
%   It is required, and there is no sensible default: it is the compute
%   budget for the whole run, since a mesh may grow by up to 4^RMAX
%   patches, and it also sizes the quadforest handed to
%   TAYLORSTATE.INTACYC / INTBCYC.
%
%   OUT = FF_AMR(..., OPTS) accepts a struct with fields:
%
%     marking   ['threshold'] 'threshold' or 'dorfler'
%     amr_tol   [1e-3]        threshold on the per-patch indicator. Unlike
%                             SIGMA_AMR's EPS_SIGMA this is an absolute
%                             quantity, in the level-0 patch's parameter
%                             units, so probe FF_INDICATOR on DOM0 before
%                             choosing it.
%     theta     [0.5]         Dorfler bulk-marking fraction
%     mode      [1]           1 or 2, selecting the fundamental form
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
%   Marking is per surface.
%
%   OUT is a struct with the final DOM, QF, P2Q, SIGMA and ALPHA, plus a
%   HISTORY struct array carrying, for each level, the patch counts, the
%   indicator statistics, the number of patches marked, and the timings.
%   MAX_ETA and L2_ETA hold the fundamental-form error, not a sigma
%   difference. OUT.REASON says why refinement stopped:
%
%     'resolved'  the indicator marked nothing
%     'depth'     everything still marked is already at depth RMAX
%     'maxlevels' RMAX passes ran with patches still eligible
%
%   See also SIGMA_AMR, SURFACEMESH/ADAP_REF, SURFACEMESH/FF_INDICATOR.

% -- options ---------------------------------------------------------------
opts = [];
if ~isempty(varargin), opts = varargin{1}; end
if isstruct(opts) && isfield(opts, 'eps_sigma')
    error('FF_AMR:amr_tol', ...
        ['FF_AMR thresholds an absolute fundamental-form error, not a ' ...
         'relative sigma indicator. Use opts.amr_tol, not opts.eps_sigma.']);
end
marking   = getopt(opts, 'marking', 'threshold');
amr_tol   = getopt(opts, 'amr_tol', 1e-3);
theta     = getopt(opts, 'theta', 0.5);
mode      = getopt(opts, 'mode', 1);
refstate  = getopt(opts, 'refstate', []);
vtkbase   = getopt(opts, 'vtkbase', '');
savebase  = getopt(opts, 'savebase', false);
verbose   = getopt(opts, 'verbose', true);

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

% -- refine -----------------------------------------------------------------
hist = cell(1, ns);
for s = 1:ns
    [~, ~, ~, hist{s}] = surfacemesh.adap_ref(base{s}, amr_tol, rmax, ...
        mode, marking=marking, theta=theta);
end
nlev = max(cellfun(@numel, hist));

% -- solve on every level ---------------------------------------------------
history = struct('level', {}, 'npatches', {}, 'npts', {}, 'nmarked', {}, ...
    'max_eta', {}, 'l2_eta', {}, 'alpha', {}, 'time_s', {}, 'err_vs_B0', {});
domk = cell(1, ns); qfk = cell(1, ns); p2qk = cell(1, ns);
eta = cell(1, ns); marked = cell(1, ns);
for lev = 0:nlev-1
    % A surface whose refinement stopped early keeps its last mesh.
    for s = 1:ns
        h = hist{s}(min(lev+1, end));
        [domk{s}, qfk{s}, p2qk{s}] = surfacemesh.refine_leaves( ...
            base{s}, h.p2q, [], rmax);
        eta{s} = h.eta;
        marked{s} = h.marked;
    end
    if verbose
        fprintf('[ff_amr] level %d: %s patches\n', lev, patchstr(domk));
    end

    t0 = tic;
    [tsk, sigmak] = solve_level(domk, qfk, p2qk, domparams, zk, flux, ...
        tol, refstate);
    tk = toc(t0);

    history(end+1) = record(lev, domk, tsk, eta, marked, tk, ...
        refstate); %#ok<AGROW>
    if ~isempty(vtkbase)
        writevtk(vtkbase, lev, domk, sigmak);
    end
    if ~isequal(savebase, false)
        savelevel(savebase, lev, domk, tsk, p2qk);
    end
    if verbose
        fprintf('    max eta = %.3e, marked %d for next level%s\n', ...
            history(end).max_eta, history(end).nmarked, errstr(history(end)));
    end
end

nleft = 0;
for s = 1:ns
    nleft = nleft + sum(p2qk{s}(marked{s}, 2) < rmax);
end
if history(end).nmarked == 0
    reason = 'resolved';
elseif nleft == 0
    reason = 'depth';
else
    reason = 'maxlevels';
end
if verbose
    fprintf('[ff_amr] stopped: %s\n', reason);
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
out.opts = struct('rmax', rmax, 'marking', marking, 'amr_tol', amr_tol, ...
    'theta', theta, 'mode', mode);

end

% ==========================================================================

function v = getopt(opts, name, default)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end

function [ts, sigma] = solve_level(dom, qf, p2q, domparams, zk, flux, ...
    tol, refstate)
%SOLVE_LEVEL Build a Domain carrying the quadforest and solve on it.
D = Domain(dom, domparams, qf, p2q);
if isempty(refstate)
    ts = TaylorState(D, domparams, zk, flux, tol);
else
    ns = numel(dom);
    B0 = cell(1, ns);
    for s = 1:ns
        B0{s} = refstate.B0fun(dom{s}, s);
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
