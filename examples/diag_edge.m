function diag_edge(n, nv, outmat)
%DIAG_EDGE Is the large near-edge error on the inner torus in the solution,
%   or only in the on-surface field?
%
%   diag_xs2.m found that at n = 11, nv = 6 the inner-surface error peaks on
%   patch v-edges (2.5e-5) and the next ring in (1.1e-5), while nodes deeper
%   in each patch sit at 3.3e-6 and converge at ~h^(p-1). The plotted error
%   comes from surface_B, which builds the on-surface field from surface
%   derivatives of the densities -- least accurate at patch edges.
%
%   This solves once (convoop2_sweep.m setup: zk = 0.5, nr = 24,
%   tol = 1e-10) and evaluates the field with interior_B at points stepped off
%   selected nodes into the shell, comparing with B0 there. A point at
%   parameter (u, v) on the torus of scale s is stepped to scale s + t (inner)
%   or s - t (outer), which stays inside the shell for 0 < t < ao - ai.
%   Selected nodes: the 10 worst inner nodes, the worst inner node at least
%   two nodes from every edge, and the same two kinds on the outer surface.
%
%   Reading the table: if the error just off a bad edge node drops to the
%   deep-node level, the solution is fine there and the problem is
%   surface_B; if it stays near the on-surface value, the densities
%   themselves are wrong near the edges.
%
%   About 1 h at n = 11, nv = 6 on jupiter.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(n),      n = 11;  end
if nargin < 2 || isempty(nv),     nv = 6;  end
if nargin < 3 || isempty(outmat), outmat = fullfile(here, 'diag_edge.mat'); end

tol    = 1e-10;
zk     = 0.5;
nr     = 24;
nt     = 40;
np     = 40;
ntheta = 200;
ao     = 1.0;
ai     = 0.6;
rmin   = 2.0;
rmaj   = 2.0;
jmag   = 1.0;
nu     = 3*nv;
ts_off = [1e-1 3e-2 1e-2 3e-3 1e-3];   % scale steps off the surface

% Torus coefficients, copied from prepare_torus.m (keep in sync).
d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.20 0];
d(3,2) = 2;
d(4,3) = 0;

% prepare_torus reorders patches when nu and nv are both powers of 2
assert(~(bitand(nv, nv-1) == 0 && bitand(nu, nu-1) == 0), ...
    'diag_edge: Morton-ordered patches are not handled');

fprintf('=== diag_edge start %s  n = %d, nv = %d  -> %s\n', datestr(now), n, nv, outmat);
fprintf('  maxNumCompThreads = %d\n', maxNumCompThreads);

[dom, qnodes, qweights] = prepare_torus(n, nu, nv, n, nu, nv, ao, ai, nr, nt, np);
B0 = {reftaylorsurffun(dom{1}, n, ntheta, rmin, rmaj, jmag, zk), ...
      reftaylorsurffun(dom{2}, n, ntheta, rmin, rmaj, jmag, zk)};
flux = zeros(1,2);
for i = 1:nr*nt
    B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{1}(:,i));
    flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
end
for i = 1:nr*np
    B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{2}(:,i));
    flux(2) = flux(2) - B0eval(3)*qweights{2}(i);
end

ts = RefTaylorState(dom, [n, nu, nv], zk, flux, B0, qnodes, qweights, tol);
t0 = tic;
ts = ts.solve(true);
t_solve = toc(t0);
B = ts.surface_B();

B0max = max(vecinfnorm(B0{1}), vecinfnorm(B0{2}));
deep = false(n);
deep(3:n-2, 3:n-2) = true;
ub = linspace(0, 2*pi, nu+1);
vb = linspace(0, 2*pi, nv+1);
sc = [ao ai];
sgn = [-1 1];                       % step direction in scale: outer in, inner out
names = {'outer', 'inner'};

sel = zeros(0, 4);                  % [surface, patch, local index, deep?]
errnode = cell(1,2);
for s = 1:2
    dB = B0{s} - B{s};
    e = zeros(n, n, nu*nv);
    for k = 1:nu*nv
        e(:,:,k) = max(cat(3, abs(dB.components{1}.vals{k}), ...
            abs(dB.components{2}.vals{k}), ...
            abs(dB.components{3}.vals{k})), [], 3)/B0max;
    end
    errnode{s} = e;
    [~, order] = sort(e(:), 'descend');
    nworst = 1 + 9*(s == 2);
    [il, kp] = ind2sub([n*n, nu*nv], order(1:nworst));
    sel = [sel; s*ones(nworst,1), kp, il, zeros(nworst,1)]; %#ok<AGROW>
    edeep = e .* repmat(deep, 1, 1, nu*nv);
    [~, idp] = max(edeep(:));
    [il, kp] = ind2sub([n*n, nu*nv], idp);
    sel = [sel; s, kp, il, 1]; %#ok<AGROW>
end

nsel = size(sel, 1);
uv = zeros(2, nsel);
pts = zeros(3, nsel*numel(ts_off));
for q = 1:nsel
    kp = sel(q,2);
    ku = floor((kp-1)/nv) + 1;
    kv = mod(kp-1, nv) + 1;
    [uu, vv] = chebpts2(n, n, [ub(ku:ku+1) vb(kv:kv+1)]);
    uv(:,q) = [uu(sel(q,3)); vv(sel(q,3))];
    for it = 1:numel(ts_off)
        c = (sc(sel(q,1)) + sgn(sel(q,1))*ts_off(it))*d;
        c(3,2) = 2;
        [I, J] = ndgrid(-1:2, -1:2);
        a = (1 - I)*uv(2,q) + J*uv(1,q);
        pts(:, (q-1)*numel(ts_off) + it) = ...
            [cos(uv(1,q))*sum(c.*cos(a), 'all'); ...
             sin(uv(1,q))*sum(c.*cos(a), 'all'); ...
             sum(c.*sin(a), 'all')];
    end
end

Bpts = ts.interior_B(pts);
B0pts = zeros(3, size(pts,2));
for i = 1:size(pts,2)
    B0pts(:,i) = reftaylor(ntheta, rmin, rmaj, jmag, zk, pts(:,i));
end
eoff = reshape(max(abs(Bpts - B0pts), [], 1)/B0max, numel(ts_off), nsel);

fprintf('\nsolve %.0f s.  err on surface: outer %.3e, inner %.3e\n', t_solve, ...
    max(errnode{1}, [], 'all'), max(errnode{2}, [], 'all'));
fprintf('\nsurface  kind   patch  local  u       v       on-surf  | off by scale step t = %s\n', ...
    num2str(ts_off, ' %-9.0e'));
for q = 1:nsel
    kp = sel(q,2); il = sel(q,3);
    [ir, ic] = ind2sub([n n], il);
    kind = 'worst';
    if sel(q,4), kind = 'deep'; end
    edge = '';
    if ir == 1 || ir == n, edge = [edge 'v']; end
    if ic == 1 || ic == n, edge = [edge 'u']; end
    fprintf('%-6s   %-5s  %4d  %2d,%2d%-2s %.4f  %.4f  %.3e |%s\n', ...
        names{sel(q,1)}, kind, kp, ir, ic, edge, uv(1,q), uv(2,q), ...
        errnode{sel(q,1)}(il + (kp-1)*n*n), num2str(eoff(:,q).', ' %.3e'));
end

save(outmat, 'sel', 'uv', 'pts', 'ts_off', 'eoff', 'errnode', 'Bpts', ...
    'B0pts', 't_solve');
fprintf('\nDONE_SENTINEL %s\n', datestr(now));
end

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf), ...
         norm(f.components{2}, inf), ...
         norm(f.components{3}, inf)]);
end
