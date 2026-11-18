function diag_xs2(n, nvs, outmat)
%DIAG_XS2 Test whether the surface-2 flux cross-section causes the erratic
%   inner-surface error of the two-surface sweep (convoop2_sweep.m).
%
%   prepare_torus builds cross-section 2 (the z-flux, fluxnodes{2}) from the
%   v = 0 curves of both tori but places it at z = 0, while those curves lie at
%   z = -0.16*s*sin(u). Its inner edge therefore sits up to 0.017 inside the
%   inner torus (outside the domain) and its outer edge stops up to 0.028 short
%   of the outer wall. The solver evaluates its own field at these nodes to
%   impose the flux constraint, so the constraint is not exactly met by B0.
%
%   For each nv this solves twice, with the convoop2_sweep.m setup (zk = 0.5,
%   nr = 24, tol = 1e-10):
%     'orig'   cross-section 2 as built by prepare_torus
%     'fixed'  same plane z = 0, radial segments from the inner torus to the
%              outer torus at their exact z = 0 crossings
%   and reports per-surface errors, where the inner-surface error peaks, and
%   the mismatch between the solved field and B0 on cross-section 2.
%   'orig' should reproduce convdata_twosurface_highp.csv.
%
%   About 1 h per solve at n = 11 on jupiter. Results are saved after every
%   solve.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(n),      n = 11;      end
if nargin < 2 || isempty(nvs),    nvs = [6 7]; end
if nargin < 3 || isempty(outmat), outmat = fullfile(here, 'diag_xs2.mat'); end

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

% Torus coefficients, copied from prepare_torus.m (keep in sync).
d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.20 0];
d(3,2) = 2;
d(4,3) = 0;
[I, J] = ndgrid(-1:2, -1:2);

% Exact z = 0 crossing near v = 0 of each torus: R(u) for the fixed section.
up = 2*pi*(0:np-1)/np;
Rcross = zeros(2, np);            % row 1 outer (ao), row 2 inner (ai)
sc = [ao ai];
for is = 1:2
    c = sc(is)*d;
    c(3,2) = 2;
    for jp = 1:np
        v = 0;
        for it = 1:50
            a = (1 - I)*v + J*up(jp);
            v = v - sum(c.*sin(a), 'all')/sum(c.*(1 - I).*cos(a), 'all');
        end
        a = (1 - I)*v + J*up(jp);
        assert(abs(sum(c.*sin(a), 'all')) < 1e-14);
        Rcross(is, jp) = sum(c.*cos(a), 'all');
    end
end
[rnodes, rwts] = chebpts(nr, [0 1], 1);
xs2fix = zeros(3, nr*np);
w2fix  = zeros(1, nr*np);
for i = 1:nr
    for jp = 1:np
        ij = (i-1)*np + jp;
        R = (1 - rnodes(i))*Rcross(2,jp) + rnodes(i)*Rcross(1,jp);
        xs2fix(:,ij) = [R*cos(up(jp)); R*sin(up(jp)); 0];
        w2fix(ij) = (2*pi/np)*rwts(i)*R*(Rcross(1,jp) - Rcross(2,jp));
    end
end

res = struct([]);
fprintf('=== diag_xs2 start %s  n = %d  -> %s\n', datestr(now), n, outmat);
fprintf('  maxNumCompThreads = %d\n', maxNumCompThreads);

for nv = nvs
    nu = 3*nv;
    for mode = {'orig', 'fixed'}
        mode = mode{1}; %#ok<FXSET>
        fprintf('\n\t========\n\tn = %d, nv = %d, %s\n', n, nv, mode);

        [dom, qnodes, qweights] = prepare_torus(n, nu, nv, n, nu, nv, ...
            ao, ai, nr, nt, np);
        if strcmp(mode, 'fixed')
            qnodes{2} = xs2fix;
            qweights{2} = w2fix;
        end
        B0 = {reftaylorsurffun(dom{1}, n, ntheta, rmin, rmaj, jmag, zk), ...
              reftaylorsurffun(dom{2}, n, ntheta, rmin, rmaj, jmag, zk)};

        flux = zeros(1,2);
        for i = 1:nr*nt
            B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{1}(:,i));
            flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
        end
        B0xs2 = zeros(3, nr*np);
        for i = 1:nr*np
            B0xs2(:,i) = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{2}(:,i));
            flux(2) = flux(2) - B0xs2(3,i)*qweights{2}(i);
        end

        ts = RefTaylorState(dom, [n, nu, nv], zk, flux, B0, qnodes, ...
            qweights, tol);
        t0 = tic;
        ts = ts.solve(true);
        t_solve = toc(t0);
        B = ts.surface_B();

        % per-node error, max over components, both surfaces
        B0max = max(vecinfnorm(B0{1}), vecinfnorm(B0{2}));
        errnode = cell(1,2);
        for s = 1:2
            dB = B0{s} - B{s};
            npatch = numel(dom{s}.x);
            e = zeros(n, n, npatch);
            for k = 1:npatch
                e(:,:,k) = max(cat(3, abs(dB.components{1}.vals{k}), ...
                    abs(dB.components{2}.vals{k}), ...
                    abs(dB.components{3}.vals{k})), [], 3)/B0max;
            end
            errnode{s} = e;
        end
        err_s = [max(errnode{1}, [], 'all'), max(errnode{2}, [], 'all')];

        % where the inner-surface error peaks (no Morton reordering for these nu, nv)
        [~, imax] = max(errnode{2}(:));
        [il, kp] = ind2sub([n*n, nu*nv], imax);
        ku = floor((kp-1)/nv) + 1;
        kv = mod(kp-1, nv) + 1;
        ub = linspace(0, 2*pi, nu+1);
        vb = linspace(0, 2*pi, nv+1);
        [uu, vv] = chebpts2(n, n, [ub(ku:ku+1) vb(kv:kv+1)]);
        xyz = [dom{2}.x{kp}(il); dom{2}.y{kp}(il); dom{2}.z{kp}(il)];

        % solved field vs B0 on cross-section 2
        Bxs2 = ts.interior_B(qnodes{2});
        disc = reshape(abs(Bxs2(3,:) - B0xs2(3,:)), np, nr);
        % The solver enforces sum((Bxs2 - B0xs2)_z * w) = 0 exactly, so the
        % informative number is the part of it from nodes inside the inner
        % torus (outside the domain), which the in-domain nodes must cancel.
        inside = hypot(qnodes{2}(1,:), qnodes{2}(2,:)) < repmat(Rcross(2,:), 1, nr);
        fluxmis = abs(sum((Bxs2(3,inside) - B0xs2(3,inside)).*qweights{2}(inside))) ...
            /abs(flux(2));

        fprintf(['n=%d nv=%d %-5s err=%.6e  outer %.6e  inner %.6e  ' ...
                 'alpha = [%s]  solve %.0f s\n'], n, nv, mode, max(err_s), ...
            err_s(1), err_s(2), num2str(ts.alpha.', '%.10e '), t_solve);
        fprintf('  inner max at u = %.4f, v = %.4f (patch %d,%d), xyz = [%.4f %.4f %.4f]\n', ...
            uu(il), vv(il), ku, kv, xyz);
        fprintf('  XS2 |Bz - B0z|, max over u, radial index 1..5: %s\n', ...
            num2str(max(disc(:,1:5), [], 1), ' %.2e'));
        fprintf('  XS2 |Bz - B0z|, max over u, radial index %d..%d: %s\n', ...
            nr-4, nr, num2str(max(disc(:,nr-4:nr), [], 1), ' %.2e'));
        fprintf('  XS2 nodes inside the inner torus: %d; their flux mismatch / flux: %.3e\n', ...
            nnz(inside), fluxmis);

        res(end+1).n = n; %#ok<AGROW>
        res(end).nv = nv;
        res(end).mode = mode;
        res(end).err_s = err_s;
        res(end).alpha = ts.alpha;
        res(end).flux = flux;
        res(end).fluxmis = fluxmis;
        res(end).ninside = nnz(inside);
        res(end).disc = disc;
        res(end).errnode = errnode;
        res(end).inner_max_uv = [uu(il) vv(il)];
        res(end).t_solve = t_solve;
        save(outmat, 'res', 'Rcross', 'xs2fix', 'w2fix');
    end
end

fprintf('\nDONE_SENTINEL %s\n', datestr(now));
end

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf), ...
         norm(f.components{2}, inf), ...
         norm(f.components{3}, inf)]);
end
