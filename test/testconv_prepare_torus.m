%% Known-good baseline: RefTaylorState convergence on prepare_torus
%
% Runs the *same* RefTaylorState pipeline used by conv_square_torus.m, but
% on the smooth prepare_torus geometry (as in examples/convoop.m), at patch
% counts matched to my_square_torus's m=6 and m=8 (nu*nv = 432 and 768),
% for multiple polynomial orders n.
%
% This gives a direct numeric baseline for "what does the surface-B relative
% error trend look like for this code path when the geometry/reference
% field are well-matched." conv_square_torus.m currently shows the error
% *increasing* with resolution (0.604 -> 0.669 -> 0.722 at m=6,8,10, before
% the orientation fix); if this test shows a clearly *decreasing* trend at
% matched patch counts, that confirms the my_square_torus pipeline has a
% real bug beyond the orientation/sign issue already fixed, and the next
% step should compare intermediate quantities (dfunc, wfunc, sigma, m0)
% between the two geometries directly.

clear
addpath('../examples')

tol = 1e-4;
zk = 0;
ns = [5 7]; % polynomial order + 1

ntheta = 1e3;
jmag = 1.0;
rmin = 2.0;
rmaj = 2.0;

nr = 12;
nt = 100;

nvs = [12 16]; % nu = 3*nv -> nu*nv = 432, 768 (matches m=6,8 for my_square_torus)

lerr    = zeros(numel(ns), numel(nvs));
lsigma  = zeros(numel(ns), numel(nvs));
npatches = zeros(numel(ns), numel(nvs));
dof     = zeros(numel(ns), numel(nvs)); % n*sqrt(nu*nv), matching conv_square_torus.m x-axis

for i_n = 1:numel(ns)
    n = ns(i_n);
    for i_nv = 1:numel(nvs)
        nv = nvs(i_nv);
        nu = 3*nv;
        domparams = [n, nu, nv];

        [dom, qnodes, qweights] = prepare_torus(n, nu, nv, nr, nt);
        dom = dom{1};
        qnodes = qnodes{1};
        qweights = qweights{1};

        B0 = reftaylorsurffun(dom, n, ntheta, rmin, rmaj, jmag, zk);

        flux = 0;
        for j = 1:nr*nt
            B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes(:,j));
            flux = flux + B0eval(2)*qweights(j);
        end

        B0c = {B0};
        qnodesc = {qnodes};
        qweightsc = {qweights};
        domc = {dom};
        ts = RefTaylorState(domc, domparams, zk, flux, B0c, qnodesc, qweightsc, tol);
        ts = ts.solve(true);
        B = ts.surface_B();

        lerr(i_n, i_nv)    = vecinfnorm(B{1} - B0) / vecinfnorm(B0);
        lsigma(i_n, i_nv)  = norm(ts.sigma{1}, inf);
        npatches(i_n, i_nv) = nu*nv;
        dof(i_n, i_nv)      = sqrt(n^2 * nu * nv);
        fprintf('prepare_torus: n=%d, nu=%d, nv=%d (npatches=%d): alpha = % .6e %+.6ei, surface-B rel. error = %.6e, |sigma|_inf = %.6e\n', ...
            n, nu, nv, nu*nv, real(ts.alpha), imag(ts.alpha), lerr(i_n, i_nv), lsigma(i_n, i_nv));
    end
end

fprintf(['\ntestconv_prepare_torus: compare this error trend to ' ...
    'conv_square_torus.m''s (0.604, 0.669, 0.722 at m=6,8,10 before the ' ...
    'orientation fix). If this baseline decreases with resolution while ' ...
    'my_square_torus does not, the bug is specific to the square-torus ' ...
    'geometry/pipeline.\n']);

colors = lines(numel(ns));

figure(1); clf;
for i_n = 1:numel(ns)
    loglog(dof(i_n,:), lerr(i_n,:), 'o-', 'Color', colors(i_n,:), ...
        'DisplayName', sprintf('$n = %d$', ns(i_n)));
    hold on;
end
xlabel('$n \sqrt{n_u n_v}$', 'Interpreter', 'latex');
ylabel('rel. $L^\infty$ surface-$B$ error', 'Interpreter', 'latex');
title('prepare\_torus: surface-B convergence', 'Interpreter', 'none');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;

figure(2); clf;
for i_n = 1:numel(ns)
    loglog(dof(i_n,:), lsigma(i_n,:), 's--', 'Color', colors(i_n,:), ...
        'DisplayName', sprintf('$n = %d$', ns(i_n)));
    hold on;
end
xlabel('$n \sqrt{n_u n_v}$', 'Interpreter', 'latex');
ylabel('$\|\sigma\|_\infty$', 'Interpreter', 'latex');
title('prepare\_torus: $\|\sigma\|_\infty$', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf), ...
    norm(f.components{2}, inf), ...
    norm(f.components{3}, inf)]);
end
