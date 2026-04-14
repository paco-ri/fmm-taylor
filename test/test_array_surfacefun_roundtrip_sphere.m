% test_array_surfacefun_roundtrip_sphere.m
%
% Generate a surfacemesh.sphere, define a spherical harmonic surfacefun on
% it, convert to array via surfacefun_to_array, convert back via
% array_to_surfacefun, and verify the result matches the original.
% Runs the checks for nsph = 3, 4, 5.

nref = 1;
nsph_vals = 3:6;
err_lb_list  = zeros(size(nsph_vals));
err_slp_list = zeros(size(nsph_vals));

% initialize u and S_u_fun to empty surfacefuns
u = surfacefun();
S_u_fun = surfacefun();
dom = surfacemesh();

for k = 1:numel(nsph_vals)
    nsph = nsph_vals(k);
    fprintf('\n=== nsph = %d ===\n', nsph);
    [err_lb_list(k), err_slp_list(k), u, S_u_fun, dom] = run_checks(nsph, nref);
end

figure;
semilogy(nsph_vals, err_lb_list,  '-o', 'DisplayName', 'Laplace-Beltrami error');
hold on;
semilogy(nsph_vals, err_slp_list, '-s', 'DisplayName', 'Single-layer potential error');
xlabel('nsph');
% ylabel('Relative error');
ylabel('Error');
title('Convergence vs nsph');
legend('Location', 'best');
grid on;

% plot u and S_u_fun for the last case
figure(2)
plot(u)
hold on
plot(dom)
colorbar

figure(3)
plot(S_u_fun)
colorbar

% For error isolation, create a surfacefun just involving S_u_fun.vals{1} defined on dom{1}
p = 1;
dom_single = surfacemesh({dom.x{p}}, {dom.y{p}}, {dom.z{p}});
u_single = surfacefun({u.vals{p}}, dom_single);
figure(4)
plot(u_single)
colorbar

S_u_fun_single = surfacefun({S_u_fun.vals{p}}, dom_single);
figure(5)
plot(S_u_fun_single)
colorbar

% -------------------------------------------------------------------------
function [err_lb, err_slp, u, S_u_fun, dom] = run_checks(nsph, nref)
% RUN_CHECKS  Build a sphere mesh with the given nsph/nref, define a
%             spherical harmonic surfacefun, and verify:
%               1. Laplace-Beltrami eigenvalue recovery
%               2. Single-layer potential eigenvalue

dom = surfacemesh.sphere(nsph, nref);
S   = surfer.surfacemesh_to_surfer(dom);

% --- define a spherical harmonic surfacefun ---
l = 1;
m = 0;
sph  = spherefun.sphharm(l, m);
sol  = surfacefun(@(x,y,z) sph(x,y,z), dom);
lambda = -l * (l + 1);
f = lambda * sol;

% --- Laplace-Beltrami eigenvalue check ---
pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
L.rankdef = true;
u = L.solve();
err_lb = norm(u - sol); %  / norm(sol);
fprintf('Laplace-Beltrami eigenvalue error: %.4e\n', err_lb);

% --- Apply single-layer potential to u and check eigenvalue ---
opts = [];
opts.format = 'rsc';
tic;
eps_slp = 1e-9;
Q = lap3d.dirichlet.get_quadrature_correction(S, eps_slp, [1.0,0], S, opts);
opts.precomp_quadrature = Q;
toc;
fprintf('Quadrature correction computed in %.2f seconds.\n', toc);
u_arr = surfacefun_to_array(u, dom, S);

S_u = lap3d.dirichlet.eval(S, u_arr, S, eps_slp, [1.0,0], opts);
S_u_fun = array_to_surfacefun(S_u, dom, S);
ratio = norm(S_u_fun) / norm(u);
fprintf('Ratio between S_u_fun and u: %.4e (should be close to 1/(2*l+1) = %.4e)\n', ratio, 1/(2*l+1));
% err_slp = norm(S_u_fun - (1 / (2 * l + 1)) * u); % / norm(u);
err_slp = norm(S_u_fun - ratio * u);
fprintf('Single-layer potential error: %.4e\n', err_slp);
end
