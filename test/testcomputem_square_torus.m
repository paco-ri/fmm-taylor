%% Final assembly (compute_m / surface_B) with "exact" alpha
%
% Isolates RefTaylorState.compute_m / TaylorState.surface_B from
% compute_sigma_alpha (RefTaylorState.m, ~line 228) by reusing the alpha
% computed at a HIGH resolution (m = 8) as a stand-in for the "exact" value,
% and plugging it into a LOW-resolution (m = 6) solve in place of that
% solve's own alpha.
%
% If using the high-res alpha at low resolution noticeably reduces the
% surface-B error (relative to the low-res solve's own alpha), the
% non-convergence seen in examples/conv_square_torus.m is dominated by
% alpha/flux bookkeeping (compute_sigma_alpha, see
% testflux_square_torus.m / testgmres_square_torus.m). If the error barely
% changes, the bottleneck is downstream, in compute_m/surface_B itself
% (near-quadrature for the final layer potentials, see
% testquadcorr_square_torus.m).

clear
addpath('../examples')

tol = 1e-4;
zk = 0;
n = 5;
ntheta = 1e3;
jmag = 1.0;
rmin = 1.0;
rmaj = 1.0;
nt = 12;

%% High-resolution solve (m = 8) -- source of "exact" alpha
m_hi = 8;
[ts_hi, B0_hi, ~] = build_and_solve(n, m_hi, zk, tol, ntheta, jmag, rmin, rmaj, nt);
alpha_hi = ts_hi.alpha;
B_hi = ts_hi.surface_B();
err_hi = vecinfnorm(B_hi{1} - B0_hi{1}) / vecinfnorm(B0_hi{1});
fprintf('m = %d: alpha = % .6e %+.6ei, surface-B rel. error = %.6e\n', ...
    m_hi, real(alpha_hi), imag(alpha_hi), err_hi);

%% Low-resolution solve (m = 6) -- own alpha
m_lo = 6;
[ts_lo, B0_lo, dom_lo] = build_and_solve(n, m_lo, zk, tol, ntheta, jmag, rmin, rmaj, nt);
alpha_lo = ts_lo.alpha;
B_lo = ts_lo.surface_B();
err_lo = vecinfnorm(B_lo{1} - B0_lo{1}) / vecinfnorm(B0_lo{1});
fprintf('m = %d: alpha = % .6e %+.6ei, surface-B rel. error = %.6e (own alpha)\n', ...
    m_lo, real(alpha_lo), imag(alpha_lo), err_lo);

%% Low-resolution assembly with the high-resolution ("exact") alpha
dfunc = solveforD(ts_lo, false);
wfunc = solveforW(ts_lo, false);

ts_lo2 = ts_lo;
ts_lo2.alpha = alpha_hi;
ts_lo2.sigma{1} = 1i*alpha_hi.*dfunc{1} - wfunc{1};
ts_lo2 = ts_lo2.compute_m();
B_lo2 = ts_lo2.surface_B();
err_lo2 = vecinfnorm(B_lo2{1} - B0_lo{1}) / vecinfnorm(B0_lo{1});
fprintf('m = %d: surface-B rel. error = %.6e (alpha from m = %d)\n', ...
    m_lo, err_lo2, m_hi);

fprintf(['\ntestcomputem_square_torus: if err_lo2 << err_lo (and closer to ' ...
    'err_hi), alpha/flux bookkeeping (compute_sigma_alpha) explains most ' ...
    'of the non-convergence; if err_lo2 ~ err_lo, look at compute_m / ' ...
    'surface_B and the near-quadrature corrections instead.\n']);

function [ts, B0, dom] = build_and_solve(n, m, zk, tol, ntheta, jmag, rmin, rmaj, nt)
nu = 3*m;
nv = 4*m;
dom = my_square_torus(n, nu, nv);
domparams = [n, nu, nv];

B0 = reftaylorsurffun(dom, n, nu, nv, ntheta, rmin, rmaj, jmag, zk);

[qnodes, qweights] = square_flux_quad(nt);
flux = 0;
for i = 1:nt*nt
    B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes(:,i));
    flux = flux + B0eval(2)*qweights(i);
end

B0c = {B0};
qnodesc = {qnodes};
qweightsc = {qweights};
domc = {dom};
ts = RefTaylorState(domc, domparams, zk, flux, B0c, qnodesc, qweightsc, tol);
ts = ts.solve(true);
B0 = B0c;
end

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf), ...
    norm(f.components{2}, inf), ...
    norm(f.components{3}, inf)]);
end
