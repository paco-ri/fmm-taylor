%% Flux bookkeeping sign check on the (reparametrized) my_square_torus
%
% Isolates fluxalpha = mtxfluxalphanontaylor(domain, xs_nodes, xs_weights,
% mH, zk, ...) -- the term \oint curl S0[mH] . yhat dA over the
% cross-section disk used in compute_sigma_alpha (RefTaylorState.m,
% nsurf==1 branch) -- from the rest of the solve, and checks its sign
% against:
%   * the A-cycle circulation of mH itself (TaylorState.intacyc, see
%     testmH_square_torus.m), via Stokes' theorem
%     (curl S0[mH] flux through the disk = circulation of S0[mH] around
%     its boundary, the A-cycle),
%   * the prescribed flux of the reference field B0 through the same
%     cross-section disk (square_flux_quad.m + reftaylor), which is the
%     quantity fluxalpha is ultimately compared against in
%     obj.alpha = -1i*(flux - fluxsigmaW)/(-fluxsigmaD + fluxalpha).
%
% A sign mismatch here (relative to circulartorus) would directly explain
% a wrong/non-converging alpha after the my_square_torus u/v swap.

clear
addpath('../examples')

zk = 0;
epstaylor = 1e-7;
epslh = 1e-7;
nt = 8;

%% Square torus
n = 5;
m = 6;
nu = 3*m;
nv = 4*m;
dom = my_square_torus(n, nu, nv);
domparams = [n, nu, nv];
domain = Domain(dom, domparams);

[qnodes, qweights] = square_flux_quad(nt);

fluxalpha_sq = TaylorState.mtxfluxalphanontaylor(domain, qnodes, qweights, ...
    domain.mH{1}, zk, epstaylor, epslh);
acyc_mH_sq = TaylorState.intacyc(domain.mH{1}, n, nu, nv);

fprintf('square torus:  fluxalpha = % .6e %+.6ei\n', real(fluxalpha_sq), imag(fluxalpha_sq));
fprintf('square torus:  mH A-cycle circulation = % .6e %+.6ei\n', ...
    real(acyc_mH_sq), imag(acyc_mH_sq));

% Reference flux of B0 through the same cross-section disk
ntheta = 1e3;
jmag = 1.0;
rmin = 1.0;
rmaj = 1.0;
flux_sq = 0;
for i = 1:nt*nt
    B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes(:,i));
    flux_sq = flux_sq + B0eval(2)*qweights(i);
end
fprintf('square torus:  reference flux through XS disk = % .6e\n', flux_sq);

%% Circular torus baseline (known-good, convoop.m geometry)
n2 = 5;
nv2 = 8;
nu2 = nv2*3;
ao = 1.0;
ai = 0.6;
nr = 16; ntfine = 40; np = 40;
domc = prepare_torus(n2,nu2,nv2,nr,ntfine);
domc = domc{1}; % prepare_torus always returns a cell of surfacemeshes
domainc = Domain(domc, [n2, nu2, nv2]);

[qnodesc, qweightsc] = square_flux_quad(nt); % same disk shape, generic check

fluxalpha_c = TaylorState.mtxfluxalphanontaylor(domainc, qnodesc, qweightsc, ...
    domainc.mH{1}, zk, epstaylor, epslh);
acyc_mH_c = TaylorState.intacyc(domainc.mH{1}, n2, nu2, nv2);

fprintf('circular torus: fluxalpha = % .6e %+.6ei\n', real(fluxalpha_c), imag(fluxalpha_c));
fprintf('circular torus: mH A-cycle circulation = % .6e %+.6ei\n', ...
    real(acyc_mH_c), imag(acyc_mH_c));

fprintf(['\ntestflux_square_torus: compare the *signs* of fluxalpha vs. ' ...
    'the mH A-cycle circulation for each geometry. If the square-torus ' ...
    'sign relationship is opposite the circular-torus one, the u/v swap ' ...
    'in my_square_torus.m introduced a cycle-orientation flip that ' ...
    'should be corrected (e.g. by reversing the v-parametrization of ' ...
    'the cross-section loop in evalSquareTorus).\n']);
