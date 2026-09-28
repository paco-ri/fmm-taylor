%% Does flipping the sign of mH fix surface-B convergence?
%
% Tests 1 and 5 found that the A-cycle loop traced by Domain.acycquad /
% TaylorState.intacyc on my_square_torus is oriented oppositely to the
% +yhat-normal cross-section disk used by square_flux_quad /
% mtxfluxalphanontaylor: the A-cycle circulation of mH is positive
% (~+2.0, testmH_square_torus.m) while fluxalpha (curl S0[mH] flux through
% the +yhat disk) is negative (~-0.24, testflux_square_torus.m) and has the
% opposite sign from the reference flux of B0 (~+0.039).
%
% obj.alpha = -1i*(flux - fluxsigmaW)/(-fluxsigmaD + fluxalpha) is linear in
% 1/fluxalpha (and hence sensitive to its sign), and obj.sigma{1}, obj.m
% both depend on alpha. Test 6 found alpha converges to a *stable* value
% (~0.228i at both m=6 and m=8) but the resulting surface-B error does not
% shrink (and slightly grows) -- consistent with alpha having converged to
% a self-consistent but *wrong* value due to this sign mismatch.
%
% This test re-runs the m=6 and m=8 solves with domain.mH{1} negated
% immediately after Domain construction (before compute_sigma_alpha uses it
% for fluxalpha and for mH itself in compute_m), and reports the resulting
% surface-B errors. If they shrink with resolution (unlike the un-flipped
% errors ~0.60 and ~0.67 from testcomputem_square_torus.m), the sign
% mismatch identified in Tests 1/5 is the root cause, and the fix belongs in
% Domain.compute_mH (or equivalently, in the orientation convention used by
% mtxfluxalphanontaylor/mtxfluxsigmanontaylor's XS quadrature).

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

for m = [6 8]
    nu = 3*m;
    nv = 4*m;
    dom = my_square_torus(n, nu, nv);
    domparams = [n, nu, nv];

    B0 = reftaylorsurffun(dom, n, ntheta, rmin, rmaj, jmag, zk);

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

    % Flip the sign of mH before it is used by compute_sigma_alpha/compute_m
    ts.domain.mH{1} = -ts.domain.mH{1};

    ts = ts.solve(true);
    B = ts.surface_B();
    err = vecinfnorm(B{1} - B0) / vecinfnorm(B0);
    fprintf('m = %d: alpha = % .6e %+.6ei, surface-B rel. error (mH negated) = %.6e\n', ...
        m, real(ts.alpha), imag(ts.alpha), err);
end

fprintf(['\ntestmHsign_square_torus: compare these errors to the ' ...
    'un-flipped errors from testcomputem_square_torus.m (~0.604 at m=6, ' ...
    '~0.669 at m=8). If negating mH gives errors that shrink with m, the ' ...
    'A-cycle/fluxalpha sign mismatch from Tests 1 and 5 is the root cause.\n']);

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf), ...
    norm(f.components{2}, inf), ...
    norm(f.components{3}, inf)]);
end
