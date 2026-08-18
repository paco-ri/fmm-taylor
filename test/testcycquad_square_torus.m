%% Cycle quadrature test on the (reparametrized) my_square_torus
%
% my_square_torus now produces a single ring of nu*nv patches with the
% same (u,v) = (toroidal angle, cross-section position) convention as
% prepare_torus. This test checks that, on this geometry:
%
%   * Domain.acycquad / TaylorState.intacyc trace the closed
%     cross-section loop at phi = 0 ("A-cycle"), and
%   * Domain.bcycquad / TaylorState.intbcyc trace the closed toroidal
%     loop at a fixed cross-section point ("B-cycle"),
%
% by checking the circulation of fields with a known closed-loop
% circulation of 2*pi, and that the *quad/*cyc routines agree with each
% other.

clear

for m = [6 8]
    n  = 6;
    nu = 3*m;  % toroidal patch count
    nv = 4*m;  % cross-section patch count (divisible by 4)

    dom = my_square_torus(n, nu, nv);
    domparams = [n, nu, nv];

    %% A-cycle: cross-section loop at phi = 0 (x > 0, y = 0 half-plane)
    % Planar vortex centered inside the rectangular cross-section
    % [rinner,router] x [zb,zt] = [0.5,1] x [-0.25,0.25]. Its circulation
    % around any loop enclosing (xc,zc) once has magnitude 2*pi.
    xc = 0.75;
    zc = 0.0;
    f = surfacefunv(@(x,y,z) -(z-zc)./((x-xc).^2+(z-zc).^2), ...
                     @(x,y,z) 0*x, ...
                     @(x,y,z)  (x-xc)./((x-xc).^2+(z-zc).^2), dom);

    integrala = TaylorState.intacyc(f, n, nu, nv);
    erra = abs(abs(integrala) - 2*pi);
    fprintf('m = %2d: A-cycle circulation = % .12e  (|err| vs 2*pi = %.3e)\n', ...
        m, integrala, erra);
    assert(erra < 1e-8, 'A-cycle circulation magnitude does not match 2*pi');

    % Cross-check against Domain.acycquad directly (used for obj.aquad).
    [xa, xva, wa] = Domain.acycquad(dom, domparams);
    assert(max(abs(xa(:,2))) < 1e-10, ...
        'Domain.acycquad points are not at phi = 0 (y ~= 0)');
    assert(all(xa(:,1) > 0), ...
        'Domain.acycquad points are not on the x > 0 half-plane');
    fa = [-(xa(:,3)-zc)./((xa(:,1)-xc).^2+(xa(:,3)-zc).^2), ...
          zeros(size(xa,1),1), ...
          (xa(:,1)-xc)./((xa(:,1)-xc).^2+(xa(:,3)-zc).^2)];
    acyc_via_quad = sum(wa .* sum(fa.*xva, 2));
    assert(abs(acyc_via_quad - integrala) < 1e-10, ...
        'Domain.acycquad and TaylorState.intacyc disagree');

    % NOTE on sign convention: the cross-section loop traced here goes
    % outer wall (z: zb->zt) -> top (r: router->rinner) -> inner wall
    % (z: zt->zb) -> bottom (r: rinner->router), which is
    % counterclockwise in the (z,x) plane, i.e. the positive orientation
    % for a disk with normal +yhat (since zhat x xhat = yhat). This is
    % the same disk/normal convention used by square_flux_quad.m. The
    % vortex field above is set up so that a counterclockwise loop in the
    % (x,z) plane gives +2*pi; since our loop is counterclockwise in
    % (z,x) (= clockwise in (x,z)), expect integrala ~= -2*pi. Test 5
    % (testflux_square_torus.m) checks this sign against
    % mtxfluxalphanontaylor's convention for curl S0[mH].

    %% B-cycle: toroidal loop at a fixed cross-section point
    % Azimuthal field with circulation 2*pi around any loop that winds
    % once around the z-axis.
    g = surfacefunv(@(x,y,z) -y./(x.^2+y.^2), ...
                     @(x,y,z)  x./(x.^2+y.^2), ...
                     @(x,y,z) 0*z, dom);

    integralb = TaylorState.intbcyc(g, n, nu, nv);
    errb = abs(abs(integralb) - 2*pi);
    fprintf('m = %2d: B-cycle circulation = % .12e  (|err| vs 2*pi = %.3e)\n', ...
        m, integralb, errb);
    assert(errb < 1e-8, 'B-cycle circulation magnitude does not match 2*pi');

    % Cross-check against Domain.bcycquad directly (used for obj.bquad).
    [xb, xub, wb] = Domain.bcycquad(dom, domparams);
    gb = [-xb(:,2)./(xb(:,1).^2+xb(:,2).^2), ...
           xb(:,1)./(xb(:,1).^2+xb(:,2).^2), ...
           zeros(size(xb,1),1)];
    bcyc_via_quad = sum(wb .* sum(gb.*xub, 2));
    assert(abs(bcyc_via_quad - integralb) < 1e-10, ...
        'Domain.bcycquad and TaylorState.intbcyc disagree');

    % The B-cycle should sit at a single (r,z) cross-section point and
    % sweep the full toroidal angle.
    rb = sqrt(xb(:,1).^2 + xb(:,2).^2);
    assert(max(rb) - min(rb) < 1e-10, ...
        'B-cycle points are not at a fixed cross-section radius');
    assert(max(abs(xb(:,3) - xb(1,3))) < 1e-10, ...
        'B-cycle points are not at a fixed z');
end

fprintf('testcycquad_square_torus: all checks passed.\n');
