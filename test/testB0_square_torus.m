%% Reference field B0 tangency check on my_square_torus
%
% conv_square_torus.m builds the reference field B0 via
% reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,zk), with rmin=rmaj=1.0
% -- the field of a planar current ring of radius rmin centered a distance
% rmaj away. For zk=0, B0 is curl-free/div-free everywhere off the ring, and
% is the *exact* harmonic field for any torus that happens to be one of its
% flux surfaces (i.e. a surface to which B0 is everywhere tangent, B0.vn=0).
%
% relative error vecinfnorm(B-B0)/vecinfnorm(B0) only converges to 0 if B0
% IS (very close to) the true solution on dom -- i.e. only if dom is
% (approximately) a flux surface of this particular ring field. If B0 is
% not tangent to dom (B0.vn = O(|B0|), not O(1) smaller), then B0 is simply
% the wrong reference field for dom, and (B - B0) converges to a fixed
% nonzero continuum difference rather than to 0, regardless of resolution.
%
% This test computes max|B0.vn| / max|B0| on my_square_torus with the
% rmin=rmaj=1.0 used by conv_square_torus.m, and on a circulartorus baseline
% with rmin/rmaj chosen so the torus is (close to) a flux surface, for
% comparison.

clear
addpath('../examples')

zk = 0;
ntheta = 1e3;
jmag = 1.0;

%% my_square_torus with conv_square_torus.m's rmin=rmaj=1.0
n = 5;
m = 6;
nu = 3*m;
nv = 4*m;
rmin = 1.0;
rmaj = 1.0;
dom = my_square_torus(n, nu, nv);
domain = Domain(dom, [n, nu, nv]);
B0 = reftaylorsurffun(dom, n, nu, nv, ntheta, rmin, rmaj, jmag, zk);

B0dotvn = B0.components{1}.*domain.vn{1}.components{1} ...
        + B0.components{2}.*domain.vn{1}.components{2} ...
        + B0.components{3}.*domain.vn{1}.components{3};

maxB0n = max([norm(B0.components{1}, inf), ...
               norm(B0.components{2}, inf), ...
               norm(B0.components{3}, inf)]);
maxB0dotvn = norm(B0dotvn, inf);

fprintf('square torus (rmin=%.2f, rmaj=%.2f): max|B0| = %.6e, max|B0.vn| = %.6e, ratio = %.6f\n', ...
    rmin, rmaj, maxB0n, maxB0dotvn, maxB0dotvn/maxB0n);

%% prepare_torus baseline, using the same (rmin,rmaj)=(2.0,2.0) as convoop.m
% (the known-converging example) for its reference field.
circ_rmin = 2.0;
circ_rmaj = 2.0;
circ_nu = 24;
circ_nv = 8;
circdom = prepare_torus(n, circ_nu, circ_nv, 12, 100);
circdom = circdom{1}; % prepare_torus always returns a cell of surfacemeshes
circdomain = Domain(circdom, [n, circ_nu, circ_nv]);
B0c = reftaylorsurffun(circdom, n, circ_nu, circ_nv, ntheta, circ_rmin, circ_rmaj, jmag, zk);

B0cdotvn = B0c.components{1}.*circdomain.vn{1}.components{1} ...
         + B0c.components{2}.*circdomain.vn{1}.components{2} ...
         + B0c.components{3}.*circdomain.vn{1}.components{3};

maxB0cn = max([norm(B0c.components{1}, inf), ...
                norm(B0c.components{2}, inf), ...
                norm(B0c.components{3}, inf)]);
maxB0cdotvn = norm(B0cdotvn, inf);

fprintf('prepare_torus (rmin=%.2f, rmaj=%.2f, as in convoop.m): max|B0| = %.6e, max|B0.vn| = %.6e, ratio = %.6f\n', ...
    circ_rmin, circ_rmaj, maxB0cn, maxB0cdotvn, maxB0cdotvn/maxB0cn);

fprintf(['\ntestB0_square_torus: if the square-torus ratio is O(1) while the ' ...
    'circulartorus ratio is small, B0 (rmin=rmaj=1.0) is not (close to) ' ...
    'tangent to my_square_torus, i.e. it is the wrong reference field for ' ...
    'this geometry -- (B - B0) would then converge to a fixed nonzero ' ...
    'continuum difference rather than to 0, independent of mesh resolution. ' ...
    'In that case, conv_square_torus.m needs rmin/rmaj values for which ' ...
    'my_square_torus is (close to) a flux surface of reftaylor''s ring field, ' ...
    'or B0 needs to be replaced by a reference field actually tangent to ' ...
    'my_square_torus.\n']);
