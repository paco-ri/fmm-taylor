%% Is the surface-B error concentrated at the 4 sharp edges?
%
% testconv_prepare_torus.m confirms the RefTaylorState pipeline itself
% gives a small, decreasing surface-B error on a smooth torus (2.6e-3 ->
% 1.5e-3). conv_square_torus.m shows the opposite trend on my_square_torus
% (increasing, ~0.60 -> 0.67 -> 0.72), even after the A-cycle orientation
% fix -- so the bug is specific to my_square_torus's geometry/pipeline.
%
% testquadcorr_square_torus.m found curl S0[mH] (a fixed test field) is
% O(1) and roughly UNIFORM across edge vs interior patches (ratio ~1.08).
% This test instead looks at the *actual solve output*: B{1}-B0, sigma{1},
% and m{1}, split into patches touching one of the 4 sharp edges vs
% interior patches, at m=6 and m=8. If the error (and/or sigma/m
% themselves) is much larger on edge patches and that gap *grows* with
% resolution, the sharp edges are the bottleneck (likely via the
% near-quadrature correction's smoothness assumptions across non-C1 patch
% boundaries); if edge and interior are comparable, the bug is global.

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
    ts = ts.solve(true);
    B = ts.surface_B();

    Bdiff = B{1} - B0;

    edge_kv = unique([1, m, m+1, 2*m, 2*m+1, 3*m, 3*m+1, nv]);

    [Bdiff_edge, Bdiff_int] = patchnorms(Bdiff, nu, nv, edge_kv);
    [sigma_edge, sigma_int] = patchnorms(ts.sigma{1}, nu, nv, edge_kv);
    [m_edge, m_int] = patchnorms(ts.m{1}, nu, nv, edge_kv);

    fprintf('m = %2d: rms |B-B0|     edge = %.6e, interior = %.6e (ratio %.3f)\n', ...
        m, Bdiff_edge, Bdiff_int, Bdiff_edge/Bdiff_int);
    fprintf('m = %2d: rms |sigma|    edge = %.6e, interior = %.6e (ratio %.3f)\n', ...
        m, sigma_edge, sigma_int, sigma_edge/sigma_int);
    fprintf('m = %2d: rms |m|        edge = %.6e, interior = %.6e (ratio %.3f)\n', ...
        m, m_edge, m_int, m_edge/m_int);
end

fprintf(['\ntestBerror_edge_square_torus: compare the edge/interior ratios ' ...
    'and their absolute values across m=6 and m=8. A ratio that is large ' ...
    'and growing with resolution implicates the sharp edges (near-' ...
    'quadrature correction across non-C1 patch boundaries); a ratio near 1 ' ...
    'that holds at both resolutions points to a global bug instead.\n']);

function [edge_rms, interior_rms] = patchnorms(f, nu, nv, edge_kv)
% f: a surfacefun or surfacefunv on a my_square_torus single-ring mesh.
if isa(f, 'surfacefunv')
    ncomp = numel(f.components);
else
    ncomp = 1;
end

edge_sum = 0;
edge_npts = 0;
interior_sum = 0;
interior_npts = 0;
for ku = 1:nu
    for kv = 1:nv
        k = (ku-1)*nv + kv;
        if ncomp == 1
            vals = f.vals{k}(:);
        else
            vals = [];
            for c = 1:ncomp
                vals = [vals; f.components{c}.vals{k}(:)]; %#ok<AGROW>
            end
        end
        s = sum(abs(vals).^2);
        if ismember(kv, edge_kv)
            edge_sum = edge_sum + s;
            edge_npts = edge_npts + numel(vals);
        else
            interior_sum = interior_sum + s;
            interior_npts = interior_npts + numel(vals);
        end
    end
end
edge_rms = sqrt(edge_sum/edge_npts);
interior_rms = sqrt(interior_sum/interior_npts);
end
