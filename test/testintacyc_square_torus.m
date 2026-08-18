clear

% A-cycle integration test on a square torus.
% The field is a planar vortex in the x-z plane centered inside the
% cross-section hole. Its circulation around one poloidal loop is |2*pi|.

n = 8;
nu = 4;
nv = 3*nu;

dom = my_square_torus(n, nu, nv);

xc = 0.75;
zc = 0.0;

f = surfacefunv(@(x,y,z) -(z - zc)./((x - xc).^2 + (z - zc).^2), ...
    @(x,y,z) 0*x, ...
    @(x,y,z)  (x - xc)./((x - xc).^2 + (z - zc).^2), dom);

% [integrala, quadpts, quadwts] = TaylorState.intacyc(f, n, nu, nv);
% [integrala, quadpts, quadwts] = TaylorState.intbcyc(f, n, nu, nv);
acyc_patch_list = 1:nv:4*nu*nv;
% acyc_patch_list = 1:nv:nu*nv;
% acyc_patch_list = acyc_patch_list + 3*nu*nv;
rotate = false(4*nu, 1);
for i = 1:nu
    rotate(nu+i) = true;
    rotate(2*nu+i) = true;
end
[integrala, quadpts, quadwts] = temp_intacyc(f, acyc_patch_list, rotate);
surfacemesh_to_vtk(dom, "square_torus_acyc.vtk", "Points", quadpts);

expected = 2*pi;
tol = 1e-8;
err = abs(abs(integrala) - expected);

if (err < tol)
    fprintf("square_torus intacyc test passed!\terror = %.6e\tintegrala = %.16e\n", ...
        err, integrala);
else
    error("square_torus intacyc test failed! error = %.6e, integrala = %.16e, expected magnitude = %.16e", ...
        err, integrala, expected);
end

function [integrala, quadpts, quadwts] = temp_intacyc(f, patch_list, rotate)
dom = f.domain;
num_patches = length(patch_list);
patch_order = size(dom(1).x{1}, 1);
nalloc = num_patches * patch_order;
xderiv = zeros(nalloc, 3);
x = zeros(nalloc, 3);
avals = zeros(nalloc, 3);
awts = zeros(nalloc, 1);

j = 1;
for i = patch_list
    idx_start = (j-1)*patch_order + 1;
    idx_end = j*patch_order;
    if rotate(j)
        x(idx_start:idx_end, :) = [dom.x{i}(:,1).'; dom.y{i}(:,1).'; dom.z{i}(:,1).'].';
        xderiv(idx_start:idx_end, :) = -1.*[dom.xv{i}(:,1).'; dom.yv{i}(:,1).'; dom.zv{i}(:,1).'].';
        avals(idx_start:idx_end, 1) = f.components{1}.vals{i}(:,1);
        avals(idx_start:idx_end, 2) = f.components{2}.vals{i}(:,1);
        avals(idx_start:idx_end, 3) = f.components{3}.vals{i}(:,1);
    else
        x(idx_start:idx_end, :) = [dom.x{i}(1,:); dom.y{i}(1,:); dom.z{i}(1,:)].';
        xderiv(idx_start:idx_end, :) = [dom.xu{i}(1,:); dom.yu{i}(1,:); dom.zu{i}(1,:)].';
        avals(idx_start:idx_end, 1) = f.components{1}.vals{i}(1,:);
        avals(idx_start:idx_end, 2) = f.components{2}.vals{i}(1,:);
        avals(idx_start:idx_end, 3) = f.components{3}.vals{i}(1,:);
    end
    [~, awts(idx_start:idx_end)] = chebpts(patch_order);
    j = j + 1;
end

axu = dot(conj(avals), xderiv, 2);
integrala = dot(awts, axu);
quadpts = x;
quadwts = awts;

end