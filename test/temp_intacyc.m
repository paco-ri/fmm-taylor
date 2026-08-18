function [integrala, quadpts, quadwts] = temp_intacyc(f, n, nu, nv, varargin)
%TEMP_INTACYC Temporary A-cycle integrator for testing.
%   This wrapper matches the TaylorState.intacyc call signature and uses
%   the same patch selection and orientation logic as
%   testintacyc_square_torus.m.

% Derive the A-cycle patch list and rotation pattern from n, nu, nv.
patch_list = 1:nv:4*nu*nv;
rotate = false(4*nu, 1);
for i = 1:nu
    rotate(nu+i) = true;
    rotate(2*nu+i) = true;
end

dom = f.domain;
num_patches = length(patch_list);
patch_order = n;
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
