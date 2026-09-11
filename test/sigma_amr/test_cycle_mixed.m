% Refine some of the patches that the A-cycle passes through, and check whether all 
% its quadrature points still lie at one toroidal angle.
fprintf('=== test_cycle_mixed start ===\n');

n = 6; nv = 4; nu = 3*nv; rmin = 2.0; rmaj = 4.0;
dom0 = circulartorus(n,nu,nv,rmin,rmaj);
np0 = length(dom0.x);
p0 = [(1:np0).', zeros(np0,2)];
rmax = 3;
phi   = @(P) atan2(P(:,2), P(:,1));
theta = @(P) atan2(P(:,3), sqrt(P(:,1).^2+P(:,2).^2) - rmaj);

cases = {1, [1 2], 1:4, [1 5 9]};
names = {'only patch 1', 'patches 1-2', 'all acyc patches 1-4', 'bcyc patches'};

for c = 1:numel(cases)
    [domR, qfR, p2qR] = surfacemesh.refine_leaves(dom0, p0, cases{c}, rmax);
    fR = surfacefunv(@(x,y,z) 0*x+1, @(x,y,z) 0*x, @(x,y,z) 0*x, domR);
    [~, qa, wa] = TaylorState.intacyc(fR, n, nu, nv, qfR, p2qR);
    [~, qb, wb] = TaylorState.intbcyc(fR, n, nu, nv, qfR, p2qR);
    qa = qa(wa~=0,:); qb = qb(wb~=0,:);
    ua = uniquetol(phi(qa), 1e-8);
    ub = uniquetol(theta(qb), 1e-8);
    fprintf('%-24s npat=%3d | acyc %2d pts, %d distinct phi %s | bcyc %2d pts, %d distinct theta\n', ...
        names{c}, length(domR.x), size(qa,1), numel(ua), ...
        mat2str(round(ua(:).',6)), size(qb,1), numel(ub));
end
fprintf('DONE_SENTINEL\n');
