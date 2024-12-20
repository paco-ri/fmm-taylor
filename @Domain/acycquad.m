function [x,xv,w] = acycquad(dom,domparams)
%ACYCQUAD computes A-cycle quadrature for a domain
%   Arguments
%     dom = domain
%     domparams = n,nu,nv of domain

n = domparams(1);
nv = domparams(3);
nnodes = n*nv;

xv = zeros(nnodes,3);
x = zeros(nnodes,3);
w = zeros(nnodes,1);

for i = 1:nv
    xv_patch = [dom.xv{i}(:,1).'; dom.yv{i}(:,1).'; dom.zv{i}(:,1).'].';
    x_patch = [dom.x{i}(:,1).'; dom.y{i}(:,1).'; dom.z{i}(:,1).'].';

    idx_range = (i-1)*n+1:i*n;
    xv(idx_range, :) = xv_patch;
    x(idx_range, :) = x_patch;
    [~, quad_weights] = chebpts(n);
    w(idx_range) = quad_weights;
end

end