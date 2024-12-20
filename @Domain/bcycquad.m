function [x,xu,w] = bcycquad(dom,domparams)
%ACYCQUAD computes A-cycle quadrature for a domain
%   Arguments
%     dom = domain
%     domparams = n,nu,nv of domain

n = domparams(1);
nu = domparams(2);
nv = domparams(3);
nnodes = n*nu;

xu = zeros(nnodes,3);
x = zeros(nnodes,3);
w = zeros(nnodes,1);

for i = 1:nu
    idx = (i-1)*nv+1;
    xu_patch = [dom.xu{idx}(1,:); dom.yu{idx}(1,:); dom.zu{idx}(1,:)].';
    x_patch = [dom.x{idx}(1,:); dom.y{idx}(1,:); dom.z{idx}(1,:)].';

    idx_range = (i-1)*n+1:i*n;
    xu(idx_range, :) = xu_patch;
    x(idx_range, :) = x_patch;
    [~, quad_weights] = chebpts(n);
    w(idx_range) = quad_weights;
end

end