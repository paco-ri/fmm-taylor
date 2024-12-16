function [x,xv,w] = acycquad(dom,mode,varargin)
%ACYCQUAD computes A-cycle quadrature for a domain
%   Arguments
%     dom = domain
%     mode = G-L on dom nodes if 0, trapezoidal off dom nodes else
%   Optional arguments
%     n,nu,nv = if mode == 0, n,nu,nv of domain
%     nnodes = if mode ~= 0, number of nodes

if mode == 0
    n = varargin{1};
    nu = varargin{2};
    nv = varargin{3};
    nnodes = n*nv;
else
    nnodes = varargin{1};
end

xv = zeros(nnodes,3);
x = zeros(nnodes,3);
w = zeros(nnodes,1);

if mode == 0

    for i = 1:nv
        xv_patch = [dom.xv{i}(:,1).', dom.yv{i}(:,1).', dom.zv{i}(:,1).'];
        x_patch = [dom.x{i}(:,1).', dom.y{i}(:,1).', dom.z{i}(:,1).'];

        idx_range = (i-1)*n+1:i*n;
        xv(idx_range, :) = xv_patch;
        x(idx_range, :) = x_patch;
        [~, quad_weights] = chebpts(n);
        w(idx_range) = quad_weights;
    end

else

    for i = 0:nnodes-1
        % WIP
    end

end

end