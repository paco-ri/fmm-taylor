function [qnodes, qweights] = square_flux_quad(nt, router, rinner, zb, zt)
%SQUARE_FLUX_QUAD Tensor-product Chebyshev quadrature for the square-torus cross-section.
%   This is the square-torus analogue of ellipsefluxquad.m, but it uses a
%   single quadrature set on the rectangular cross-section with the toroidal
%   angle fixed.
%
%   Arguments:
%     nt     : number of nodes per side
%     router : outer radius of the square cross-section (default 1)
%     rinner : inner radius of the square cross-section (default 0.5)
%     zb     : bottom z value (default -0.25)
%     zt     : top z value (default 0.25)
%
%   Returns:
%     qnodes   : [double(3,nt*nt)] quadrature nodes
%     qweights : [double(1,nt*nt)] quadrature weights

if nargin < 1 || isempty(nt)
    nt = 12;
end
if nargin < 2 || isempty(router)
    router = 1;
end
if nargin < 3 || isempty(rinner)
    rinner = 0.5 * router;
end
if nargin < 4 || isempty(zb)
    zb = -0.25 * router;
end
if nargin < 5 || isempty(zt)
    zt = 0.25 * router;
end

[rnodes, rwts] = chebpts(nt, [rinner router], 1);
[znodes, zwts] = chebpts(nt, [zb zt], 1);

qnodes = zeros(3, nt*nt);
qweights = zeros(1, nt*nt);

for i = 1:nt
    for j = 1:nt
        ij = (i-1)*nt + j;
        rr = rnodes(i);
        zz = znodes(j);
        qnodes(:,ij) = [rr; 0; zz];
        qweights(ij) = rwts(i) * zwts(j);
    end
end

end