function B = ringfield(rings, zk, X)
%RINGFIELD Reference Taylor-state field of several current rings.
%   B = RINGFIELD(RINGS, ZK, X) is B0 summed over the struct array RINGS at
%   the 3 x N points X, or a surfacefunv when X is a surfacemesh. Each ring
%   has fields c (centre), Q (rotation; Q(:,1) is the ring normal), a
%   (radius), ntheta (quadrature points) and jmag, and gives reftaylor's
%   field rotated by Q and centred at c. Unlike reftaylor at ZK = 0, the
%   ZK = 0 field here is the ZK -> 0 limit.
if isa(X, 'surfacemesh')
    dom = X;
    xc = dom.x; yc = dom.y; zc = dom.z;
    sz = [size(xc{1}), numel(xc)];
    B = ringfield(rings, zk, [reshape(cat(3, xc{:}), 1, []); ...
        reshape(cat(3, yc{:}), 1, []); reshape(cat(3, zc{:}), 1, [])]);
    comp = cell(1, 3);
    for c = 1:3
        comp{c} = reshape(num2cell(reshape(B(c,:), sz), [1 2]), [], 1);
    end
    B = surfacefunv(comp{:}, dom);
    return
end
B = zeros(size(X));
for r = rings
    th = 2*pi*(0:r.ntheta-1)/r.ntheta;
    S = r.c + r.Q*[zeros(size(th)); r.a*cos(th); r.a*sin(th)];
    J = r.jmag*r.a*r.Q*[zeros(size(th)); -sin(th); cos(th)];
    for j = 1:r.ntheta
        d = X - S(:,j);
        dist = vecnorm(d);
        G = exp(1i*zk*dist)./dist;
        gradG = d.*(G.*(1i*zk*dist - 1)./dist.^2);
        gxJ = [gradG(2,:)*J(3,j) - gradG(3,:)*J(2,j); ...
               gradG(3,:)*J(1,j) - gradG(1,:)*J(3,j); ...
               gradG(1,:)*J(2,j) - gradG(2,:)*J(1,j)];
        B = B + (zk*G.*J(:,j) + gxJ)*r.a/(2*r.ntheta);
    end
end
end
