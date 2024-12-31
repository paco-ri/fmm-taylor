function [qnodes, qweights] = ellipsefluxquad(nr, nt, a, a0, b)
%TORUSFLUXQUAD Computes quadrature for cross-section of
%   twisted_ellipse_torus
%   Gauss-Legendre in r, periodic trapezoidal in theta
% 
%   Arguments:
%     nr: [int] number of r nodes
%     nt: [int] number of theta nodes
%   Returns:
%     qnodes: [double(3,nr*nt)] quadrature nodes
%     qweights: [double(1,nr*nt)] quadrature weights

[rnodes, rweights] = chebpts(nr,[0 1],1);
qnodes = zeros([3 nr*nt]);
qweights = zeros([1 nr*nt]);
for i = 1:nr
    rr = rnodes(i);
    wr = rweights(i);
    for j = 1:nt
        ij = (i-1)*nt+j;
        tt = 2*pi*(j-1)/nt;
        [g1, ~, g2] = evalTorus(0,tt,a,a0,b);
        [dg1, ~, dg2] = dvEvalTorus(0,tt,a,b);
        g1 = g1 - a0; 
        qnodes(:,ij) = rr*[g1; 0; g2] + [a0; 0; 0];
        qweights(1,ij) = (2*pi/nt) * wr*rr*(g1*dg2-g2*dg1);
    end
end

end

function [x, y, z] = evalTorus(u, v, a, a0, b)

r = a*cos(v).*cos(u) - b*sin(v).*sin(u) + a0;

x = r.*cos(u);
y = r.*sin(u);
z = a*cos(v).*sin(u) + b*cos(u).*sin(v);
 
end

function [x, y, z] = dvEvalTorus(u, v, a, b)

drdv = -a*sin(v).*cos(u) - b*cos(v).*sin(u);

x = drdv.*cos(u);
y = drdv.*sin(u);
z = -a*sin(v).*sin(u) + b*cos(u).*cos(v);
 
end