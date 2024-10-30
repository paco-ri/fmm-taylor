function [qnodes, qweights] = toroidalfluxquadaxi(nr,nt,ro,ao,ri,ai)
%TOROIDALFLUXQUAD Computes quadrature for toroidal cross-section
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
        [go1, ~, go2] = evalTorus(0,tt,ro,ao);
        [gi1, ~, gi2] = evalTorus(0,tt,ri,ai);
        [dgo1, ~, dgo2] = dvEvalTorus(0,tt,ao);
        [dgi1, ~, dgi2] = dvEvalTorus(0,tt,ai);
        qnodes(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[go1; 0; go2];
        qweights(1,ij) = (2*pi/nt) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

end

function [x, y, z] = evalTorus(u, v, rmaj, rmin)

x = (rmaj + rmin*cos(v)).*cos(u);
y = (rmaj + rmin*cos(v)).*sin(u);
z = rmin*sin(v);
 
end

function [x, y, z] = dvEvalTorus(u, v, rmin)

x = -rmin*sin(v).*cos(u);
y = -rmin*sin(v).*sin(u);
z = rmin*cos(v);

end

