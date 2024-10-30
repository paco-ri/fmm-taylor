function [qnodes, qweights] = poloidalfluxquadaxi(nr,np,ro,ao,ri,ai)
%POLOIDALFLUXQUAD Computes quadrature for poloidal cross-section of
%   toroidal_shell
%   Gauss-Legendre in r, periodic trapezoidal in phi
% 
%   Arguments:
%     nr: [int] number of r nodes
%     np: [int] number of phi nodes
%   Returns:
%     qnodes: [double(3,nr*np)] quadrature nodes
%     qweights: [double(1,nr*np)] quadrature weights

[rnodes, rweights] = chebpts(nr,[0 1],1);
qnodes = zeros([3 nr*np]);
qweights = zeros([1 nr*np]);
for i = 1:nr
    rr = rnodes(i);
    wr = rweights(i);
    for j = 1:np
        ij = (i-1)*np+j;
        pp = 2*pi*(j-1)/np;
        [go1, go2] = evalTorus(pp,0,ro,ao);
        [gi1, gi2] = evalTorus(pp,0,ri,ai);
        [dgo1, dgo2] = duEvalTorus(pp,0,ro,ao);
        [dgi1, dgi2] = duEvalTorus(pp,0,ri,ai);
        qnodes(:,ij) = (1-rr)*[gi1; gi2; 0] + rr*[go1; go2; 0];
        qweights(1,ij) = (2*pi/np) ...
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

function [x, y, z] = duEvalTorus(u, v, rmaj, rmin)

x = -(rmaj + rmin*cos(v)).*sin(u);
y = (rmaj + rmin*cos(v)).*cos(u);
z = 0.*u;

end