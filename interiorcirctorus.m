function [interior,interiorwts] = interiorcirctorus(nphi,ntheta,nr,rmin,rmaj)
%INTERIORCIRCTORUS Generate points and quadrature weights in the interior
%   of the circular torus
%   Arguments:
%     nphi: [int] number of points in phi-dir.
%     ntheta: [int] number of points in theta-dir.
%     nr: [int] number of points in r-dir.
%     rmin: [double] minor radius of torus
%     rmaj: [double] major radius of torus

phis = 2*pi.*(0:nphi-1)./nphi;
thetas = 2*pi.*(0:ntheta-1)./ntheta;
[rs, rwts] = legpts(nr, [0 rmin]); % skip zero
interior = zeros([3 (nr*ntheta+1)*nphi]);

% points on axis
for i = 1:nphi
    interior(1,i) = rmaj*cos(phis(i));
    interior(2,i) = rmaj*sin(phis(i));
    interior(3,i) = 0;
end
interiorwts = zeros([1 (nr*ntheta+1)*nphi]);
for i = 1:nphi
    for j = 1:ntheta
        for k = 1:nr
            interior(1,(i-1)*ntheta*nr+(j-1)*nr+k+nphi) = ...
                (rmaj+rs(k)*cos(thetas(j)))*cos(phis(i));
            interior(2,(i-1)*ntheta*nr+(j-1)*nr+k+nphi) = ...
                (rmaj+rs(k)*cos(thetas(j)))*sin(phis(i));
            interior(3,(i-1)*ntheta*nr+(j-1)*nr+k+nphi) = ...
                rs(k)*sin(thetas(j));
            interiorwts(1,(i-1)*ntheta*nr+(j-1)*nr+k+nphi) = ...
                4*pi*pi*rwts(k)/nphi/ntheta;
        end
    end
end
end