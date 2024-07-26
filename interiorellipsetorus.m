function [interior,interiorwts] = interiorellipsetorus(nphi,ntheta,nr,a,a0,b)
%INTERIORELLIPSETORUS Generate points and quadrature weights in the 
% interior of twisted_ellipse_torus
%   Arguments:
%     nphi: [int] number of points in phi-dir.
%     ntheta: [int] number of points in theta-dir.
%     nr: [int] number of points in r-dir.

phis = 2*pi.*(0:nphi-1)./nphi;
thetas = 2*pi.*(0:ntheta-1)./ntheta;
[rs, rwts] = legpts(nr, [0 1]); % skip zero 
interior = zeros([3 (nr*ntheta+1)*nphi]);

% points on axis
for i = 1:nphi
    interior(1,i) = a0*cos(phis(i));
    interior(2,i) = a0*sin(phis(i));
    interior(3,i) = 0;
end
interiorwts = zeros([1 (nr*ntheta+1)*nphi]);
ind = nphi+1;
for i = 1:nphi
    phi = phis(i);
    for j = 1:ntheta
        theta = thetas(j);
        for k = 1:nr
            r = rs(k);
            [interior(1,ind), interior(2,ind), interior(3,ind)] = ...
                evalTorus(phi,theta,a,a0,b,r);
            interiorwts(1,ind) = 4*pi*pi*rwts(k)/nphi/ntheta;
            ind = ind+1;
        end
    end
end
% plot3(interior(1,:),interior(2,:),interior(3,:),'.')
end

function [x, y, z] = evalTorus(u, v, a, a0, b, varargin)

if nargin == 6
    scale = varargin{1};
else
    scale = 1;
end

r = scale.*(a*cos(v).*cos(u) - b*sin(v).*sin(u)) + a0;

x = r.*cos(u);
y = r.*sin(u);
z = scale.*(a*cos(v).*sin(u) + b*cos(u).*sin(v));
 
end