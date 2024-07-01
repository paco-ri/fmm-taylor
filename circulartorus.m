function dom = circulartorus(n, nu, nv, rmin, rmaj)

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = nu;
end

if ( nargin < 4 )
    rmin = 1.0;
end

if ( nargin < 5 )
    rmaj = 2.0;
end

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

% npatches = nu*nv;
% npols = n*n;
% srcvals = zeros([12 npols*npatches]);
% xint = legpts(n);
% xcheb = chebpts(n,[-1,1]);
% eval = barymat(xint,xcheb);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, rmin, rmaj);
        
        % [uuu, vvv] = meshgrid(legpts(n,ubreaks(ku:ku+1)), ...
        %     legpts(n,vbreaks(kv:kv+1)));
        % [x2, y2, z2] = evalTorus(uuu, vvv, rmin, rmaj);
        % [dxdu, dydu, dzdu] = evalTorusDv(uuu, vvv, rmin, rmaj);
        % [dxdv, dydv, dzdv] = evalTorusDu(uuu, vvv, rmin, rmaj);
        % [nx, ny, nz] = evalTorusNormal(uuu, vvv, rmin, rmaj);
        % x2 = x2';
        % y2 = y2';
        % z2 = z2';
        % dxdu = dxdu';
        % dydu = dydu';
        % dzdu = dzdu';
        % dxdv = dxdv';
        % dydv = dydv';
        % dzdv = dzdv';
        % nx = nx';
        % ny = ny';
        % nz = nz';
        % kstart = (k-1)*npols+1;
        % kend = k*npols;
        % srcvals(1,kstart:kend) = x2(:);
        % srcvals(2,kstart:kend) = y2(:);
        % srcvals(3,kstart:kend) = z2(:);
        % srcvals(4,kstart:kend) = dxdu(:);
        % srcvals(5,kstart:kend) = dydu(:);
        % srcvals(6,kstart:kend) = dzdu(:);
        % srcvals(7,kstart:kend) = dxdv(:);
        % srcvals(8,kstart:kend) = dydv(:);
        % srcvals(9,kstart:kend) = dzdv(:);
        % srcvals(10,kstart:kend) = nx(:);
        % srcvals(11,kstart:kend) = ny(:);
        % srcvals(12,kstart:kend) = nz(:);

        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

% iptype = 11;
% S = surfer(npatches,n-1,srcvals,iptype);

end

function [x, y, z] = evalTorus(u, v, rmin, rmaj)

x = (rmaj + rmin*cos(u)).*cos(v);
y = (rmaj + rmin*cos(u)).*sin(v);
z = rmin*sin(u);
% x = (rmaj + rmin*cos(v)).*cos(u);
% y = (rmaj + rmin*cos(v)).*sin(u);
% z = rmin*sin(v);
 
end

function [x, y, z] = evalTorusDv(u, v, rmin, rmaj)

x = -(rmaj + rmin*cos(u)).*sin(v);
y = (rmaj + rmin*cos(u)).*cos(v);
z = 0;

end

function [x, y, z] = evalTorusDu(u, v, rmin, rmaj)

x = -rmin*sin(u).*cos(v);
y = -rmin*sin(u).*sin(v);
z = rmin*cos(u);
 
end

function [x, y, z] = evalTorusNormal(u, v, rmin, rmaj)

x = cos(v).*cos(u);
y = sin(v).*cos(u);
z = sin(u);
 
end