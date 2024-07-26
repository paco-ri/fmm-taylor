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
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

% iptype = 11;
% S = surfer(npatches,n-1,srcvals,iptype);

end

function [x, y, z] = evalTorus(u, v, rmin, rmaj)

x = (rmaj + rmin*cos(v)).*cos(u);
y = (rmaj + rmin*cos(v)).*sin(u);
z = rmin*sin(v);
 
end