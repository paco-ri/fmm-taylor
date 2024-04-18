function dom = circulartorus(n, nu, nv)

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = nu;
end

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalTorus(u, v)

rmaj = 5.0;
rmin = 2.0;
x = (rmaj + rmin*cos(u)).*cos(v);
y = (rmaj + rmin*cos(u)).*sin(v);
z = rmin*sin(u);
 
end