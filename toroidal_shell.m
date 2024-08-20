function dom = toroidal_shell(ro, ao, bo, ri, ai, bi, n, nu, nv)

if ( nargin < 8 )
    nu = 8;
end

if ( nargin < 9 )
    nv = nu;
end

x = cell(2*nu*nv, 1);
y = cell(2*nu*nv, 1);
z = cell(2*nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, ro, ao, bo);
        k = k+1;
    end
end
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, ri, ai, bi);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalTorus(u, v, r, a, b)

R = a*cos(v).*cos(u) - b*sin(v).*sin(u) + r;
x = R.*cos(u);
y = R.*sin(u);
z = a*cos(v).*sin(u) + b*cos(u).*sin(v);

end