function dom = twisted_ellipse_torus(a, a0, b, n, nu, nv)

if ( nargin < 5 )
    nu = 8;
end

if ( nargin < 6 )
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
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, a, a0, b);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalTorus(u, v, a, a0, b)

r = a*cos(v).*cos(u) - b*sin(v).*sin(u) + a0;

x = r.*cos(u);
y = r.*sin(u);
z = a*cos(v).*sin(u) + b*cos(u).*sin(v);
 
end
