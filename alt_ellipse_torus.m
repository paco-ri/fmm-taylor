ro = 2.0;
ao = 0.7;
bo = 1.0;
ri = 2.0;
ai = 0.3;
bi = 0.55;

n = 5;
nv = 4;
nu = 8;

dom = twist(ro, ao, bo, n, nu, nv);
plot(dom)
hold on

nphi = 20;
options = optimset('Display','off');
for j = 0:nphi-1
    phi = 3*j*pi/nphi;
    f = @(u) evalTorusZ(u,phi,ao,bo);
    theta = fsolve(f,-2/3*phi,options);
    [x,y,z] = evalTorus(theta,phi,ro,ao,bo);
    plot3(x,y,z,'o')
end

function dom = twist(r, a, b, n, nu, nv)

if ( nargin < 6 )
    nu = 8;
end

if ( nargin < 7 )
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
        [x{k}, y{k}, z{k}] = evalTorus(uu, vv, r, a, b);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalTorus(u, v, r, a, b)

uu = 3/2.*u;
R = a*cos(v).*cos(uu) - b*sin(v).*sin(uu) + r;
x = R.*cos(u);
y = R.*sin(u);
z = a*cos(v).*sin(uu) + b*cos(uu).*sin(v);

end

function z = evalTorusZ(u, v, a, b)

uu = 3/2.*u;
z = a*cos(v).*sin(uu) + b*cos(uu).*sin(v);
end