function dom = my_square_torus(n, nu, nv)
%MY_SQUARE_TORUS  Genus-one surface with a rectangular cross-section.
%   dom = MY_SQUARE_TORUS(n, nu, nv) builds a single ring of nu*nv patches,
%   following the same (u,v) convention as prepare_torus: u is the toroidal
%   angle (period 2*pi, nu patches) and v parametrizes one trip around the
%   rectangular cross-section (nv patches, must be divisible by 4: bottom
%   annulus, rinner wall, top annulus, router wall, in that order). This
%   traces the cross-section loop clockwise in the (x,z) plane at phi=0,
%   matching the +yhat-normal orientation assumed by square_flux_quad.m.

if ( nargin < 2 )
    nu = 8;
end

if ( nargin < 3 )
    nv = nu;
end

if ( mod(nv, 4) ~= 0 )
    error('my_square_torus:nv', 'nv must be divisible by 4.');
end

router = 1;
rinner = 0.5;
zb = -0.25;
zt = 0.25;

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);

ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 4, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalSquareTorus(uu, vv, router, rinner, zb, zt);
        k = k+1;
    end
end

dom = surfacemesh(x, y, z);

end

function [x, y, z] = evalSquareTorus(u, v, router, rinner, zb, zt)

phi = u;

r = zeros(size(v));
z = zeros(size(v));

% Face 1, v in [0,1]: bottom annulus (z = zb), r from router to rinner
m = (v >= 0) & (v <= 1);
s = v(m);
r(m) = router + s.*(rinner - router);
z(m) = zb;

% Face 2, v in [1,2]: inner wall (r = rinner), z from zb to zt
m = (v > 1) & (v <= 2);
s = v(m) - 1;
r(m) = rinner;
z(m) = zb + s.*(zt - zb);

% Face 3, v in [2,3]: top annulus (z = zt), r from rinner to router
m = (v > 2) & (v <= 3);
s = v(m) - 2;
r(m) = rinner + s.*(router - rinner);
z(m) = zt;

% Face 4, v in [3,4]: outer wall (r = router), z from zt to zb
m = (v > 3) & (v <= 4);
s = v(m) - 3;
r(m) = router;
z(m) = zt + s.*(zb - zt);

x = r .* cos(phi);
y = r .* sin(phi);

end
