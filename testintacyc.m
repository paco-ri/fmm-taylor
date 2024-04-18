n = 8; 
nu = 4;
nv = nu*3;
dom = circulartorus(n,nu,nv);

rmaj = 5.0;
rmin = 2.0;
% cos phi = x/(rmaj + rmin cos theta), sin phi = y/(rmaj + rmin cos theta)
costheta = @(x,y,z) sign(x.^2+y.^2-rmaj^2).*sqrt(1 - z.^2./(rmin^2));
sintheta = @(x,y,z) z/rmin;
cosphi = @(x,y,z) x./(rmaj + rmin*costheta(x,y,z));
sinphi = @(x,y,z) y./(rmaj + rmin*costheta(x,y,z));
f = surfacefunv(@(x,y,z) -sintheta(x,y,z).*cosphi(x,y,z), ...
    @(x,y,z) -sintheta(x,y,z).*sinphi(x,y,z), ...
    @(x,y,z) costheta(x,y,z), dom);

plot(norm(f))
hold on
plot(dom)
colorbar

integ = intacyc(f,n,nu,nv);
disp(integ)
