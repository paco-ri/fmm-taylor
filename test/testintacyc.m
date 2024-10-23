addpath('../')

n = 4; 
nv = 10;
nu = nv*3;
rmaj = 4.0;
rmin = 2.0;
dom = circulartorus(n,nu,nv,rmin,rmaj);

% cos phi = x/(rmaj + rmin cos theta), sin phi = y/(rmaj + rmin cos theta)
costheta = @(x,y,z) sign(x.^2+y.^2-rmaj^2).*sqrt(1 - z.^2./(rmin^2));
sintheta = @(x,y,z) z/rmin;
cosphi = @(x,y,z) x./(rmaj + rmin*costheta(x,y,z));
sinphi = @(x,y,z) y./(rmaj + rmin*costheta(x,y,z));
f = surfacefunv(@(x,y,z) -sintheta(x,y,z).*cosphi(x,y,z), ...
    @(x,y,z) -sintheta(x,y,z).*sinphi(x,y,z), ...
    @(x,y,z) costheta(x,y,z), dom);


% plot(norm(f))
% hold on
% plot(dom)
% colorbar
% 
% g = surfacefunv(@(x,y,z) -costheta(x,y,z).^1.*sintheta(x,y,z).*cosphi(x,y,z), ...
%     @(x,y,z) -costheta(x,y,z).^1.*sintheta(x,y,z).*sinphi(x,y,z), ...
%     @(x,y,z) costheta(x,y,z).^1.*costheta(x,y,z), dom);

g = surfacefunv(@(x,y,z) -costheta(x,y,z).*sinphi(x,y,z), ...
    @(x,y,z) costheta(x,y,z).*cosphi(x,y,z), ...
    @(x,y,z) 1.*z, dom);

S = surfer.surfacemesh_to_surfer(dom);
G = surfacefun_to_array(g,dom,S);
gg = array_to_surfacefun(G,dom,S);
diff = g - gg;

% plot(norm(diff))
% hold on
% colorbar

% h = f + 2i.*f + g - 1i.*g;
h = surfacefunv(@(x,y,z) 1.*x, @(x,y,z) 1i.*x+x.^3, @(x,y,z) 1.*z, dom);
% integ = TaylorState.intacyc(h,n,nv);
integ = TaylorState.intbcyc(h,n,nu);
disp(integ)
