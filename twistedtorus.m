clear

n = 20;
a = 1;
a0 = 3;
b = .4;
nconst = 1;
nu = nconst*8;
nv = nconst*24;

dom = twisted_ellipse_torus(a, a0, b, n, nu, nv);

% [xt, yt, zt] = xthetaoversqrtdetG(a, a0, b, n, nu, nv);
% [xp, yp, zp] = xphioversqrtdetG(a, a0, b, n, nu, nv);
eu = surfacefunv(dom.xu, dom.yu, dom.zu, dom);
ev = surfacefunv(dom.xv, dom.yv, dom.zv, dom); 
G = surfacefun(dom.J, dom); 
guu = surfacefun(dom.E, dom); 
guv = surfacefun(dom.F, dom); 

vn = cross(eu, ev); 

% f = surfacefunv(xp, yp, zp, dom);
% f = surfacefunv(xt, yt, zt, dom); 
f1 = eu./(guu.*sqrt(G));
% f2 = (-guv.*eu + guu.*ev)./(guu.*sqrt(G)); 
divf = divergence(f1);

plot(divf)
hold on
colormap('turbo')
colorbar

maxEst(divf)

function g = detG(a, a0, b, u, v)

r = a*cos(v);
z = b*sin(v);
rp = -a*sin(v);
zp = b*cos(v);

g = .5*(((3 + cos(2*u)).*r.^2 + 4*cos(u).*r.*(a0 - sin(u).*z) + ...
    2*(a0 - sin(u).*z).^2).*rp.^2 + 4*r.*z.*rp.*zp + ...
    (2*a0^2 + 2*cos(u).^2.*r.^2 - 4*a0*sin(u).*z - ...
    (-3 + cos(2*u)).*z.^2 + 4*cos(u).*r.*(a0 - sin(u).*z)).*zp.^2);

end

function [x, y, z] = xthetaoversqrtdetG(a, a0, b, n, nu, nv)

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        x{k} = (-a*sin(vv).*cos(uu) - b*cos(vv).*sin(uu)).*cos(uu);
        y{k} = (-a*sin(vv).*cos(uu) - b*cos(vv).*sin(uu)).*sin(uu);
        z{k} = -a*sin(vv).*sin(uu) + b*cos(vv).*cos(uu);
        sdG = sqrt(detG(a, a0, b, uu, vv));
        x{k} = x{k}./sdG;
        y{k} = y{k}./sdG;
        z{k} = z{k}./sdG;
        k = k+1;
    end
end

end

function [x, y, z] = xphioversqrtdetG(a, a0, b, n, nu, nv)

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        
        rr = a*cos(vv);
        zz = b*sin(vv);
        
        c = rr.*cos(uu) - zz.*sin(uu) + a0;
        cphi = -rr.*sin(uu) - zz.*cos(uu);

        x{k} = cphi.*cos(uu) - c.*sin(uu);
        y{k} = cphi.*sin(uu) + c.*cos(uu);
        z{k} = rr.*cos(uu) - zz.*sin(uu);
        sdG = sqrt(detG(a, a0, b, uu, vv));
        x{k} = x{k}./sdG;
        y{k} = y{k}./sdG;
        z{k} = z{k}./sdG;
        k = k+1;
    end
end

end