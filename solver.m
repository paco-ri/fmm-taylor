clear

% Get torus geometry from surface-hps
rmin = 2;
rmaj = 10;
nu = 8;
nv = 8;
n = 20;

x = cell(nu*nv,1);
y = cell(nu*nv,1);
z = cell(nu*nv,1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        x{k} = (rmin*cos(uu)+rmaj).*cos(vv);
        y{k} = (rmin*cos(uu)+rmaj).*sin(vv);
        z{k} = rmin*sin(uu);
        k = k+1;
    end
end

dom = surfacemesh(x,y,z);
plot(dom)