function [dom, fluxnodes, fluxwts] = alt_stellarator(n, nu, nv, scale)

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
        [x{k}, y{k}, z{k}] = evalStellarator(uu, vv, scale);
        k = k+1;
    end
end

% Try to reorder to create better merge indices:
if ( bitand(nv, nv-1) == 0 && bitand(nu, nu-1) == 0 )
    ordering = morton(nv, nu);
    ordering = ordering(:);
    x(ordering) = x;
    y(ordering) = y;
    z(ordering) = z;
end

dom = surfacemesh(x, y, z);

fluxnodes = zeros(3,n*nv);
for kv = 1:nv
    v = chebpts(n, vbreaks(kv:kv+1));
    idx = n*(kv-1)+1:n*kv;
    [fluxnodes(1,idx), fluxnodes(2,idx), fluxnodes(3,idx)] = evalStellarator(0, v, scale);
end

fluxwts = zeros(1,nv);

end

function [x, y, z] = evalStellarator(v, u, scale)

Q = 3;
R = 1;
r0 = 0;
z0 = 0;
d = [0.15  0.09  0.00  0.00  0.00;
     0.00  1.00  0.03 -0.01  0.00;
     0.08  4.00 -0.01 -0.02  0.00;
     0.01 -0.28 -0.28  0.03  0.02;
     0.00  0.09 -0.03  0.06  0.00;
     0.01 -0.02  0.02  0.00 -0.02];
d = scale.*d;
d(3,2) = 4.00;

rz = zeros(size(u));
for j = -1:4
    for k = -1:3
        jj = j+2;
        kk = k+2;
        rz = rz + d(jj,kk)*exp(-1i*j*u+1i*k*Q*v);
    end
end
rz = exp(1i*u).*rz;
r1 = real(rz);
z1 = imag(rz);
r = r0 + R*(r1-r0);
z = z0 + R*(z1-z0);

[x, y] = pol2cart(v, r);

end