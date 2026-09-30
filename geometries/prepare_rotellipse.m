function [doms, fluxnodes, fluxwts] = prepare_rotellipse(n, nu, nv, nr, nt)
%PREPARE_ROTELLIPSE One-surface rotating-ellipse stellarator.
%   Same surface as surfacefun-plus's surfacemesh.torus, with u toroidal and
%   v poloidal (swapped there) so normals match prepare_torus. The flux
%   cross-section is the u = pi half-plane (y = 0, x < 0), where the section is
%   an ellipse; nodes fan out from R = 4.458, about which it is star-shaped.

xscenter = -4.458;

x = cell(nu*nv, 1);
y = cell(nu*nv, 1);
z = cell(nu*nv, 1);
ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);

k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        [x{k}, y{k}, z{k}] = evalRotEllipse(uu, vv);
        k = k+1;
    end
end

if ( bitand(nv, nv-1) == 0 && bitand(nu, nu-1) == 0 )
    ordering = morton(nv, nu);
    ordering = ordering(:);
    x(ordering) = x;
    y(ordering) = y;
    z(ordering) = z;
end

doms = {surfacemesh(x, y, z)};

[rnodes, rwts] = chebpts(nr,[0 1],1);
fluxnodes = {zeros(3,nt*nr)};
fluxwts = {zeros(1,nt*nr)};
for i = 1:nr
    rr = rnodes(i);
    wr = rwts(i);
    for j = 1:nt
        ij = (i-1)*nt+j;
        tt = 2*pi*(j-1)/nt;
        [gi1, ~, gi2] = evalRotEllipse(pi, tt);
        [dgi1, ~, dgi2] = dvEvalRotEllipse(pi, tt);
        fluxnodes{1}(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[xscenter; 0; 0];
        fluxwts{1}(1,ij) = (2*pi/nt) ...
            * wr*((xscenter-gi1)*(1-rr)*dgi2 + gi2*(1-rr)*dgi1);
    end
end

end

function [x, y, z] = evalRotEllipse(u, v)

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));
for i = -1:2
    for j = -1:2
        ph = (1-i)*v + j*u;
        x = x + d(i+2,j+2)*cos(u).*cos(ph);
        y = y + d(i+2,j+2)*sin(u).*cos(ph);
        z = z + d(i+2,j+2)*sin(ph);
    end
end

end

function [x, y, z] = dvEvalRotEllipse(u, v)

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));
for i = -1:2
    for j = -1:2
        ph = (1-i)*v + j*u;
        x = x - (1-i)*d(i+2,j+2)*cos(u).*sin(ph);
        y = y - (1-i)*d(i+2,j+2)*sin(u).*sin(ph);
        z = z + (1-i)*d(i+2,j+2)*cos(ph);
    end
end

end
