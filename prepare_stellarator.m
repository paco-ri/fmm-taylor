function [doms, fluxnodes, fluxwts] = prepare_stellarator(n, nu, nv, varargin)

if nargin == 5
    scales = 1;
    ns = n;
    nus = nu;
    nvs = nv;
    nr = varargin{1};
    nt = varargin{2};
    nsurf = 1;
elseif nargin == 11
    ns = [n varargin{1}];
    nus = [nu varargin{2}];
    nvs = [nv varargin{3}];
    scales = sort([varargin{4} varargin{5}]);
    nr = varargin{6};
    nt = varargin{7};
    np = varargin{8};
    nsurf = 2;
else
    disp('invalid number of args (4 or 8)')
end
    
doms = cell(1,nsurf);
for i = 1:nsurf
    n = ns(i);
    nu = nus(i);
    nv = nvs(i);
    rr = scales(i);
    
    x = cell(nu*nv, 1);
    y = cell(nu*nv, 1);
    z = cell(nu*nv, 1);
    ubreaks = linspace(0, 2*pi, nu+1);
    vbreaks = linspace(0, 2*pi, nv+1);
    
    k = 1;
    for ku = 1:nu
        for kv = 1:nv
            [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
            [x{k}, y{k}, z{k}] = evalStellarator(uu, vv, rr);
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
    
    doms{i} = surfacemesh(x, y, z);

end

fluxnodes = cell(1,nsurf);
fluxwts = cell(1,nsurf);

[rnodes, rwts] = chebpts(nr,[0 1],1);
fluxnodes{1} = zeros(3,nt*nr);
fluxwts{1} = zeros(1,nt*nr);
if nsurf == 1
    scales = [0 1];
end
for i = 1:nr
    rr = rnodes(i);
    wr = rwts(i);
    for j = 1:nt
        ij = (i-1)*nt+j;
        tt = 2*pi*(j-1)/nt;
        [go1, ~, go2] = evalStellarator(0,tt,scales(2));
        [gi1, ~, gi2] = evalStellarator(0,tt,scales(1));
        [dgo1, ~, dgo2] = dvEvalStellarator(0,tt,scales(2));
        [dgi1, ~, dgi2] = dvEvalStellarator(0,tt,scales(1));
        fluxnodes{1}(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[go1; 0; go2];
        fluxwts{1}(1,ij) = (2*pi/nt) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

if nsurf == 2
    fluxnodes{2} = zeros([3 nr*np]);
    fluxwts{2} = zeros([1 nr*np]);
    for i = 1:nr
        rr = rnodes(i);
        wr = rwts(i);
        for j = 1:np
            ij = (i-1)*np+j;
            pp = 2*pi*(j-1)/np;
            [go1, go2] = evalStellarator(pp,0,scales(2));
            [gi1, gi2] = evalStellarator(pp,0,scales(1));
            [dgo1, dgo2] = duEvalStellarator(pp,0,scales(2));
            [dgi1, dgi2] = duEvalStellarator(pp,0,scales(1));
            fluxnodes{2}(:,ij) = (1-rr)*[gi1; gi2; 0] + rr*[go1; go2; 0];
            fluxwts{2}(1,ij) = (2*pi/np) ...
                * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
                - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
        end
    end
end

end

function [x, y, z] = evalStellarator(u, v, scale)

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

rz = zeros(size(v));
for j = -1:4
    for k = -1:3
        jj = j+2;
        kk = k+2;
        rz = rz + d(jj,kk)*exp(-1i*j*v+1i*k*Q*u);
    end
end
rz = exp(1i*v).*rz;
r1 = real(rz);
z1 = imag(rz);
r = r0 + R*(r1-r0);
z = z0 + R*(z1-z0);

[x, y] = pol2cart(u, r);

end

function [x, y, z] = dvEvalStellarator(u, v, scale)

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

rz = zeros(size(v));
for j = -1:4
    for k = -1:3
        jj = j+2;
        kk = k+2;
        rz = rz + 1i*(1-j)*d(jj,kk)*exp(1i*(1-j)*v+1i*k*Q*u);
    end
end
r1 = real(rz);
z1 = imag(rz);
r = r0 + R*(r1-r0);
z = z0 + R*(z1-z0);

[x, y] = pol2cart(u, r);

end

function [x, y] = dthpol2cart(th, r)
[xx, yy] = pol2cart(th, r);
x = -yy;
y = xx;
end

function [x, y] = drpol2cart(th, r)
[xx, yy] = pol2cart(th, r);
x = xx/r;
y = yy/r;
end

function [x, y, z] = duEvalStellarator(u, v, scale)

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

rz = zeros(size(v));
durz = zeros(size(v));
for j = -1:4
    for k = -1:3
        jj = j+2;
        kk = k+2;
        rz = rz + d(jj,kk)*exp(-1i*j*v+1i*k*Q*u);
        durz = durz + 1i*k*Q*d(jj,kk)*exp(-1i*j*v+1i*k*Q*u);
    end
end
rz = exp(1i*v).*rz;
r1 = real(rz);
% z1 = imag(rz);
r = r0 + R*(r1-r0);
durz = exp(1i*v).*durz;
dur1 = real(durz);
duz1 = imag(durz);
% z = z0 + R*(z1-z0);
z = R*duz1;

% [x, y] = pol2cart(u, r);
[x, y] = dthpol2cart(u, r);
[x1, y1] = drpol2cart(u, r);
x = x + x1*dur1;
y = y + y1*dur1;

end