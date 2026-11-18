function [doms, fluxnodes, fluxwts] = prepare_torus(n, nu, nv, varargin)

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
    scales = sort([varargin{4} varargin{5}], 'descend');
    nr = varargin{6};
    nt = varargin{7};
    np = varargin{8};
    nsurf = 2;
else
    disp('invalid number of args (4 or 8)')
end

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];
d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25  -.20 0];
d(3,2)=2.;
d(4,3)=0;%-.2; % critical mode
    
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
            [x{k}, y{k}, z{k}] = evalTorus(uu, vv, rr, d);
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
        [go1, ~, go2] = evalTorus(0,tt,scales(1),d);
        [gi1, ~, gi2] = evalTorus(0,tt,scales(2),d);
        [dgo1, ~, dgo2] = dvEvalTorus(0,tt,scales(1),d);
        [dgi1, ~, dgi2] = dvEvalTorus(0,tt,scales(2),d);
        fluxnodes{1}(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[go1; 0; go2];
        fluxwts{1}(1,ij) = (2*pi/nt) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

if nsurf == 2
    % Cross-section 2 is the annulus in the plane z = 0 between the tori's
    % outboard z = 0 crossings. Those lie near v = 0 but not on it: the v = 0
    % curve sits at z = -0.16*scale*sin(u).
    pp = 2*pi*(0:np-1)/np;
    Rx = zeros(2, np);              % row 1 outer, row 2 inner
    for is = 1:2
        for j = 1:np
            v = 0;
            for it = 1:50
                [~, ~, zz] = evalTorus(pp(j),v,scales(is),d);
                [~, ~, dzz] = dvEvalTorus(pp(j),v,scales(is),d);
                v = v - zz/dzz;
            end
            [xx, yy, zz] = evalTorus(pp(j),v,scales(is),d);
            assert(abs(zz) < 1e-13, 'prepare_torus: z = 0 crossing not found');
            Rx(is,j) = hypot(xx, yy);
        end
    end
    fluxnodes{2} = zeros([3 nr*np]);
    fluxwts{2} = zeros([1 nr*np]);
    for i = 1:nr
        for j = 1:np
            ij = (i-1)*np+j;
            R = (1-rnodes(i))*Rx(2,j) + rnodes(i)*Rx(1,j);
            fluxnodes{2}(:,ij) = [R*cos(pp(j)); R*sin(pp(j)); 0];
            fluxwts{2}(1,ij) = (2*pi/np)*rwts(i)*R*(Rx(1,j) - Rx(2,j));
        end
    end
end

end

function [x, y, z] = evalTorus(u, v, scale, d)

d = scale.*d;
d(3,2) = 2.0;

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));

for i = -1:2
    for j = -1:2
        ii = i+2;
        jj = j+2;
        x = x + d(ii,jj)*cos(u).*cos((1-i)*v+j*u);
        y = y + d(ii,jj)*sin(u).*cos((1-i)*v+j*u);
        z = z + d(ii,jj)*sin((1-i)*v+j*u);
    end
end
 
end

function [x, y, z] = dvEvalTorus(u, v, scale, d)

d = scale.*d;
d(3,2) = 2.0;

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));

for i = -1:2
    for j = -1:2
        ii = i+2;
        jj = j+2;
        x = x - d(ii,jj)*(1-i)*cos(u).*sin((1-i)*v+j*u);
        y = y - d(ii,jj)*(1-i)*sin(u).*sin((1-i)*v+j*u);
        z = z + d(ii,jj)*(1-i)*cos((1-i)*v+j*u);
    end
end

end