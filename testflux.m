% 01/02/24 ref Taylor state flux test

% define surface 
n = 4; % polynomial order 
nu = 6; 
nv = nu*3;
% dom = surfacemesh.torus(n, nu, nv);
dom = circulartorus(n,nu,nv);

% compute reference Taylor state
rmaj = 5.0;
rmin = 3.0;
ntheta = 1e2;%1e3;
jmag = 1.0;
lambda = 1.0;

B0 = surfacefunv(dom);
B0x = cell(nu*nv,1);
B0y = cell(nu*nv,1);
B0z = cell(nu*nv,1);
for i = 1:nu*nv
    B0x{i} = zeros(n);
    B0y{i} = zeros(n);
    B0z{i} = zeros(n);
    for j = 1:n
        for k = 1:n
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,...
                [dom.x{i}(j,k) dom.y{i}(j,k) dom.z{i}(j,k)]);
            B0x{i}(j,k) = B0eval(1);
            B0y{i}(j,k) = B0eval(2);
            B0z{i}(j,k) = B0eval(3);
        end
    end
end
B0.components{1} = surfacefun(B0x,dom);
B0.components{2} = surfacefun(B0y,dom);
B0.components{3} = surfacefun(B0z,dom);

% integrate on a circle
nr = 100; % number of disc pts in radial dir
nt = 100; % in angular dir
surfint = 0;
domrmaj = 5.0;
domrmin = 2.0;
shift = [domrmaj 0 0];
B0evals = zeros([nr nt]);
rs = zeros([nr 1]);
ts = zeros([nt 1]);
for i = 1:nr
    r = (2*i+1)*domrmin/(2*nr);
    rs(i) = r;
    for j = 1:nt
        t = 2*pi*(2*j+1)/(2*nt);
        ts(j) = t;
        targpt = [r*cos(t) 0 r*sin(t)];
        B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targpt+shift);
        B0evals(i,j) = -B0eval(2);
        surfint = surfint - B0eval(2)*r*domrmin*2*pi/(nr*nt);
    end
end

disp(intacyc(B0,n,nu,nv));
disp(surfint);
