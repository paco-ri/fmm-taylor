% --- Tolerances ---
tol = 1e-4;

% --- Beltrami parameter ---
zk = 1.0;

% --- Geometry parameters ---
ns = 5;%[5 7 9]; % polynomial order + 1
nvs = [6 8];% 10]; % number of patches in poloidal direction

lerr = zeros(4,size(ns,2)*size(nvs,2));
lind = 1;

for n = ns
    for nv = nvs
        fprintf('\t========\n\tn = %d, nv = %d \n', n, nv)
        nu = nv*3; % number of patches in toroidal direction
        % a = 2.0; % minor radius, horiz. axis
        % a0 = 5.0; % major radius
        % b = 3.0; % minor radius vert. axis
        % dom = twisted_ellipse_torus(a,a0,b,n,nu,nv); % surfacemesh
        % domparams = [n, nu, nv];
        dom = alt_stellarator(n,nu,nv);
        domparams = [n,nu,nv];

        % --- Compute B0 ---
        ntheta = 1e3;
        rmin = 5.0;
        rmaj = 5.0;
        jmag = 1.0;
        B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,zk);

        % --- Compute XS quad. and flux ---
        % Specific for twisted ellipse geometry
        nr = 12;
        nt = 100;
        [qnodes, qweights] = ellipsefluxquad(nr,nt,a,a0,b);
        flux = 0;
        for i = 1:nr*nt
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes(:,i));
            flux = flux + B0eval(2)*qweights(i);
        end
        
        B0 = {B0};
        qnodes = {qnodes};
        qweights = {qweights};
        ts = RefTaylorState(dom,domparams,zk,flux,B0,qnodes,qweights,tol);
        ts = ts.solve(true);
        B = ts.surface_B();

        lerr(1,lind) = n; % number of points on each patch
        lerr(2,lind) = nv; % geom. param. 
        lerr(3,lind) = nu; % geom. param. 
        lerr(4,lind) = vecinfnorm(B0{1}-B{1})/vecinfnorm(B0{1}); 

        lind = lind + 1;
    end
    tol = tol*1e-2;
end

loglog(sqrt(lerr(1,1:size(nvs,2)).^2.*lerr(2,1:size(nvs,2)).* ...
    lerr(3,1:size(nvs,2))), lerr(4,1:size(nvs,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nvs,2)+1:(i+1)*size(nvs,2);
    loglog(sqrt(lerr(1,ind).^2.*lerr(2,ind).* ...
        lerr(3,ind)), lerr(4,ind), 'o-')
end

h1 = 50:90;
h2 = 80:120;
h3 = 100:150;
loglog(h1,1e1.*h1.^(-1),'--')
loglog(h2,1e4.*h2.^(-3),'--')
loglog(h3,1e7.*h3.^(-5),'--')

legend('p=4','p=6','p=8','O(h^{-1})','O(h^{-3})','O(h^{-5})')

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf) ...
        norm(f.components{2}, inf), ...
        norm(f.components{3}, inf)]);
end
