%%%% CONVQAS %%%%

% --- Tolerances ---
tol = 1e-4;

% --- Beltrami parameter ---
zk = 0;

% --- Geometry parameters ---
ns = [5 7 9]; % polynomial order + 1
nvs = [4 6 8]; % number of patches in poloidal direction
fid = fopen('qas.txt');
qasmodes = textscan(fid, '%f %f %f');
fclose(fid);
fid = fopen('qascoils.txt');
coilmodes = textscan(fid, '%f %f %f');
fclose(fid);
modes = {coilmodes, qasmodes};

% --- Reference Taylor state parameters ---
rref = 2.5;

lerr = zeros(4,size(ns,2)*size(nvs,2));
lind = 1;

for n = ns
    for nv = nvs
        fprintf('\t========\n\tn = %d, nv = %d \n', n, nv)
        nu = nv*3; % number of patches in toroidal direction

        % --- Get domain, compute XS quad. and flux ---
        nr = 16;
        nt = 40;
        np = 40;
        [dom, qnodes, qweights] = prepare2([n n],[nu nu],[nv nv],nr,nt,np,modes);
        domparams = [n, nu, nv];
        domo = dom{1};
        domi = dom{2};

        % --- Compute B0 ---
        ntheta = 1e3;
        jmag = 1.0;
        B0o = reftaylorsurffun(domo,n,nu,nv,ntheta,rref,rref,jmag,zk);
        B0i = reftaylorsurffun(domi,n,nu,nv,ntheta,rref,rref,jmag,zk);
        
        
        flux = zeros(1,2);
        for i = 1:nr*nt
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes{1}(:,i));
            flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
        end
        for i = 1:nr*np
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes{2}(:,i));
            flux(2) = flux(2) - B0eval(3)*qweights{2}(i);
        end
        
        B0 = {B0o,B0i};
        tols = [tol .01*tol .01*tol tol tol]; 
        ts = RefTaylorState(dom,domparams,zk,flux,B0,qnodes,qweights,tols);
        ts = ts.solve(true);
        B = ts.surface_B();

        lerr(1,lind) = n; % number of points on each patch
        lerr(2,lind) = nv; % geom. param. 
        lerr(3,lind) = nu; % geom. param. 
        lerr(4,lind) = max(vecinfnorm(B0{1}-B{1}),vecinfnorm(B0{2}-B{2}))...
            /max(vecinfnorm(B0{1}),vecinfnorm(B0{2})); 

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

h1 = 34:68;
h2 = 50:90;
h3 = 62:124;
loglog(h1,5e1.*h1.^(-2),'--')
loglog(h2,9e4.*h2.^(-4),'--')
loglog(h3,1e8.*h3.^(-6),'--')

legend('p=poly. order=4','p=6','p=8','O(h^{-2})','O(h^{-4})','O(h^{-6})','location','southwest')
xlabel('h, panel size')
ylabel('||B-B_0||_\infty/||B_0||_\infty')