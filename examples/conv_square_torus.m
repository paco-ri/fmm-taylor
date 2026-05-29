% --- Tolerances ---
tol = 1e-4;

% --- Beltrami parameter ---
zk = 0;

% --- Geometry parameters ---
ns = [5]; % 7 9]; % polynomial order + 1
nvs = 3 .* [6 8];% 10 12]; % number of patches in

lerr = zeros(4,size(ns,2)*size(nvs,2));
lind = 1;

for n = ns
    for nv = nvs
        fprintf('\t========\n\tn = %d, nv = %d \n', n, nv)
        nu = nv / 3; % idivide(nv, int32(3)); % number of patches in toroidal direction
        dom = surfacemesh.sharptorus(n, nu, nv);
        domparams = [n, nu, nv];

        % --- Compute B0 ---
        ntheta = 1e3;
        jmag = 1.0;
        rmin = 1.0;
        rmaj = 1.0;
        B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,zk);

        % --- Compute XS quad. and flux ---
        nt = 12;
        [qnodes, qweights] = square_flux_quad(nt);
        flux = 0;
        for i = 1:nt*nt
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes(:,i));
            flux = flux + B0eval(2)*qweights(i);
        end
        
        B0 = {B0};
        qnodes = {qnodes};
        qweights = {qweights};
        dom = {dom};
        ts = RefTaylorState(dom,domparams,zk,flux,B0,qnodes,qweights,tol);
        ts = ts.solve(true);
    t1 = tic;
        B = ts.surface_B();
    t2 = toc;
        fprintf('get surface B: %f s\n',t2)
        lerr(:,lind) = [n; nv; nu; norm(B{1}-B0{1},inf)/norm(B0{1},inf)];
        lind = lind + 1;
    end
end