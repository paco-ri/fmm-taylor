%% testing the effect of the tolerance

% Notes:
% - for n = 9, nv = 7, ideal tolerance is around 1e-5
% - as the tolerance decreases, the number of iterations for the Aw = B0.n 
%   solve increases beyond the expected value of .5*(number of iterations
%   for the A*d = Phi solve)

n = 9;
nv = 7;
nu = nv*3;
zk = 0;

tols = [1e-4 1e-5 1e-6 1e-7];
ntol = size(tols,2);

err = zeros(4,ntol);
ind = 1;

tss = cell(1,ntol);

for tol = tols
    fprintf('----------\nn=9, nv=7, tol=%f \n----------\n',tol)

    ao = 1.0; % outer minor radius
    ai = 0.6; % inner minor radius

    nr = 16;
    nt = 40;
    np = 40;
    [dom, qnodes, qweights] = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);

    ntheta = 1e3;
    rmin = 2.0;
    rmaj = 2.0;
    jmag = 1.0;
    B0o = reftaylorsurffun(dom{1},n,nu,nv,ntheta,rmin,rmaj,jmag,zk);
    B0i = reftaylorsurffun(dom{2},n,nu,nv,ntheta,rmin,rmaj,jmag,zk);
        
    flux = zeros(1,2);
    for i = 1:nr*nt
        B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes{1}(:,i));
        flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
    end
    for i = 1:nr*np
        [B0eval, B0srcpts] = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes{2}(:,i));
        flux(2) = flux(2) - B0eval(3)*qweights{2}(i);
    end
        
    B0 = {B0o,B0i};
    domparams = [n nu nv];
    ts = RefTaylorState(dom,domparams,zk,flux,B0,qnodes,qweights,tol);
    ts = ts.solve(true);
    B = ts.surface_B();

    err(1,ind) = n;
    err(2,ind) = nv;
    err(3,ind) = nu;
    err(4,ind) = max(vecinfnorm(B0{1}-B{1}),vecinfnorm(B0{2}-B{2}))...
            /max(vecinfnorm(B0{1}),vecinfnorm(B0{2}));

    tss{ind} = ts;
    ind = ind + 1;
end