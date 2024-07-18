ns = 5; % 4:2:10;
nus = 4:2:8; % 4:4:16;
lerr = zeros([15 size(ns,2)*size(nus,2)]);
lind = 1;

% wavenumber 
zk = 1.0 + 0.0i; 
lambda = real(zk); 

whichgeom = 1; % circular torus

% domain
domrmin = 1.0;
if whichgeom == 1
    domrmaj = 2.0;
else
    domrmaj = 4.5; % for surfacemesh.torus
end

% interior points at which B is computed for convergence analysis
nintphi = 32;
ninttheta = 32;
nintr = 64;

% quadrature options
eps = 1e-7;
opts_quad = [];
opts_quad.format='rsc';

for n = ns
for nu = nus

% define surface 
nv = nu*3;
fprintf('n = %d, nu = %d, nv = %d\n',n,nu,nv)

if whichgeom == 1
    dom = circulartorus(n,nu,nv,domrmin,domrmaj);
else
    dom = surfacemesh.torus(n, nu, nv);
end

vn = normal(dom);

% get harmonic surface vector field 
sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
phihat = surfacefunv(@(x,y,z) -sinphi(x,y,z), ...
                     @(x,y,z) cosphi(x,y,z), ...
                     @(x,y,z) 0.*z, dom);
tauhat = cross(vn, phihat); 
[u, v, vH, curlfree, divfree] = hodge(tauhat);
mH = vH + 1i.*cross(vn,vH);

% harmonic surface vector field in axisymmetric case -- no L-B solve
sintheta = @(x,y,z) z./domrmin; 
costheta = @(x,y,z) (sqrt(x.^2 + y.^2) - domrmaj)./domrmin;
% theta routines are wrong for surfacemesh.torus
tauhat2 = surfacefunv(@(x,y,z) -sintheta(x,y,z).*cosphi(x,y,z), ...
                     @(x,y,z) -sintheta(x,y,z).*sinphi(x,y,z), ...
                     @(x,y,z) costheta(x,y,z), dom);
overr = surfacefun(@(x,y,z) 1./sqrt(x.^2 + y.^2), dom);
mH2 = tauhat2.*overr - 1i.*phihat.*overr;
vn2 = surfacefunv(@(x,y,z) costheta(x,y,z).*cosphi(x,y,z), ...
                 @(x,y,z) costheta(x,y,z).*sinphi(x,y,z), ...
                 @(x,y,z) sintheta(x,y,z), dom); 
mH = mH2;
vn = vn2;

% compute reference Taylor state
if whichgeom == 1
    rmaj = 2.0; % dist. of center of current ring from origin
    rmin = 2.0; % radius of current ring 
else
    rmaj = 4.5; % for surfacemesh.torus
    rmin = 3.0; % for surfacemesh.torus
end
ntheta = 1e3; % number of disc. points
jmag = 1.0; % current magnitude 
B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,lambda);
nB0 = dot(vn,B0);

% convert surfacemesh dom to surfer
S = surfer.surfacemesh_to_surfer(dom);

% compute near quadrature correction for taylor.static/dynamic routines
t1 = tic;
if zk == 0
    Q = taylor.static.get_quadrature_correction(S,eps,S,opts_quad);
else
    Q = taylor.dynamic.get_quadrature_correction(S,zk,eps,S,opts_quad);
end
t2 = toc(t1);
fprintf('on-surface quadrature: %f s\n', t2)

% do GMRES to solve A11*W = nB0
b = surfacefun_to_array(nB0,dom,S); 
t1 = tic;
[W, flag1, relres1, iter1] = gmres(@(s) A(s,dom,S,zk,eps,Q),b,[],eps,50);
t2 = toc(t1);
fprintf('GMRES for A11*W = n.B0: %f s / %d iter. = %f s\n', ...
    t2, iter1(2), t2/iter1(2))
wfunc = array_to_surfacefun(W,dom,S);

% do GMRES to solve A11*D = A12
Balpha = mtxBalpha(S,dom,mH,zk,eps,S,Q);
b = surfacefun_to_array(Balpha,dom,S);
t1 = tic;
[D, flag2, relres2, iter2] = gmres(@(s) A(s,dom,S,zk,eps,Q),b,[],eps,50);
t2 = toc(t1);
fprintf('GMRES for A11*D = A12: %f s / %d iter. = %f s\n', ...
    t2, iter2(2), t2/iter2(2))
dfunc = array_to_surfacefun(D,dom,S);

% === GMRES CHECK ===

bcheck = mtxBsigma(S,dom,dfunc,zk,eps,S,Q);
bcheck = surfacefun_to_array(bcheck,dom,S);   
berr = bcheck - b;

% ===================

% post-GMRES processing

% using *nontaylor routines
% get cross-section quadrature
nr = 80; % # of nodes in r-dir. (80 for circ)
nt = 150; % # of nodes in theta-dir. (150 for circ)
if whichgeom == 1
    shift = [domrmaj; 0; 0];
    [rr,wr] = legpts(nr,[0 domrmin]);
    tt = 2*pi.*(0:nt-1)./nt;
    wt = 2*pi/nt;
    
    % simultaneously compute B0 flux
    qnodes = zeros([3 nr*nt]);
    qweights = zeros([1 nr*nt]);
    flux = 0;
    for i = 1:nr
        for j = 1:nt
            ii = nt*(i-1) + j;
            targpt = [rr(i)*cos(tt(j)); 0; rr(i)*sin(tt(j))];
            targpt = targpt + shift;
            qnodes(:,ii) = targpt;
            qweights(ii) = rr(i)*wr(i)*wt;
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targpt);
            flux = flux + B0eval(2)*qweights(ii);
        end
    end
else
    [qnodes, qweights] = torusfluxquad(nr, nt);
    flux = 0;
    for i = 1:nr*nt
        B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,qnodes(:,i));
        flux = flux + B0eval(2)*qweights(i);
    end
end
fprintf('B0 flux = %f\n',flux)

targinfoflux = [];
targinfoflux.r = qnodes;
t1 = tic;
if zk == 0
    Qflux = taylor.static.get_quadrature_correction(S,eps,targinfoflux,opts_quad);
else
    Qflux = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfoflux,opts_quad);
end
t2 = toc(t1);
fprintf('flux quadrature: %f s\n',t2)

t1 = tic;
fluxsigmaD = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,dfunc,zk,eps,Qflux);
fluxsigmaW = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,wfunc,zk,eps,Qflux);
fluxalpha = mtxfluxalphanontaylor(S,dom,qnodes,qweights,mH,zk,eps,Qflux);

% compute alpha and sigma
alpha = -1i*(flux - fluxsigmaW)/(-fluxsigmaD + fluxalpha);
sigma = 1i*alpha.*dfunc - wfunc;
t2 = toc(t1);
fprintf('alpha and sigma: %f s\n', t2)

m0 = debyem0(sigma,zk);
m0err = div(m0) - 1i*zk.*sigma;
fprintf('m0err = %f\n', norm(m0err))

% compute B on-surface
m = m0 + alpha.*mH; 
nxm = cross(vn,m);

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';
if zk == 0
    gradSsigma = taylor.static.eval_gradS0(S,sigmavals,eps,S,Q);
else
    gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,eps,S,Q);
end
mvals = surfacefun_to_array(m,dom,S);
mvals = mvals.';
if zk == 0
    curlSm = taylor.static.eval_curlS0(S,mvals,eps,S,Q);
else
    curlSm = taylor.dynamic.eval_curlSk(S,zk,mvals,eps,S,Q);
end

vnvals = surfacefun_to_array(vn,dom,S);
vnvals = vnvals.';
nxmvals = surfacefun_to_array(nxm,dom,S);
nxmvals = nxmvals.';
B = -vnvals.*sigmavals./2 + 1i.*nxmvals./2 - gradSsigma + 1i.*curlSm;
if zk ~= 0
    Qhelm = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[1.0 0],S);
    opts_helm = [];
    opts_helm.precomp_quadrature = Qhelm;
    opts_helm.format = 'rsc';
    Sm = zeros(size(mvals));
    for j = 1:3
        Sm(j,:) = helm3d.dirichlet.eval(S,mvals(j,:),S,eps,zk,[1.0 0],opts_helm);
    end
    B = B + 1i*zk.*Sm;
end
B = array_to_surfacefun(B.',dom,S);

% compute B throughout the interior
[interior, interiorwts] = interiorcirctorus(nintphi,ninttheta,nintr, ...
    domrmin,domrmaj);
targinfoint = [];
targinfoint.r = interior;

t1 = tic;
if zk == 0
    Qint = taylor.static.get_quadrature_correction(S,eps,targinfoint,opts_quad);
else
    Qint = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfoint,opts_quad);
end
t2 = toc(t1);
fprintf('interior quadrature: %f s\n', t2)

t1 = tic;
if zk == 0
    gradSsigmaint = taylor.static.eval_gradS0(S,sigmavals,eps,targinfoint,Qint);
    curlSmint = taylor.static.eval_curlS0(S,mvals,eps,targinfoint,Qint);
else
    gradSsigmaint = taylor.dynamic.eval_gradSk(S,zk,sigmavals,eps,targinfoint,Qint);
    curlSmint = taylor.dynamic.eval_curlSk(S,zk,mvals,eps,targinfoint,Qint);
end
Bint = -gradSsigmaint + 1i.*curlSmint;
if zk ~= 0
    Qhelmint = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[1.0 0],targinfoint);
    opts_helmint = [];
    opts_helmint.precomp_quadrature = Qhelmint;
    opts_helmint.format = 'rsc';
    Smint = zeros(size(interior));
    for j = 1:3
        Smint(j,:) = helm3d.dirichlet.eval(S,mvals(j,:),targinfoint, ...
            eps,zk,[1.0 0],opts_helmint);
    end
    Bint = Bint + 1i*zk.*Smint;
end
t2 = toc(t1); 
fprintf('interior B: %f s\n', t2)

B0int = zeros(size(interior));
for j = 1:size(interior,2)
    B0int(:,j) = reftaylor(ntheta,rmin,rmaj,jmag,lambda,interior(:,j));
end

% n.B/B0 plots
figure(1)
% subplot(1,3,1)
plot(dot(vn,B0-B))
colorbar

% subplot(1,3,2)
% plot(dot(vn,B0))
% colorbar
% 
% subplot(1,3,3)
% plot(dot(vn,B))
% colorbar

figure(2)
% subplot(1,3,1)
plot(norm(B0-B))
colorbar
% subplot(1,3,2)
% plot(norm(B0))
% colorbar
% subplot(1,3,3)
% plot(norm(B))
% colorbar

lerr(1,lind) = n;
lerr(2,lind) = nu;
lerr(3,lind) = nv;
lerr(4,lind) = max(abs(Bint-B0int),[],'all')/max(abs(B0int),[],'all');
lerr(5,lind) = max(real(Bint-B0int),[],'all')/max(real(B0int),[],'all');
lerr(6,lind) = dot(dot(Bint-B0int,Bint-B0int),interiorwts) ...
    /dot(dot(B0int,B0int),interiorwts);
lerr(7,lind) = dot(dot(real(Bint-B0int),real(Bint-B0int)),interiorwts) ...
    /dot(dot(real(B0int),real(B0int)),interiorwts);
lerr(8,lind) = dot(sum(abs(Bint-B0int),1),interiorwts) ...
    /dot(sum(abs(B0int),1),interiorwts);
lerr(9,lind) = dot(sum(abs(real(Bint-B0int)),1),interiorwts) ...
    /dot(sum(abs(real(B0int)),1),interiorwts);
lerr(10,lind) = norm(dot(vn,B0-B),inf)/norm(dot(vn,B0),inf);
lerr(11,lind) = norm(dot(vn,B0-B),2)/norm(dot(vn,B0),2);
lerr(12,lind) = norm(dot(vn,B0-B),1)/norm(dot(vn,B0),1);
% lerr(13,lind) = norm(norm(B0-B),inf)/norm(norm(B0),inf);
numer = max([norm(B0.components{1}-B.components{1}, inf) ...
    norm(B0.components{2}-B.components{2}, inf), ...
    norm(B0.components{3}-B.components{3}, inf)]);
denom = max([norm(B0.components{1}, inf) ...
    norm(B0.components{2}, inf), ...
    norm(B0.components{3}, inf)]);
lerr(13,lind) = numer/denom;
lerr(14,lind) = norm(norm(B0-B),2)/norm(norm(B0),2);
lerr(15,lind) = norm(norm(B0-B),1)/norm(norm(B0),1);
lind = lind+1;

fprintf('==============\n')

end
eps = eps*1e-2;
end

% B flux
nr = 41; % # of nodes in r-dir.
nt = 301; % 23; % # of nodes in theta-dir.
if whichgeom == 1
    shift = [domrmaj; 0; 0];
    [rr,wr] = legpts(nr,[0 domrmin]);
    tt = 2*pi.*(0:nt-1)./nt;
    wt = 2*pi/nt;
    
    qnodes = zeros([3 nr*nt]);
    qweights = zeros([1 nr*nt]);
    for i = 1:nr
        for j = 1:nt
            ii = nt*(i-1) + j;
            targpt = [rr(i)*cos(tt(j)); 0; rr(i)*sin(tt(j))];
            targpt = targpt + shift;
            qnodes(:,ii) = targpt;
            qweights(ii) = rr(i)*wr(i)*wt;
        end
    end
else
    [qnodes, qweights] = torusfluxquad(nr, nt);
end

targinfoflux = [];
targinfoflux.r = qnodes;
t1 = tic;
if zk == 0
    Qflux = taylor.static.get_quadrature_correction(S,eps,targinfoflux,opts_quad);
else
    Qflux = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfoflux,opts_quad);
end
t2 = toc(t1);
fprintf('flux quadrature: %f s\n',t2)

Bfluxsigma = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,sigma,zk,eps,Qflux);
BfluxmH = mtxfluxalphanontaylor(S,dom,qnodes,qweights,mH,zk,eps,Qflux);
Bflux = -Bfluxsigma + 1i*alpha*BfluxmH;
fprintf('B flux = %f\n', Bflux)
fprintf('flux difference = %f\n', Bflux - flux)

% function in solve for D
function A = A(s,dom,S,zk,eps,Q)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = mtxBsigma(S,dom,sigma,zk,eps,S,Q);
A = surfacefun_to_array(Bsigma,dom,S);
end

