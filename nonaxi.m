ns = 6:2:10;
nus = 4:2:12; %4:4:16;
lerr = zeros([15 size(ns,2)*size(nus,2)]);
lind = 1;

whichgeom = 1; % circular torus

% domain
domrmin = 1.0;
if whichgeom == 1
    domrmaj = 2.0;
else
    domrmaj = 4.5; % for surfacemesh.torus
end

% interior points at which B is computed for convergence analysis
nintphis = 30;
intphis = 2*pi.*(0:nintphis-1)./nintphis;
nintthetas = 30;
intthetas = 2*pi.*(0:nintthetas-1)./nintthetas;
nintrs = 16;
[interior, interiorwts] = interiorcirctorus(nintphis,nintthetas,nintrs, ...
    domrmin,domrmaj);
targinfoint = [];
targinfoint.r = interior;

% quadrature options
eps = 1e-6;
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
% wrong for surfacemesh.torus
% sintheta = @(x,y,z) z./domrmin; 
% costheta = @(x,y,z) (sqrt(x.^2 + y.^2) - domrmaj)./domrmin;
phihat = surfacefunv(@(x,y,z) -sinphi(x,y,z), ...
                     @(x,y,z) cosphi(x,y,z), ...
                     @(x,y,z) 0.*z, dom);
tauhat = cross(vn, phihat); 
[u, v, vH, curlfree, divfree] = hodge(tauhat);
mH = vH + 1i.*cross(vn,vH);

% harmonic surface vector field in axisymmetric case -- no L-B solve
% tauhat2 = surfacefunv(@(x,y,z) -sintheta(x,y,z).*cosphi(x,y,z), ...
%                      @(x,y,z) -sintheta(x,y,z).*sinphi(x,y,z), ...
%                      @(x,y,z) costheta(x,y,z), dom);
% overr = surfacefun(@(x,y,z) 1./sqrt(x.^2 + y.^2), dom);
% mH2 = tauhat2.*overr - 1i.*phihat.*overr;
% vn2 = surfacefunv(@(x,y,z) costheta(x,y,z).*cosphi(x,y,z), ...
%                  @(x,y,z) costheta(x,y,z).*sinphi(x,y,z), ...
%                  @(x,y,z) sintheta(x,y,z), dom); 
% mH = mH2;
% vn = vn2;

% test properties of mH
% fprintf('surface div of mH = %f\n',integral2(div(mH))/surfacearea(dom));
% fprintf('surface div of n x mH = %f\n',integral2(div(cross(vn,mH)))/surfacearea(dom));
% fprintf('n . mH = %f\n',integral2(dot(vn,mH))/surfacearea(dom));
% fprintf('n x mH + i mH = %f\n',integral2(norm(cross(vn,mH)+1i.*mH))/surfacearea(dom));
% 
% fprintf('surface div of mH2 = %f\n',integral2(div(mH2))/surfacearea(dom));
% fprintf('surface div of n2 x mH2 = %f\n',integral2(div(cross(vn2,mH2)))/surfacearea(dom));
% fprintf('n2 . mH2 = %f\n',integral2(dot(vn2,mH2))/surfacearea(dom));
% fprintf('n2 x mH2 + i mH2 = %f\n',integral2(norm(cross(vn2,mH2)+1i.*mH2))/surfacearea(dom));

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
lambda = 0; 
B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,lambda);
nB0 = dot(vn,B0);

% convert surfacemesh dom to surfer
S = surfer.surfacemesh_to_surfer(dom);

% compute near quadrature correction for taylor.static routines
t1 = tic;
Q = taylor.static.get_quadrature_correction(S,eps,S,opts_quad);
t2 = toc(t1);
fprintf('on-surface quadrature: %f s\n', t2)

t1 = tic;
Qint = taylor.static.get_quadrature_correction(S,eps,targinfoint,opts_quad);
t2 = toc(t1);
fprintf('interior quadrature: %f s\n', t2)

% do GMRES to solve A11*W = nB0
b = surfacefun_to_array(nB0,dom,S); 
t1 = tic;
[W, flag1, relres1, iter1] = gmres(@(s) A(s,dom,S,eps,Q),b,[],eps,50);
t2 = toc(t1);
fprintf('GMRES for A11*W = n.B0: %f s / %d iter. = %f s\n', ...
    t2, iter1(2), t2/iter1(2))
wfunc = array_to_surfacefun(W,dom,S);

% do GMRES to solve A11*D = A12
Balphar = mtxBalpha(S,dom,real(mH),eps,S,Q);
Balphai = mtxBalpha(S,dom,imag(mH),eps,S,Q);
br = surfacefun_to_array(Balphar,dom,S);
bi = surfacefun_to_array(Balphai,dom,S);
t1 = tic;
[Dr, flag2, relres2, iter2] = gmres(@(s) A(s,dom,S,eps,Q),br,[],eps,50);
[Di, flag3, relres3, iter3] = gmres(@(s) A(s,dom,S,eps,Q),bi,[],eps,50);
t2 = toc(t1);
fprintf('GMRES for A11*D = A12: %f s / %d iter. = %f s\n', ...
    t2, iter2(2)+iter3(2), t2/(iter2(2)+iter3(2)))
D = Dr + 1i*Di;
dfunc = array_to_surfacefun(D,dom,S);

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
opts_quad = [];
opts_quad.format = 'rsc';
t1 = tic;
Qflux = taylor.static.get_quadrature_correction(S,eps,targinfoflux,opts_quad);
t2 = toc(t1);
fprintf('flux quadrature: %f s\n',t2)

t1 = tic;
fluxsigmaDr = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,real(dfunc),eps,Qflux);
fluxsigmaDi = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,imag(dfunc),eps,Qflux);
fluxsigmaD = fluxsigmaDr + 1i*fluxsigmaDi;
fluxsigmaW = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,wfunc,eps,Qflux);
fluxalphar = mtxfluxalphanontaylor(S,dom,qnodes,qweights,real(mH),eps,Qflux);
fluxalphai = mtxfluxalphanontaylor(S,dom,qnodes,qweights,imag(mH),eps,Qflux);
fluxalpha = fluxalphar + 1i*fluxalphai;

% compute alpha and sigma
alpha = -1i*(flux - fluxsigmaW)/(-fluxsigmaD + fluxalpha);
sigma = 1i*alpha.*dfunc - wfunc;

% compute B
m = alpha.*mH; % will later involve call to debyem0
t2 = toc(t1);
fprintf('alpha and sigma: %f s\n', t2)

% interior point at which is curl is estimated with finite differences
% h = 1e-1;
% xx = 2.1*cos(5*pi/4);%pi/2-.1);
% yy = 2.1*sin(5*pi/4);%pi/2-.1);
% zz = .1;
% interior = [xx xx+h xx-h xx   xx   xx   xx;
%             yy yy   yy   yy+h yy-h yy   yy;
%             zz zz   zz   zz   zz   zz+h zz-h];

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';
gradS0sigmar = taylor.static.eval_gradS0(S,real(sigmavals),eps,S,Q);
gradS0sigmai = taylor.static.eval_gradS0(S,imag(sigmavals),eps,S,Q);
gradS0sigma = gradS0sigmar + 1i.*gradS0sigmai;
mvals = surfacefun_to_array(m,dom,S);
mvals = mvals.';
curlS0mr = taylor.static.eval_curlS0(S,real(mvals),eps,S,Q);
curlS0mi = taylor.static.eval_curlS0(S,imag(mvals),eps,S,Q);
curlS0m = curlS0mr + 1i*curlS0mi;

% for interior B
t1 = tic;
gradS0sigmarint = taylor.static.eval_gradS0(S,real(sigmavals),eps,targinfoint,Qint);
gradS0sigmaiint = taylor.static.eval_gradS0(S,imag(sigmavals),eps,targinfoint,Qint);
gradS0sigmaint = gradS0sigmarint + 1i.*gradS0sigmaiint;
curlS0mrint = taylor.static.eval_curlS0(S,real(mvals),eps,targinfoint,Qint);
curlS0miint = taylor.static.eval_curlS0(S,imag(mvals),eps,targinfoint,Qint);
curlS0mint = curlS0mrint + 1i*curlS0miint;
t2 = toc(t1); 
fprintf('interior B: %f s\n', t2)

vnvals = surfacefun_to_array(vn,dom,S);
vnvals = vnvals.';
B = -vnvals.*sigmavals./2 + mvals./2 - gradS0sigma + 1i.*curlS0m;
% if eval routs are average of interior and exterior limits
% B = vnvals.*sigmavals./2 - mvals./2 - gradS0sigma + 1i.*curlS0m;
% if eval routs already compute exterior limit
% B = -vnvals.*sigmavals + mvals - gradS0sigma + 1i.*curlS0m;
% no identity terms
% B = -gradS0sigma + 1i.*curlS0m;
gradS0sigmafun = array_to_surfacefun(gradS0sigma.',dom,S);
curlS0mfun = array_to_surfacefun(curlS0m.',dom,S);
B = array_to_surfacefun(B.',dom,S);
Bint = -gradS0sigmaint + 1i.*curlS0mint;
B0int = zeros(size(interior));
for j = 1:size(interior,2)
    B0int(:,j) = reftaylor(ntheta,rmin,rmaj,jmag,lambda,interior(:,j));
end

% n.B/B0 plots
% figure(1)
% subplot(1,3,1)
% plot(dot(vn,B0-B))
% colorbar
% 
% subplot(1,3,2)
% plot(dot(vn,B0))
% colorbar
% 
% subplot(1,3,3)
% plot(dot(vn,B))
% colorbar
% 
% figure(2)
% subplot(1,3,1)
% plot(norm(B0-B))
% colorbar
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
lerr(13,lind) = norm(B0-B,inf)/norm(B0,inf);
lerr(14,lind) = norm(norm(B0-B),2)/norm(norm(B0),2);
lerr(15,lind) = norm(norm(B0-B),1)/norm(norm(B0),1);
lind = lind+1;

fprintf('==============\n')

end
eps = eps*1e-2;
end

% estimate curl with finite differences
% DxBint = (Bint(3,5)-Bint(3,4)-Bint(2,7)+Bint(2,6))/(2*h);
% DyBint = (Bint(1,7)-Bint(1,6)-Bint(3,3)+Bint(3,2))/(2*h);
% DzBint = (Bint(2,3)-Bint(2,2)-Bint(1,5)+Bint(1,4))/(2*h);
% disp(DxBint)
% disp(DyBint)
% disp(DzBint)

% disp(norm(B0.components{1}-B.components{1},inf)./norm(B0.components{1},inf))
% disp(norm(B0.components{2}-B.components{2},inf)./norm(B0.components{2},inf))
% disp(norm(B0.components{3}-B.components{3},inf)./norm(B0.components{3},inf))

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
opts_quad = [];
opts_quad.format = 'rsc';
t1 = tic;
Qflux = taylor.static.get_quadrature_correction(S,eps,targinfoflux,opts_quad);
t2 = toc(t1);
fprintf('flux quadrature: %f s\n',t2)

Bfluxsigmar = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,real(sigma),eps,Qflux);
Bfluxsigmai = mtxfluxsigmanontaylor(S,dom,qnodes,qweights,imag(sigma),eps,Qflux);
Bfluxmr = mtxfluxalphanontaylor(S,dom,qnodes,qweights,real(m),eps,Qflux);
Bfluxmi = mtxfluxalphanontaylor(S,dom,qnodes,qweights,imag(m),eps,Qflux);
Bflux = -Bfluxsigmar - 1i.*Bfluxsigmai + 1i*(Bfluxmr + 1i.*Bfluxmi);
fprintf('B flux = %f\n', Bflux)
fprintf('flux difference = %f\n', Bflux - flux)

% function in solve for D
function A = A(s,dom,S,eps,Q)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = mtxBsigma(S,dom,sigma,eps,S,Q);
A = surfacefun_to_array(Bsigma,dom,S);
end

