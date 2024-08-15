% Running the code without doing error analysis

% ===== Tolerances and quadrature options =====
macheps = eps;
eps = 1e-7;
epsqc = eps;
epsqclh = eps;
epsgmres = eps;
optsqc = [];
optsqc.format = 'rsc';
optsqclh = [];
optsqclh.format = 'rsc';

% ===== Beltrami parameter =====
zk = 0;%complex(0.5); 

% ===== Set geometry =====

% --- Geometry parameters ---
n = 5; % polynomial order + 1
nv = 7; % number of patches in poloidal direction
nu = nv*3; % number of patches in toroidal direction
a = 2.0; % minor radius, horiz. axis
a0 = 5.0; % major radius
b = 3.0; % minor radius vert. axis

% --- surfacemesh and surfer construction ---
dom = twisted_ellipse_torus(a,a0,b,n,nu,nv); % surfacemesh
S = surfer.surfacemesh_to_surfer(dom); % surfer
vn = normal(dom); % normal vector

% --- Specify flux ---
flux = 1.0;

% ===== Solve for harmonic surface vector field(s) =====
mH = getmH(dom);

% ===== Get near quadrature corrections =====
Q = getqc(S,zk,epsqc,optsqc);
Qlh = getqclh(S,zk,epsqclh,optsqclh);

% ===== Do GMRES to solve A11*D = A12 =====
[D,dfunc] = solveforD(S,dom,mH,zk,epsgmres,Q);

% ===== Compute sigma and alpha =====
fluxsigmaD = mtxfluxsigma(S,dom,n,nu,nv,dfunc,zk,eps,S,Q,Qlh,optsqclh);
fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,mH,zk,eps,S,Q,Qlh,optsqclh);
alpha = flux/(1i*fluxsigmaD + fluxalpha);
sigma = 1i*alpha.*dfunc;

% ===== Compute B on-surface =====
m0 = debyem0(sigma,zk);
m = m0 + alpha.*mH;
B = surfaceB(sigma,m,S,dom,vn,zk,eps,Q,optsqc,optsqclh);

% ===== Plot =====
% plot(dot(vn,B))
% colorbar
disp(norm(dot(vn,B),inf))

% ============================
% ===== Helper functions =====
% ============================

% Get harmonic surface vector field
function mH = getmH(dom)
sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
phihat = surfacefunv(@(x,y,z) -sinphi(x,y,z), ...
                     @(x,y,z) cosphi(x,y,z), ...
                     @(x,y,z) 0.*z, dom);
vn = normal(dom);
tauhat = cross(vn, phihat); 
[~, ~, vH] = hodge(tauhat);
mH = vH + 1i.*cross(vn,vH);
end

% Get near quadrature correction for taylor routines
function Q = getqc(S,lambda,epsqc,optsqc)
if abs(lambda) < 2.5e-16
    Q = taylor.static.get_quadrature_correction(S,epsqc,S,optsqc);
else
    Q = taylor.dynamic.get_quadrature_correction(S,lambda,epsqc,S,optsqc);
end
end

% Get near quadrature correction for Laplace/Helmholtz potentials
function Qlh = getqclh(S,lambda,epsqclh,optsqclh)
if abs(lambda) < 2.5e-16
    Qlh = lap3d.dirichlet.get_quadrature_correction(S,epsqclh, ...
        [1.0,0.0],S,optsqclh);
else
    Qlh = helm3d.dirichlet.get_quadrature_correction(S,epsqclh,lambda, ...
        [1.0,0.0],S,optsqclh);
end
end

% "Matrix" A in GMRES 
function A = A(s,dom,S,zk,eps,Q)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = mtxBsigma(S,dom,sigma,zk,eps,S,Q);
A = surfacefun_to_array(Bsigma,dom,S);
end

% Solve A11*D = A12
function [D,dfunc] = solveforD(S,dom,mH,lambda,eps,Q)
Balpha = mtxBalpha(S,dom,mH,lambda,eps,S,Q);
b = surfacefun_to_array(Balpha,dom,S);
t1 = tic;
[D,~,~,iter] = gmres(@(s) A(s,dom,S,lambda,eps,Q),b,[],eps,50);
t2 = toc(t1);
fprintf('GMRES for A11*D = A12: %f s / %d iter. = %f s\n', ...
    t2, iter(2), t2/iter(2))
dfunc = array_to_surfacefun(D,dom,S);
end

% Compute B on-surface
function B = surfaceB(sigma,m,S,dom,vn,zk,eps,Q,optsqc,optsqclh)
sigmavals = surfacefun_to_array(sigma,dom,S);
mvals = surfacefun_to_array(m,dom,S);
if abs(zk) < 2.5e-16
    gradSksigma = taylor.static.eval_gradS0(S,sigmavals.',eps,S,Q,optsqc);
    gradSksigma = array_to_surfacefun(gradSksigma.',dom,S);
    curlSkm = taylor.static.eval_curlS0(S,mvals.',eps,S,Q,optsqc);
    curlSkm = array_to_surfacefun(curlSkm.',dom,S);

    B = -sigma.*vn./2 + m./2 - gradSksigma + 1i.*curlSkm;
else
    gradSksigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals.',eps, ...
        S,Q,optsqc);
    gradSksigma = array_to_surfacefun(gradSksigma.',dom,S);
    curlSkm = taylor.dynamic.eval_curlSk(S,zk,mvals.',eps,S,Q,optsqc);
    curlSkm = array_to_surfacefun(curlSkm.',dom,S);

    dpars = [1.0, 0.0];
    Qlh = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,dpars, ...
        S,optsqclh);

    optslh = [];
    optslh.format = 'rsc';
    optslh.precomp_quadrature = Qlh;
    Skm1 = helm3d.dirichlet.eval(S,mvals(:,1),S,eps,zk,dpars,optslh);
    Skm2 = helm3d.dirichlet.eval(S,mvals(:,2),S,eps,zk,dpars,optslh);
    Skm3 = helm3d.dirichlet.eval(S,mvals(:,3),S,eps,zk,dpars,optslh);
    Skm = [Skm1 Skm2 Skm3];
    Skm = array_to_surfacefun(Skm,dom,S);

    B = -sigma.*vn./2 + m./2 + 1i.*zk.*Skm - gradSksigma + 1i.*curlSkm;
end
end