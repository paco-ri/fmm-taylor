clear

% define surface 
n = 4; % polynomial order 
nu = 6; 
nv = nu*3;
% dom = surfacemesh.torus(n, nu, nv);
dom = circulartorus(n,nu,nv);

% get harmonic surface vector field 
x0 = 0.2; y0 = 0.2; z0 = 0.2;
denom = @(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2);
V = cross([0 1 1], surfacefunv(@(x,y,z) (x-x0)./denom(x,y,z), ...
                               @(x,y,z) (y-y0)./denom(x,y,z), ...
                               @(x,y,z) (z-z0)./denom(x,y,z), dom));

vn = normal(dom);
f = -cross(vn, vn, V);

[u, v, mH, curlfree, divfree] = hodge(f);

% set wavenumber and zpars for single-layer potential eval.
lambda = 1e-8;
zpars = complex([lambda, 1.0, 0.0]);

% compute reference Taylor state
rmaj = 5.0;
rmin = 3.0;
ntheta = 1e3;
jmag = 1.0;

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

nB0 = dot(vn,B0);

% convert surfacemesh dom to surfer
[S,domsurfer] = surfer.surfacemesh_to_surfer(dom);

% compute near quadrature correction
eps = 1e-8;
opts_quad = [];
opts_quad.format='rsc';
Q = helm3d.dirichlet.get_quadrature_correction(S,zpars,eps,S,opts_quad);

% do GMRES to solve A11*W = nB0
b = surfacefun_to_array(nB0,dom,S); 
b = -b;
disp('start gmres 1')
[W, flag1, relres1, iter1] = gmres(@(s) A(s,dom,S,zpars,eps,Q),b,[],1e-12,50);
disp('end gmres 1')
wfunc = array_to_surfacefun(W,dom,S);

% do GMRES to solve A11*D = A12
Balpha = mtxBalpha(S,dom,zpars,mH,eps,Q);
b = surfacefun_to_array(Balpha,dom,S);
b = -b; 
disp('start gmres 2')
[D, flag2, relres2, iter2] = gmres(@(s) A(s,dom,S,zpars,eps,Q),b,[],1e-12,50);
disp('end gmres 2')
dfunc = array_to_surfacefun(D,dom,S);

% compute flux 
flux = intacyc(B0,n,nu,nv)/lambda;
fluxsigmaD = mtxfluxsigma(S,dom,n,nu,nv,zpars,dfunc,eps,Q);
fluxsigmaW = mtxfluxsigma(S,dom,n,nu,nv,zpars,wfunc,eps,Q);
fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,zpars,mH,eps,Q);

% compute alpha and sigma
% alpha = flux/(fluxsigmaD + fluxalpha);
% alpha = (flux - fluxsigmaW)/(fluxsigmaD + fluxalpha); % <-- OG
alpha = (flux + fluxsigmaW)/(fluxsigmaD + fluxalpha);
% sigma = alpha*dfunc;
% sigma = wfunc - alpha*dfunc; % <-- OG
sigma = alpha*dfunc - wfunc;

% compute B
m0 = debyem0(sigma,lambda);
m = m0 + mH./(1/alpha); 

mvals = surfacefun_to_array(m,dom,S);
Skm1 = helm3d.dirichlet.eval(S,zpars,mvals(:,1),eps,S,Q);
Skm2 = helm3d.dirichlet.eval(S,zpars,mvals(:,2),eps,S,Q);
Skm3 = helm3d.dirichlet.eval(S,zpars,mvals(:,3),eps,S,Q);
Skm = [Skm1 Skm2 Skm3];
Skm = array_to_surfacefun(Skm,dom,S);

sigmavals = surfacefun_to_array(sigma,dom,S);
Sksigma = helm3d.dirichlet.eval(S,zpars,sigmavals,eps,S,Q);
Sksigma = array_to_surfacefun(Sksigma,dom,S);

% NEED TO OVERRIDE times()
B = surfacefunv(dom);
nxm = cross(vn,m);
gradSksigma = grad(Sksigma);
curlSkm = curl(Skm);
for k = 1:3
    B.components{k} = vn.components{k}*sigma/(-2) + 1i*nxm.components{k}/2 ...
        + 1i*lambda*Skm.components{k} - gradSksigma.components{k} ...
        + 1i*curlSkm.components{k};
end

plot(norm(B0 - B)./norm(B0))
hold on
colorbar

% function in solve for D
function A = A(s,dom,S,zpars,eps,Q)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = mtxBsigma(S,dom,zpars,sigma,eps,Q);
A = surfacefun_to_array(Bsigma,dom,S);
end

