% define surface 
% n = 4; % polynomial order 
% nu = 6; 
% nv = nu*3;

% small problem
n = 3;
nu = 3;
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

[u, v, vH, curlfree, divfree] = hodge(f);
mH = vH + cross(vn,vH)./(-1i);

% set dpars for single-layer potential eval.
dpars = [1.0, 0.0];

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
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,0,...
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
eps = 1e-4;
opts_quad = [];
opts_quad.format='rsc';
Q = taylor.static.get_quadrature_correction(S,eps,S,opts_quad);

% do GMRES to solve A11*W = nB0
b = surfacefun_to_array(nB0,dom,S); 
disp('start gmres 1')
[W, flag1, relres1, iter1] = gmres(@(s) A(s,dom,S,eps),b,[],eps,50);
disp('end gmres 1')
wfunc = array_to_surfacefun(W,dom,S);

% do GMRES to solve A11*D = A12
Balpha = mtxBalpha(S,dom,mH,eps,S,Q);
b = surfacefun_to_array(Balpha,dom,S);
disp('start gmres 2')
[D, flag2, relres2, iter2] = gmres(@(s) A(s,dom,S,eps),b,[],eps,50);
disp('end gmres 2')
dfunc = array_to_surfacefun(D,dom,S);

Qlap = lap3d.dirichlet.get_quadrature_correction(S,dpars,eps,S,opts_quad);

% compute flux (see (27))
j0 = cross(vn,B0);
j0vals = surfacefun_to_array(j0,dom,S);
S0j01 = lap3d.dirichlet.eval(S,dpars,j0vals(:,1),eps,S,Qlap);
S0j02 = lap3d.dirichlet.eval(S,dpars,j0vals(:,2),eps,S,Qlap);
S0j03 = lap3d.dirichlet.eval(S,dpars,j0vals(:,3),eps,S,Qlap);
S0j0 = [S0j01 S0j02 S0j03];
S0j0 = array_to_surfacefun(S0j0,dom,S);
flux = intacyc(S0j0,n,nu,nv);

fluxsigmaD = mtxfluxsigma(S,dom,n,nu,nv,dfunc,eps,S,Q);
fluxsigmaW = mtxfluxsigma(S,dom,n,nu,nv,wfunc,eps,S,Q);
fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,mH,eps,S,Q);

% compute alpha and sigma
% alpha = flux/(fluxsigmaD + fluxalpha);
alpha = -1i*(flux - fluxsigmaW)/(-fluxsigmaD + fluxalpha); % <-- OG
% alpha = (flux + fluxsigmaW)/(fluxsigmaD + fluxalpha);
% sigma = alpha*dfunc;
% sigma = wfunc - alpha*dfunc; % <-- OG
sigma = 1i*alpha*dfunc - wfunc;

% compute B
m0 = 0; % will later be call to debyem0
m = m0 + mH./(1/alpha); 

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';
gradS0sigma = taylor.static.eval_gradS0(S,sigmavals,eps,S,Q);

mvals = surfacefun_to_array(m,dom,S);
mvals = mvals.';
curlS0m = taylor.static.eval_curlS0(S,mvals,eps,S,Q);

% NEED TO OVERRIDE times()
vnvals = surfacefun_to_array(vn,dom,S);
vnvals = vnvals.';
B = -vnvals.*sigmavals./2 + mvals./2 - gradS0sigma + 1i.*curlS0m;
B = array_to_surfacefun(B.',dom,S);

plot(norm(B0 - B)./norm(B0))
hold on
colorbar

% function in solve for D
function A = A(s,dom,S,eps)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = mtxBsigma(S,dom,sigma,eps);
A = surfacefun_to_array(Bsigma,dom,S);
end

