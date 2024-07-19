% test \nabla_\Gamma \cdot m = i \lambda \sigma

ns = 4:2:12;
nus = 8:4:20;
numns = size(ns,2);
numnus = size(nus,2);
errs = zeros([4 numns*numnus]);

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

ii = 1;
for n = ns
for nu = nus

% define surface 
nv = nu*3;
% fprintf('n = %d, nu = %d, nv = %d\n',n,nu,nv)

if whichgeom == 1
    dom = circulartorus(n,nu,nv,domrmin,domrmaj);
else
    dom = surfacemesh.torus(n, nu, nv);
end

vn = normal(dom);

% get harmonic surface vector field 
sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
sigma = surfacefun(@(x,y,z) sinphi(x,y,z)+cosphi(x,y,z)+z,dom);
sigmao = resample(sigma,n*2);

m0 = debyem0(sigma,lambda);
m0err = div(m0) - 1i*lambda.*sigma;

m0o = debyem0(sigmao,lambda);
m0oerr = div(m0o) - 1i*lambda.*sigmao;

errs(1,ii) = n;
errs(2,ii) = nu;
errs(3,ii) = norm(m0err);
errs(4,ii) = norm(m0oerr);

ii = ii + 1;
end
end

figure(1)
semilogy(sqrt(errs(1,1:numnus).*errs(2,1:numnus)), errs(3,1:numnus), 'o-')
hold on
for j = 1:(numns-1)
    jj = (numnus*j+1):(numnus*(j+1));
    semilogy(sqrt(errs(1,jj).*errs(2,jj)), errs(3,jj), 'o-')
end

figure(2)
semilogy(sqrt(errs(1,1:numnus).*errs(2,1:numnus)), errs(4,1:numnus), 'o-')
hold on
for j = 1:(numns-1)
    jj = (numnus*j+1):(numnus*(j+1));
    semilogy(sqrt(errs(1,jj).*errs(2,jj)), errs(4,jj), 'o-')
end

% nsph = 8;
% nref = 3;
% dom = surfacemesh.sphere(nsph,nref);
% S = surfer.surfacemesh_to_surfer(dom);
% [srcvals,~,~,~,~,wts] = extract_arrays(S);
% vn = normal(dom);
% 
% ndeg = 3; 
% mdeg = 1;
% 
% f = spherefun.sphharm(ndeg,mdeg);
% sigma = surfacefun(@(x,y,z) f(x,y,z),dom);
% sigvals = surfacefun_to_array(sigma,dom,S);
% 
% lambda = 1.0;
% m0 = debyem0(sigma,lambda);
% norm(div(m0)-1i*lambda.*sigma)
