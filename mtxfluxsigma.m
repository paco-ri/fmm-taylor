function fluxsigma = mtxfluxsigma(S,dom,n,nu,nv,zpars,sigma,eps,Q)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
%   Detailed explanation goes here

lambda = zpars(1);
m0 = debyem0(sigma,lambda);
m0vals = surfacefun_to_array(m0,dom,S);

% single layer potential evals
Skm01 = helm3d.dirichlet.eval(S,zpars,m0vals(:,1),eps,S,Q);
Skm02 = helm3d.dirichlet.eval(S,zpars,m0vals(:,2),eps,S,Q);
Skm03 = helm3d.dirichlet.eval(S,zpars,m0vals(:,3),eps,S,Q);
Skm0 = [Skm01 Skm02 Skm03];
Skm0 = array_to_surfacefun(Skm0,dom,S);

vn = normal(dom);
nxm0 = cross(vn,m0);

% derivatives of layer potentials
curlSkm0 = curl(Skm0);

integrand = surfacefunv(dom);
for j = 1:3
    integrand.components{j} = 1i*(Skm0.components{j} ...
        + (nxm0.components{j}/2 + curlSkm0.components{j})/lambda);
end
fluxsigma = intacyc(integrand,n,nu,nv);

end