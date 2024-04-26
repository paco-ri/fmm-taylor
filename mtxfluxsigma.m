function fluxsigma = mtxfluxsigma(S,dom,n,nu,nv,zpars,sigma,eps,Q)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
%   Detailed explanation goes here

lambda = zpars(1);
m0 = debyem0(sigma,lambda);
m0vals = surfacefun_to_array(m0,dom,S);
m0vals = m0vals';

% single layer potential evals
S0m01 = helm3d.dirichlet.eval(S,zpars,m0vals(1,:),eps,S,Q);
S0m02 = helm3d.dirichlet.eval(S,zpars,m0vals(2,:),eps,S,Q);
S0m03 = helm3d.dirichlet.eval(S,zpars,m0vals(3,:),eps,S,Q);
S0m0 = [S0m01 S0m02 S0m03];
S0m0 = array_to_surfacefun(S0m0,dom,S);

vn = normal(dom);
nxm0 = cross(vn,m0);

% derivatives of layer potentials
[curlS0m0,~] = virtualcasing.evalgradcurlS0(S,m0vals,m0vals(1,:),eps,S,Q);
curlS0m0 = array_to_surfacefun(curlS0m0',dom,S); % note transpose

integrand = surfacefunv(dom);
for j = 1:3
    integrand.components{j} = 1i*(S0m0.components{j} ...
        + (nxm0.components{j}/2 + curlS0m0.components{j})/lambda);
end
fluxsigma = intacyc(integrand,n,nu,nv);

end