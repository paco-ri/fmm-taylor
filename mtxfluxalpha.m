function fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,zpars,mH,eps,Q)
%UNTITLED3 compute alpha coefficient in flux condition
%   Detailed explanation goes here

lambda = zpars(1);
mHvals = surfacefun_to_array(mH,dom,S);

% single layer potential evals
SkmH1 = helm3d.dirichlet.eval(S,zpars,mHvals(:,1),eps,S,Q);
SkmH2 = helm3d.dirichlet.eval(S,zpars,mHvals(:,2),eps,S,Q);
SkmH3 = helm3d.dirichlet.eval(S,zpars,mHvals(:,3),eps,S,Q);
SkmH = [SkmH1 SkmH2 SkmH3];
SkmH = array_to_surfacefun(SkmH,dom,S);

vn = normal(dom);
nxmH = cross(vn,mH);

% derivative of a layer potential
curlSkmH = curl(SkmH);

integrand = surfacefunv(dom);
for j = 1:3
    integrand.components{j} = 1i*(SkmH.components{j} ...
        + (nxmH.components{j}/2 + curlSkmH.components{j})/lambda);
end
fluxalpha = intacyc(integrand,n,nu,nv);

end