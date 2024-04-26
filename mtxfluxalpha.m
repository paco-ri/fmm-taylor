function fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,zpars,mH,eps,Q)
%MTXFLUXALPHA compute alpha coefficient in flux condition
%   Detailed explanation goes here

lambda = zpars(1);
mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

% single layer potential evals
S0mH1 = helm3d.dirichlet.eval(S,zpars,mHvals(1,:),eps,S,Q);
S0mH2 = helm3d.dirichlet.eval(S,zpars,mHvals(2,:),eps,S,Q);
S0mH3 = helm3d.dirichlet.eval(S,zpars,mHvals(3,:),eps,S,Q);
S0mH = [S0mH1 S0mH2 S0mH3];
S0mH = array_to_surfacefun(S0mH,dom,S);

vn = normal(dom);
nxmH = cross(vn,mH);

% derivative of layer potential
[curlS0mH,~] = virtualcasing.evalgradcurlS0(S,mHvals,mHvals(1,:),eps,S,Q);
curlS0mH = array_to_surfacefun(curlS0mH',dom,S);

integrand = surfacefunv(dom);
for j = 1:3
    integrand.components{j} = 1i*(S0mH.components{j} ...
        + (nxmH.components{j}/2 + curlS0mH.components{j})/lambda);
end
fluxalpha = intacyc(integrand,n,nu,nv);

end