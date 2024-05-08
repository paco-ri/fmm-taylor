function fluxsigma = mtxfluxsigma(S,dom,n,nu,nv,sigma,eps,Q)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
%   Detailed explanation goes here

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals';

% derivatives of layer potentials
dummyvals = repmat(sigmavals,3,1);
[~,gradS0sigma] = virtualcasing.evalgradcurlS0(S,dummyvals,sigma,eps,S,Q);
gradS0sigma = array_to_surfacefun(gradS0sigma',dom,S); % note transpose
vn = normal(dom);
nxgradS0sigma = cross(vn,gradS0sigma);
nxvals = surfacefun_to_array(nxgradS0sigma,dom,S);

% single layer potential evals
zpars = complex([1.0,0.0]);
S0nx1 = lap3d.dirichlet.eval(S,zpars,nxvals(1,:),eps,S,Q);
S0nx2 = lap3d.dirichlet.eval(S,zpars,nxvals(2,:),eps,S,Q);
S0nx3 = lap3d.dirichlet.eval(S,zpars,nxvals(3,:),eps,S,Q);
S0nx = [S0nx1 S0nx2 S0nx3];
S0nx = array_to_surfacefun(S0nx,dom,S);

fluxsigma = intacyc(S0nx,n,nu,nv);
fluxsigma = -fluxsigma;

end