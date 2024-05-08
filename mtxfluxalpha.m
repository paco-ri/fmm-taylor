function fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,mH,eps,Q)
%MTXFLUXALPHA compute alpha coefficient in flux condition
%   Detailed explanation goes here

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

% derivative of layer potential
[curlS0mH,~] = virtualcasing.evalgradcurlS0(S,mHvals,mHvals(1,:),eps,S,Q);
curlS0mH = array_to_surfacefun(curlS0mH',dom,S);
vn = normal(dom);
nxcurlS0mH = curl(vn,curlS0mH);
nxminusmH = nxcurlS0mH - mH/2;
nxminusmHvals = surfacefun_to_array(nxminusmH,dom,S);

% single layer potential evals
zpars = complex([1.0,0.0]);
S0nx1 = lap3d.dirichlet.eval(S,zpars,nxminusmHvals(1,:),eps,S,Q);
S0nx2 = lap3d.dirichlet.eval(S,zpars,nxminusmHvals(2,:),eps,S,Q);
S0nx3 = lap3d.dirichlet.eval(S,zpars,nxminusmHvals(3,:),eps,S,Q);
S0nx = [S0nx1 S0nx2 S0nx3];
S0nx = array_to_surfacefun(S0nx,dom,S);

fluxalpha = intacyc(S0nx,n,nu,nv);
fluxalpha = 1i*fluxalpha;

end