function Balpha = mtxBalpha(S,dom,zpars,mH,eps,Q)
%MTXBALPHA compute alpha coefficient for surface magnetic field
%   Detailed explanation goes here

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

% single layer potential evals
S0mH1 = helm3d.dirichlet.eval(S,zpars,mHvals(1,:),eps,S,Q);
S0mH2 = helm3d.dirichlet.eval(S,zpars,mHvals(2,:),eps,S,Q);
S0mH3 = helm3d.dirichlet.eval(S,zpars,mHvals(3,:),eps,S,Q);
S0mH = [S0mH1 S0mH2 S0mH3];
S0mH = array_to_surfacefun(S0mH,dom,S);

% derivatives of layer potentials 
n = normal(dom);
nS0mH = dot(n,S0mH);
[curlS0mH,~] = virtualcasing.evalgradcurlS0(S,mHvals,mHvals(1,:),eps,S,Q);
curlS0mH = array_to_surfacefun(curlS0mH',dom,S);
ncurlS0mH = dot(n,curlS0mH);

% construct Balpha
lambda = zpars(1);
Balpha = 1i*(lambda*nS0mH + ncurlS0mH);

end