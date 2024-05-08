function Balpha = mtxBalpha(S,dom,mH,eps,Q)
%MTXBALPHA compute alpha coefficient for surface magnetic field
%   Detailed explanation goes here

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

% derivatives of layer potentials 
n = normal(dom);
[curlS0mH,~] = virtualcasing.evalgradcurlS0(S,mHvals,mHvals(1,:),eps,S,Q);
curlS0mH = array_to_surfacefun(curlS0mH',dom,S);
ncurlS0mH = dot(n,curlS0mH);

% construct Balpha
Balpha = 1i*ncurlS0mH;

end