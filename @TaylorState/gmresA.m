function A = gmresA(s,dom,S,zk,epstaylor,epslh,opts,optslh)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = TaylorState.mtxBsigma(S,dom,sigma,zk,epstaylor,epslh,S,opts,optslh);
A = surfacefun_to_array(Bsigma,dom,S);
end