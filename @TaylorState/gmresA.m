function A = gmresA(s,dom,io,S,zk,epstaylor,epslh,opts,optslh)
sigma = array_to_surfacefun(s,dom,S);
Bsigma = TaylorState.mtxBsigma(S,dom,io,sigma,zk,epstaylor,epslh,S,opts,optslh);
A = surfacefun_to_array(Bsigma,dom,S);
end