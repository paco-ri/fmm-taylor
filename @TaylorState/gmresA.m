function A = gmresA(s,dom,S,L,zk,epstaylor,epslh,opts,optslh)

% cell arrays: S, dom, opts, optslh
Bsigma = TaylorState.mtxBsigma(S,dom,L,s,zk,epstaylor, ...
    epslh,S,opts,optslh);

if length(dom) == 1
    % A = surfacefun_to_array(Bsigma{1},dom{1},S{1});
    A = surfacefun_to_array(Bsigma{1,1},dom{1},S{1});
else
    Ao = surfacefun_to_array(Bsigma{1},dom{1},S{1});
    Ai = surfacefun_to_array(Bsigma{2},dom{2},S{2});
    A = [Ao; Ai];
    % A11 = surfacefun_to_array(Bsigma{1,1},dom{1},S{1});
    % A12 = surfacefun_to_array(Bsigma{1,2},dom{1},S{1});
    % A21 = surfacefun_to_array(Bsigma{2,1},dom{2},S{2});
    % A22 = surfacefun_to_array(Bsigma{1,1},dom{2},S{2});
    % A = [A11 A12; A21 A22];
    % A = A(:,i);
end

end